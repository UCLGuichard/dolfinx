// Copyright (C) 2006-2017 Anders Logg and Garth N. Wells
//
// This file is part of DOLFINX (https://www.fenicsproject.org)
//
// SPDX-License-Identifier:    LGPL-3.0-or-later

#include "TopologyComputation.h"
#include "Connectivity.h"
#include "Mesh.h"
#include "MeshEntity.h"
#include "MeshIterator.h"
#include "Topology.h"
#include "cell_types.h"
#include <Eigen/Dense>
#include <algorithm>
#include <boost/unordered_map.hpp>
#include <cstdint>
#include <dolfinx/common/MPI.h>
#include <dolfinx/common/Timer.h>
#include <dolfinx/common/log.h>
#include <dolfinx/common/utils.h>
#include <memory>
#include <numeric>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using namespace dolfinx;
using namespace dolfinx::mesh;

namespace
{

// Takes an Eigen::Array and obtains the sort permutation to reorder the
// rows in ascending order. Each row must be sorted beforehand.
template <typename T>
std::vector<int> sort_by_perm(
    const Eigen::Ref<const Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic,
                                        Eigen::RowMajor>>& arr_data)
{
  // Sort an Eigen::Array by creating a permutation vector
  std::vector<int> index(arr_data.rows());
  std::iota(index.begin(), index.end(), 0);
  const int cols = arr_data.cols();

  // Lambda with capture for sort comparison
  const auto cmp = [&arr_data, &cols](int a, int b) {
    const T* row_a = arr_data.row(a).data();
    const T* row_b = arr_data.row(b).data();
    return std::lexicographical_compare(row_a, row_a + cols, row_b,
                                        row_b + cols);
  };

  std::sort(index.begin(), index.end(), cmp);
  return index;
}
//-----------------------------------------------------------------------------
std::tuple<std::shared_ptr<Connectivity>, std::shared_ptr<Connectivity>,
           std::int32_t>
compute_entities_by_key_matching_new(const Mesh& mesh, int dim)
{
  if (dim == 0)
  {
    throw std::runtime_error(
        "Cannot create vertices for topology. Should already exist.");
  }

  // Get mesh topology and connectivity
  const Topology& topology = mesh.topology();
  const int tdim = topology.dim();

  // Check if entities have already been computed
  if (topology.connectivity(dim, 0))
  {
    // Check that we have cell-entity connectivity
    if (!topology.connectivity(tdim, dim))
      throw std::runtime_error("Missing cell-entity connectivity");

    return {nullptr, nullptr, topology.size(dim)};
  }

  // Start timer
  common::Timer timer("Compute entities of dim = " + std::to_string(dim));

  // Initialize local array of entities
  const std::int8_t num_entities
      = mesh::cell_num_entities(mesh.cell_type(), dim);
  const int num_vertices
      = mesh::num_cell_vertices(mesh::cell_entity_type(mesh.cell_type(), dim));

  // Create map from cell vertices to entity vertices
  Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> e_vertices
      = mesh::get_entity_vertices(mesh.cell_type(), dim);

  // List of vertices for each entity in each cell.
  Eigen::Array<std::int32_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
      entity_list(mesh.num_entities(tdim) * num_entities, num_vertices);

  int k = 0;
  for (auto& c : MeshRange(mesh, tdim, MeshRangeType::ALL))
  {
    // Get vertices from cell
    const std::int32_t* vertices = c.entities(0);
    assert(vertices);

    // Iterate over entities of cell
    for (int i = 0; i < num_entities; ++i)
    {
      // Get entity vertices
      for (int j = 0; j < num_vertices; ++j)
        entity_list(k, j) = vertices[e_vertices(i, j)];

      ++k;
    }
  }
  assert(k == entity_list.rows());

  std::vector<std::int32_t> entity_index(entity_list.rows());
  std::int32_t entity_count = 0;

  // Copy list and sort vertices of each entity into order
  Eigen::Array<std::int32_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
      entity_list_sorted = entity_list;
  for (int i = 0; i < entity_list_sorted.rows(); ++i)
    std::sort(entity_list_sorted.row(i).data(),
              entity_list_sorted.row(i).data() + num_vertices);

  // Sort the list and label (first pass)
  std::vector<std::int32_t> sort_order
      = sort_by_perm<std::int32_t>(entity_list_sorted);
  std::int32_t last = sort_order[0];
  entity_index[last] = 0;
  for (std::size_t i = 1; i < sort_order.size(); ++i)
  {
    std::int32_t j = sort_order[i];
    if ((entity_list_sorted.row(j) != entity_list_sorted.row(last)).any())
      ++entity_count;
    entity_index[j] = entity_count;
    last = j;
  }
  ++entity_count;

  // FIXME: Need to find ghosts, so we can put at end of range

  // Get a single row in entity list for each entity
  std::vector<std::int32_t> unique_row(entity_count);
  for (int i = 0; i < entity_list.rows(); ++i)
    unique_row[entity_index[i]] = i;

  int mpi_size = dolfinx::MPI::size(mesh.mpi_comm());
  std::vector<std::vector<std::int64_t>> send_entities(mpi_size);
  std::vector<std::vector<std::int64_t>> send_index(mpi_size);
  std::vector<std::vector<std::int64_t>> recv_entities(mpi_size);

  // Get all "possibly shared" entities, based on vertex sharing
  // Send to other processes, and see if we get the same back
  const std::map<std::int32_t, std::set<std::int32_t>>& shared_vertices
      = mesh.topology().shared_entities(0);
  const std::vector<std::int64_t>& global_vertex_indices
      = mesh.topology().global_indices(0);

  std::stringstream s;

  std::map<int, int> procs;
  for (int i : unique_row)
  {
    procs.clear();
    for (int j = 0; j < num_vertices; ++j)
    {
      const int v = entity_list_sorted(i, j);
      const auto it = shared_vertices.find(v);
      if (it != shared_vertices.end())
        for (std::int32_t p : it->second)
          ++procs[p];
    }
    for (const auto& q : procs)
      if (q.second == num_vertices)
      {
        const int p = q.first;

        s << "Entity " << entity_index[i] << " may be shared with process " << p
          << "\n";
        for (int j = 0; j < num_vertices; ++j)
        {
          std::int64_t vglobal = global_vertex_indices[entity_list(i, j)];
          send_entities[p].push_back(vglobal);
        }
        send_index[p].push_back(entity_index[i]);
        std::sort(send_entities[p].end() - num_vertices,
                  send_entities[p].end());
      }
  }

  dolfinx::MPI::all_to_all(mesh.mpi_comm(), send_entities, recv_entities);

  // Compare received with sent for each process
  for (std::size_t p = 0; p < send_entities.size(); ++p)
  {
    const std::vector<std::int64_t>& sendp = send_entities[p];
    const std::vector<std::int64_t>& recvp = recv_entities[p];
    std::set<std::vector<std::int64_t>> recv_set;
    for (std::size_t i = 0; i < recvp.size() / num_vertices; ++i)
    {
      recv_set.insert(
          std::vector<std::int64_t>(recvp.begin() + i * num_vertices,
                                    recvp.begin() + (i + 1) * num_vertices));
    }
    for (std::size_t i = 0; i < sendp.size() / num_vertices; ++i)
    {
      const std::vector<std::int64_t> b(sendp.begin() + i * num_vertices,
                                        sendp.begin() + (i + 1) * num_vertices);
      if (recv_set.find(b) == recv_set.end())
      {
        s << "Sent, but did not receive back - " << send_index[p][i]
          << " not shared\n";
      }
      else
      {
        s << "entity " << send_index[p][i] << " shared with " << p << "\n";
        // shared with process p - decide owner
      }
    }
  }

  std::cout << s.str() << "\n";

  Eigen::Array<std::int32_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
      connectivity_ce(mesh.num_entities(tdim), num_entities);
  std::copy(entity_index.begin(), entity_index.end(), connectivity_ce.data());

  // Cell-entity connectivity
  auto ce = std::make_shared<Connectivity>(connectivity_ce);

  // Entity-vertex connectivity
  Eigen::Array<std::int32_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
      connectivity_ev(entity_count, num_vertices);
  for (int i = 0; i < entity_list.rows(); ++i)
    connectivity_ev.row(entity_index[i]) = entity_list.row(i);

  auto ev = std::make_shared<Connectivity>(connectivity_ev);

  return {ce, ev, entity_count};
}
//-----------------------------------------------------------------------------
// Compute connectivity from transpose
Connectivity compute_from_transpose(const Mesh& mesh, int d0, int d1)
{
  // The transpose is computed in three steps:
  //
  //   1. Iterate over entities of dimension d1 and count the number
  //      of connections for each entity of dimension d0
  //
  //   2. Allocate memory / prepare data structures
  //
  //   3. Iterate again over entities of dimension d1 and add connections
  //      for each entity of dimension d0

  LOG(INFO) << "Computing mesh connectivity " << d0 << " - " << d1
            << "from transpose.";

  // Get mesh topology and connectivity
  const Topology& topology = mesh.topology();

  // Need connectivity d1 - d0
  if (!topology.connectivity(d1, d0))
    throw std::runtime_error("Missing required connectivity d1-d0.");

  // Compute number of connections for each e0
  std::vector<std::int32_t> num_connections(topology.size(d0), 0);
  for (auto& e1 : MeshRange(mesh, d1, MeshRangeType::ALL))
    for (auto& e0 : EntityRange(e1, d0))
      num_connections[e0.index()]++;

  // Compute offsets
  std::vector<std::int32_t> offsets(num_connections.size() + 1, 0);
  std::partial_sum(num_connections.begin(), num_connections.end(),
                   offsets.begin() + 1);

  std::vector<std::int32_t> counter(num_connections.size(), 0);
  std::vector<std::int32_t> connections(offsets.back());
  for (auto& e1 : MeshRange(mesh, d1, MeshRangeType::ALL))
    for (auto& e0 : EntityRange(e1, d0))
      connections[offsets[e0.index()] + counter[e0.index()]++] = e1.index();

  return Connectivity(connections, offsets);
}
//-----------------------------------------------------------------------------
// Direct lookup of entity from vertices in a map
Connectivity compute_from_map(const Mesh& mesh, int d0, int d1)
{
  assert(d1 > 0);
  assert(d0 > d1);

  // Get the type of entity d0
  mesh::CellType cell_type = mesh::cell_entity_type(mesh.cell_type(), d0);

  // Make a map from the sorted d1 entity vertices to the d1 entity index
  boost::unordered_map<std::vector<std::int32_t>, std::int32_t> entity_to_index;
  entity_to_index.reserve(mesh.num_entities(d1));

  const std::size_t num_verts_d1
      = mesh::num_cell_vertices(mesh::cell_entity_type(mesh.cell_type(), d1));

  std::vector<std::int32_t> key(num_verts_d1);
  for (auto& e : MeshRange(mesh, d1, MeshRangeType::ALL))
  {
    std::partial_sort_copy(e.entities(0), e.entities(0) + num_verts_d1,
                           key.begin(), key.end());
    entity_to_index.insert({key, e.index()});
  }

  Eigen::Array<std::int32_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
      connections(mesh.num_entities(d0),
                  mesh::cell_num_entities(cell_type, d1));

  // Search for d1 entities of d0 in map, and recover index
  std::vector<std::int32_t> entities;
  const Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
      e_vertices_ref = mesh::get_entity_vertices(cell_type, d1);
  Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> keys
      = e_vertices_ref;
  for (auto& e : MeshRange(mesh, d0, MeshRangeType::ALL))
  {
    entities.clear();
    const std::int32_t* e0 = e.entities(0);
    for (Eigen::Index i = 0; i < e_vertices_ref.rows(); ++i)
      for (Eigen::Index j = 0; j < e_vertices_ref.cols(); ++j)
        keys(i, j) = e0[e_vertices_ref(i, j)];
    for (Eigen::Index i = 0; i < keys.rows(); ++i)
    {
      std::partial_sort_copy(keys.row(i).data(),
                             keys.row(i).data() + keys.row(i).cols(),
                             key.begin(), key.end());
      const auto it = entity_to_index.find(key);
      assert(it != entity_to_index.end());
      entities.push_back(it->second);
    }
    for (std::size_t k = 0; k < entities.size(); ++k)
      connections(e.index(), k) = entities[k];
  }

  return Connectivity(connections);
}
} // namespace

//-----------------------------------------------------------------------------
void TopologyComputation::compute_entities(Mesh& mesh, int dim)
{
  LOG(INFO) << "Computing mesh entities of dimension " << dim;

  // Check if entities have already been computed
  Topology& topology = mesh.topology();

  // Vertices must always exist
  if (dim == 0)
    return;

  if (topology.connectivity(dim, 0))
  {
    // Make sure we really have the connectivity
    if (!topology.connectivity(topology.dim(), dim))
    {
      throw std::runtime_error(
          "Cannot compute topological entities. Entities of topological "
          "dimension "
          + std::to_string(dim) + " exist but connectivity is missing.");
    }
    return;
  }

  std::tuple<std::shared_ptr<Connectivity>, std::shared_ptr<Connectivity>,
             std::int32_t>
      data = compute_entities_by_key_matching_new(mesh, dim);

  // Set cell-entity connectivity
  if (std::get<0>(data))
    topology.set_connectivity(std::get<0>(data), topology.dim(), dim);

  // Set entity-vertex connectivity
  if (std::get<1>(data))
    topology.set_connectivity(std::get<1>(data), dim, 0);

  // Initialise ghost entity offset
  topology.init_ghost(dim, std::get<2>(data));
}
//-----------------------------------------------------------------------------
void TopologyComputation::compute_connectivity(Mesh& mesh, int d0, int d1)
{
  // This is where all the logic takes place to find a strategy for
  // the connectivity computation. For any given pair (d0, d1), the
  // connectivity is computed by suitably combining the following
  // basic building blocks:
  //
  //   1. compute_entities():     d  - 0  from dim - 0
  //   2. compute_transpose():    d0 - d1 from d1 - d0
  //   3. compute_intersection(): d0 - d1 from d0 - d' - d1
  //   4. compute_from_map():     d0 - d1 from d1 - 0 and d0 - 0
  // Each of these functions assume a set of preconditions that we
  // need to satisfy.

  LOG(INFO) << "Requesting connectivity " << d0 << " - " << d1;

  // Get mesh topology and connectivity
  Topology& topology = mesh.topology();

  // Return if connectivity has already been computed
  if (topology.connectivity(d0, d1))
    return;

  // Make sure entities exist
  if (d0 > 0)
    assert(topology.connectivity(d0, 0));
  if (d1 > 0)
    assert(topology.connectivity(d1, 0));

  // Start timer
  common::Timer timer("Compute connectivity " + std::to_string(d0) + "-"
                      + std::to_string(d1));

  // Decide how to compute the connectivity
  if (d0 == d1)
  {
    // For d0-d1, use indentity connecticity
    std::vector<std::vector<std::size_t>> connectivity_dd(
        topology.size(d0), std::vector<std::size_t>(1));
    for (auto& e : MeshRange(mesh, d0, MeshRangeType::ALL))
      connectivity_dd[e.index()][0] = e.index();
    auto connectivity = std::make_shared<Connectivity>(connectivity_dd);
    topology.set_connectivity(connectivity, d0, d1);
  }
  else if (d0 < d1)
  {
    // Compute connectivity d1 - d0 and take transpose
    compute_connectivity(mesh, d1, d0);
    auto c
        = std::make_shared<Connectivity>(compute_from_transpose(mesh, d0, d1));
    topology.set_connectivity(c, d0, d1);
  }
  else if (d0 > d1)
  {
    // Compute by mapping vertices from a lower dimension entity to
    // those of a higher dimension entity
    auto c = std::make_shared<Connectivity>(compute_from_map(mesh, d0, d1));
    topology.set_connectivity(c, d0, d1);
  }
  else
    throw std::runtime_error("Entity dimension error when computing topology.");
}
//--------------------------------------------------------------------------
