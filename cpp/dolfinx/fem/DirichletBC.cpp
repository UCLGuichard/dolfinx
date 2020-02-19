// Copyright (C) 2007-2020 Anders Logg and Garth N. Wells
//
// This file is part of DOLFINX (https://www.fenicsproject.org)
//
// SPDX-License-Identifier:    LGPL-3.0-or-later

#include "DirichletBC.h"
#include "DofMap.h"
#include "FiniteElement.h"
#include <algorithm>
#include <array>
#include <dolfinx/common/IndexMap.h>
#include <dolfinx/fem/CoordinateElement.h>
#include <dolfinx/function/Function.h>
#include <dolfinx/function/FunctionSpace.h>
#include <dolfinx/mesh/Mesh.h>
#include <dolfinx/mesh/MeshEntity.h>
#include <dolfinx/mesh/MeshIterator.h>
#include <dolfinx/mesh/cell_types.h>
#include <map>
#include <numeric>
#include <utility>

using namespace dolfinx;
using namespace dolfinx::fem;

namespace
{
//-----------------------------------------------------------------------------
/// Find DOFs on this processes that are constrained by a Dirichlet
/// condition detected by another process
///
/// @param[in] map The IndexMap with the dof layout
/// @param[in] dofs_local The IndexMap with the dof layout
/// @return List of local dofs with boundary conditions applied but
///   detected by other processes. It may contain duplicate entries.
std::vector<std::int32_t>
get_remote_bcs1(const common::IndexMap& map,
                const std::vector<std::int32_t>& dofs_local)
{
  // Get number of processes in neighbourhood
  MPI_Comm comm = map.mpi_comm_neighborhood();
  int num_neighbours(-1), outdegree(-2), weighted(-1);
  MPI_Dist_graph_neighbors_count(comm, &num_neighbours, &outdegree, &weighted);
  assert(num_neighbours == outdegree);

  // Return early if there are no neighbours
  if (num_neighbours == 0)
    return std::vector<std::int32_t>();

  std::vector<int> neighbours(num_neighbours), neighbours1(num_neighbours),
      weights(num_neighbours), weights1(num_neighbours);
  MPI_Dist_graph_neighbors(comm, num_neighbours, neighbours.data(),
                           weights.data(), num_neighbours, neighbours1.data(),
                           weights1.data());

  std::vector<std::int32_t> send_sizes(num_neighbours, 0);
  std::vector<std::int32_t> shared_dofs(dofs_local.size());
  std::map<int, std::set<int>> shared_indices = map.compute_shared_indices();
  for (std::size_t i = 0; i < dofs_local.size(); i++)
  {
    auto it = shared_indices.find(dofs_local[i]);
    if (it != shared_indices.end())
    {
      std::set<int>& sharing_proc = it->second;
      for (int p : sharing_proc)
      {
        auto iter = std::find(neighbours.begin(), neighbours.end(), p);
        assert(iter != neighbours.end());
        const int np = std::distance(neighbours.begin(), iter);
        int pos
            = std::accumulate(send_sizes.begin(), send_sizes.begin() + np, 0);
        shared_dofs.insert(shared_dofs.begin() + pos, dofs_local[i]);
        send_sizes[np]++;
      }
    }
  }

  // Figure out how many entries to receive from each neighbour
  std::vector<int> recv_sizes(num_neighbours);
  MPI_Neighbor_alltoall(send_sizes.data(), 1, MPI::mpi_type<int>(),
                        recv_sizes.data(), 1, MPI::mpi_type<int>(), comm);

  std::vector<int> send_disp(num_neighbours + 1, 0);
  std::vector<int> recv_disp(num_neighbours + 1, 0);
  std::partial_sum(send_sizes.begin(), send_sizes.end(), send_disp.begin() + 1);
  shared_dofs.resize(send_disp.back());
  const std::vector<std::int64_t> dofs_global
      = map.local_to_global(shared_dofs, false);
  std::partial_sum(recv_sizes.begin(), recv_sizes.end(), recv_disp.begin() + 1);

  //  May have repeated shared indices with different processes
  std::vector<std::int64_t> recv_dofs(recv_disp.back());
  MPI_Neighbor_alltoallv(
      dofs_global.data(), send_sizes.data(), send_disp.data(), MPI_INT64_T,
      recv_dofs.data(), recv_sizes.data(), recv_disp.data(), MPI_INT64_T, comm);

  // Build vector of local dof indicies that have been marked by another
  // process
  std::vector<std::int32_t> dofs = map.global_to_local(recv_dofs, false);
  dofs.erase(std::remove(dofs.begin(), dofs.end(), -1), dofs.end());

  return dofs;
}
//-----------------------------------------------------------------------------

/// Find DOFs on this processes that are constrained by a Dirichlet
/// condition detected by another process
///
/// @param[in] map0 The IndexMap with the dof layout
/// @param[in] map1 The IndexMap with the dof layout
/// @param[in] dofs_local The IndexMap with the dof layout
/// @return List of local dofs with boundary conditions applied but
///   detected by other processes. It may contain duplicate entries.
std::vector<std::array<std::int32_t, 2>>
get_remote_bcs2(const common::IndexMap& map0, const common::IndexMap& map1,
                const std::vector<std::array<std::int32_t, 2>>& dofs_local)
{

  std::array<MPI_Request, 2> requests;
  // Get number of processes in neighbourhood
  MPI_Comm comm = map0.mpi_comm_neighborhood();
  int num_neighbours(-1), outdegree(-2), weighted(-1);
  MPI_Dist_graph_neighbors_count(comm, &num_neighbours, &outdegree, &weighted);
  assert(num_neighbours == outdegree);

  // Return early if there are no neighbours
  if (num_neighbours == 0)
    return std::vector<std::array<std::int32_t, 2>>();

  std::vector<int> neighbours(num_neighbours), neighbours1(num_neighbours),
      weights(num_neighbours), weights1(num_neighbours);
  MPI_Dist_graph_neighbors(comm, num_neighbours, neighbours.data(),
                           weights.data(), num_neighbours, neighbours1.data(),
                           weights1.data());

  std::array<std::reference_wrapper<const common::IndexMap>, 2> maps{
      {map0, map1}};
  std::array<std::vector<std::int32_t>, 2> dofs_local_i;
  std::array<std::vector<std::int32_t>, 2> dofs_i;
  std::array<std::vector<std::int64_t>, 2> recv_dofs_i;

  for (std::size_t i = 0; i < 2; i++)
  {
    std::vector<std::int32_t>& dofs_i = dofs_local_i[i];
    dofs_i.reserve(dofs_local.size());
    for (std::size_t j = 0; j < dofs_local.size(); ++j)
      dofs_i[j] = dofs_local[j][i];
  }

  for (std::size_t i = 0; i < 2; i++)
  {
    const common::IndexMap& map = maps[i];
    std::vector<std::int32_t> send_sizes(num_neighbours, 0);
    std::vector<std::int32_t> shared_dofs(dofs_local.size());
    std::map<int, std::set<int>> shared_indices = map.compute_shared_indices();

    for (std::size_t j = 0; j < dofs_local_i[i].size(); j++)
    {
      auto it = shared_indices.find(dofs_local_i[i][j]);
      if (it != shared_indices.end())
      {
        std::set<int>& sharing_proc = it->second;
        for (int p : sharing_proc)
        {
          auto iter = std::find(neighbours.begin(), neighbours.end(), p);
          assert(iter != neighbours.end());
          const int np = std::distance(neighbours.begin(), iter);
          int pos
              = std::accumulate(send_sizes.begin(), send_sizes.begin() + np, 0);
          shared_dofs.insert(shared_dofs.begin() + pos, dofs_local_i[i][j]);
          send_sizes[np]++;
        }
      }
    }
    // Figure out how many entries to receive from each neighbour
    std::vector<int> recv_sizes(num_neighbours);
    MPI_Neighbor_alltoall(send_sizes.data(), 1, MPI::mpi_type<int>(),
                          recv_sizes.data(), 1, MPI::mpi_type<int>(), comm);

    std::vector<int> send_disp(num_neighbours + 1, 0);
    std::vector<int> recv_disp(num_neighbours + 1, 0);
    std::partial_sum(send_sizes.begin(), send_sizes.end(),
                     send_disp.begin() + 1);
    shared_dofs.resize(send_disp.back());
    const std::vector<std::int64_t> dofs_global
        = map.local_to_global(shared_dofs, false);
    std::partial_sum(recv_sizes.begin(), recv_sizes.end(),
                     recv_disp.begin() + 1);

    //  May have repeated shared indices with different processes
    std::vector<std::int64_t>& recv_dofs = recv_dofs_i[i];
    recv_dofs.reserve(recv_disp.back());
    MPI_Ineighbor_alltoallv(dofs_global.data(), send_sizes.data(),
                            send_disp.data(), MPI_INT64_T, recv_dofs.data(),
                            recv_sizes.data(), recv_disp.data(), MPI_INT64_T,
                            comm, &requests[i]);
  }

  // Build vector of local dof indicies that have been marked by another
  // process
  for (std::size_t i = 0; i < 2; i++)
  {
    int index;
    MPI_Waitany(2, requests.data(), &index, MPI_STATUS_IGNORE);
    const common::IndexMap& map = maps[index];
    dofs_i[index] = map.global_to_local(recv_dofs_i[index], false);
    dofs_i[index].erase(
        std::remove(dofs_i[index].begin(), dofs_i[index].end(), -1),
        dofs_i[index].end());
  }

  std::vector<std::array<std::int32_t, 2>> dofs;
  dofs.reserve(dofs_i[0].size());
  for (std::size_t i = 0; i < dofs_i[0].size(); ++i)
    dofs.push_back({dofs_i[0][i], dofs_i[1][i]});

  return dofs;
}
//-----------------------------------------------------------------------------
Eigen::Array<std::int32_t, Eigen::Dynamic, 2> _locate_dofs_topological(
    const std::vector<std::reference_wrapper<function::FunctionSpace>>& V,
    const int dim, const Eigen::Ref<const Eigen::ArrayXi>& entities,
    bool remote)
{
  const function::FunctionSpace& V0 = V.at(0).get();
  const function::FunctionSpace& V1 = V.at(1).get();

  // Get mesh
  assert(V0.mesh());
  assert(V1.mesh());
  if (V0.mesh() != V1.mesh())
    throw std::runtime_error("Meshes are not the same.");
  const mesh::Mesh& mesh = *V0.mesh();
  const std::size_t tdim = mesh.topology().dim();

  assert(V0.element());
  assert(V1.element());
  if (!V0.has_element(*V1.element()))
  {
    throw std::runtime_error("Function spaces must have the same elements or "
                             "one be a subelement of another.");
  }

  // Get dofmaps
  assert(V0.dofmap());
  assert(V1.dofmap());
  const DofMap& dofmap0 = *V0.dofmap();
  const DofMap& dofmap1 = *V1.dofmap();

  // Initialise entity-cell connectivity
  mesh.create_entities(tdim);
  mesh.create_connectivity(dim, tdim);

  // Allocate space
  assert(dofmap0.element_dof_layout);
  const int num_entity_dofs
      = dofmap0.element_dof_layout->num_entity_closure_dofs(dim);

  // Build vector local dofs for each cell facet
  std::vector<Eigen::Array<int, Eigen::Dynamic, 1>> entity_dofs;
  for (int i = 0; i < mesh::cell_num_entities(mesh.topology().cell_type(), dim);
       ++i)
  {
    entity_dofs.push_back(
        dofmap0.element_dof_layout->entity_closure_dofs(dim, i));
  }

  // Iterate over marked facets
  std::vector<std::array<std::int32_t, 2>> bc_dofs;
  for (Eigen::Index e = 0; e < entities.rows(); ++e)
  {
    // Create facet and attached cell
    const mesh::MeshEntity entity(mesh, dim, entities[e]);
    const std::size_t cell_index = entity.entities(tdim)[0];
    const mesh::MeshEntity cell(mesh, tdim, cell_index);

    // Get cell dofmap
    auto cell_dofs0 = dofmap0.cell_dofs(cell.index());
    auto cell_dofs1 = dofmap1.cell_dofs(cell.index());

    // Loop over facet dofs
    const int entity_local_index = cell.index(entity);
    for (int i = 0; i < num_entity_dofs; ++i)
    {
      const int index = entity_dofs[entity_local_index][i];
      bc_dofs.push_back(
          {(std::int32_t)cell_dofs0[index], (std::int32_t)cell_dofs1[index]});
    }
  }

  // TODO: is removing duplicates at this point worth the effort?
  // Remove duplicates
  std::sort(bc_dofs.begin(), bc_dofs.end());
  bc_dofs.erase(std::unique(bc_dofs.begin(), bc_dofs.end()), bc_dofs.end());

  if (remote)
  {
    // Get bc dof indices (local) in (V, Vg) spaces on this process that
    // were found by other processes, e.g. a vertex dof on this process
    // that has no connected facets on the boundary.
    const std::vector<std::array<std::int32_t, 2>> dofs_remote
        = get_remote_bcs2(*V0.dofmap()->index_map, *V1.dofmap()->index_map,
                          bc_dofs);

    // Add received bc indices to dofs_local
    bc_dofs.insert(bc_dofs.end(), dofs_remote.begin(), dofs_remote.end());

    // Remove duplicates
    std::sort(bc_dofs.begin(), bc_dofs.end());
    bc_dofs.erase(std::unique(bc_dofs.begin(), bc_dofs.end()), bc_dofs.end());
  }

  Eigen::Array<std::int32_t, Eigen::Dynamic, 2> dofs(bc_dofs.size(), 2);
  for (std::size_t i = 0; i < bc_dofs.size(); ++i)
  {
    dofs(i, 0) = bc_dofs[i][0];
    dofs(i, 1) = bc_dofs[i][1];
  }

  return dofs;
}
//-----------------------------------------------------------------------------

/// TODO: Add doc
Eigen::Array<std::int32_t, Eigen::Dynamic, 1>
_locate_dofs_topological(const function::FunctionSpace& V, const int entity_dim,
                         const Eigen::Ref<const Eigen::ArrayXi>& entities,
                         bool remote)
{
  assert(V.dofmap());
  const DofMap& dofmap = *V.dofmap();
  assert(V.mesh());
  mesh::Mesh mesh = *V.mesh();

  const int tdim = mesh.topology().dim();

  // Initialise entity-cell connectivity
  mesh.create_entities(tdim);
  mesh.create_connectivity(entity_dim, tdim);

  // Prepare an element - local dof layout for dofs on entities of the
  // entity_dim
  const int num_cell_entities
      = mesh::cell_num_entities(mesh.topology().cell_type(), entity_dim);
  std::vector<Eigen::Array<int, Eigen::Dynamic, 1>> entity_dofs;
  for (int i = 0; i < num_cell_entities; ++i)
  {
    entity_dofs.push_back(
        dofmap.element_dof_layout->entity_closure_dofs(entity_dim, i));
  }

  const int num_entity_closure_dofs
      = dofmap.element_dof_layout->num_entity_closure_dofs(entity_dim);
  std::vector<std::int32_t> dofs;
  for (Eigen::Index i = 0; i < entities.rows(); ++i)
  {
    // Create entity and attached cell
    const mesh::MeshEntity entity(mesh, entity_dim, entities[i]);
    const int cell_index = entity.entities(tdim)[0];
    const mesh::MeshEntity cell(mesh, tdim, cell_index);

    // Get cell dofmap
    auto cell_dofs = dofmap.cell_dofs(cell_index);

    // Loop over entity dofs
    const int entity_local_index = cell.index(entity);
    for (int j = 0; j < num_entity_closure_dofs; j++)
    {
      const int index = entity_dofs[entity_local_index][j];
      dofs.push_back(cell_dofs[index]);
    }
  }

  // TODO: is removing duplicates at this point worth the effort?
  // Remove duplicates
  std::sort(dofs.begin(), dofs.end());
  dofs.erase(std::unique(dofs.begin(), dofs.end()), dofs.end());

  if (remote)
  {
    const std::vector<std::int32_t> dofs_remote
        = get_remote_bcs1(*V.dofmap()->index_map, dofs);

    // Add received bc indices to dofs_local
    dofs.insert(dofs.end(), dofs_remote.begin(), dofs_remote.end());

    // Remove duplicates
    std::sort(dofs.begin(), dofs.end());
    dofs.erase(std::unique(dofs.begin(), dofs.end()), dofs.end());
  }

  // Copy to Eigen array
  Eigen::Array<std::int32_t, Eigen::Dynamic, 1> _dofs
      = Eigen::Map<Eigen::Array<std::int32_t, Eigen::Dynamic, 1>>(dofs.data(),
                                                                  dofs.size());

  return _dofs;
}
//-----------------------------------------------------------------------------
Eigen::Array<std::int32_t, Eigen::Dynamic, 2> _locate_dofs_geometrical(
    const std::vector<std::reference_wrapper<function::FunctionSpace>>& V,
    marking_function marker)
{
  // FIXME: Calling V.tabulate_dof_coordinates() is very expensive,
  // especially when we usually want the boundary dofs only. Add
  // interface that computes dofs coordinates only for specified cell.

  // Get function spaces
  const function::FunctionSpace& V0 = V.at(0).get();
  const function::FunctionSpace& V1 = V.at(1).get();

  // Get mesh
  assert(V0.mesh());
  assert(V1.mesh());
  if (V0.mesh() != V1.mesh())
    throw std::runtime_error("Meshes are not the same.");
  const mesh::Mesh& mesh = *V1.mesh();
  const std::size_t tdim = mesh.topology().dim();

  assert(V0.element());
  assert(V1.element());
  if (!V0.has_element(*V1.element()))
  {
    throw std::runtime_error("Function spaces must have the same elements or "
                             "one be a subelement of another.");
  }

  // Compute dof coordinates
  const Eigen::Array<double, 3, Eigen::Dynamic, Eigen::RowMajor> dof_coordinates
      = V1.tabulate_dof_coordinates().transpose();

  // Evaluate marker for each dof coordinate
  const Eigen::Array<bool, Eigen::Dynamic, 1> marked_dofs
      = marker(dof_coordinates);

  // Get dofmaps
  assert(V0.dofmap());
  assert(V1.dofmap());
  const DofMap& dofmap0 = *V0.dofmap();
  const DofMap& dofmap1 = *V1.dofmap();

  // Iterate over cells
  const mesh::Topology& topology = mesh.topology();
  std::vector<std::array<std::int32_t, 2>> bc_dofs;
  for (int c = 0; c < topology.connectivity(tdim, 0)->num_nodes(); ++c)
  {
    // Get cell dofmap
    auto cell_dofs0 = dofmap0.cell_dofs(c);
    auto cell_dofs1 = dofmap1.cell_dofs(c);

    // Loop over cell dofs and add to bc_dofs if marked.
    for (Eigen::Index i = 0; i < cell_dofs1.rows(); ++i)
    {
      if (marked_dofs[cell_dofs1[i]])
      {
        bc_dofs.push_back(
            {(std::int32_t)cell_dofs0[i], (std::int32_t)cell_dofs1[i]});
      }
    }
  }

  // Remove duplicates
  std::sort(bc_dofs.begin(), bc_dofs.end());
  bc_dofs.erase(std::unique(bc_dofs.begin(), bc_dofs.end()), bc_dofs.end());

  // Copy to Eigen array
  Eigen::Array<std::int32_t, Eigen::Dynamic, 2> dofs(bc_dofs.size(), 2);
  for (std::size_t i = 0; i < bc_dofs.size(); ++i)
  {
    dofs(i, 0) = bc_dofs[i][0];
    dofs(i, 1) = bc_dofs[i][1];
  }

  return dofs;
}
//-----------------------------------------------------------------------------
Eigen::Array<std::int32_t, Eigen::Dynamic, 1>
_locate_dofs_geometrical(const function::FunctionSpace& V,
                         marking_function marker)
{
  // FIXME: Calling V.tabulate_dof_coordinates() is very expensive,
  // especially when we usually want the boundary dofs only. Add
  // interface that computes dofs coordinates only for specified cell.

  // Compute dof coordinates
  const Eigen::Array<double, 3, Eigen::Dynamic, Eigen::RowMajor> dof_coordinates
      = V.tabulate_dof_coordinates().transpose();

  // Compute marker for each dof coordinate
  const Eigen::Array<bool, Eigen::Dynamic, 1> marked_dofs
      = marker(dof_coordinates);

  std::vector<std::int32_t> dofs;
  dofs.reserve(marked_dofs.count());
  for (Eigen::Index i = 0; i < marked_dofs.rows(); ++i)
  {
    if (marked_dofs[i])
      dofs.push_back(i);
  }

  return Eigen::Map<Eigen::Array<std::int32_t, Eigen::Dynamic, 1>>(dofs.data(),
                                                                   dofs.size());
}
} // namespace

//-----------------------------------------------------------------------------
Eigen::Array<std::int32_t, Eigen::Dynamic, Eigen::Dynamic>
fem::locate_dofs_topological(
    const std::vector<std::reference_wrapper<function::FunctionSpace>>& V,
    const int dim, const Eigen::Ref<const Eigen::ArrayXi>& entities,
    bool remote)
{
  if (V.size() == 2)
    return _locate_dofs_topological(V, dim, entities, remote);
  else if (V.size() == 1)
    return _locate_dofs_topological(V[0].get(), dim, entities, remote);
  else
    throw std::runtime_error("Expected only 1 or 2 function spaces.");
}
//-----------------------------------------------------------------------------
Eigen::Array<std::int32_t, Eigen::Dynamic, Eigen::Dynamic>
fem::locate_dofs_geometrical(
    const std::vector<std::reference_wrapper<function::FunctionSpace>>& V,
    marking_function marker)
{
  if (V.size() == 2)
    return _locate_dofs_geometrical(V, marker);
  else if (V.size() == 1)
    return _locate_dofs_geometrical(V[0].get(), marker);
  else
    throw std::runtime_error("Expected only 1 or 2 function spaces.");
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
DirichletBC::DirichletBC(
    std::shared_ptr<const function::Function> g,
    const Eigen::Ref<const Eigen::Array<std::int32_t, Eigen::Dynamic, 1>>&
        V_dofs)
    : _function_space(g->function_space()), _g(g), _dofs(V_dofs.rows(), 2)
{
  // Stack indices as columns, fits column-major _dofs layout
  _dofs.col(0) = V_dofs;
  _dofs.col(1) = V_dofs;

  const int owned_size = _function_space->dofmap()->index_map->block_size
                         * _function_space->dofmap()->index_map->size_local();
  auto it = std::lower_bound(_dofs.col(0).data(),
                             _dofs.col(0).data() + _dofs.rows(), owned_size);
  _owned_indices = std::distance(_dofs.col(0).data(), it);
}
//-----------------------------------------------------------------------------
DirichletBC::DirichletBC(
    std::shared_ptr<const function::Function> g,
    const Eigen::Ref<const Eigen::Array<std::int32_t, Eigen::Dynamic, 2>>&
        V_g_dofs,
    std::shared_ptr<const function::FunctionSpace> V)
    : _function_space(V), _g(g), _dofs(V_g_dofs)
{
  const int owned_size = _function_space->dofmap()->index_map->block_size
                         * _function_space->dofmap()->index_map->size_local();
  auto it = std::lower_bound(_dofs.col(0).data(),
                             _dofs.col(0).data() + _dofs.rows(), owned_size);
  _owned_indices = std::distance(_dofs.col(0).data(), it);
}
//-----------------------------------------------------------------------------
std::shared_ptr<const function::FunctionSpace>
DirichletBC::function_space() const
{
  return _function_space;
}
//-----------------------------------------------------------------------------
std::shared_ptr<const function::Function> DirichletBC::value() const
{
  return _g;
}
//-----------------------------------------------------------------------------
const Eigen::Array<std::int32_t, Eigen::Dynamic, 2>& DirichletBC::dofs() const
{
  return _dofs;
}
//-----------------------------------------------------------------------------
const Eigen::Ref<const Eigen::Array<std::int32_t, Eigen::Dynamic, 2>>
DirichletBC::dofs_owned() const
{
  return _dofs.block<Eigen::Dynamic, 2>(0, 0, _owned_indices, 2);
}
// -----------------------------------------------------------------------------
void DirichletBC::set(
    Eigen::Ref<Eigen::Matrix<PetscScalar, Eigen::Dynamic, 1>> x,
    double scale) const
{
  // FIXME: This one excludes ghosts. Need to straighten out.
  assert(_g);
  la::VecReadWrapper g(_g->vector().vec(), false);
  for (Eigen::Index i = 0; i < _dofs.rows(); ++i)
  {
    if (_dofs(i, 0) < x.rows())
      x[_dofs(i, 0)] = scale * g.x[_dofs(i, 1)];
  }
}
//-----------------------------------------------------------------------------
void DirichletBC::set(
    Eigen::Ref<Eigen::Matrix<PetscScalar, Eigen::Dynamic, 1>> x,
    const Eigen::Ref<const Eigen::Matrix<PetscScalar, Eigen::Dynamic, 1>>& x0,
    double scale) const
{
  // FIXME: This one excludes ghosts. Need to straighten out.
  assert(_g);
  assert(x.rows() <= x0.rows());
  la::VecReadWrapper g(_g->vector().vec(), false);
  for (Eigen::Index i = 0; i < _dofs.rows(); ++i)
  {
    if (_dofs(i, 0) < x.rows())
      x[_dofs(i, 0)] = scale * (g.x[_dofs(i, 1)] - x0[_dofs(i, 0)]);
  }
}
//-----------------------------------------------------------------------------
void DirichletBC::dof_values(
    Eigen::Ref<Eigen::Matrix<PetscScalar, Eigen::Dynamic, 1>> values) const
{
  assert(_g);
  la::VecReadWrapper g(_g->vector().vec());
  for (Eigen::Index i = 0; i < _dofs.rows(); ++i)
    values[_dofs(i, 0)] = g.x[_dofs(i, 1)];
}
//-----------------------------------------------------------------------------
void DirichletBC::mark_dofs(std::vector<bool>& markers) const
{
  for (Eigen::Index i = 0; i < _dofs.rows(); ++i)
  {
    assert(_dofs(i, 0) < (std::int32_t)markers.size());
    markers[_dofs(i, 0)] = true;
  }
}
//-----------------------------------------------------------------------------
