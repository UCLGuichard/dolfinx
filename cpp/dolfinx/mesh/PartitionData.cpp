// Copyright (C) 2018 Garth N. Wells
//
// This file is part of DOLFINX (https://www.fenicsproject.org)
//
// SPDX-License-Identifier:    LGPL-3.0-or-later

#include "PartitionData.h"
#include <algorithm>
#include <set>
#include <sstream>

using namespace dolfinx;
using namespace dolfinx::mesh;

//-----------------------------------------------------------------------------
PartitionData::PartitionData(
    const std::vector<int>& cell_partition,
    const std::map<std::int64_t, std::vector<int>>& ghost_procs)
    : _offset(1)

{
  for (std::size_t i = 0; i < cell_partition.size(); ++i)
  {
    auto it = ghost_procs.find(i);
    if (it == ghost_procs.end())
      _dest_processes.push_back(cell_partition[i]);
    else
    {
      _dest_processes.insert(_dest_processes.end(), it->second.begin(),
                             it->second.end());
    }
    _offset.push_back(_dest_processes.size());
  }
}
//-----------------------------------------------------------------------------
PartitionData::PartitionData(
    const std::pair<std::vector<int>, std::map<std::int64_t, std::vector<int>>>&
        data)
    : PartitionData(data.first, data.second)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
std::int32_t PartitionData::num_procs(std::int32_t i) const
{
  return _offset[i + 1] - _offset[i];
}
//-----------------------------------------------------------------------------
const std::int32_t* PartitionData::procs(std::int32_t i) const
{
  return _dest_processes.data() + _offset[i];
}
//-----------------------------------------------------------------------------
std::int32_t PartitionData::size() const { return _offset.size() - 1; }
//-----------------------------------------------------------------------------
std::int32_t PartitionData::num_ghosts() const
{
  return _dest_processes.size() - _offset.size() + 1;
}
//-----------------------------------------------------------------------------
MPI_Comm PartitionData::neighbour_comm(MPI_Comm mpi_comm) const
{
  const int num_processes = MPI::size(mpi_comm);
  const int mpi_rank = MPI::rank(mpi_comm);

  std::int32_t partition_size = _offset.size() - 1;

  std::vector<std::set<std::int32_t>> edges_per_proc(num_processes);
  for (std::int32_t i = 0; i < partition_size; i++)
  {
    std::int32_t num_procs = _offset[i + 1] - _offset[i];
    std::vector<std::int32_t> procs(_dest_processes.data() + _offset[i],
                                    _dest_processes.data() + _offset[i]
                                        + num_procs);
    for (std::int32_t j = 0; j < num_procs; j++)
      edges_per_proc[procs[j]].insert(procs.begin(), procs.end());
  }

  std::vector<std::vector<std::int32_t>> send_buffer(num_processes);
  std::vector<std::vector<std::int32_t>> received_buffer(num_processes);

  for (std::int32_t i = 0; i < num_processes; i++)
    send_buffer[i].assign(edges_per_proc[i].begin(), edges_per_proc[i].end());

  dolfinx::MPI::all_to_all(mpi_comm, send_buffer, received_buffer);

  std::set<std::int32_t> neighbors_set(send_buffer[mpi_rank].begin(),
                                       send_buffer[mpi_rank].end());
  for (std::int32_t i = 0; i < num_processes; i++)
    neighbors_set.insert(received_buffer[i].begin(), received_buffer[i].end());

  std::vector<std::int32_t> neighbors(neighbors_set.begin(),
                                      neighbors_set.end());

  MPI_Comm neighbour_comm;

  // Create neighbourhood communicator. No communication is needed to
  // build the graph with complete adjacency information
  MPI_Dist_graph_create_adjacent(mpi_comm, neighbors.size(), neighbors.data(),
                                 MPI_UNWEIGHTED, neighbors.size(),
                                 neighbors.data(), MPI_UNWEIGHTED,
                                 MPI_INFO_NULL, false, &neighbour_comm);

#if DEBUG
  {
    int indegree(-1), outdegree(-2), weighted(-1);
    MPI_Dist_graph_neighbors_count(neighbour_comm, &indegree, &outdegree,
                                   &weighted);
    assert(indegree == outdegree);
    std::vector<int> neighbours(indegree), neighbours1(indegree),
        weights(indegree), weights1(indegree);

    MPI_Dist_graph_neighbors(neighbour_comm, indegree, neighbours.data(),
                             weights.data(), outdegree, neighbours1.data(),
                             weights1.data());
    assert(neighbours == neighbours1);
  }
#endif

  return neighbour_comm;
}
