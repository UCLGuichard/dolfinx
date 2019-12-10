// Copyright (C) 2007-2019 Garth N. Wells
//
// This file is part of DOLFINX (https://www.fenicsproject.org)
//
// SPDX-License-Identifier:    LGPL-3.0-or-later

#include "SparsityPatternBuilder.h"
#include <dolfinx/common/IndexMap.h>
#include <dolfinx/common/MPI.h>
#include <dolfinx/fem/DofMap.h>
#include <dolfinx/fem/MultiPointConstraint.h>
#include <dolfinx/la/SparsityPattern.h>
#include <dolfinx/mesh/Mesh.h>
#include <dolfinx/mesh/MeshEntity.h>
#include <dolfinx/mesh/MeshIterator.h>

using namespace dolfinx;
using namespace dolfinx::fem;

//-----------------------------------------------------------------------------
void SparsityPatternBuilder::cells(
    la::SparsityPattern& pattern, const mesh::Mesh& mesh,
    const std::array<const fem::DofMap*, 2> dofmaps)
{
  assert(dofmaps[0]);
  assert(dofmaps[1]);
  const int D = mesh.topology().dim();
  for (auto& cell : mesh::MeshRange(mesh, D))
  {
    pattern.insert_local(dofmaps[0]->cell_dofs(cell.index()),
                         dofmaps[1]->cell_dofs(cell.index()));
  }
}
//-----------------------------------------------------------------------------
void SparsityPatternBuilder::interior_facets(
    la::SparsityPattern& pattern, const mesh::Mesh& mesh,
    const std::array<const fem::DofMap*, 2> dofmaps)
{
  assert(dofmaps[0]);
  assert(dofmaps[1]);

  const std::size_t D = mesh.topology().dim();
  mesh.create_entities(D - 1);
  mesh.create_connectivity(D - 1, D);

  // Array to store macro-dofs, if required (for interior facets)
  std::array<Eigen::Array<PetscInt, Eigen::Dynamic, 1>, 2> macro_dofs;
  std::shared_ptr<const mesh::Connectivity> connectivity
      = mesh.topology().connectivity(D - 1, D);
  if (!connectivity)
    throw std::runtime_error("Facet-cell connectivity has not been computed.");

  for (auto& facet : mesh::MeshRange(mesh, D - 1))
  {
    // Continue if facet is exterior facet
    if (connectivity->size_global(facet.index()) == 1)
      continue;

    // FIXME: sort out ghosting

    // Get cells incident with facet
    assert(connectivity->size(facet.index()) == 2);
    const mesh::MeshEntity cell0(mesh, D, facet.entities(D)[0]);
    const mesh::MeshEntity cell1(mesh, D, facet.entities(D)[1]);

    // Tabulate dofs for each dimension on macro element
    for (std::size_t i = 0; i < 2; i++)
    {
      auto cell_dofs0 = dofmaps[i]->cell_dofs(cell0.index());
      auto cell_dofs1 = dofmaps[i]->cell_dofs(cell1.index());
      macro_dofs[i].resize(cell_dofs0.size() + cell_dofs1.size());
      std::copy(cell_dofs0.data(), cell_dofs0.data() + cell_dofs0.size(),
                macro_dofs[i].data());
      std::copy(cell_dofs1.data(), cell_dofs1.data() + cell_dofs1.size(),
                macro_dofs[i].data() + cell_dofs0.size());
    }

    pattern.insert_local(macro_dofs[0], macro_dofs[1]);
  }
}
//-----------------------------------------------------------------------------
void SparsityPatternBuilder::exterior_facets(
    la::SparsityPattern& pattern, const mesh::Mesh& mesh,
    const std::array<const fem::DofMap*, 2> dofmaps)
{
  const std::size_t D = mesh.topology().dim();
  mesh.create_entities(D - 1);
  mesh.create_connectivity(D - 1, D);

  std::shared_ptr<const mesh::Connectivity> connectivity
      = mesh.topology().connectivity(D - 1, D);
  if (!connectivity)
    throw std::runtime_error("Facet-cell connectivity has not been computed.");
  for (auto& facet : mesh::MeshRange(mesh, D - 1))
  {
    // Skip interior facets
    if (connectivity->size_global(facet.index()) > 1)
      continue;

    // FIXME: sort out ghosting

    assert(connectivity->size(facet.index()) == 1);
    mesh::MeshEntity cell(mesh, D, facet.entities(D)[0]);
    pattern.insert_local(dofmaps[0]->cell_dofs(cell.index()),
                         dofmaps[1]->cell_dofs(cell.index()));
  }
}
//-----------------------------------------------------------------------------
void SparsityPatternBuilder::MultiPointConstraint(
    la::SparsityPattern& pattern, const mesh::Mesh& mesh,
    const std::array<const fem::DofMap*, 2> dofmaps,
    fem::MultiPointConstraint& mpc)
{
  assert(dofmaps[0]);
  assert(dofmaps[1]);
  const std::unordered_map<std::size_t, std::size_t> pairs
      = mpc.slave_to_master();

  const int D = mesh.topology().dim();
  for (auto& cell : mesh::MeshRange(mesh, D))
  {

    // Master slave sparsity pattern
    std::array<Eigen::Array<PetscInt, Eigen::Dynamic, 1>, 2> master_slave_dofs;
    // Dofs previously owned by slave dof
    std::array<Eigen::Array<PetscInt, Eigen::Dynamic, 1>, 2> new_master_dofs;

    for (std::size_t i = 0; i < 2; i++)
    {
      auto cell_dof_list = dofmaps[i]->cell_dofs(cell.index());
      new_master_dofs[i].resize(cell_dof_list.size());
      std::copy(cell_dof_list.data(),
                cell_dof_list.data() + cell_dof_list.size(),
                new_master_dofs[i].data());
      for (auto it = pairs.begin(); it != pairs.end(); ++it)
      {
        for (std::size_t j = 0; j < unsigned(cell_dof_list.size()); ++j)
        {

          if (it->first == unsigned(cell_dof_list[j]))
          {
            new_master_dofs[i](j) = it->second;
            master_slave_dofs[i].conservativeResize(master_slave_dofs[i].size()
                                                    + 2);
            master_slave_dofs[i].row(master_slave_dofs[i].rows() - 1)
                = it->first;
            master_slave_dofs[i].row(master_slave_dofs[i].rows() - 2)
                = it->second;
          }
        }
      }
    }
    pattern.insert_local(new_master_dofs[0], new_master_dofs[1]);
    pattern.insert_local(master_slave_dofs[0], master_slave_dofs[1]);
  }
}
//-----------------------------------------------------------------------------