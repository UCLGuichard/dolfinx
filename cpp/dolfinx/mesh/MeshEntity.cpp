// Copyright (C) 2006-2011 Anders Logg
//
// This file is part of DOLFINX (https://www.fenicsproject.org)
//
// SPDX-License-Identifier:    LGPL-3.0-or-later

#include "MeshEntity.h"
#include "Geometry.h"
#include "Mesh.h"
#include "MeshEntity.h"
#include "MeshIterator.h"
#include "Topology.h"
#include <dolfinx/common/log.h>

using namespace dolfinx;
using namespace dolfinx::mesh;

//-----------------------------------------------------------------------------
int MeshEntity::index(const MeshEntity& entity) const
{
  // Must be in the same mesh to be incident
  if (_mesh != entity._mesh)
    throw std::runtime_error("Mesh entity is defined on a different mesh");

  // Get list of entities for given topological dimension
  const std::int32_t* entities = _mesh->topology()
                                     .connectivity(_dim, entity._dim)
                                     ->connections(_local_index);
  // Check if any entity matches
  for (int i = 0; i < num_entities(entity._dim); ++i)
    if (entities[i] == entity._local_index)
      return i;

  // Entity was not found
  throw std::runtime_error("Mesh entity was not found");

  return -1;
}
//-----------------------------------------------------------------------------
int MeshEntity::facet_permutation(const MeshEntity& entity) const
{
  // FIXME: cache this somewhere to avoid computing it every time its needed
  // AIM:
  // return facet_perm[entity._dim][index(entity)]

  // Must be in the same mesh to be incident
  if (_mesh != entity._mesh)
    throw std::runtime_error("Mesh entity is defined on a different mesh");

  // If the entity is a point, no permutation is required
  if (entity._dim == 0)
    return 0;

  // If the entity is an interval, it should be oriented pointing from the
  // lowest numbered vertex to the highest numbered vertex
  if (entity._dim == 1)
  {
    const std::int32_t* vertices = entity.entities(0);
    const int e_vertices[2] = {get_vertex_local_index(vertices[0]),
                               get_vertex_local_index(vertices[1])};
    // Return the number of reflections
    return e_vertices[1] < e_vertices[0];
  }

  // Triangles and quadrilaterals
  if (entity._dim == 2)
  {
    if (entity.num_entities(0) == 3)
    { // triangle
      // Find the index of the lowest numbered vertex
      int num_min = -1;
      const std::int32_t* vertices = entity.entities(0);
      const int e_vertices[3] = {get_vertex_local_index(vertices[0]),
                                 get_vertex_local_index(vertices[1]),
                                 get_vertex_local_index(vertices[2])};
      for (int v = 0; v < 3; ++v)
        if (num_min == -1 || e_vertices[v] < e_vertices[num_min])
          num_min = v;
      // pre is the number of the next vertex clockwise from the lowest numbered
      // vertex
      const int pre = num_min == 0 ? e_vertices[entity.num_entities(0) - 1]
                                   : e_vertices[num_min - 1];
      // post is the number of the next vertex anticlockwise from the lowest
      // numbered vertex
      const int post = num_min == entity.num_entities(0) - 1
                           ? e_vertices[0]
                           : e_vertices[num_min + 1];
      // Orient that triangle so the the lowest numbered vertex is the origin,
      // and the next vertex anticlockwise from the lowest has a lower number
      // than the next vertex clockwise. Return N such that N % 2 is the number
      // of reflections and N / 2 is the number of rotations
      return 2 * num_min + (post > pre);
    }
    if (entity.num_entities(0) == 4)
    { // quadrilateral
      // Find the index of the lowest numbered vertex
      int num_min = -1;
      const std::int32_t* vertices = entity.entities(0);
      const int e_vertices[4] = {get_vertex_local_index(vertices[0]),
                                 get_vertex_local_index(vertices[1]),
                                 get_vertex_local_index(vertices[2]),
                                 get_vertex_local_index(vertices[3])};
      for (int v = 0; v < 4; ++v)
        if (num_min == -1 || e_vertices[v] < e_vertices[num_min])
          num_min = v;
      // rots is the number of rotations to get the lowest numbered vertex to
      // the origin
      int rots = num_min;
      // pre is the (local) number of the next vertex clockwise from the lowest
      // numbered vertex
      int pre = 2;
      // post is the (local) number of the next vertex anticlockwise from the
      // lowest numbered vertex
      int post = 1;
      // The tensor product ordering of quads must be taken into account
      if (num_min == 1)
      {
        pre = 0;
        post = 3;
      }
      else if (num_min == 2)
      {
        pre = 3;
        post = 0;
        rots = 3;
      }
      else if (num_min == 3)
      {
        pre = 1;
        post = 2;
        rots = 2;
      }
      // Orient that quad so the the lowest numbered vertex is the origin, and
      // the next vertex anticlockwise from the lowest has a lower number than
      // the next vertex clockwise. Return N such that N % 2 is the number of
      // reflections and N / 2 is the number of rotations
      return 2 * rots + (e_vertices[post] > e_vertices[pre]);
    }
  }
  LOG(WARNING) << "No facet permutation was found for a facet. Integrals "
                  "containing jumps may be incorrect.";
  return 0;
}
//-----------------------------------------------------------------------------
std::string MeshEntity::str(bool verbose) const
{
  if (verbose)
    LOG(WARNING) << "Verbose output for MeshEntityIterator not implemented.";

  std::stringstream s;
  s << "<Mesh entity " << index() << " of topological dimension " << dim()
    << ">";
  return s.str();
}
//-----------------------------------------------------------------------------
