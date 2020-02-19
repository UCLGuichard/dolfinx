// Copyright (C) 2009 Anders Logg
//
// This file is part of DOLFINX (https://www.fenicsproject.org)
//
// SPDX-License-Identifier:    LGPL-3.0-or-later

#include "utils.h"
#include <cstdlib>
#include <sstream>

//-----------------------------------------------------------------------------
std::string dolfinx::common::indent(std::string block)
{
  std::string indentation("  ");
  std::stringstream s;

  s << indentation;
  for (std::size_t i = 0; i < block.size(); ++i)
  {
    s << block[i];
    if (block[i] == '\n' && i < block.size() - 1)
      s << indentation;
  }

  return s.str();
}
//-----------------------------------------------------------------------------
// std::vector<std::int64_t> dolfinx::common::transform_to_contiguous(
//     MPI_Comm comm, const std::vector<std::int64_t>& indices)
// {
//   // Get max index across all processes
//   auto it = std::max_element(indices.begin(), indices.end());
//   std::int64_t max = (it != indices.end()) ? *it : -1;
//   max = dolfinx::MPI::max(comm, max);
//   assert(max >= 0);

//   // Divide (0, max) equally across processes
//   const int size = MPI::size(comm);
//   std::vector<std::int64_t> ranges(size);
//   MPI::all_gather(comm, (std::int64_t)points.rows(), ranges);
//   // for (std::size_t i = 1; i < ranges.size(); ++i)
//   //   ranges[i] += ranges[i - 1];
//   // ranges.insert(ranges.begin(), 0);


//   return std::vector<std::int64_t>();
// }
//-----------------------------------------------------------------------------
