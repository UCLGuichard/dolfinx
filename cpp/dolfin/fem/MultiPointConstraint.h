// Copyright (C) 2019 Jorgen S. Dokken
//
// This file is part of DOLFIN (https://www.fenicsproject.org)
//
// SPDX-License-Identifier:    LGPL-3.0-or-later

#pragma once

#include <dolfin/function/FunctionSpace.h>

namespace dolfin
{
namespace fem
{

  /// This class provides the interface for setting multi-point constraints.
  ///
  ///   u_i = u_j,
  ///
  /// where u_i and u_j denotes the i-th and j-th global degree of freedom in the corresponding
  /// function space.
  /// A MultiPointBC is specified by the function space (trial space), and a vector of pairs, connecting
  /// master and slave nodes with a linear dependency.

  class MultiPointConstraint
  {

  public:
  /// Create multipointconstraint
  ///
  /// @param[in] V The functionspace on which the multipoint constraint
  /// condition is applied
  /// @param[in] slave_to_master Mapping specifying the relationship between a slave and a master node.
	MultiPointConstraint(std::shared_ptr<const function::FunctionSpace> V,  std::unordered_map<std::size_t, std::size_t> slave_to_master);

	/// Slave to master map
	std::unordered_map<std::size_t, std::size_t> slave_to_master() const;

  private:
	std::shared_ptr<const function::FunctionSpace> _function_space;
	const std::unordered_map<std::size_t, std::size_t> _slave_to_master;
   };

}
}