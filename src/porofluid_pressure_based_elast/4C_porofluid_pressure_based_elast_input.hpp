// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_POROFLUID_PRESSURE_BASED_ELAST_INPUT_HPP
#define FOUR_C_POROFLUID_PRESSURE_BASED_ELAST_INPUT_HPP


#include "4C_config.hpp"

#include "4C_io_input_spec.hpp"

#include <map>

FOUR_C_NAMESPACE_OPEN

// forward declaration
namespace Core::Conditions
{
  class ConditionDefinition;
}  // namespace Core::Conditions

namespace PoroPressureBased
{
  /// type of coupling strategy for porofluid-elasticity problems
  enum class SolutionSchemePorofluidElast
  {
    undefined,
    twoway_partitioned,
    twoway_monolithic
  };

  //! relaxation methods for partitioned coupling
  enum class RelaxationMethods
  {
    none,
    constant,
    aitken
  };

  /// set the valid parameters for porofluid-elasticity problems
  void set_valid_parameters_porofluid_elast(std::map<std::string, Core::IO::InputSpec>& list);

}  // namespace PoroPressureBased


FOUR_C_NAMESPACE_CLOSE

#endif
