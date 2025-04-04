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
  /// Type of coupling strategy for POROMULTIPHASE problems
  enum SolutionSchemeOverFields
  {
    solscheme_undefined,
    solscheme_twoway_partitioned,
    solscheme_twoway_monolithic
  };

  //! relaxation methods for partitioned coupling
  enum RelaxationMethods
  {
    relaxation_none,
    relaxation_constant,
    relaxation_aitken
  };

  /// set the valid parameters
  void set_valid_parameters_porofluid_elast(std::map<std::string, Core::IO::InputSpec>& list);

}  // namespace PoroPressureBased


FOUR_C_NAMESPACE_CLOSE

#endif
