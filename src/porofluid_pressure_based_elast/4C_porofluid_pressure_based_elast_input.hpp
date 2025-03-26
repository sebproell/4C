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

namespace POROMULTIPHASE
{
  /// Type of coupling strategy for POROMULTIPHASE problems
  enum SolutionSchemeOverFields
  {
    solscheme_undefined,
    solscheme_twoway_partitioned,
    solscheme_twoway_monolithic
  };

  /// type of finite difference check
  enum FdCheck
  {
    fdcheck_none,
    fdcheck_global
  };

  /// type of norm to be calculated
  enum VectorNorm
  {
    norm_undefined,
    norm_l1,         //!< L1/linear norm
    norm_l1_scaled,  //!< L1/linear norm scaled by length of vector
    norm_l2,         //!< L2/Euclidean norm
    norm_rms,        //!< root mean square (RMS) norm
    norm_inf         //!< Maximum/infinity norm
  };

  //! relaxation methods for partitioned coupling
  enum RelaxationMethods
  {
    relaxation_none,
    relaxation_constant,
    relaxation_aitken
  };

  /// set the POROMULTIPHASE parameters
  void set_valid_parameters(std::map<std::string, Core::IO::InputSpec>& list);

}  // namespace POROMULTIPHASE


FOUR_C_NAMESPACE_CLOSE

#endif
