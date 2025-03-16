// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_INPAR_POROELAST_HPP
#define FOUR_C_INPAR_POROELAST_HPP

#include "4C_config.hpp"

#include "4C_io_input_spec.hpp"
#include "4C_utils_exceptions.hpp"

#include <map>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
namespace Inpar
{
  namespace PoroElast
  {
    /// Type of coupling strategy for poroelasticity problems
    enum SolutionSchemeOverFields
    {
      //    OneWay,
      //   SequStagg,
      //   IterStagg,
      Partitioned,
      Monolithic,
      Monolithic_structuresplit,
      Monolithic_fluidsplit,
      Monolithic_nopenetrationsplit,
      Monolithic_meshtying
    };

    /// flag for control, in which equation the transient terms are included
    enum TransientEquationsOfPoroFluid
    {
      transient_none,
      transient_momentum_only,
      transient_continuity_only,
      transient_all
    };

    /// @name Solution technique and related

    /// type of norm to check for convergence
    enum ConvNorm
    {
      convnorm_undefined,
      convnorm_abs_global,       ///< absolute norm of global solution vectors
      convnorm_abs_singlefields  ///< absolute norm of single field solution vectors
      //    convnorm_rel_global,       ///< absolute norm of global solution vectors
      //    convnorm_rel_singlefields  ///< absolute norm of single field solution vectors
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

    /// type of norm to check for convergence
    enum BinaryOp
    {
      bop_undefined,
      bop_and,  ///<  and
      bop_or    ///<  or
    };

    /// type of initial field for poroelasticity problem
    enum InitialField
    {
      //  initfield_zero_field,
      initfield_field_by_function
      //   initfield_field_by_condition
    };

    //! map enum term to std::string
    static inline std::string vector_norm_string(const enum VectorNorm norm  //!< input enum term
    )
    {
      switch (norm)
      {
        case Inpar::PoroElast::norm_l1:
          return "L1";
          break;
        case Inpar::PoroElast::norm_l1_scaled:
          return "L1_scaled";
          break;
        case Inpar::PoroElast::norm_l2:
          return "L2";
          break;
        case Inpar::PoroElast::norm_rms:
          return "Rms";
          break;
        case Inpar::PoroElast::norm_inf:
          return "Inf";
          break;
        default:
          FOUR_C_THROW("Cannot make std::string to vector norm {}", norm);
          return "";
      }
    }

    //@}


    /// set the poroelast parameters
    void set_valid_parameters(std::map<std::string, Core::IO::InputSpec>& list);

  }  // namespace PoroElast

}  // namespace Inpar

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
