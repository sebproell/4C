// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_INPAR_TSI_HPP
#define FOUR_C_INPAR_TSI_HPP

#include "4C_config.hpp"

#include "4C_io_input_spec.hpp"
#include "4C_utils_exceptions.hpp"

#include <map>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |                                                           dano 11/09 |
 *----------------------------------------------------------------------*/
namespace Inpar
{
  namespace TSI
  {
    //! Type of coupling strategy for TSI problems
    enum SolutionSchemeOverFields
    {
      OneWay,
      SequStagg,
      IterStagg,
      IterStaggAitken,
      IterStaggAitkenIrons,
      IterStaggFixedRel,
      Monolithic
    };

    //! @name Solution technique and related stuff

    //! type of norm to check for convergence
    enum ConvNorm
    {
      convnorm_abs,  //!< absolute norm
      convnorm_rel,  //!< relative norm of TSI problem with initial TSI rhs
      convnorm_mix   //!< mixed absolute-relative norm
    };

    //! type of norm to check for convergence
    enum BinaryOp
    {
      bop_and,               //!< and
      bop_or,                //!< or
      bop_coupl_or_single,   //!< either TSI problem or single field problems converged
      bop_coupl_and_single,  //!< either TSI problem or single field problems converged
      bop_and_single,        //!< and in single field problems
      bop_or_single          //!< or in single field problems
    };

    //! type of solution techniques
    enum NlnSolTech
    {
      soltech_newtonfull,  //!< full Newton-Raphson iteration
      soltech_ptc,         //!< pseudo transient continuation nonlinear iteration
    };

    //! type of line-search strategy
    enum LineSearch
    {
      LS_none = 0,   //!< no line search
      LS_structure,  //!< line-search based on structural residual
      LS_thermo,     //!< line-search based on thermal residual
      LS_or,         //!< line-search based on structural or thermal residual
      LS_and,        //!< line-search based on structural and thermal residual
    };

    //! Map solution technique enum to std::string
    static inline std::string nln_sol_tech_string(const enum NlnSolTech name  //!< enum to convert
    )
    {
      switch (name)
      {
        case soltech_newtonfull:
          return "fullnewton";
          break;
        case soltech_ptc:
          return "ptc";
          break;
        default:
          FOUR_C_THROW("Cannot make std::string for solution technique {}", name);
          return "";
      }
    }

    //@}

    //! @name General
    //@{

    //! type of vector norm used for error/residual vectors
    enum VectorNorm
    {
      norm_vague = 0,  //!< undetermined norm
      norm_l1,         //!< L1/linear norm
      norm_l1_scaled,  //!< L1/linear norm scaled by length of vector
      norm_l2,         //!< L2/Euclidean norm
      norm_rms,        //!< root mean square (RMS) norm
      norm_inf         //!< Maximum/infinity norm
    };

    //! map enum term to std::string
    static inline std::string vector_norm_string(const enum VectorNorm norm  //!< input enum term
    )
    {
      switch (norm)
      {
        case Inpar::TSI::norm_vague:
          return "Vague";
          break;
        case Inpar::TSI::norm_l1:
          return "L1";
          break;
        case Inpar::TSI::norm_l1_scaled:
          return "L1_scaled";
          break;
        case Inpar::TSI::norm_l2:
          return "L2";
          break;
        case Inpar::TSI::norm_rms:
          return "Rms";
          break;
        case Inpar::TSI::norm_inf:
          return "Inf";
          break;
        default:
          FOUR_C_THROW("Cannot make std::string to vector norm {}", norm);
          return "";
      }
    }

    //! Method used to calculate plastic dissipation
    enum DissipationMode
    {
      pl_multiplier,  //!< Dissipation = yield stress times plastic multiplier
      pl_flow,        //!< Dissipation = Mandel stress : sym(L^p)
      Taylor_Quinney  //!< Dissipation based on Taylor Quinney factor
    };

    //@}

    /// set the tsi parameters
    void set_valid_parameters(std::map<std::string, Core::IO::InputSpec>& list);


  }  // namespace TSI

}  // namespace Inpar

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
