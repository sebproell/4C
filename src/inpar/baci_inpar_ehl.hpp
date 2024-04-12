/*--------------------------------------------------------------------------*/
/*! \file

\brief Elastohydrodynamic lubrication (lubrication structure interaction) parameters

\level 3

*/
/*--------------------------------------------------------------------------*/
#ifndef FOUR_C_INPAR_EHL_HPP
#define FOUR_C_INPAR_EHL_HPP


#include "baci_config.hpp"

#include "baci_utils_exceptions.hpp"
#include "baci_utils_parameter_list.hpp"

BACI_NAMESPACE_OPEN

// forward declaration
namespace INPUT
{
  class ConditionDefinition;
}

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
namespace INPAR
{
  namespace EHL
  {
    /// Type of coupling strategy for EHL problems
    enum SolutionSchemeOverFields
    {
      ehl_IterStagg,
      ehl_Monolithic
    };

    /// Type of coupling strategy between the two fields of the EHL problems
    enum FieldCoupling
    {
      coupling_none,
      coupling_matching
    };

    //! type of norm to check for convergence
    enum ConvNorm
    {
      convnorm_abs,  //!< absolute norm
      convnorm_rel,  //!< relative norm of EHL problem with inital EHL rhs
      convnorm_mix   //!< mixed absolute-relative norm
    };

    //! type of norm to check for convergence
    enum BinaryOp
    {
      bop_and,              //!< and
      bop_or,               //!< or
      bop_coupl_or_singl,   //!< either EHL problem or single field problems converged
      bop_coupl_and_singl,  //!< either EHL problem or single field problems converged
      bop_and_singl,        //!< and in single field problems
      bop_or_singl          //!< or in single field problems
    };

    //! type of solution techniques
    enum NlnSolTech
    {
      soltech_newtonfull  //!< full Newton-Raphson iteration
    };

    //! Map solution technique enum to std::string
    static inline std::string NlnSolTechString(const enum NlnSolTech name  //!< enum to convert
    )
    {
      switch (name)
      {
        case soltech_newtonfull:
          return "fullnewton";
          break;
        default:
          dserror("Cannot make std::string for solution technique %d", name);
          return "";
      }
    }

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


    /// set the ehl parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific ehl conditions
    void SetValidConditions(std::vector<Teuchos::RCP<INPUT::ConditionDefinition>>& condlist);

  }  // namespace EHL

}  // namespace INPAR


BACI_NAMESPACE_CLOSE

#endif