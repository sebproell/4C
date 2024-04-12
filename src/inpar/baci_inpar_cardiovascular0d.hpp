/*----------------------------------------------------------------------*/
/*! \file

\brief Input parameters for 0d cardiovascular-structure coupling

\level 2

*/

/*----------------------------------------------------------------------*/

#ifndef FOUR_C_INPAR_CARDIOVASCULAR0D_HPP
#define FOUR_C_INPAR_CARDIOVASCULAR0D_HPP

#include "baci_config.hpp"

#include "baci_utils_parameter_list.hpp"

BACI_NAMESPACE_OPEN

// forward declaration
namespace INPUT
{
  class ConditionDefinition;
}

namespace INPAR
{
  namespace CARDIOVASCULAR0D
  {
    /// possible 0D cardiovascular-structural solvers
    enum Cardvasc0DSolveAlgo
    {
      cardvasc0dsolve_direct,  ///< build monolithic 0D cardiovascular-structural system
      cardvasc0dsolve_simple,  ///< use simple preconditioner for iterative solve
      cardvasc0dsolve_AMGnxn
    };

    enum Cardvasc0DAtriumModel
    {
      atr_prescribed,
      atr_elastance_0d,
      atr_structure_3d
    };

    enum Cardvasc0DVentricleModel
    {
      ventr_prescribed,
      ventr_elastance_0d,
      ventr_structure_3d
    };

    enum Cardvasc0DRespiratoryModel
    {
      resp_none,
      resp_standard
    };

    /// set the 0Dcardiovascular parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific 0Dcardiovascular conditions
    void SetValidConditions(std::vector<Teuchos::RCP<INPUT::ConditionDefinition>>& condlist);

  }  // namespace CARDIOVASCULAR0D
}  // namespace INPAR
/*----------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif