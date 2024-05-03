/*----------------------------------------------------------------------*/
/*! \file
 \brief  base class for partitioned porous multiphase flow through elastic medium problems

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROMULTIPHASE_PARTITIONED_HPP
#define FOUR_C_POROMULTIPHASE_PARTITIONED_HPP

#include "4C_config.hpp"

#include "4C_poromultiphase_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace POROMULTIPHASE
{
  //! Base class of all solid-scatra algorithms
  class PoroMultiPhasePartitioned : public PoroMultiPhaseBase
  {
   public:
    /// create using a Epetra_Comm
    PoroMultiPhasePartitioned(const Epetra_Comm& comm,
        const Teuchos::ParameterList& globaltimeparams);  // Problem builder


  };  // PoroMultiPhasePartitioned


}  // namespace POROMULTIPHASE


FOUR_C_NAMESPACE_CLOSE

#endif