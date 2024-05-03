/*! \file
\brief Nullspace computation for a Discretization
\level 0
*/

#ifndef FOUR_C_LIB_DISCRET_NULLSPACE_HPP
#define FOUR_C_LIB_DISCRET_NULLSPACE_HPP

#include "4C_config.hpp"

#include <Epetra_MultiVector.h>
#include <Teuchos_RCPDecl.hpp>

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;
}

namespace DRT
{
  /*!
   \brief Calculate the nullspace based on a given discretization

  The nullspace is build by looping over all nodes of a discretization and stored
          in the respective variable.

     \param dis (in): discretization
     \param numdf (in): number of degrees of freedom
     \param dimns (in): nullspace dimension
     \param map (in): nullspace map
      */
  Teuchos::RCP<Epetra_MultiVector> ComputeNullSpace(const DRT::Discretization& dis, const int numdf,
      const int dimns, const Teuchos::RCP<Epetra_Map> dofmap);
}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif