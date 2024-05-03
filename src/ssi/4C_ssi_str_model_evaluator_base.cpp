/*----------------------------------------------------------------------*/
/*! \file
\brief structural model evaluator for scalar-structure interaction

\level 2

    */
/*----------------------------------------------------------------------*/
#include "4C_ssi_str_model_evaluator_base.hpp"

#include "4C_adapter_str_ssiwrapper.hpp"
#include "4C_comm_exporter.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_discretization_fem_general_utils_gauss_point_postprocess.hpp"
#include "4C_io.hpp"
#include "4C_lib_utils_gid_vector.hpp"
#include "4C_structure_new_model_evaluator_data.hpp"
#include "4C_structure_new_timint_basedataglobalstate.hpp"

#include <Epetra_IntVector.h>
#include <Epetra_Vector.h>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BaseSSI::DetermineStressStrain()
{
  // extract raw data for element-wise stresses
  const std::vector<char>& stressdata = EvalData().StressData();

  // initialize map for element-wise stresses
  const auto stresses =
      Teuchos::rcp(new std::map<int, Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>>);

  // initialize position pointer
  std::vector<char>::size_type position(0);

  // loop over all row elements
  for (int i = 0; i < Discret().ElementRowMap()->NumMyElements(); ++i)
  {
    // initialize matrix for stresses associated with current element
    const auto stresses_ele = Teuchos::rcp(new CORE::LINALG::SerialDenseMatrix);

    // extract stresses
    CORE::COMM::ParObject::ExtractfromPack(position, stressdata, *stresses_ele);

    // store stresses
    (*stresses)[Discret().ElementRowMap()->GID(i)] = stresses_ele;
  }

  // export map to column format
  CORE::COMM::Exporter exporter(
      *Discret().ElementRowMap(), *Discret().ElementColMap(), Discret().Comm());
  exporter.Export(*stresses);

  // prepare nodal stress vectors
  Epetra_MultiVector nodal_stresses_source(*Discret().NodeRowMap(), 6);

  Discret().Evaluate(
      [&](DRT::Element& ele)
      {
        CORE::FE::ExtrapolateGaussPointQuantityToNodes(
            ele, *stresses->at(ele.Id()), Discret(), nodal_stresses_source);
      });

  const auto* nodegids = Discret().NodeRowMap();
  for (int i = 0; i < nodegids->NumMyElements(); ++i)
  {
    const int nodegid = nodegids->GID(i);

    // extract lid of node as multi-vector is sorted according to the node ids
    const DRT::Node* const node = Discret().gNode(nodegid);
    const int nodelid = Discret().NodeRowMap()->LID(nodegid);

    // extract dof lid of first degree of freedom associated with current node in second nodeset
    const int dofgid_epetra = Discret().Dof(2, node, 0);
    const int doflid_epetra = mechanical_stress_state_->Map().LID(dofgid_epetra);
    if (doflid_epetra < 0) FOUR_C_THROW("Local ID not found in epetra vector!");

    (*mechanical_stress_state_)[doflid_epetra] = (*nodal_stresses_source(0))[nodelid];
    (*mechanical_stress_state_)[doflid_epetra + 1] = (*nodal_stresses_source(1))[nodelid];
    (*mechanical_stress_state_)[doflid_epetra + 2] = (*nodal_stresses_source(2))[nodelid];
    (*mechanical_stress_state_)[doflid_epetra + 3] = (*nodal_stresses_source(3))[nodelid];
    (*mechanical_stress_state_)[doflid_epetra + 4] = (*nodal_stresses_source(4))[nodelid];
    (*mechanical_stress_state_)[doflid_epetra + 5] = (*nodal_stresses_source(5))[nodelid];
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> STR::MODELEVALUATOR::BaseSSI::GetBlockDofRowMapPtr() const
{
  CheckInitSetup();
  return GState().DofRowMap();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BaseSSI::Setup()
{
  // check initialization
  CheckInit();

  if (Discret().NumDofSets() - 1 == 2)
    mechanical_stress_state_ = Teuchos::rcp(new Epetra_Vector(*Discret().DofRowMap(2), true));

  // set flag
  issetup_ = true;
}
FOUR_C_NAMESPACE_CLOSE