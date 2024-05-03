/*-----------------------------------------------------------*/
/*! \file

\brief Base class for modelevaluators in partitioned algorithms.


\level 3
*/
/*-----------------------------------------------------------*/


#ifndef FOUR_C_STRUCTURE_NEW_MODEL_EVALUATOR_MULTIPHYSICS_HPP
#define FOUR_C_STRUCTURE_NEW_MODEL_EVALUATOR_MULTIPHYSICS_HPP


#include "4C_config.hpp"

#include "4C_structure_new_model_evaluator_generic.hpp"

// forward declaration
class Epetra_Vector;
class Epetra_Map;
namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos

FOUR_C_NAMESPACE_OPEN

namespace CORE::LINALG
{
  class SparseOperator;
}  // namespace CORE::LINALG
namespace IO
{
  class DiscretizationWriter;
  class DiscretizationReader;
}  // namespace IO
namespace DRT
{
  class Discretization;
}  // namespace DRT
namespace STR
{
  class Integrator;
  namespace TIMINT
  {
    class BaseDataGlobalState;
    class BaseDataIO;
    class Base;
  }  // namespace TIMINT
  namespace MODELEVALUATOR
  {
    class Data;


    //! supported multiphysic problems
    enum MultiphysicType
    {
      mt_none = 0,  //!< none specific default value
      mt_fsi = 1,   //!< multiphysics type fluid-structure-interaction
      mt_ssi = 2    //!< multiphysics type structure-scalar-interaction
    };              // MultiphysicType


    /*! \brief This is the base class for all multiphysics models.
     *
     *  This class summarizes the functionality which all model multiphysics model
     *  evaluators share.
     *
     *  \date 11/16
     *  \author rauch */
    class Multiphysics : public Generic
    {
     public:
      //! constructor
      Multiphysics();


      //! initialize the class variables
      void Init(const Teuchos::RCP<STR::MODELEVALUATOR::Data>& eval_data_ptr,
          const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr,
          const Teuchos::RCP<STR::TIMINT::BaseDataIO>& gio_ptr,
          const Teuchos::RCP<STR::Integrator>& int_ptr,
          const Teuchos::RCP<const STR::TIMINT::Base>& timint_ptr, const int& dof_offset) override;

      //! setup class variables
      void Setup() override;

      //! set the active model type wrapped in this class.
      //! only active model type is evaluated.
      //! e.g. mt_fsi in case fluid-structure interaction is to be evaluated
      void SetActiveModelType(enum STR::MODELEVALUATOR::MultiphysicType mtype)
      {
        active_mt_ = mtype;
      };

      void CheckActiveModelType() const
      {
        if (active_mt_ == mt_none) FOUR_C_THROW("No active model evaluator set for Multiphysics");
      };

      //! @name Functions which are derived from the base generic class
      //! @{
      //! [derived]
      INPAR::STR::ModelType Type() const override { return INPAR::STR::model_partitioned_coupling; }

      //! reset class variables (without jacobian) [derived]
      void Reset(const Epetra_Vector& x) override;

      //! [derived]
      bool EvaluateForce() override;

      //! [derived]
      bool EvaluateStiff() override;

      //! [derived] not needed in partitioned scheme
      bool EvaluateForceStiff() override;

      //! derived
      void PreEvaluate() override{};

      //! derived
      void PostEvaluate() override{};

      //! derived
      bool AssembleForce(Epetra_Vector& f, const double& timefac_np) const override;

      //! Assemble the jacobian at \f$t_{n+1}\f$ not needed in partitioned scheme
      bool AssembleJacobian(
          CORE::LINALG::SparseOperator& jac, const double& timefac_np) const override;

      //! [derived]
      void WriteRestart(
          IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const override{};

      //! [derived]
      void ReadRestart(IO::DiscretizationReader& ioreader) override{};

      //! [derived]
      void Predict(const INPAR::STR::PredEnum& pred_type) override{};

      //! derived
      void RunPreComputeX(const Epetra_Vector& xold, Epetra_Vector& dir_mutable,
          const NOX::NLN::Group& curr_grp) override{};

      //! recover condensed Lagrange multipliers
      void RunPostComputeX(const Epetra_Vector& xold, const Epetra_Vector& dir,
          const Epetra_Vector& xnew) override{};

      //! derived
      void RunPostIterate(const ::NOX::Solver::Generic& solver) override{};

      //! [derived]
      void UpdateStepState(const double& timefac_n) override;

      //! [derived]
      void UpdateStepElement() override{};

      //! [derived]
      void DetermineStressStrain() override{};

      //! [derived]
      void DetermineEnergy() override{};

      //! [derived]
      void DetermineOptionalQuantity() override{};

      //! [derived]
      void OutputStepState(IO::DiscretizationWriter& iowriter) const override{};

      //! derived
      void ResetStepState() override{};

      //! [derived]
      void PostOutput() override{};

      //! @name Accessors to model specific things
      //! @{

      //! Returns a pointer to the model specific dof row map
      Teuchos::RCP<const Epetra_Map> GetBlockDofRowMapPtr() const override
      {
        return Teuchos::null;
      };

      //! Returns a pointer to the current model solution vector (usually the Lagrange multiplier
      //! vector)
      Teuchos::RCP<const Epetra_Vector> GetCurrentSolutionPtr() const override
      {
        return Teuchos::null;
      };

      //! Returns a pointer to the model solution vector of the last time step (usually the Lagrange
      //! multiplier vector)
      Teuchos::RCP<const Epetra_Vector> GetLastTimeStepSolutionPtr() const override
      {
        return Teuchos::null;
      };

      //! @}

     protected:
      //! map containing the model evaluators of the sub modules
      std::map<enum STR::MODELEVALUATOR::MultiphysicType,
          Teuchos::RCP<STR::MODELEVALUATOR::Generic>>
          me_map_;

      //! currently active model evaluator type
      STR::MODELEVALUATOR::MultiphysicType active_mt_;

      //! return reference to map containing the model evaluators
      std::map<enum STR::MODELEVALUATOR::MultiphysicType,
          Teuchos::RCP<STR::MODELEVALUATOR::Generic>>&
      GetModelEvalutaorMap()
      {
        return me_map_;
      };

     public:
      //! return RCP to model evaluator of specific MultiphysicType
      Teuchos::RCP<STR::MODELEVALUATOR::Generic> GetModelEvaluatorFromMap(
          enum STR::MODELEVALUATOR::MultiphysicType mtype) const
      {
        return me_map_.at(mtype);
      }


    };  // class Multiphysics

  }  // namespace MODELEVALUATOR
}  // namespace STR

FOUR_C_NAMESPACE_CLOSE

#endif