/*----------------------------------------------------------------------*/
/*! \file

 \brief porous medium algorithm with matrix split for condensation of
      no-penetration constraint

\level 2

 *----------------------------------------------------------------------*/


#ifndef FOUR_C_POROELAST_MONOLITHICSPLIT_NOPENETRATION_HPP
#define FOUR_C_POROELAST_MONOLITHICSPLIT_NOPENETRATION_HPP

#include "baci_config.hpp"

#include "baci_poroelast_monolithicsplit.hpp"

BACI_NAMESPACE_OPEN

namespace CORE::LINALG
{
  class MatrixRowTransform;
  class MatrixColTransform;
  class MatrixRowColTransform;
}  // namespace CORE::LINALG

namespace ADAPTER
{
  class CouplingNonLinMortar;
}

namespace POROELAST
{
  //! monolithic structure split for condensing DOFs, when using the brinkman-equation
  class MonolithicSplitNoPenetration : public MonolithicSplit
  {
   public:
    //! create using a Epetra_Comm
    explicit MonolithicSplitNoPenetration(const Epetra_Comm& comm,
        const Teuchos::ParameterList& timeparams,
        Teuchos::RCP<CORE::LINALG::MapExtractor> porosity_splitter);

    //! Setup the monolithic system
    void SetupSystem() override;

    //! setup composed right hand side from field solvers
    void SetupRHS(bool firstcall = false) override;

    //! setup composed system matrix from field solvers
    void SetupSystemMatrix(CORE::LINALG::BlockSparseMatrixBase& mat) override;

    //! take current results for converged and save for next time step
    void Update() override;

    //! read restart data
    void ReadRestart(const int step) override;

    //! contains text to PrintNewtonIter
    void PrintNewtonIterTextStream(std::ostringstream& oss) override;

    //! contains header to PrintNewtonIter
    void PrintNewtonIterHeaderStream(std::ostringstream& oss) override;


   protected:
    //! @name Apply current field state to system

    //! Evaluate mechanical-fluid system matrix
    void ApplyStrCouplMatrix(
        Teuchos::RCP<CORE::LINALG::SparseOperator> k_sf  //!< mechanical-fluid stiffness matrix
        ) override;

    //! Evaluate fluid-mechanical system matrix
    void ApplyFluidCouplMatrix(
        Teuchos::RCP<CORE::LINALG::SparseOperator> k_fs  //!< fluid-mechanical tangent matrix
        ) override;

    //! recover Lagrange multiplier \f$\lambda_\Gamma\f$ at the interface at the end of each
    //! iteration step (i.e. condensed forces onto the structure) needed for rhs in next newton step
    void RecoverLagrangeMultiplierAfterNewtonStep(Teuchos::RCP<const Epetra_Vector> x) override;

    //! recover Lagrange multiplier \f$\lambda_\Gamma\f$ at the interface at the end of each time
    //! step (i.e. condensed forces onto the structure) needed for rhs in next time step
    void RecoverLagrangeMultiplierAfterTimeStep() override;

    //! output
    void Output(bool forced_writerestart = false) override;

    //! setup of coupling object and system matrices
    void SetupCouplingAndMatrices() override;

    //! start a new time step
    void PrepareTimeStep() override;

    //! convergence check for Newton solver
    void BuildConvergenceNorms() override;

    //!@}
   private:
    //! build block vector from field vectors, e.g. rhs, increment vector
    void SetupVector(Epetra_Vector& f,         //!< vector of length of all dofs
        Teuchos::RCP<const Epetra_Vector> sv,  //!< vector containing only structural dofs
        Teuchos::RCP<const Epetra_Vector> fv   //!< vector containing only fluid dofs
        ) override;

    //! @name Global matrices and vectors

    Teuchos::RCP<CORE::LINALG::SparseOperator> k_struct_;
    Teuchos::RCP<CORE::LINALG::SparseOperator> k_fluid_;

    Teuchos::RCP<CORE::LINALG::SparseMatrix> k_lambda_;
    Teuchos::RCP<CORE::LINALG::SparseMatrix> k_D_;
    Teuchos::RCP<CORE::LINALG::SparseMatrix> k_invD_;

    Teuchos::RCP<CORE::LINALG::SparseMatrix> k_Dn_;

    Teuchos::RCP<CORE::LINALG::SparseMatrix> k_lambdainvD_;

    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> k_porodisp_;
    Teuchos::RCP<CORE::LINALG::SparseOperator> k_porofluid_;

    Teuchos::RCP<Epetra_Vector> nopenetration_rhs_;

    //! transform object for k_D matrix \f$D\f$
    Teuchos::RCP<CORE::LINALG::MatrixColTransform> k_D_transform_;
    //! transform object for k_D matrix \f$D\f$
    Teuchos::RCP<CORE::LINALG::MatrixRowTransform> k_invD_transform_;

    //! transform object for linearization of k_D matrix \f$D\f$
    Teuchos::RCP<CORE::LINALG::MatrixColTransform> k_DLin_transform_;

    //! Lagrange multiplier \f$\lambda_\Gamma^{n+1}\f$ at the interface (ie condensed forces onto
    //! the structure) evaluated at actual iteration step \f$t_{n+1}\f$ but needed for next
    //! iteration step
    Teuchos::RCP<Epetra_Vector> lambdanp_;

    //!@}

    //! @name Some quantities to recover the Langrange multiplier at the end of each iteration step

    //! block \f$F_{\Gamma I,i+1}\f$ of fluid matrix at current iteration \f$i+1\f$
    Teuchos::RCP<const CORE::LINALG::SparseMatrix> fgicur_;

    //! block \f$S_{\Gamma\Gamma,i+1}\f$ of fluid matrix at current iteration \f$i+1\f$
    Teuchos::RCP<const CORE::LINALG::SparseMatrix> fggcur_;

    //! block \f$Cfs_{\Gamma\Gamma,i+1}\f$ of fs-coupling matrix at current iteration \f$i+1\f$
    Teuchos::RCP<const CORE::LINALG::SparseMatrix> cfsgicur_;

    //! block \f$Cfs_{\Gamma\Gamma,i+1}\f$ of fs-coupling matrix at current iteration \f$i+1\f$
    Teuchos::RCP<const CORE::LINALG::SparseMatrix> cfsggcur_;

    Teuchos::RCP<Epetra_Vector> rhs_fgcur_;

    //!@}

    //! @name Iterative solution technique

    //!< norm of no-penetration constraint
    double normrhs_nopenetration_;
    //!@}

    Teuchos::RCP<ADAPTER::CouplingNonLinMortar> mortar_adapter_;
  };

}  // namespace POROELAST


BACI_NAMESPACE_CLOSE

#endif