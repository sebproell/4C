/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_scatra_monolithic_twoway.cpp

 \brief two-way coupled monolithic algorithm for scalar transport within multiphase porous medium

   \level 3

   \maintainer  Johannes Kremheller
                kremheller@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/




#include "poromultiphase_scatra_monolithic_twoway.H"
#include <Teuchos_TimeMonitor.hpp>


#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_poromultiphase/poromultiphase_base.H"
#include "../drt_adapter/ad_porofluidmultiphase_wrapper.H"
#include "../drt_adapter/ad_str_structure.H"

#include "../drt_scatra_ele/scatra_ele_action.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_assemblestrategy.H"

#include "../drt_adapter/adapter_scatra_base_algorithm.H"
#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_io/io_control.H"
#include "../linalg/linalg_solver.H"
#include "../drt_inpar/inpar_solver.H"

#include "../linalg/linalg_utils.H"

/*----------------------------------------------------------------------*
 | constructor                                              vuong 08/16 |
 *----------------------------------------------------------------------*/
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraMonolithicTwoWay::PoroMultiPhaseScaTraMonolithicTwoWay(
    const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams):
    PoroMultiPhaseScaTraMonolithic(comm, globaltimeparams),
    directsolve_(true),
    ittolinc_(0.0),
    ittolres_(0.0),
    itmax_(0),
    itmin_(1),
    itnum_(0),
    blockrowdofmap_(Teuchos::null),
    tolinc_(0.0),
    tolfres_(0.0),
    tolinc_struct_(0.0),
    tolfres_struct_(0.0),
    tolinc_fluid_(0.0),
    tolfres_fluid_(0.0),
    tolinc_scatra_(0.0),
    tolfres_scatra_(0.0),
    normrhs_(0.0),
    norminc_(0.0),
    normrhsfluid_(0.0),
    normincfluid_(0.0),
    normrhsstruct_(0.0),
    normincstruct_(0.0),
    normrhsscatra_(0.0),
    normincscatra_(0.0),
    vectornormfres_(INPAR::POROMULTIPHASESCATRA::norm_undefined),
    vectornorminc_(INPAR::POROMULTIPHASESCATRA::norm_undefined),
    timernewton_(comm),
    dtsolve_(0.0),
    dtele_(0.0),
    fdcheck_(INPAR::POROMULTIPHASESCATRA::FDCheck::fdcheck_none)
{

}

/*----------------------------------------------------------------------*
 | initialization                                            vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraMonolithicTwoWay::Init(
    const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& algoparams,
    const Teuchos::ParameterList& poroparams,
    const Teuchos::ParameterList& structparams,
    const Teuchos::ParameterList& fluidparams,
    const Teuchos::ParameterList& scatraparams,
    const std::string& struct_disname,
    const std::string& fluid_disname,
    const std::string& scatra_disname,
    bool isale,
    int nds_disp,
    int nds_vel,
    int nds_solidpressure,
    int ndsporofluid_scatra)
{
  //call base class
  POROMULTIPHASESCATRA::PoroMultiPhaseScaTraMonolithic::Init(
      globaltimeparams,
      algoparams,
      poroparams,
      structparams,
      fluidparams,
      scatraparams,
      struct_disname,
      fluid_disname,
      scatra_disname,
      isale,
      nds_disp,
      nds_vel,
      nds_solidpressure,
      ndsporofluid_scatra);

  // read input variables
  itmax_ = algoparams.get<int>("ITEMAX");
  ittolinc_ = algoparams.get<double>("TOLINC_GLOBAL");
  ittolres_ = algoparams.get<double>("TOLRES_GLOBAL");

  blockrowdofmap_ = Teuchos::rcp(new LINALG::MultiMapExtractor);

  fdcheck_ = DRT::INPUT::IntegralValue<INPAR::POROMULTIPHASESCATRA::FDCheck>(algoparams,"FDCHECK");

}

/*----------------------------------------------------------------------*
 | setup the fully monolithic fluid-structure-scatra system (called in  |
 | poromultiphase_scatra_dyn.cpp)                      kremheller 03/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraMonolithicTwoWay::SetupSystem()
{
  //setup the poro subsystem first
  PoroField()->SetupSystem();

  // create combined map
  std::vector<Teuchos::RCP<const Epetra_Map> > vecSpaces;

  {
    vecSpaces.push_back(PoroField()->StructDofRowMap());
    vecSpaces.push_back(PoroField()->FluidDofRowMap());
    const Epetra_Map* dofrowmapscatra = (ScatraAlgo()->ScaTraField()->Discretization())->DofRowMap(0);
    vecSpaces.push_back(Teuchos::rcp(dofrowmapscatra,false));
  }

  if (vecSpaces[0]->NumGlobalElements() == 0)
    dserror("No poro structure equation. Panic.");
  if (vecSpaces[1]->NumGlobalElements() == 0)
    dserror("No poro fluid equation. Panic.");
  if (vecSpaces[2]->NumGlobalElements()==0)
    dserror("No scatra equation. Panic.");

  // full fluid-structure-scatra-map
  fullmap_ = LINALG::MultiMapExtractor::MergeMaps(vecSpaces);

  // full Poromultiphase-elasticity-blockmap
  blockrowdofmap_->Setup(*fullmap_, vecSpaces);

  //-----------------------------------build map of global dofs with DBC
  BuildCombinedDBCMap();
  // -------------------------------------------------------------

  // initialize Poroscatra-systemmatrix_
  systemmatrix_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
                                      *Extractor(),
                                      *Extractor(),
                                      81,
                                      false,
                                      true));

  //! structure-scatra coupling matrix k_pss_ --> equal to zero so far
  //! fluid-scatra coupling matrix
  k_pfs_ = Teuchos::rcp(new LINALG::SparseMatrix(
                        *(PoroField()->FluidDofRowMap()),
                        //*(FluidField()->DofRowMap()),
                        81, true, true));

  //! scatra-structure coupling matrix
  k_sps_ = Teuchos::rcp(new LINALG::SparseMatrix(
                      *(ScatraAlgo()->ScaTraField()->Discretization()->DofRowMap()), 81, true, true));
  //! scatra-fluid coupling matrix
  k_spf_ = Teuchos::rcp(new LINALG::SparseMatrix(
                        *(ScatraAlgo()->ScaTraField()->Discretization()->DofRowMap()),
                        //*(FluidField()->DofRowMap()),
                        81, true, true));

  return;
}

/*-----------------------------------------------------------------------/
|  build the combined dbcmap                           kremheller 06/17  |
/-----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraMonolithicTwoWay::BuildCombinedDBCMap()
{
  // Combined DBC map of poromultielast-problem
  const Teuchos::RCP<const Epetra_Map> porocondmap =
      PoroField()->CombinedDBCMap();
  const Teuchos::RCP<const Epetra_Map> scatracondmap =
      ScatraAlgo()->ScaTraField()->DirichMaps()->CondMap();
  combinedDBCMap_ = LINALG::MergeMap(porocondmap, scatracondmap, false);

  return;
}

/*----------------------------------------------------------------------*
 | setup the solver if necessary                        kremheller 03/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraMonolithicTwoWay::SetupSolver()
{
  //  solver
  // create a linear solver
  // get dynamic section of poroelasticity
  const Teuchos::ParameterList& poromultscatradyn = DRT::Problem::Instance()->PoroMultiPhaseScatraDynamicParams();
  // get the solver number used for linear poroelasticity solver
  const int linsolvernumber = poromultscatradyn.get<int>("LINEAR_SOLVER");
  // check if the poroelasticity solver has a valid solver number
  if (linsolvernumber == (-1))
    dserror("no linear solver defined for poromultiphaseflow with scatra coupling.\n"
        " Please set LINEAR_SOLVER in POROMULTIPHASESCATRA DYNAMIC to a valid number!");
  const Teuchos::ParameterList& solverparams =
    DRT::Problem::Instance()->SolverParams(linsolvernumber);
  const int solvertype = DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(
    solverparams, "SOLVER");

  directsolve_ = (   solvertype == INPAR::SOLVER::umfpack
                  or solvertype == INPAR::SOLVER::superlu
                  or solvertype == INPAR::SOLVER::amesos_klu_nonsym);

  if (directsolve_)
  {
    solver_ = Teuchos::rcp(new LINALG::Solver( solverparams,
                                      Comm(),
                                      DRT::Problem::Instance()->ErrorFile()->Handle())
                 );
  }
  else
    dserror("only direct solvers work so far");


  vectornormfres_ = DRT::INPUT::IntegralValue<INPAR::POROMULTIPHASESCATRA::VectorNorm>(
      poromultscatradyn, "VECTORNORM_RESF");
  vectornorminc_ = DRT::INPUT::IntegralValue<INPAR::POROMULTIPHASESCATRA::VectorNorm>(
      poromultscatradyn, "VECTORNORM_INC");
}

/*----------------------------------------------------------------------*
 |                                                    kremheller 06/17  |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraMonolithicTwoWay::TimeStep()
{
  // Prepare stuff
  SetupNewton();
  PrintHeader();

  // Evaluate
  Evaluate(iterinc_);

  // Newton-Loop
  while ( (not Converged() and itnum_ < itmax_) or (itnum_ < itmin_) )
  {
    // increment number of iteration
    itnum_++;

    // Solve
    LinearSolve();
    solver_->ResetTolerance();

    // Build Convergence Norms
    BuildConvergenceNorms();

    // Evaluate
    if(not Converged())
    {
      Evaluate(iterinc_);
      // perform FD Check of monolithic system matrix
      if(fdcheck_ == INPAR::POROMULTIPHASESCATRA::fdcheck_global)
        PoroMultiPhaseScaTraFDCheck();
    }
    else
    {
      // convergence check is based on residual(phi_i) < tol and phi_i+1 - phi_i < tol
      // in this function we update phi_i+1 as phi_i+1 = phi_i + iterinc for all fields
      // even though we have not evaluated the residual of phi_i+1 it will still be more exact than
      // the one at phi_i
      UpdateFieldsAfterConvergence();
    }

    // print output
    NewtonOutput();
  }

  // Error-Check
  NewtonErrorCheck();

  return;
}

/*----------------------------------------------------------------------*
 | Evaluate (build global Matrix and RHS)            kremheller 06/17   |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraMonolithicTwoWay::Evaluate(Teuchos::RCP<const Epetra_Vector> x)
{
  TEUCHOS_FUNC_TIME_MONITOR("POROMULTIPHASESCATRA::PoroMultiPhaseScaTraMonolithicTwoWay::Evaluate");

  // reset timer
  timernewton_.ResetStartTime();
  // *********** time measurement ***********
  double dtcpu = timernewton_.WallTime();
  // *********** time measurement ***********

  // displacement, fluid variable and scatra variable incremental vector
  Teuchos::RCP<const Epetra_Vector> porostructinc;
  Teuchos::RCP<const Epetra_Vector> porofluidinc;
  Teuchos::RCP<const Epetra_Vector> scatrainc;
  ExtractFieldVectors(x, porostructinc, porofluidinc, scatrainc);

  // (1) Newton update of the scatra field
  ScatraAlgo()->ScaTraField()->UpdateIter(scatrainc);

  // (2) set scatra solution on fluid field
  SetScatraSolution();

  // (3) access poro problem to build poro-poro block
  PoroField()->Evaluate(porostructinc,porofluidinc,itnum_==0);

  // (4) set fluid and structure solution on scatra field
  SetPoroSolution();

  // (5) access ScaTra problem to build scatra-scatra block
  ScatraAlgo()->ScaTraField()->PrepareLinearSolve();

  // (6) Build the monolithic system matrix
  SetupSystemMatrix();

  // check whether we have a sanely filled tangent matrix
  if (not systemmatrix_->Filled())
  {
    dserror("Effective tangent matrix must be filled here");
  }

  // (7) Build the monolithic system vector
  SetupRHS();

  // *********** time measurement ***********
  dtele_ = timernewton_.WallTime() - dtcpu;
  // *********** time measurement ***********

}

/*----------------------------------------------------------------------*
 | setup monolithic system matrix                      kremheller 06/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraMonolithicTwoWay::SetupSystemMatrix()
{
  // set loma block matrix to zero
  systemmatrix_->Zero();

  //----------------------------------------------------------------------
  // 1st diagonal block (upper left): poro weighting - poro solution
  // has dimensions ((ndim+n_phases)*n_nodes)x((ndim+n_phases)*n_nodes)
  //----------------------------------------------------------------------
  // get matrix block
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> mat_pp = PoroField()->BlockSystemMatrix();

  // uncomplete matrix block (appears to be required in certain cases)
  mat_pp->UnComplete();

  // assign matrix block
  systemmatrix_->Assign(0,0,LINALG::View,mat_pp->Matrix(0,0));
  systemmatrix_->Assign(0,1,LINALG::View,mat_pp->Matrix(0,1));
  systemmatrix_->Assign(1,0,LINALG::View,mat_pp->Matrix(1,0));
  systemmatrix_->Assign(1,1,LINALG::View,mat_pp->Matrix(1,1));

  //----------------------------------------------------------------------
  // 2nd diagonal block (lower right): scatra weighting - scatra solution
  // has dimensions (n_species*n_nodes)x(n_species*n_nodes)
  //----------------------------------------------------------------------
  // get matrix block
  Teuchos::RCP<LINALG::SparseMatrix> mat_ss = ScatraAlgo()->ScaTraField()->SystemMatrix();

  // uncomplete matrix block (appears to be required in certain cases)
  mat_ss->UnComplete();

  // assign matrix block
  systemmatrix_->Assign(2,2,LINALG::View,*mat_ss);

  // complete scatra block matrix
  systemmatrix_->Complete();

  //----------------------------------------------------------------------
  // 1st off-diagonal block k_ps (upper right): poro weighting - scatra solution
  // has dimensions ((ndim+n_phases)*n_nodes)x(n_species*n_nodes)
  // so far no coupling of structure with scatra --> k_pss_ = 0
  // --> dimensions (n_phases*n_nodes)x(n_species*n_nodes)
  //----------------------------------------------------------------------

  // create empty matrix
  Teuchos::RCP<LINALG::SparseMatrix> k_pfs = PoroFluidScatraCouplingMatrix();

  // call the porofluid-elements and calculate the off-diagonal scatra matrix block
  ApplyPoroFluidScatraCouplMatrix(k_pfs);

  // apply DBC's also on off-diagonal fluid-scatra coupling block (main-diagonal blocks have already been set, either in
  // poromultielast_monolithic.cpp or in the respective evalute calls)
  k_pfs->ApplyDirichlet(*PoroField()->FluidField()->GetDBCMapExtractor()->CondMap(),false);

  // uncomplete matrix block (appears to be required in certain cases)
  //k_pss_->UnComplete();
  k_pfs->UnComplete();

  // assign matrix block
  //systemmatrix_->Assign(0,2,LINALG::View,*(k_pss_)); --> zero
  systemmatrix_->Assign(1,2,LINALG::View,*(k_pfs));

  //----------------------------------------------------------------------
  // 2nd off-diagonal block k_sp (lower left): scatra weighting - poro solution
  // has dimensions (n_species*n_nodes)x((ndim+n_phases)*n_nodes)
  //----------------------------------------------------------------------

  // create empty matrix
  Teuchos::RCP<LINALG::SparseMatrix> k_sps = ScatraStructCouplingMatrix();

  // call the scatra-elements and calculate the off-diagonal structure matrix block
  ApplyScatraStructCouplMatrix(k_sps);

  // apply DBC's also on off-diagonal scatra-structure coupling block (main-diagonal blocks have already been set, either in
  // poromultielast_monolithic.cpp or in the respective evalute calls)
  k_sps->ApplyDirichlet(*ScatraAlgo()->ScaTraField()->DirichMaps()->CondMap(),false);

  // create empty matrix
  Teuchos::RCP<LINALG::SparseMatrix> k_spf = ScatraPoroFluidCouplingMatrix();

  // call the scatra-elements and calculate the off-diagonal structure matrix block
  ApplyScatraPoroFluidCouplMatrix(k_spf);

  // apply DBC's also on off-diagonal scatra-fluid coupling block (main-diagonal blocks have already been set, either in
  // poromultielast_monolithic.cpp or in the respective evalute calls)
  k_spf->ApplyDirichlet(*ScatraAlgo()->ScaTraField()->DirichMaps()->CondMap(),false);

  // uncomplete matrix block (appears to be required in certain cases)
  k_sps->UnComplete();
  k_spf->UnComplete();

  // assign matrix block
  systemmatrix_->Assign(2,0,LINALG::View,*(k_sps));
  systemmatrix_->Assign(2,1,LINALG::View,*(k_spf));

  // complete block matrix
  systemmatrix_->Complete();

  // Debug: matlab output of system matrix
  bool matlab = false;
  if (matlab)
  {
    std::string filename = "../o/mymatrix.dat";
    LINALG::PrintBlockMatrixInMatlabFormat(filename, *systemmatrix_);
    dserror("exit");
  }
}

/*----------------------------------------------------------------------*
 | get porofluid-scatra coupling sparse matrix         kremheller 06/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> POROMULTIPHASESCATRA::PoroMultiPhaseScaTraMonolithicTwoWay::PoroFluidScatraCouplingMatrix()
{
  Teuchos::RCP<LINALG::SparseMatrix> sparse = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_pfs_);
  if(sparse==Teuchos::null)
    dserror("cast to LINALG::SparseMatrix failed!");

  return sparse;
}

/*----------------------------------------------------------------------*
 | get scatra-structure coupling sparse matrix         kremheller 07/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> POROMULTIPHASESCATRA::PoroMultiPhaseScaTraMonolithicTwoWay::ScatraStructCouplingMatrix()
{
  Teuchos::RCP<LINALG::SparseMatrix> sparse = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_sps_);
  if(sparse==Teuchos::null)
    dserror("cast to LINALG::SparseMatrix failed!");

  return sparse;
}

/*----------------------------------------------------------------------*
 | get scatra-fluid coupling sparse matrix             kremheller 07/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> POROMULTIPHASESCATRA::PoroMultiPhaseScaTraMonolithicTwoWay::ScatraPoroFluidCouplingMatrix()
{
  Teuchos::RCP<LINALG::SparseMatrix> sparse = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_spf_);
  if(sparse==Teuchos::null)
    dserror("cast to LINALG::SparseMatrix failed!");

  return sparse;
}

/*----------------------------------------------------------------------*
 | evaluate porofluid-scatra coupling sparse matrix    kremheller 06/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraMonolithicTwoWay::ApplyPoroFluidScatraCouplMatrix(
    Teuchos::RCP< LINALG::SparseOperator> k_pfs //!< off-diagonal tangent matrix term
  )
{

  //reset
  k_pfs->Zero();
  //evaluate
  PoroField()->FluidField()->AssembleFluidScatraCouplingMat(k_pfs);
  //complete
  k_pfs->Complete(ScatraAlgo()->ScaTraField()->SystemMatrix()->RangeMap(), PoroField()->FluidField()->SystemMatrix()->RangeMap());

  return;
}

/*----------------------------------------------------------------------*
 | evaluate scatra-structure coupling sparse matrix    kremheller 07/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraMonolithicTwoWay::ApplyScatraStructCouplMatrix(
    Teuchos::RCP< LINALG::SparseOperator> k_sps //!< off-diagonal tangent matrix term
  )
{
  // create the parameters for the discretization
  Teuchos::ParameterList sparams_struct;

  k_sps->Zero();

  sparams_struct.set<int>("action", SCATRA::calc_scatra_mono_odblock_mesh);
  // other parameters that might be needed by the elements
  sparams_struct.set("delta time", Dt());
  sparams_struct.set("total time", Time());

  // provide element parameter list with numbers of dofsets associated with displacement and velocity dofs on scatra discretization
  sparams_struct.set<int>("ndsdisp",ScatraAlgo()->ScaTraField()->NdsDisp());
  sparams_struct.set<int>("ndsvel",ScatraAlgo()->ScaTraField()->NdsVel());
  sparams_struct.set<int>("ndspres",2);

  ScatraAlgo()->ScaTraField()->Discretization()->ClearState();
  ScatraAlgo()->ScaTraField()->Discretization()->SetState(0,"hist",ScatraAlgo()->ScaTraField()->Hist());
  ScatraAlgo()->ScaTraField()->Discretization()->SetState(0,"phinp",ScatraAlgo()->ScaTraField()->Phinp());

  // build specific assemble strategy for mechanical-fluid system matrix
  // from the point of view of StructureField:
  // structdofset = 0, fluiddofset = 1
  DRT::AssembleStrategy scatrastrategy_struct(
      0,               // scatradofset for row
      1,               // structuredofset for column
      k_sps,          // scatra-structure coupling matrix
      Teuchos::null ,
      Teuchos::null ,
      Teuchos::null,
      Teuchos::null
  );

  ScatraAlgo()->ScaTraField()->Discretization()->Evaluate(sparams_struct, scatrastrategy_struct);

  //complete
  k_sps->Complete(PoroField()->StructureField()->SystemMatrix()->RangeMap(), ScatraAlgo()->ScaTraField()->SystemMatrix()->RangeMap());

  ScatraAlgo()->ScaTraField()->Discretization()->ClearState();

  return;
}

/*----------------------------------------------------------------------*
 | evaluate scatra-porofluid coupling sparse matrix    kremheller 07/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraMonolithicTwoWay::ApplyScatraPoroFluidCouplMatrix(
    Teuchos::RCP< LINALG::SparseOperator> k_spf //!< off-diagonal tangent matrix term
  )
{
  // create the parameters for the discretization
  Teuchos::ParameterList sparams_fluid;

  k_spf->Zero();

  sparams_fluid.set<int>("action", SCATRA::calc_scatra_mono_odblock_fluid);
  // other parameters that might be needed by the elements
  sparams_fluid.set("delta time", Dt());
  sparams_fluid.set("total time", Time());

  // provide element parameter list with numbers of dofsets associated with displacement and velocity dofs on scatra discretization
  sparams_fluid.set<int>("ndsdisp",ScatraAlgo()->ScaTraField()->NdsDisp());
  sparams_fluid.set<int>("ndsvel",ScatraAlgo()->ScaTraField()->NdsVel());
  sparams_fluid.set<int>("ndspres",2);

  ScatraAlgo()->ScaTraField()->Discretization()->ClearState();
  ScatraAlgo()->ScaTraField()->Discretization()->SetState(0,"hist",ScatraAlgo()->ScaTraField()->Hist());
  ScatraAlgo()->ScaTraField()->Discretization()->SetState(0,"phinp",ScatraAlgo()->ScaTraField()->Phinp());


  // build specific assemble strategy for mechanical-fluid system matrix
  // from the point of view of StructureField:
  // structdofset = 0, fluiddofset = 1
  DRT::AssembleStrategy scatrastrategy_fluid(
      0,               // scatradofset for row
      2,               // fluiddofset for column
      k_spf,          // scatra-structure coupling matrix
      Teuchos::null ,
      Teuchos::null ,
      Teuchos::null,
      Teuchos::null
  );

  ScatraAlgo()->ScaTraField()->Discretization()->Evaluate(sparams_fluid, scatrastrategy_fluid);

  //complete
  k_spf->Complete(PoroField()->FluidField()->SystemMatrix()->RangeMap(), ScatraAlgo()->ScaTraField()->SystemMatrix()->RangeMap());

  ScatraAlgo()->ScaTraField()->Discretization()->ClearState();

  return;
}

/*-----------------------------------------------------------------------------*
 | update fields after convergence as phi_i+1=phi_i+iterinc   kremheller 07/17 |
 *-----------------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraMonolithicTwoWay::UpdateFieldsAfterConvergence()
{
  // displacement, fluid variable and scatra variable incremental vector
  Teuchos::RCP<const Epetra_Vector> porostructinc;
  Teuchos::RCP<const Epetra_Vector> porofluidinc;
  Teuchos::RCP<const Epetra_Vector> scatrainc;
  ExtractFieldVectors(iterinc_, porostructinc, porofluidinc, scatrainc);

  // update ScaTra field
  ScatraAlgo()->ScaTraField()->UpdateIter(scatrainc);

  // update structure and fluid field
  PoroField()->UpdateFieldsAfterConvergence(porostructinc, porofluidinc);
}

/*----------------------------------------------------------------------*
 | setup RHS (like fsimon)                             kremheller 06/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraMonolithicTwoWay::SetupRHS()
{
  // create full monolithic rhs vector
  if(rhs_==Teuchos::null)
    rhs_ = Teuchos::rcp(new Epetra_Vector(*DofRowMap(), true));

  //note: rhs of fluid-structure system already setup in evaluate call

  // fill the Poroelasticity rhs vector rhs_ with the single field rhss
  SetupVector(*rhs_, PoroField()->RHS(), ScatraAlgo()->ScaTraField()->Residual());
}

/*----------------------------------------------------------------------*
 | setup vector of the poro and scatra field           kremheller 06/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraMonolithicTwoWay::SetupVector(Epetra_Vector &f,
                                        Teuchos::RCP<const Epetra_Vector> pv,
                                        Teuchos::RCP<const Epetra_Vector> sv)
{
  // extract dofs of the two fields
  // and put the poro/scatra field vector into the global vector f
  // noticing the block number

//  Teuchos::RCP<const Epetra_Vector> psx;
//  Teuchos::RCP<const Epetra_Vector> pfx;

  Extractor()->InsertVector(*(PoroField()->Extractor()->ExtractVector(pv,0)), 0, f);
  Extractor()->InsertVector(*(PoroField()->Extractor()->ExtractVector(pv,1)), 1, f);
  Extractor()->InsertVector(*sv, 2, f);
}

/*----------------------------------------------------------------------*
 | extract field vectors for calling Evaluate() of the  kremheller 03/17|
 | single fields                                                        |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraMonolithicTwoWay::ExtractFieldVectors(
    Teuchos::RCP<const Epetra_Vector> x,
    Teuchos::RCP<const Epetra_Vector>& stx,
    Teuchos::RCP<const Epetra_Vector>& flx,
    Teuchos::RCP<const Epetra_Vector>& scx)
{
  TEUCHOS_FUNC_TIME_MONITOR("POROMULTIPHASESCATRA::PoroMultiPhaseScaTraMonolithicTwoWay::ExtractFieldVectors");

  // process structure unknowns of the first field
  stx = Extractor()->ExtractVector(x, 0);

  // process fluid unknowns of the second field
  flx = Extractor()->ExtractVector(x, 1);

  // process scatra unknowns of the third field
  scx = Extractor()->ExtractVector(x, 2);
}

/*----------------------------------------------------------------------*
 | Solve linear Poromultiphase-elasticity system     kremheller 06/17   |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraMonolithicTwoWay::LinearSolve()
{
  // reset timer
  timernewton_.ResetStartTime();
  // *********** time measurement ***********
  double dtcpu = timernewton_.WallTime();
  // *********** time measurement ***********

  /*if (solveradapttol_ and (iter_ > 1))
  {
    double worst = normrhs_;
    double wanted = tolfres_;
    solver_->AdaptTolerance(wanted, worst, solveradaptolbetter_);
  }*/
  iterinc_->PutScalar(0.0);  // Useful? depends on solver and more

  // equilibrate global system of equations if necessary
  //EquilibrateSystem(systemmatrix_,rhs_);

  if(directsolve_)
  {
    // merge blockmatrix to SparseMatrix
    Teuchos::RCP<LINALG::SparseMatrix> sparse = systemmatrix_->Merge();

    // standard solver call
    // system is ready to solve since Dirichlet Boundary conditions have been applied in SetupSystemMatrix
    // or Evaluate
    solver_->Solve( sparse->EpetraOperator(),
                    iterinc_, rhs_,
                    true,
                    itnum_==1
                    );
  }
  else
    dserror("only direct solvers work so far");

  // *********** time measurement ***********
  dtsolve_ = timernewton_.WallTime() - dtcpu;
  // *********** time measurement ***********

  return;
}

/*----------------------------------------------------------------------*
 | simple convergence check                            kremheller 06/17 |
 *----------------------------------------------------------------------*/
bool POROMULTIPHASESCATRA::PoroMultiPhaseScaTraMonolithicTwoWay::Converged()
{

  return (norminc_ < ittolinc_ && normincfluid_ < ittolinc_ && normincstruct_ < ittolinc_ && normincscatra_ < ittolinc_ &&
      normrhs_ < ittolres_ && normrhsfluid_ < ittolres_ && normrhsstruct_ < ittolres_ && normrhsscatra_ < ittolres_);
}

/*----------------------------------------------------------------------*
 | Build necessary norms                               kremheller 06/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraMonolithicTwoWay::BuildConvergenceNorms()
{
  //------------------------------------------------------------ build residual force norms
  normrhs_ = UTILS::CalculateVectorNorm(vectornormfres_, rhs_);
  Teuchos::RCP<const Epetra_Vector> rhs_st;
  Teuchos::RCP<const Epetra_Vector> rhs_fl;
  Teuchos::RCP<const Epetra_Vector> rhs_sc;

  // get structure and fluid RHS
  ExtractFieldVectors(rhs_,rhs_st,rhs_fl,rhs_sc);

  // build also norms for structure, fluid and scatra
  normrhsstruct_ = UTILS::CalculateVectorNorm(vectornormfres_, rhs_st);
  normrhsfluid_ = UTILS::CalculateVectorNorm(vectornormfres_, rhs_fl);
  normrhsscatra_ = UTILS::CalculateVectorNorm(vectornormfres_, rhs_sc);

  //------------------------------------------------------------- build residual increment norms
  norminc_ = UTILS::CalculateVectorNorm(vectornorminc_, iterinc_);

  // displacement and fluid velocity & pressure incremental vector
  Teuchos::RCP<const Epetra_Vector> iterincst;
  Teuchos::RCP<const Epetra_Vector> iterincfl;
  Teuchos::RCP<const Epetra_Vector> iterincsc;

  // get structure and fluid increment
  ExtractFieldVectors(iterinc_,iterincst,iterincfl,iterincsc);

  // build also norms for fluid and structure
  normincstruct_ = UTILS::CalculateVectorNorm(vectornorminc_, iterincst);
  normincfluid_ = UTILS::CalculateVectorNorm(vectornorminc_, iterincfl);
  normincscatra_ = UTILS::CalculateVectorNorm(vectornorminc_, iterincsc);

  // build the total solution vector to build increment norm
  Teuchos::RCP<Epetra_Vector> sol_vec = Teuchos::rcp(new Epetra_Vector(*DofRowMap(), true));
  Extractor()->InsertVector(PoroField()->StructDispnp(), 0, sol_vec);
  Extractor()->InsertVector(PoroField()->FluidPhinp(), 1, sol_vec);
  Extractor()->InsertVector(ScatraAlgo()->ScaTraField()->Phinp(), 2, sol_vec);

  double dispnorm = UTILS::CalculateVectorNorm(vectornorminc_, PoroField()->StructureField()->Dispnp());
  double fluidnorm = UTILS::CalculateVectorNorm(vectornorminc_, PoroField()->FluidField()->Phinp());
  double scatranorm = UTILS::CalculateVectorNorm(vectornorminc_, ScatraAlgo()->ScaTraField()->Phinp());
  double totalnorm = UTILS::CalculateVectorNorm(vectornorminc_, sol_vec);

  // take care of very small norms
  if(dispnorm < 1.0e-6) dispnorm = 1.0;
  if(fluidnorm < 1.0e-6) fluidnorm = 1.0;
  if(scatranorm < 1.0e-6) scatranorm = 1.0;
  if(totalnorm < 1.0e-6) totalnorm = 1.0;

  // build relative increment norm
  normincstruct_/=dispnorm;
  normincfluid_/=fluidnorm;
  normincscatra_/=scatranorm;
  norminc_/=totalnorm;

  return;
}

/*----------------------------------------------------------------------*
 | Setup Newton-Raphson iteration                    kremheller 06/17   |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraMonolithicTwoWay::SetupNewton()
{

  // initialise equilibrium loop and norms
  itnum_ = 0;
  normrhs_ = 0.0;
  norminc_ = 0.0;
  normrhsfluid_ = 0.0;
  normincfluid_ = 0.0;
  normrhsstruct_ = 0.0;
  normincstruct_ = 0.0;
  normrhsscatra_ = 0.0;
  normincscatra_ = 0.0;
  tolinc_ = 0.0;
  tolfres_ = 0.0;
  tolinc_struct_ = 0.0;
  tolfres_struct_ = 0.0;
  tolinc_fluid_ = 0.0;
  tolfres_fluid_ = 0.0;
  tolinc_scatra_ = 0.0;
  tolfres_scatra_ = 0.0;

  // incremental solution vector with length of all dofs
  if(iterinc_==Teuchos::null) iterinc_ = LINALG::CreateVector(*DofRowMap(), true);
  else iterinc_->PutScalar(0.0);

  // a zero vector of full length
  if(zeros_==Teuchos::null) zeros_ = LINALG::CreateVector(*DofRowMap(), true);
  else zeros_->PutScalar(0.0);

  //AitkenReset();

  return;
}

/*----------------------------------------------------------------------*
 | Newton Output (adapted form tsi)                    kremheller 06/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraMonolithicTwoWay::NewtonOutput()
{

  // print the incremental based convergence check to the screen
  if (Comm().MyPID()==0 )
  {
    if(itnum_==1)
      printf("+--------------+---------------+---------------+---------------+---------------+---------------+\n");
    printf("|-  step/max  -|-  total-inc  -|-  fluid-inc  -|-  displ-inc  -|-  scatra-inc  -|-  norm-rhs  -| (ts =%10.3E,",
        dtsolve_);
    printf("\n");
    printf("|   %3d/%3d    |  %10.3E   |  %10.3E   |  %10.3E   |  %10.3E    |  %10.3E  |  te =%10.3E)",
         itnum_,itmax_,norminc_,normincfluid_,normincstruct_,normincscatra_,normrhs_,dtele_);
    printf("\n");
    printf("+--------------+---------------+---------------+---------------+----------------+--------------+\n");
  }

  return;
}

/*----------------------------------------------------------------------*
 | Error-Check and final output                        kremheller 06/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraMonolithicTwoWay::NewtonErrorCheck()
{

  // build the maximum value of the residuals and increments
  const double maxinc = std::max(norminc_,std::max(normincfluid_,std::max(normincscatra_,normincstruct_)));
  const double maxres = std::max(normrhs_,std::max(normrhsfluid_,std::max(normrhsscatra_,normrhsstruct_)));

  // print the incremental based convergence check to the screen
  if(Converged()) //norminc_ < ittolinc_ && normrhs_ < ittolinc_ && normincfluid_ < ittolinc_ && normincstruct_ < ittolinc_
  {
    if (Comm().MyPID()==0 )
    {
      printf("|  Monolithic iteration loop converged after iteration %3d/%3d !                               |\n", itnum_,itmax_);
      printf("|  Quantity           [norm]:                 TOL                                              |\n");
      printf("|  Max. rel. increment [%3s]:  %10.3E  < %10.3E                                        |\n",
          VectorNormString(vectornorminc_).c_str(),maxinc, ittolinc_);
      printf("|  Maximum    residual [%3s]:  %10.3E  < %10.3E                                        |\n",
          VectorNormString(vectornormfres_).c_str(),maxres, ittolres_);
      printf("+--------------+---------------+---------------+---------------+----------------+--------------+\n");
      printf("\n");
    }
  }
  else
  {
    if ((Comm().MyPID()==0) )
    {
      printf("|     >>>>>> not converged in %3d steps!                                                       |\n", itmax_);
      printf("|  Quantity           [norm]:                 TOL                                              |\n");
      printf("|  Max. rel. increment [%3s]:  %10.3E    %10.3E                                        |\n",
          VectorNormString(vectornorminc_).c_str(),maxinc, ittolinc_);
      printf("|  Maximum    residual [%3s]:  %10.3E    %10.3E                                        |\n",
          VectorNormString(vectornormfres_).c_str(),maxres, ittolres_);
      printf("+--------------+---------------+---------------+---------------+----------------+--------------+\n");
      printf("\n");
      printf("\n");
    }
    dserror("The monolithic solver did not converge in ITEMAX steps!");
  }


  return;
}

/*----------------------------------------------------------------------*
 | get the dofrowmap                                   kremheller 06/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> POROMULTIPHASESCATRA::PoroMultiPhaseScaTraMonolithicTwoWay::DofRowMap()
{
  return blockrowdofmap_->FullMap();
}

/*----------------------------------------------------------------------*
 | Print Header                                      kremheller 06/17   |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraMonolithicTwoWay::PrintHeader()
{

  if (Comm().MyPID()==0)
  {
    std::cout<<"+----------------------------------------------------------------------------------------------+" << std::endl;
    std::cout<<"| MONOLITHIC POROMULTIPHASE-SCATRA SOLVER                                                      |" << std::endl;
    std::cout<<"| STEP: " << std::setw(5) << std::setprecision(4) << std::scientific << Step() << "/"
        << std::setw(5) << std::setprecision(4) << std::scientific << NStep() << ", Time: "
        << std::setw(11) << std::setprecision(4) << std::scientific << Time() << "/"
        << std::setw(11) << std::setprecision(4) << std::scientific << MaxTime() << ", Dt: "
        << std::setw(11) << std::setprecision(4) << std::scientific << Dt() <<
        "                            |"<< std::endl;
  }
}

/*----------------------------------------------------------------------*
 |  check tangent stiffness matrix via finite differences      vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraMonolithicTwoWay::PoroMultiPhaseScaTraFDCheck()
{
  std::cout << "\n******************finite difference check***************" << std::endl;

  int dof_struct = (PoroField()->StructureField()->DofRowMap()->NumGlobalElements());
  int dof_fluid = (PoroField()->FluidField()->DofRowMap()->NumGlobalElements());
  int dof_scatra = (ScatraAlgo()->ScaTraField()->DofRowMap()->NumGlobalElements());

  std::cout << "structure field has " << dof_struct << " DOFs" << std::endl;
  std::cout << "fluid field has " << dof_fluid << " DOFs" << std::endl;
  std::cout << "scatra field has " << dof_scatra << " DOFs" << std::endl;

  Teuchos::RCP<Epetra_Vector> iterinc = Teuchos::null;
  iterinc = LINALG::CreateVector(*DofRowMap(), true);

  const int dofs = iterinc->GlobalLength();
  std::cout << "in total " << dofs << " DOFs" << std::endl;
  const double delta = 1e-8;

  iterinc->PutScalar(0.0);

  iterinc->ReplaceGlobalValue(0, 0, delta);

  Teuchos::RCP<Epetra_CrsMatrix> stiff_approx = Teuchos::null;
  stiff_approx = LINALG::CreateMatrix(*DofRowMap(), 81);

  Teuchos::RCP<Epetra_Vector> rhs_old = Teuchos::rcp(new Epetra_Vector(*DofRowMap(),
      true));
  rhs_old->Update(1.0, *rhs_, 0.0);
  Teuchos::RCP<Epetra_Vector> rhs_copy = Teuchos::rcp(new Epetra_Vector(*DofRowMap(),
      true));

  Teuchos::RCP<LINALG::SparseMatrix> sparse = systemmatrix_->Merge();
  Teuchos::RCP<LINALG::SparseMatrix> sparse_copy = Teuchos::rcp(
      new LINALG::SparseMatrix(sparse->EpetraMatrix(),LINALG::Copy));

  if (false)
  {
    std::cout << "iterinc_" << std::endl << *iterinc_ << std::endl;
    std::cout << "iterinc" << std::endl << *iterinc << std::endl;
    //std::cout << "meshdisp: " << std::endl << *(PoroField()->FluidField()-> ->Dispnp());
    std::cout << "disp: " << std::endl << *(PoroField()->StructureField()->Dispnp());
    //std::cout << "fluid vel" << std::endl << *(PoroField()->FluidField()->Velnp());
    //std::cout << "fluid acc" << std::endl << *(PoroField()->FluidField()->Accnp());
    //std::cout << "gridvel fluid" << std::endl << *(PoroField()->FluidField()->GridVel());
    std::cout << "gridvel struct" << std::endl << *(PoroField()->StructureField()->Velnp());
  }

  const int zeilennr = -1;
  const int spaltenr = -1;
  for (int i = 0; i < dofs; ++i)
  {
    if (CombinedDBCMap()->MyGID(i))
    {
      iterinc->ReplaceGlobalValue(i, 0, 0.0);
    }

    if (i == spaltenr)
      std::cout << "\n******************" << spaltenr + 1
          << ". Spalte!!***************" << std::endl;

    Evaluate(iterinc);
    SetupRHS();

    rhs_copy->Update(1.0, *rhs_, 0.0);

    iterinc_->PutScalar(0.0); // Useful? depends on solver and more
    LINALG::ApplyDirichlettoSystem(sparse_copy, iterinc_, rhs_copy,
        Teuchos::null, zeros_, *CombinedDBCMap());


    if (i == spaltenr)
    {

      std::cout << "rhs_: " << (*rhs_copy)[zeilennr] << std::endl;
      std::cout << "rhs_old: " << (*rhs_old)[zeilennr] << std::endl;
    }

    rhs_copy->Update(-1.0, *rhs_old, 1.0);
    rhs_copy->Scale(-1.0 / delta);

    int* index = &i;
    for (int j = 0; j < dofs; ++j)
    {
      double value = (*rhs_copy)[j];
      stiff_approx->InsertGlobalValues(j, 1, &value, index);

      if ((j == zeilennr) and (i == spaltenr))
      {
        std::cout << "\n******************" << zeilennr + 1
            << ". Zeile!!***************" << std::endl;
        std::cout << "iterinc_" << std::endl << *iterinc_ << std::endl;
        std::cout << "iterinc" << std::endl << *iterinc << std::endl;
        //std::cout << "meshdisp: " << std::endl << *(PoroField()->FluidField()->Dispnp());
        //std::cout << "meshdisp scatra: " << std::endl << *(ScaTraField()->Discretization()->GetState(ScaTraField()->NdsDisp(),"dispnp"));
        std::cout << "disp: " << std::endl << *(PoroField()->StructureField()->Dispnp());
        //std::cout << "fluid vel" << std::endl << *(PoroField()->FluidField()->Velnp());
        //std::cout << "scatra vel" << std::endl << *(ScaTraField()->Discretization()->GetState(ScaTraField()->NdsVel(),"velocity field"));
        //std::cout << "fluid acc" << std::endl << *(PoroField()->FluidField()->Accnp());
        //std::cout << "gridvel fluid" << std::endl << *(PoroField()->FluidField()->GridVel());
        std::cout << "gridvel struct" << std::endl << *(PoroField()->StructureField()->Velnp());

        std::cout << "stiff_apprx(" << zeilennr << "," << spaltenr << "): "
            << (*rhs_copy)[zeilennr] << std::endl;

        std::cout << "value(" << zeilennr << "," << spaltenr << "): " << value
            << std::endl;
        std::cout << "\n******************" << zeilennr + 1
            << ". Zeile Ende!!***************" << std::endl;
      }
    }

    if (not CombinedDBCMap()->MyGID(i))
      iterinc->ReplaceGlobalValue(i, 0, -delta);

    iterinc->ReplaceGlobalValue(i - 1, 0, 0.0);

    if (i != dofs - 1)
      iterinc->ReplaceGlobalValue(i + 1, 0, delta);

    if (i == spaltenr)
      std::cout << "\n******************" << spaltenr + 1
          << ". Spalte Ende!!***************" << std::endl;

  }

  Evaluate(iterinc);
  SetupRHS();

  stiff_approx->FillComplete();

  Teuchos::RCP<LINALG::SparseMatrix> stiff_approx_sparse = Teuchos::null;
  stiff_approx_sparse = Teuchos::rcp(new LINALG::SparseMatrix(stiff_approx,LINALG::Copy));

  stiff_approx_sparse->Add(*sparse_copy, false, -1.0, 1.0);

  Teuchos::RCP<Epetra_CrsMatrix> sparse_crs = sparse_copy->EpetraMatrix();

  Teuchos::RCP<Epetra_CrsMatrix> error_crs =
      stiff_approx_sparse->EpetraMatrix();

  error_crs->FillComplete();
  sparse_crs->FillComplete();

  bool success = true;
  double error_max_rel = 0.0;
  double error_max_abs = 0.0;
  for (int i = 0; i < dofs; ++i)
  {
    if (not CombinedDBCMap()->MyGID(i))
    {
      for (int j = 0; j < dofs; ++j)
      {
        if (not CombinedDBCMap()->MyGID(j))
        {

          double stiff_approx_ij=0.0;
          double sparse_ij=0.0;
          double error_ij=0.0;

          {
            // get error_crs entry ij
            int errornumentries;
            int errorlength = error_crs->NumGlobalEntries(i);
            std::vector<double> errorvalues(errorlength);
            std::vector<int> errorindices(errorlength);
            //int errorextractionstatus =
            error_crs->ExtractGlobalRowCopy(i,errorlength,errornumentries,&errorvalues[0],&errorindices[0]);
            for(int k = 0; k < errorlength; ++k)
            {
              if(errorindices[k] == j)
              {
                error_ij=errorvalues[k];
                break;
              }
              else
                error_ij = 0.0;
            }
          }

          // get sparse_ij entry ij
          {
            int sparsenumentries;
            int sparselength = sparse_crs->NumGlobalEntries(i);
            std::vector<double> sparsevalues(sparselength);
            std::vector<int> sparseindices(sparselength);
           // int sparseextractionstatus =
                sparse_crs->ExtractGlobalRowCopy(i,sparselength,sparsenumentries,&sparsevalues[0],&sparseindices[0]);
            for(int k = 0; k < sparselength; ++k)
            {
              if(sparseindices[k] == j)
              {
                sparse_ij=sparsevalues[k];
                break;
              }
              else
                sparse_ij=0.0;
            }
          }

          // get stiff_approx entry ij
          {
            int approxnumentries;
            int approxlength = stiff_approx->NumGlobalEntries(i);
            std::vector<double> approxvalues(approxlength);
            std::vector<int> approxindices(approxlength);
           // int approxextractionstatus =
                stiff_approx->ExtractGlobalRowCopy(i,approxlength,approxnumentries,&approxvalues[0],&approxindices[0]);
            for(int k = 0; k < approxlength; ++k)
            {
              if(approxindices[k] == j)
              {
                stiff_approx_ij=approxvalues[k];
                break;
              }
              else
                stiff_approx_ij=0.0;
            }
          }

          double error = 0.0;
          if (abs(stiff_approx_ij) > 1e-5)
            error = error_ij / (stiff_approx_ij);
          else if (abs(sparse_ij) > 1e-5)
            error = error_ij / (sparse_ij);

          if (abs(error) > abs(error_max_rel))
            error_max_rel = abs(error);
          if (abs(error_ij) > abs(error_max_abs))
            error_max_abs = abs(error_ij);

          if ((abs(error) > 1e-4))
          {
            if ((abs(error_ij) > 1e-5))
            //  if( (sparse_ij>1e-1) or (stiff_approx_ij>1e-1) )
            {
              std::cout << "finite difference check failed entry (" << i << "," << j
                  << ")! stiff: " << sparse_ij << ", approx: "
                  << stiff_approx_ij << " ,abs. error: " << error_ij
                  << " , rel. error: " << error << std::endl;

              success = false;
            }
          }
        }
      }
    }
  }

  if(success)
  {
    std::cout << "finite difference check successful, max. rel. error: "
        << error_max_rel << " , max. abs. error: "<< error_max_abs<<std::endl;
    std::cout << "******************finite difference check done***************\n\n"
        << std::endl;
  }
  else
    dserror("PoroFDCheck failed");

  return;
}