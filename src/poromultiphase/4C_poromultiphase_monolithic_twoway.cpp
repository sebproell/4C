// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_poromultiphase_monolithic_twoway.hpp"

#include "4C_adapter_art_net.hpp"
#include "4C_adapter_porofluidmultiphase_wrapper.hpp"
#include "4C_adapter_str_wrapper.hpp"
#include "4C_fem_condition_locsys.hpp"
#include "4C_fem_general_assemblestrategy.hpp"
#include "4C_fem_general_elements_paramsminimal.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_solver.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_equilibrate.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linalg_utils_sparse_algebra_print.hpp"
#include "4C_linear_solver_method.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_linear_solver_method_parameters.hpp"
#include "4C_poromultiphase_utils.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | constructor                                          kremheller 03/17 |
 *----------------------------------------------------------------------*/
POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::PoroMultiPhaseMonolithicTwoWay(
    MPI_Comm comm, const Teuchos::ParameterList& globaltimeparams)
    : PoroMultiPhaseMonolithic(comm, globaltimeparams),
      ittolinc_(0.0),
      ittolres_(0.0),
      itmax_(0),
      itmin_(1),
      itnum_(0),
      solveradaptolbetter_(0.0),
      solveradapttol_(false),
      blockrowdofmap_(nullptr),
      equilibration_(nullptr),
      equilibration_method_(Core::LinAlg::EquilibrationMethod::none),
      tolinc_(0.0),
      tolfres_(0.0),
      tolinc_struct_(0.0),
      tolfres_struct_(0.0),
      tolinc_fluid_(0.0),
      tolfres_fluid_(0.0),
      normrhs_(0.0),
      normrhsfluid_(0.0),
      normincfluid_(0.0),
      normrhsstruct_(0.0),
      normincstruct_(0.0),
      normrhsart_(0.0),
      normincart_(0.0),
      arterypressnorm_(0.0),
      maxinc_(0.0),
      maxres_(0.0),
      vectornormfres_(Inpar::POROMULTIPHASE::norm_undefined),
      vectornorminc_(Inpar::POROMULTIPHASE::norm_undefined),
      timernewton_("", true),
      dtsolve_(0.0),
      dtele_(0.0),
      fdcheck_(Inpar::POROMULTIPHASE::FdCheck::fdcheck_none)
{
}

/*----------------------------------------------------------------------*
 | initialization                                       kremheller 03/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::init(
    const Teuchos::ParameterList& globaltimeparams, const Teuchos::ParameterList& algoparams,
    const Teuchos::ParameterList& structparams, const Teuchos::ParameterList& fluidparams,
    const std::string& struct_disname, const std::string& fluid_disname, bool isale, int nds_disp,
    int nds_vel, int nds_solidpressure, int ndsporofluid_scatra,
    const std::map<int, std::set<int>>* nearbyelepairs)
{
  // call base class
  POROMULTIPHASE::PoroMultiPhaseMonolithic::init(globaltimeparams, algoparams, structparams,
      fluidparams, struct_disname, fluid_disname, isale, nds_disp, nds_vel, nds_solidpressure,
      ndsporofluid_scatra, nearbyelepairs);

  // inform user that structure field will not be solved but displacements will just be set to zero
  if (not solve_structure_) print_structure_disabled_info();

  // Get the parameters for the convergence_check
  itmax_ = algoparams.get<int>("ITEMAX");
  ittolres_ = algoparams.sublist("MONOLITHIC").get<double>("TOLRES_GLOBAL");
  ittolinc_ = algoparams.sublist("MONOLITHIC").get<double>("TOLINC_GLOBAL");

  blockrowdofmap_ = std::make_shared<Core::LinAlg::MultiMapExtractor>();

  fdcheck_ = Teuchos::getIntegralValue<Inpar::POROMULTIPHASE::FdCheck>(
      algoparams.sublist("MONOLITHIC"), "FDCHECK");

  equilibration_method_ = Teuchos::getIntegralValue<Core::LinAlg::EquilibrationMethod>(
      algoparams.sublist("MONOLITHIC"), "EQUILIBRATION");

  solveradaptolbetter_ = algoparams.sublist("MONOLITHIC").get<double>("ADAPTCONV_BETTER");
  solveradapttol_ = algoparams.sublist("MONOLITHIC").get<bool>("ADAPTCONV");
}

/*----------------------------------------------------------------------*
 | setup the system if necessary (called in poromultiphase_dyn.cpp)     |
 |                                                     kremheller 03/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::setup_system()
{
  // -------------------------------------------------------------create combined map
  setup_maps();

  // check global map extractor
  blockrowdofmap_->check_for_valid_map_extractor();

  //-----------------------------------build map of global dofs with DBC
  build_combined_dbc_map();
  // -------------------------------------------------------------

  // initialize Poromultiphase-elasticity-systemmatrix_
  systemmatrix_ =
      std::make_shared<Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>>(
          *extractor(), *extractor(), 81, false, true);

  // Initialize rhs
  rhs_ = std::make_shared<Core::LinAlg::Vector<double>>(*dof_row_map(), true);

  k_sf_ = std::make_shared<Core::LinAlg::SparseMatrix>(*(struct_dof_row_map()), 81, true, true);
  k_fs_ = std::make_shared<Core::LinAlg::SparseMatrix>(*(fluid_dof_row_map()), 81, true, true);

  // instantiate appropriate equilibration class
  auto equilibration_method =
      std::vector<Core::LinAlg::EquilibrationMethod>(1, equilibration_method_);
  equilibration_ = Core::LinAlg::build_equilibration(
      Core::LinAlg::MatrixType::block_field, equilibration_method, fullmap_);

  // structure_field: check whether we have locsys BCs, i.e. inclined structural
  //  Dirichlet BC

  std::vector<Core::Conditions::Condition*> locsysconditions(0);
  (structure_field()->discretization())->get_condition("Locsys", locsysconditions);

  // if there are inclined structural Dirichlet BC, get the structural LocSysManager
  if (locsysconditions.size()) locsysman_ = structure_field()->locsys_manager();

  return;
}

/*----------------------------------------------------------------------*
 | setup the map                                       kremheller 04/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::setup_maps()
{
  std::vector<std::shared_ptr<const Epetra_Map>> vecSpaces;

  vecSpaces.push_back(struct_dof_row_map());

  vecSpaces.push_back(fluid_dof_row_map());

  if (vecSpaces[0]->NumGlobalElements() == 0) FOUR_C_THROW("No structure equation. Panic.");
  if (vecSpaces[1]->NumGlobalElements() == 0) FOUR_C_THROW("No fluid equation. Panic.");

  // full Poromultiphase-elasticity-map
  fullmap_ = Core::LinAlg::MultiMapExtractor::merge_maps(vecSpaces);

  // full Poromultiphase-elasticity-blockmap
  blockrowdofmap_->setup(*fullmap_, vecSpaces);

  return;
}

/*----------------------------------------------------------------------*
 | Monolithic Time Step                                kremheller 03/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::time_step()
{
  // Prepare stuff
  setup_newton();
  print_header();

  // Evaluate
  evaluate(iterinc_);

  // Newton-Loop
  while ((not converged() and itnum_ < itmax_) or (itnum_ < itmin_))
  {
    // increment number of iteration
    itnum_++;

    // Solve
    linear_solve();
    solver_->reset_tolerance();

    // Build Convergence Norms
    build_convergence_norms();

    if (not converged())
    {
      // Evaluate
      evaluate(iterinc_);

      // perform FD Check of monolithic system matrix
      if (fdcheck_ == Inpar::POROMULTIPHASE::fdcheck_global) poro_fd_check();
    }
    else
    {
      // convergence check is based on residual(phi_i) < tol and phi_i+1 - phi_i < tol
      // in this function we update phi_i+1 as phi_i+1 = phi_i + iterinc for all fields
      // even though we have not evaluated the residual of phi_i+1 it will still be more exact than
      // the one at phi_i
      update_fields_after_convergence();
    }

    // print output
    newton_output();
  }

  // Error-Check
  newton_error_check();

  return;
}

/*-----------------------------------------------------------------------/
|  build the combined dbcmap                           kremheller 03/17  |
/-----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::build_combined_dbc_map()
{
  // get structure and fluid dbc maps
  const std::shared_ptr<const Epetra_Map> scondmap =
      structure_field()->get_dbc_map_extractor()->cond_map();
  const std::shared_ptr<const Epetra_Map> fcondmap =
      fluid_field()->get_dbc_map_extractor()->cond_map();
  // merge them
  combinedDBCMap_ = Core::LinAlg::merge_map(scondmap, fcondmap, false);

  return;
}

/*----------------------------------------------------------------------*
 | Evaluate (build global Matrix and RHS)            kremheller 03/17   |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::evaluate(
    std::shared_ptr<const Core::LinAlg::Vector<double>> iterinc)
{
  TEUCHOS_FUNC_TIME_MONITOR("POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::Evaluate");

  // reset timer
  timernewton_.reset();
  // *********** time measurement ***********
  double dtcpu = timernewton_.wallTime();
  // *********** time measurement ***********

  // displacement and fluid velocity & pressure incremental vector
  std::shared_ptr<const Core::LinAlg::Vector<double>> s_iterinc;
  std::shared_ptr<const Core::LinAlg::Vector<double>> f_iterinc;
  extract_field_vectors(iterinc, s_iterinc, f_iterinc);

  evaluate(s_iterinc, f_iterinc, itnum_ == 0);

  // *********** time measurement ***********
  dtele_ = timernewton_.wallTime() - dtcpu;
  // *********** time measurement ***********

  return;
}
/*----------------------------------------------------------------------*
 | Evaluate (build global Matrix and RHS, public --> allows access      |
 | from outside --> monolithic scatra-coupling)        kremheller 06/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::evaluate(
    std::shared_ptr<const Core::LinAlg::Vector<double>> sx,
    std::shared_ptr<const Core::LinAlg::Vector<double>> fx, const bool firstcall)
{
  // (1) Update fluid Field and reconstruct pressures and saturations
  fluid_field()->update_iter(fx);

  if (solve_structure_)
  {
    // (2) set fluid solution in structure field
    structure_field()->discretization()->set_state(1, "porofluid", fluid_field()->phinp());

    // (3) evaluate structure
    if (firstcall)  // first call (iterinc_ = 0) --> sx = 0
      structure_field()->evaluate();
    else  //(this call will also update displacements and velocities)
      structure_field()->evaluate(sx);

    // (4) Set structure solution on fluid field
    set_struct_solution(structure_field()->dispnp(), structure_field()->velnp());
  }
  else
  {
    // (4) Set structure solution on fluid field
    set_struct_solution(struct_zeros_, struct_zeros_);
    structure_field()->system_matrix()->zero();
    structure_field()->system_matrix()->complete(structure_field()->system_matrix()->range_map(),
        structure_field()->system_matrix()->range_map());
  }

  // (5) Evaluate the fluid
  fluid_field()->evaluate();

  // (6) Build the monolithic system matrix
  setup_system_matrix();

  // check whether we have a sanely filled tangent matrix
  if (not systemmatrix_->filled())
  {
    FOUR_C_THROW("Effective tangent matrix must be filled here");
  }

  // (7) Build the monolithic system vector
  setup_rhs();
}

/*----------------------------------------------------------------------*
 | setup system matrix of poromultiphase-elasticity   kremheller 03/17  |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::setup_system_matrix(
    Core::LinAlg::BlockSparseMatrixBase& mat)
{
  // TEUCHOS_FUNC_TIME_MONITOR("PoroElast::Monolithic::setup_system_matrix");

  // pure structural part k_ss ((ndim*n_nodes)x(ndim*n_nodes))

  // build pure structural block k_ss
  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here. Extract Jacobian matrices and put them into composite system
  // matrix W
  std::shared_ptr<Core::LinAlg::SparseMatrix> k_ss = structure_field()->system_matrix();

  if (k_ss == nullptr) FOUR_C_THROW("structure system matrix null pointer!");

  // Copy from TSI
  if (locsysman_ != nullptr)
  {
    // rotate k_ss to local coordinate system --> k_ss^{~}
    locsysman_->rotate_global_to_local(k_ss);
    // apply apply_dirichlet_with_trafo() on rotated block k_ss^{~}
    // --> if dof has an inclined DBC: blank the complete row, the '1.0' is set
    //     on diagonal of row, i.e. on diagonal of k_ss
    k_ss->apply_dirichlet_with_trafo(
        *locsysman_->trafo(), *structure_field()->get_dbc_map_extractor()->cond_map(), true);
  }  // end locsys
  // default: (locsysman_ == nullptr), i.e. NO inclined Dirichlet BC
  else
    k_ss->apply_dirichlet(*structure_field()->get_dbc_map_extractor()->cond_map(), true);

  /*----------------------------------------------------------------------*/
  // structural part k_sf ((ndim*n_nodes)x(n_phases*n_nodes))
  // build mechanical-fluid block

  // create empty matrix
  std::shared_ptr<Core::LinAlg::SparseMatrix> k_sf = struct_fluid_coupling_matrix();

  // call the element and calculate the matrix block
  apply_str_coupl_matrix(k_sf);

  // Copy from TSI
  // apply dirichlet boundary conditions properly on matrix k_sf, i.e. blank row
  // if dof is a structural DBC
  // Normally, DBC should be applied on complete systemmatrix mat, but for
  // diagonal blocks (here k_ss, k_tt) DBC are ALREADY applied in
  // prepare_system_for_newton_solve() included in evaluate(sx)
  //
  // to avoid double work, we only call ApplyDirichlet for the off-diagonal blocks,
  // here k_sf
  // k_sf is an off-diagonal block --> pass the bool diagonal==false
  // ApplyDirichlet*() expect filled matrix
  //
  // in case of inclined STR-DBC
  //   1.) transform the off-diagonal block k_sf to the local system --> k_st^{~}
  //   2.) apply apply_dirichlet_with_trafo() on rotated block k_sf^{~}
  //              --> blank the row, which has a DBC

  // to apply Multiply in LocSys, k_st has to be FillCompleted
  k_sf->complete(
      fluid_field()->system_matrix()->range_map(), structure_field()->system_matrix()->range_map());

  if (locsysman_ != nullptr)
  {
    // rotate k_st to local coordinate system --> k_st^{~}
    locsysman_->rotate_global_to_local(k_sf);
    // apply apply_dirichlet_with_trafo() on rotated block k_st^{~}
    // --> if dof has an inclined DBC: blank the complete row, the '1.0' is set
    //     on diagonal of row, i.e. on diagonal of k_ss
    k_sf->apply_dirichlet_with_trafo(
        *locsysman_->trafo(), *structure_field()->get_dbc_map_extractor()->cond_map(), false);
  }  // end locsys
  // default: (locsysman_ == nullptr), i.e. NO inclined Dirichlet BC
  else
    k_sf->apply_dirichlet(*structure_field()->get_dbc_map_extractor()->cond_map(), false);

  /*----------------------------------------------------------------------*/
  // pure fluid part k_ff ( (n_phases*n_nodes)x(n_phases*n_nodes) )

  // build pure fluid block k_ff
  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here. Extract Jacobian matrices and put them into composite system
  // matrix W
  // NOTE: DBC's have already been applied within Evaluate (prepare_system_for_newton_solve())
  std::shared_ptr<Core::LinAlg::SparseMatrix> k_ff = fluid_field()->system_matrix();

  if (k_ff == nullptr) FOUR_C_THROW("fluid system matrix null pointer!");

  /*----------------------------------------------------------------------*/
  // fluid part k_fs ( (n_phases*n_nodes)x(ndim*n_nodes) )
  // build fluid-mechanical block

  // create empty matrix
  std::shared_ptr<Core::LinAlg::SparseMatrix> k_fs = fluid_struct_coupling_matrix();

  // call the element and calculate the matrix block
  apply_fluid_coupl_matrix(k_fs);

  // apply DBC's also on off-diagonal fluid-structure coupling block
  k_fs->apply_dirichlet(*fluid_field()->get_dbc_map_extractor()->cond_map(), false);

  // uncomplete matrix block (appears to be required in certain cases (locsys+iterative solver))
  if (solve_structure_)
  {
    k_ss->un_complete();
    k_sf->un_complete();
  }
  k_fs->un_complete();
  k_ff->un_complete();

  // assign structure part to the Poroelasticity matrix
  mat.assign(0, 0, Core::LinAlg::View, *k_ss);
  // assign coupling part to the Poroelasticity matrix
  mat.assign(0, 1, Core::LinAlg::View, *k_sf);
  // assign fluid part to the poroelasticity matrix
  mat.assign(1, 1, Core::LinAlg::View, *k_ff);
  // assign coupling part to the Poroelasticity matrix
  mat.assign(1, 0, Core::LinAlg::View, *k_fs);

  /*----------------------------------------------------------------------*/
  // done. make sure all blocks are filled.
  mat.complete();

}  // setup_system_matrix

/*----------------------------------------------------------------------*
 | get fluid structure-coupling sparse matrix           kremheller 03/17 |
 *----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::SparseMatrix>
POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::fluid_struct_coupling_matrix()
{
  std::shared_ptr<Core::LinAlg::SparseMatrix> sparse =
      std::dynamic_pointer_cast<Core::LinAlg::SparseMatrix>(k_fs_);
  if (sparse == nullptr) FOUR_C_THROW("cast to Core::LinAlg::SparseMatrix failed!");

  return sparse;
}  // fluid_struct_coupling_matrix()

/*----------------------------------------------------------------------*
 | get structure fluid-coupling sparse matrix           kremheller 03/17 |
 *----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::SparseMatrix>
POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::struct_fluid_coupling_matrix()
{
  std::shared_ptr<Core::LinAlg::SparseMatrix> sparse =
      std::dynamic_pointer_cast<Core::LinAlg::SparseMatrix>(k_sf_);
  if (sparse == nullptr) FOUR_C_THROW("cast to Core::LinAlg::SparseMatrix failed!");

  return sparse;
}  // fluid_struct_coupling_matrix()

/*----------------------------------------------------------------------*
 | evaluate fluid-structural system matrix at state    kremheller 03/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::apply_fluid_coupl_matrix(
    std::shared_ptr<Core::LinAlg::SparseOperator> k_fs  //!< off-diagonal tangent matrix term
)
{
  // reset
  k_fs->zero();
  if (solve_structure_) fluid_field()->assemble_fluid_struct_coupling_mat(k_fs);
  k_fs->complete(
      structure_field()->system_matrix()->range_map(), fluid_field()->system_matrix()->range_map());

  return;
}

/*-----------------------------------------------------------------------------*
 | update fields after convergence as phi_i+1=phi_i+iterinc   kremheller 07/17 |
 *-----------------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::update_fields_after_convergence()
{
  // displacement and fluid velocity & pressure incremental vector
  std::shared_ptr<const Core::LinAlg::Vector<double>> sx;
  std::shared_ptr<const Core::LinAlg::Vector<double>> fx;
  extract_field_vectors(iterinc_, sx, fx);

  update_fields_after_convergence(sx, fx);

  return;
}

/*-----------------------------------------------------------------------------*
 | update fields after convergence as phi_i+1=phi_i+iterinc   kremheller 07/17 |
 *-----------------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::update_fields_after_convergence(
    std::shared_ptr<const Core::LinAlg::Vector<double>>& sx,
    std::shared_ptr<const Core::LinAlg::Vector<double>>& fx)
{
  // (1) Update fluid Field and reconstruct pressures and saturations
  fluid_field()->update_iter(fx);
  fluid_field()->reconstruct_pressures_and_saturations();
  fluid_field()->reconstruct_flux();

  if (solve_structure_) structure_field()->evaluate(sx);

  // (4) Set structure solution on fluid field
  set_struct_solution(structure_field()->dispnp(), structure_field()->velnp());

  return;
}

/*----------------------------------------------------------------------*
 | evaluate structural-fluid system matrix at state    kremheller 03/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::apply_str_coupl_matrix(
    std::shared_ptr<Core::LinAlg::SparseOperator> k_sf  //!< off-diagonal tangent matrix term
)
{
  k_sf->zero();

  if (solve_structure_)
  {
    // create the parameters for the discretization
    Teuchos::ParameterList sparams;

    //! pointer to the model evaluator data container
    std::shared_ptr<Core::Elements::ParamsMinimal> params =
        std::make_shared<Core::Elements::ParamsMinimal>();

    // set parameters needed for element evaluation
    params->set_action_type(Core::Elements::struct_poro_calc_fluidcoupling);
    params->set_total_time(time());
    params->set_delta_time(dt());
    // std::cout << Dt() << std::endl;

    sparams.set<std::shared_ptr<Core::Elements::ParamsInterface>>("interface", params);
    sparams.set<std::string>("action", "struct_poro_calc_fluidcoupling");
    sparams.set<double>("delta time", dt());
    sparams.set<double>("total time", time());

    structure_field()->discretization()->clear_state();
    structure_field()->discretization()->set_state(0, "displacement", structure_field()->dispnp());
    structure_field()->discretization()->set_state(0, "velocity", structure_field()->velnp());
    structure_field()->discretization()->set_state(1, "porofluid", fluid_field()->phinp());

    // build specific assemble strategy for mechanical-fluid system matrix
    // from the point of view of structure_field:
    // structdofset = 0, fluiddofset = 1
    Core::FE::AssembleStrategy structuralstrategy(0,  // structdofset for row
        1,                                            // fluiddofset for column
        k_sf,                                         // mechanical-fluid coupling matrix
        nullptr, nullptr, nullptr, nullptr);

    // evaluate the mechanical-fluid system matrix on the structural element
    structure_field()->discretization()->evaluate(sparams, structuralstrategy);

    structure_field()->discretization()->clear_state();

    // scale with time integration factor
    k_sf->scale(1.0 - structure_field()->tim_int_param());
  }

  return;
}
/*----------------------------------------------------------------------*
 | setup solver for monolithic problem                kremheller 03/17  |
 *----------------------------------------------------------------------*/
bool POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::setup_solver()
{
  //  solver
  // create a linear solver
  // get dynamic section of poroelasticity
  const Teuchos::ParameterList& poromultdyn =
      Global::Problem::instance()->poro_multi_phase_dynamic_params();
  // get the solver number used for linear poroelasticity solver
  const int linsolvernumber = poromultdyn.sublist("MONOLITHIC").get<int>("LINEAR_SOLVER");
  // check if the poroelasticity solver has a valid solver number
  if (linsolvernumber == (-1))
    FOUR_C_THROW(
        "no linear solver defined for poromultiphaseflow. Please set LINEAR_SOLVER in "
        "POROMULTIPHASE DYNAMIC to a valid number!");
  const Teuchos::ParameterList& solverparams =
      Global::Problem::instance()->solver_params(linsolvernumber);
  const auto solvertype =
      Teuchos::getIntegralValue<Core::LinearSolver::SolverType>(solverparams, "SOLVER");

  create_linear_solver(solverparams, solvertype);

  vectornormfres_ = Teuchos::getIntegralValue<Inpar::POROMULTIPHASE::VectorNorm>(
      poromultdyn.sublist("MONOLITHIC"), "VECTORNORM_RESF");
  vectornorminc_ = Teuchos::getIntegralValue<Inpar::POROMULTIPHASE::VectorNorm>(
      poromultdyn.sublist("MONOLITHIC"), "VECTORNORM_INC");

  return true;
}

/*----------------------------------------------------------------------*
 | Create linear (iterative) solver                  kremheller 08/17   |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::create_linear_solver(
    const Teuchos::ParameterList& solverparams, const Core::LinearSolver::SolverType solvertype)
{
  solver_ = std::make_shared<Core::LinAlg::Solver>(solverparams, get_comm(),
      Global::Problem::instance()->solver_params_callback(),
      Teuchos::getIntegralValue<Core::IO::Verbositylevel>(
          Global::Problem::instance()->io_params(), "VERBOSITY"));
  // no need to do the rest for direct solvers
  if (solvertype == Core::LinearSolver::SolverType::umfpack or
      solvertype == Core::LinearSolver::SolverType::superlu)
    return;

  if (solvertype != Core::LinearSolver::SolverType::belos)
  {
    std::cout << "!!!!!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!" << std::endl;
    std::cout << " Note: the BGS2x2 preconditioner now " << std::endl;
    std::cout << " uses the structural solver and fluid solver blocks" << std::endl;
    std::cout << " for building the internal inverses" << std::endl;
    std::cout << " Remove the old BGS PRECONDITIONER BLOCK entries " << std::endl;
    std::cout << " in the input files!" << std::endl;
    std::cout << "!!!!!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!" << std::endl;
    FOUR_C_THROW("Iterative solver expected");
  }
  const auto azprectype =
      Teuchos::getIntegralValue<Core::LinearSolver::PreconditionerType>(solverparams, "AZPREC");

  // plausibility check
  switch (azprectype)
  {
    case Core::LinearSolver::PreconditionerType::block_teko:
    {
      // no plausibility checks here
      // if you forget to declare an xml file you will get an error message anyway
    }
    break;
    default:
      FOUR_C_THROW("Block preconditioner expected");
      break;
  }

  // build the null spaces of the single blocks
  build_block_null_spaces(solver_);
}

/*-----------------------------------------------------------------------------------*
 | build null spaces associated with blocks of global system matrix kremheller 08/17 |
 *-----------------------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::build_block_null_spaces(
    std::shared_ptr<Core::LinAlg::Solver>& solver)
{
  Teuchos::ParameterList& structure_params = solver->params().sublist("Inverse1");
  Core::LinearSolver::Parameters::compute_solver_parameters(
      *structure_field()->discretization(), structure_params);

  Teuchos::ParameterList& fluid_params = solver->params().sublist("Inverse2");
  Core::LinearSolver::Parameters::compute_solver_parameters(
      *fluid_field()->discretization(), fluid_params);
}

/*----------------------------------------------------------------------*
 | Setup Newton-Raphson iteration                    kremheller 03/17   |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::setup_newton()
{
  // initialise equilibrium loop and norms
  itnum_ = 0;
  normrhs_ = 0.0;
  normrhsfluid_ = 0.0;
  normincfluid_ = 0.0;
  normrhsstruct_ = 0.0;
  normincstruct_ = 0.0;
  tolinc_ = 0.0;
  tolfres_ = 0.0;
  tolinc_struct_ = 0.0;
  tolfres_struct_ = 0.0;
  tolinc_fluid_ = 0.0;
  tolfres_fluid_ = 0.0;
  normrhsart_ = 0.0;
  normincart_ = 0.0;
  arterypressnorm_ = 0.0;
  maxinc_ = 0.0;
  maxres_ = 0.0;

  // incremental solution vector with length of all dofs
  if (iterinc_ == nullptr)
    iterinc_ = Core::LinAlg::create_vector(*dof_row_map(), true);
  else
    iterinc_->put_scalar(0.0);

  // a zero vector of full length
  if (zeros_ == nullptr)
    zeros_ = Core::LinAlg::create_vector(*dof_row_map(), true);
  else
    zeros_->put_scalar(0.0);

  // AitkenReset();

  return;
}

/*----------------------------------------------------------------------*
 | Print Header                                      kremheller 03/17   |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::print_header()
{
  if (Core::Communication::my_mpi_rank(get_comm()) == 0)
  {
    if (!solve_structure_) print_structure_disabled_info();
    std::cout
        << "+----------------------------------------------------------------------------------"
           "----------+"
        << std::endl;
    std::cout << "| MONOLITHIC POROMULTIPHASE SOLVER                                               "
                 "            |"
              << std::endl;
    std::cout << "| STEP: " << std::setw(5) << std::setprecision(4) << std::scientific << step()
              << "/" << std::setw(5) << std::setprecision(4) << std::scientific << n_step()
              << ", Time: " << std::setw(11) << std::setprecision(4) << std::scientific << time()
              << "/" << std::setw(11) << std::setprecision(4) << std::scientific << max_time()
              << ", Dt: " << std::setw(11) << std::setprecision(4) << std::scientific << dt()
              << "                          |" << std::endl;
  }
}

/*----------------------------------------------------------------------*
 | Build necessary norms                               kremheller 03/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::build_convergence_norms()
{
  //------------------------------------------------------------ build residual force norms
  normrhs_ = Utils::calculate_vector_norm(vectornormfres_, *rhs_);
  std::shared_ptr<const Core::LinAlg::Vector<double>> rhs_s;
  std::shared_ptr<const Core::LinAlg::Vector<double>> rhs_f;

  // get structure and fluid RHS
  extract_structure_and_fluid_vectors(rhs_, rhs_s, rhs_f);

  // build also norms for fluid and structure
  normrhsstruct_ = Utils::calculate_vector_norm(vectornormfres_, *rhs_s);
  normrhsfluid_ = Utils::calculate_vector_norm(vectornormfres_, *rhs_f);

  //------------------------------------------------------------- build residual increment norms
  // displacement and fluid velocity & pressure incremental vector
  std::shared_ptr<const Core::LinAlg::Vector<double>> iterincs;
  std::shared_ptr<const Core::LinAlg::Vector<double>> iterincf;

  // get structure and fluid increment
  extract_structure_and_fluid_vectors(iterinc_, iterincs, iterincf);

  // build also norms for fluid and structure
  normincstruct_ = Utils::calculate_vector_norm(vectornorminc_, *iterincs);
  normincfluid_ = Utils::calculate_vector_norm(vectornorminc_, *iterincf);

  double dispnorm = Utils::calculate_vector_norm(vectornorminc_, (*structure_field()->dispnp()));
  double fluidnorm = Utils::calculate_vector_norm(vectornorminc_, (*fluid_field()->phinp()));

  // take care of very small norms
  if (dispnorm < 1.0e-6) dispnorm = 1.0;
  if (fluidnorm < 1.0e-6) fluidnorm = 1.0;
  if (arterypressnorm_ < 1.0e-6) arterypressnorm_ = 1.0;

  // build relative increment norm
  normincstruct_ /= dispnorm;
  normincfluid_ /= fluidnorm;
  normincart_ /= arterypressnorm_;

  // build the maximum value of the residuals and increments
  maxinc_ = std::max({normincart_, normincfluid_, normincstruct_});
  maxres_ = std::max({normrhs_, normrhsart_, normrhsfluid_, normrhsstruct_});

  return;
}

/*----------------------------------------------------------------------*
 | Newton Output (adapted form tsi)                    kremheller 03/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::newton_output()
{
  // print the incremental based convergence check to the screen
  if (Core::Communication::my_mpi_rank(get_comm()) == 0)
  {
    if (itnum_ == 1)
      printf(
          "+--------------+--------------+--------------+--------------+--------------+"
          "-----------------+\n");
    printf(
        "|-  step/max  -|-   max-inc  -|- fluid-inc  -|-  disp-inc  -|-   1D-inc   -|- "
        "norm(tot-rhs) -| (ts =%10.3E,",
        dtsolve_);
    printf("\n");
    printf(
        "|   %3d/%3d    | %10.3E   | %10.3E   | %10.3E   | %10.3E   | %10.3E      "
        "|  te =%10.3E)",
        itnum_, itmax_, maxinc_, normincfluid_, normincstruct_, normincart_, normrhs_, dtele_);
    printf("\n");
    printf(
        "+--------------+--------------+--------------+--------------+--------------+"
        "-----------------+\n");
  }

  return;
}

/*----------------------------------------------------------------------*
 | Error-Check and final output                        kremheller 03/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::newton_error_check()
{
  // print the incremental based convergence check to the screen
  if (converged())  // norminc_ < ittol_ && normrhs_ < ittol_ && normincfluid_ < ittol_ &&
                    // normincstruct_ < ittol_
  {
    if (Core::Communication::my_mpi_rank(get_comm()) == 0)
    {
      printf(
          "|  Monolithic iteration loop converged after iteration %3d/%3d !                        "
          "     |\n",
          itnum_, itmax_);
      printf(
          "|  Quantity           [norm]:                 TOL                                       "
          "     |\n");
      printf(
          "|  Max. rel. increment [%3s]:  %10.3E  < %10.3E                                      "
          "|\n",
          vector_norm_string(vectornorminc_).c_str(), maxinc_, ittolinc_);
      printf(
          "|  Maximum    residual [%3s]:  %10.3E  < %10.3E                                      "
          "|\n",
          vector_norm_string(vectornormfres_).c_str(), maxres_, ittolres_);
      printf(
          "+--------------+--------------+--------------+--------------+--------------+"
          "-----------------+\n");
      printf("\n");
    }
  }
  else
  {
    if ((Core::Communication::my_mpi_rank(get_comm()) == 0))
    {
      printf(
          "|     >>>>>> not converged in %3d steps!                                                "
          "   |\n",
          itmax_);
      printf(
          "|  Max. rel. increment [%3s]:  %10.3E    %10.3E                                    |\n",
          vector_norm_string(vectornorminc_).c_str(), maxinc_, ittolinc_);
      printf(
          "|  Maximum    residual [%3s]:  %10.3E    %10.3E                                    |\n",
          vector_norm_string(vectornormfres_).c_str(), maxres_, ittolres_);
      printf(
          "+--------------+----------------+------------------+--------------------+---------------"
          "---+\n");
      printf("\n");
      printf("\n");
    }
    FOUR_C_THROW("The monolithic solver did not converge in ITEMAX steps!");
  }


  return;
}

/*----------------------------------------------------------------------*
 | simple convergence check                            kremheller 03/17 |
 *----------------------------------------------------------------------*/
bool POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::converged()
{
  return (normincfluid_ < ittolinc_ && normincstruct_ < ittolinc_ && normincart_ < ittolinc_ &&
          normrhs_ < ittolres_ && normrhsstruct_ < ittolres_ && normrhsfluid_ < ittolres_ &&
          normrhsart_ < ittolres_);
}

/*----------------------------------------------------------------------*
 | Solve linear Poromultiphase-elasticity system     kremheller 03/17   |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::linear_solve()
{
  // reset timer
  timernewton_.reset();
  // *********** time measurement ***********
  double dtcpu = timernewton_.wallTime();
  // *********** time measurement ***********

  Core::LinAlg::SolverParams solver_params;
  if (solveradapttol_ and (itnum_ > 1))
  {
    solver_params.nonlin_tolerance = tolfres_;
    solver_params.nonlin_residual = std::max(maxres_, maxinc_);
    solver_params.lin_tol_better = solveradaptolbetter_;
  }
  iterinc_->put_scalar(0.0);  // Useful? depends on solver and more

  // equilibrate global system of equations if necessary
  equilibration_->equilibrate_system(systemmatrix_, rhs_, blockrowdofmap_);

  // standard solver call
  // system is ready to solve since Dirichlet Boundary conditions have been applied in
  // setup_system_matrix or Evaluate

  solver_params.refactor = true;
  solver_params.reset = itnum_ == 1;
  solver_->solve(systemmatrix_->epetra_operator(), iterinc_, rhs_, solver_params);

  equilibration_->unequilibrate_increment(iterinc_);

  // *********** time measurement ***********
  dtsolve_ = timernewton_.wallTime() - dtcpu;
  // *********** time measurement ***********

  return;
}

/*----------------------------------------------------------------------*
 | get the dof row map                                 kremheller 03/17 |
 *----------------------------------------------------------------------*/
std::shared_ptr<const Epetra_Map> POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::dof_row_map()
{
  return blockrowdofmap_->full_map();
}

/*----------------------------------------------------------------------*
 | setup RHS (like fsimon)                             kremheller 03/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::setup_rhs()
{
  // get structure part
  std::shared_ptr<Core::LinAlg::Vector<double>> str_rhs = setup_structure_partof_rhs();

  // fill the Poroelasticity rhs vector rhs_ with the single field rhss
  setup_vector(*rhs_, str_rhs, fluid_field()->rhs());

}  // setup_rhs()

/*----------------------------------------------------------------------*
 | setup RHS (like fsimon)                             kremheller 05/18 |
 *----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::Vector<double>>
POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::setup_structure_partof_rhs()
{
  // Copy from TSI
  std::shared_ptr<Core::LinAlg::Vector<double>> str_rhs = struct_zeros_;
  if (solve_structure_)
    str_rhs = std::make_shared<Core::LinAlg::Vector<double>>(*structure_field()->rhs());
  if (locsysman_ != nullptr) locsysman_->rotate_global_to_local(*str_rhs);

  return str_rhs;
}

/*----------------------------------------------------------------------*
 | setup vector of the structure and fluid field         kremheller 03/17|
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::setup_vector(Core::LinAlg::Vector<double>& f,
    std::shared_ptr<const Core::LinAlg::Vector<double>> sv,
    std::shared_ptr<const Core::LinAlg::Vector<double>> fv)
{
  extractor()->insert_vector(*sv, 0, f);

  f.scale(-1);
  extractor()->insert_vector(*fv, 1, f);
}
/*----------------------------------------------------------------------*
 | extract field vectors for calling evaluate() of the  kremheller 03/17|
 | single fields                                                        |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::extract_field_vectors(
    std::shared_ptr<const Core::LinAlg::Vector<double>> x,
    std::shared_ptr<const Core::LinAlg::Vector<double>>& sx,
    std::shared_ptr<const Core::LinAlg::Vector<double>>& fx)
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::extract_field_vectors");

  // process structure unknowns of the first field
  sx = extractor()->extract_vector(*x, 0);

  // process fluid unknowns of the second field
  fx = extractor()->extract_vector(*x, 1);
}
/*----------------------------------------------------------------------*
 | extract 3D field vectors (structure and fluid)    kremheller 10/20   |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::extract_structure_and_fluid_vectors(
    std::shared_ptr<const Core::LinAlg::Vector<double>> x,
    std::shared_ptr<const Core::LinAlg::Vector<double>>& sx,
    std::shared_ptr<const Core::LinAlg::Vector<double>>& fx)
{
  PoroMultiPhaseMonolithicTwoWay::extract_field_vectors(x, sx, fx);
}

/*----------------------------------------------------------------------*
 | inform user that structure is not solved            kremheller 08/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::print_structure_disabled_info()
{
  // print out Info
  if (Core::Communication::my_mpi_rank(get_comm()) == 0)
  {
    std::cout << "\n";
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                 "++++++++++++++++++++++++++++++++\n";
    std::cout << " INFO:    STRUCTURE FIELD IS NOT SOLVED; MAKE SURE YOU HAVE CONSTRAINED ALL DOFS "
                 "IN YOUR STRUCTURE WITH A DBC\n";
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                 "++++++++++++++++++++++++++++++++\n";
  }
}

/*----------------------------------------------------------------------*
 |  check tangent stiffness matrix via finite differences     vuong 01/12 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::poro_fd_check()
{
  std::cout << "\n******************finite difference check***************" << std::endl;

  int dof_struct = (structure_field()->dof_row_map()->NumGlobalElements());
  int dof_fluid = (fluid_field()->dof_row_map()->NumGlobalElements());

  std::cout << "structure field has " << dof_struct << " DOFs" << std::endl;
  std::cout << "fluid field has " << dof_fluid << " DOFs" << std::endl;

  std::shared_ptr<Core::LinAlg::Vector<double>> iterinc = nullptr;
  std::shared_ptr<Core::LinAlg::Vector<double>> abs_iterinc = nullptr;
  iterinc = Core::LinAlg::create_vector(*dof_row_map(), true);
  abs_iterinc = Core::LinAlg::create_vector(*dof_row_map(), true);

  const int dofs = iterinc->global_length();
  std::cout << "in total " << dofs << " DOFs" << std::endl;
  const double delta = 1e-8;

  iterinc->put_scalar(0.0);

  iterinc->replace_global_value(0, 0, delta);

  abs_iterinc->update(1.0, *iterinc_, 0.0);

  std::shared_ptr<Epetra_CrsMatrix> stiff_approx = nullptr;
  stiff_approx = Core::LinAlg::create_matrix(*dof_row_map(), 81);

  Core::LinAlg::Vector<double> rhs_old(*dof_row_map(), true);
  rhs_old.update(1.0, *rhs_, 0.0);
  Core::LinAlg::Vector<double> rhs_copy(*dof_row_map(), true);

  std::shared_ptr<Core::LinAlg::SparseMatrix> sparse = systemmatrix_->merge();
  Core::LinAlg::SparseMatrix sparse_copy(sparse->epetra_matrix(), Core::LinAlg::Copy);


  const int zeilennr = -1;
  const int spaltenr = -1;
  for (int i = 0; i < dofs; ++i)
  {
    if (combined_dbc_map()->MyGID(i))
    {
      iterinc->replace_global_value(i, 0, 0.0);
    }
    abs_iterinc->update(1.0, *iterinc, 1.0);

    if (i == spaltenr)
      std::cout << "\n******************" << spaltenr + 1 << ". Spalte!!***************"
                << std::endl;

    evaluate(iterinc);

    rhs_copy.update(1.0, *rhs_, 0.0);

    iterinc_->put_scalar(0.0);  // Useful? depends on solver and more
    Core::LinAlg::apply_dirichlet_to_system(
        sparse_copy, *iterinc_, rhs_copy, *zeros_, *combined_dbc_map());
    std::shared_ptr<Epetra_CrsMatrix> test_crs = sparse_copy.epetra_matrix();
    int sparsenumentries;
    int sparselength = test_crs->NumGlobalEntries(i);
    std::vector<double> sparsevalues(sparselength);
    std::vector<int> sparseindices(sparselength);
    // int sparseextractionstatus =
    test_crs->ExtractGlobalRowCopy(
        i, sparselength, sparsenumentries, sparsevalues.data(), sparseindices.data());


    if (i == spaltenr)
    {
      std::cout << "rhs_: " << (rhs_copy)[zeilennr] << std::endl;
      std::cout << "rhs_old: " << (rhs_old)[zeilennr] << std::endl;
    }
    // rhs_copy = ( rhs_disturb - rhs_old ) . (-1)/delta with rhs_copy==rhs_disturb
    rhs_copy.update(-1.0, rhs_old, 1.0);
    rhs_copy.scale(-1.0 / delta);

    if (i == spaltenr)
    {
      std::cout << "( rhs_disturb - rhs_old )               "
                << (rhs_copy)[zeilennr] * (-1.0) * delta << std::endl;
      std::cout << "( rhs_disturb - rhs_old ) . (-1)/delta: " << (rhs_copy)[zeilennr] << std::endl;
    }
    int* index = &i;
    for (int j = 0; j < dofs; ++j)
    {
      double value = (rhs_copy)[j];
      stiff_approx->InsertGlobalValues(j, 1, &value, index);

      if ((j == zeilennr) and (i == spaltenr))
      {
        std::cout << "\n******************" << zeilennr + 1 << ". Zeile!!***************"
                  << std::endl;
        // std::cout << "disp: " << std::endl << *(structure_field()->dispnp()->get);
        // std::cout << "gridvel struct" << std::endl << *(structure_field()->velnp());

        // std::cout << "stiff_apprx(" << zeilennr << "," << spaltenr << "): " <<
        // (*rhs_copy)[zeilennr]
        //           << std::endl;

        // std::cout << "value(" << zeilennr << "," << spaltenr << "): " << value << std::endl;
        std::cout << "\n******************" << zeilennr + 1 << ". Zeile End!!***************"
                  << std::endl;
      }
    }

    if (not combined_dbc_map()->MyGID(i)) iterinc->replace_global_value(i, 0, -delta);

    iterinc->replace_global_value(i - 1, 0, 0.0);

    if (i != dofs - 1) iterinc->replace_global_value(i + 1, 0, delta);

    if (i == spaltenr)
      std::cout << "\n******************" << spaltenr + 1 << ". Spalte End!!***************"
                << std::endl;
  }

  evaluate(iterinc);

  stiff_approx->FillComplete();

  std::shared_ptr<Core::LinAlg::SparseMatrix> stiff_approx_sparse = nullptr;
  stiff_approx_sparse =
      std::make_shared<Core::LinAlg::SparseMatrix>(stiff_approx, Core::LinAlg::Copy);

  stiff_approx_sparse->add(sparse_copy, false, -1.0, 1.0);

  std::shared_ptr<Epetra_CrsMatrix> sparse_crs = sparse_copy.epetra_matrix();

  std::shared_ptr<Epetra_CrsMatrix> error_crs = stiff_approx_sparse->epetra_matrix();

  error_crs->FillComplete();
  sparse_crs->FillComplete();

  bool success = true;
  double error_max = 0.0;
  double abs_error_max = 0.0;
  for (int i = 0; i < dofs; ++i)
  {
    if (not combined_dbc_map()->MyGID(i))
    {
      for (int j = 0; j < dofs; ++j)
      {
        if (not combined_dbc_map()->MyGID(j))
        {
          double stiff_approx_ij = 0.0;
          double sparse_ij = 0.0;
          double error_ij = 0.0;

          {
            // get error_crs entry ij
            int errornumentries;
            int errorlength = error_crs->NumGlobalEntries(i);
            std::vector<double> errorvalues(errorlength);
            std::vector<int> errorindices(errorlength);
            // int errorextractionstatus =
            error_crs->ExtractGlobalRowCopy(
                i, errorlength, errornumentries, errorvalues.data(), errorindices.data());
            for (int k = 0; k < errorlength; ++k)
            {
              if (errorindices[k] == j)
              {
                error_ij = errorvalues[k];
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
            sparse_crs->ExtractGlobalRowCopy(
                i, sparselength, sparsenumentries, sparsevalues.data(), sparseindices.data());
            for (int k = 0; k < sparselength; ++k)
            {
              if (sparseindices[k] == j)
              {
                sparse_ij = sparsevalues[k];
                break;
              }
              else
                sparse_ij = 0.0;
            }
          }

          // get stiff_approx entry ij
          {
            int approxnumentries;
            int approxlength = stiff_approx->NumGlobalEntries(i);
            std::vector<double> approxvalues(approxlength);
            std::vector<int> approxindices(approxlength);
            // int approxextractionstatus =
            stiff_approx->ExtractGlobalRowCopy(
                i, approxlength, approxnumentries, approxvalues.data(), approxindices.data());
            for (int k = 0; k < approxlength; ++k)
            {
              if (approxindices[k] == j)
              {
                stiff_approx_ij = approxvalues[k];
                break;
              }
              else
                stiff_approx_ij = 0.0;
            }
          }

          double error = 0.0;
          if (abs(stiff_approx_ij) > 1e-5)
            error = error_ij / (stiff_approx_ij);
          else if (abs(sparse_ij) > 1e-5)
            error = error_ij / (sparse_ij);

          if (abs(error) > abs(error_max)) error_max = abs(error);
          if (abs(error_ij) > abs(abs_error_max)) abs_error_max = abs(error_ij);

          if ((abs(error) > 1e-4))
          {
            if ((abs(error_ij) > 1e-5))
            {
              // if( (sparse_ij>1e-1) or (stiff_approx_ij>1e-1) )
              {
                std::cout << "finite difference check failed entry (" << i << "," << j
                          << ")! stiff: " << sparse_ij << ", approx: " << stiff_approx_ij
                          << " ,abs. error: " << error_ij << " , rel. error: " << error
                          << std::endl;

                success = false;
              }
            }
          }
        }
      }
    }
  }

  if (success)
  {
    std::cout << "finite difference check successful, max. rel. error: " << error_max
              << "  (max. abs. error: " << abs_error_max << ")" << std::endl;
    std::cout << "******************finite difference check done***************\n\n" << std::endl;
  }
  else
    FOUR_C_THROW("PoroFDCheck failed in step: {}, iter: {}", step(), itnum_);

  return;
}

/*----------------------------------------------------------------------*
 | constructor                                         kremheller 05/18 |
 *----------------------------------------------------------------------*/
POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWayArteryCoupling::
    PoroMultiPhaseMonolithicTwoWayArteryCoupling(
        MPI_Comm comm, const Teuchos::ParameterList& globaltimeparams)
    : PoroMultiPhaseMonolithicTwoWay(comm, globaltimeparams)
{
  blockrowdofmap_artporo_ = std::make_shared<Core::LinAlg::MultiMapExtractor>();

  return;
}

/*----------------------------------------------------------------------*
 | setup the map                                       kremheller 04/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWayArteryCoupling::setup_maps()
{
  std::vector<std::shared_ptr<const Epetra_Map>> vecSpaces;

  vecSpaces.push_back(struct_dof_row_map());

  vecSpaces.push_back(fluid_dof_row_map());

  vecSpaces.push_back(artery_dof_row_map());

  if (vecSpaces[0]->NumGlobalElements() == 0) FOUR_C_THROW("No structure equation. Panic.");
  if (vecSpaces[1]->NumGlobalElements() == 0) FOUR_C_THROW("No fluid equation. Panic.");
  if (vecSpaces[2]->NumGlobalElements() == 0) FOUR_C_THROW("No fluid equation. Panic.");

  // full Poromultiphase-elasticity-map
  fullmap_ = Core::LinAlg::MultiMapExtractor::merge_maps(vecSpaces);

  // full Poromultiphase-elasticity-blockmap
  blockrowdofmap_->setup(*fullmap_, vecSpaces);

  // full map of artery and poromulti DOFs
  fullmap_artporo_ = Core::LinAlg::MultiMapExtractor::merge_maps({vecSpaces[1], vecSpaces[2]});

  // full artery-poromulti-blockmap
  blockrowdofmap_artporo_->setup(*fullmap_artporo_, {vecSpaces[1], vecSpaces[2]});

  return;
}

/*----------------------------------------------------------------------*
 | extract field vectors for calling evaluate() of the  kremheller 04/18|
 | single fields                                                        |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWayArteryCoupling::build_convergence_norms()
{
  std::shared_ptr<const Core::LinAlg::Vector<double>> arteryrhs =
      extractor()->extract_vector(*rhs_, 2);
  std::shared_ptr<const Core::LinAlg::Vector<double>> arteryinc =
      extractor()->extract_vector(*iterinc_, 2);

  // build also norms for artery
  normrhsart_ = Utils::calculate_vector_norm(vectornormfres_, *arteryrhs);
  normincart_ = Utils::calculate_vector_norm(vectornorminc_, *arteryinc);
  arterypressnorm_ = Utils::calculate_vector_norm(
      vectornorminc_, (*fluid_field()->art_net_tim_int()->pressurenp()));

  // call base class
  PoroMultiPhaseMonolithicTwoWay::build_convergence_norms();

  return;
}
/*----------------------------------------------------------------------*
 | extract field vectors for calling evaluate() of the  kremheller 04/18|
 | single fields                                                        |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWayArteryCoupling::extract_field_vectors(
    std::shared_ptr<const Core::LinAlg::Vector<double>> x,
    std::shared_ptr<const Core::LinAlg::Vector<double>>& sx,
    std::shared_ptr<const Core::LinAlg::Vector<double>>& fx)
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWayArteryCoupling::extract_field_vectors");

  // process structure unknowns of the first field
  sx = extractor()->extract_vector(*x, 0);

  // process artery and porofluid unknowns
  std::shared_ptr<const Core::LinAlg::Vector<double>> porofluid =
      extractor()->extract_vector(*x, 1);
  std::shared_ptr<const Core::LinAlg::Vector<double>> artery = extractor()->extract_vector(*x, 2);

  std::shared_ptr<Core::LinAlg::Vector<double>> dummy =
      std::make_shared<Core::LinAlg::Vector<double>>(*fullmap_artporo_);

  blockrowdofmap_artporo_->insert_vector(*porofluid, 0, *dummy);
  blockrowdofmap_artporo_->insert_vector(*artery, 1, *dummy);

  fx = dummy;

  return;
}

/*----------------------------------------------------------------------*
 | setup system matrix of poromultiphase-elasticity with artery         |
 | coupling                                           kremheller 05/18  |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWayArteryCoupling::setup_system_matrix(
    Core::LinAlg::BlockSparseMatrixBase& mat)
{
  PoroMultiPhaseMonolithicTwoWay::setup_system_matrix(mat);

  // pure artery part
  mat.assign(2, 2, Core::LinAlg::View, artery_porofluid_sysmat()->matrix(1, 1));
  // artery-porofluid part
  mat.assign(2, 1, Core::LinAlg::View, artery_porofluid_sysmat()->matrix(1, 0));
  // porofluid-artery part
  mat.assign(1, 2, Core::LinAlg::View, artery_porofluid_sysmat()->matrix(0, 1));

  return;
}

/*----------------------------------------------------------------------*
 | setup rhs of poromultiphase-elasticity with artery coupling          |
 |                                                    kremheller 05/18  |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWayArteryCoupling::setup_rhs()
{
  // get structure part
  std::shared_ptr<Core::LinAlg::Vector<double>> str_rhs = setup_structure_partof_rhs();

  // insert and scale
  extractor()->insert_vector(*str_rhs, 0, *rhs_);
  rhs_->scale(-1.0);

  // insert artery part and porofluid part
  extractor()->insert_vector(
      *(blockrowdofmap_artporo_->extract_vector(*fluid_field()->artery_porofluid_rhs(), 0)), 1,
      *rhs_);
  extractor()->insert_vector(
      *(blockrowdofmap_artporo_->extract_vector(*fluid_field()->artery_porofluid_rhs(), 1)), 2,
      *rhs_);

  return;
}

/*-----------------------------------------------------------------------/
|  build the combined dbcmap                           kremheller 05/18  |
/-----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWayArteryCoupling::build_combined_dbc_map()
{
  PoroMultiPhaseMonolithicTwoWay::build_combined_dbc_map();

  const std::shared_ptr<const Epetra_Map> artcondmap =
      fluid_field()->art_net_tim_int()->get_dbc_map_extractor()->cond_map();

  // merge them
  combinedDBCMap_ = Core::LinAlg::merge_map(combinedDBCMap_, artcondmap, false);

  return;
}
/*----------------------------------------------------------------------------*
 | build null space for artery block of global system matrix kremheller 05/18 |
 *--------------------------------------------------------------------------- */
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWayArteryCoupling::build_artery_block_null_space(
    std::shared_ptr<Core::LinAlg::Solver>& solver, const int& arteryblocknum)
{
  // equip smoother for fluid matrix block with empty parameter sublists to trigger null space
  // computation
  Teuchos::ParameterList& blocksmootherparams3 =
      solver->params().sublist("Inverse" + std::to_string(arteryblocknum));
  blocksmootherparams3.sublist("Belos Parameters");
  blocksmootherparams3.sublist("MueLu Parameters");

  // build null space of complete discretization
  fluid_field()->art_net_tim_int()->discretization()->compute_null_space_if_necessary(
      blocksmootherparams3);
  // fix the null space if some DOFs are condensed out
  Core::LinearSolver::Parameters::fix_null_space("Artery",
      *(fluid_field()->art_net_tim_int()->discretization()->dof_row_map(0)),
      *(fluid_field()->artery_dof_row_map()), blocksmootherparams3);

  return;
}

/*----------------------------------------------------------------------*
 | Create linear (iterative) solver                    kremheller 05/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWayArteryCoupling::create_linear_solver(
    const Teuchos::ParameterList& solverparams, const Core::LinearSolver::SolverType solvertype)
{
  PoroMultiPhaseMonolithicTwoWay::create_linear_solver(solverparams, solvertype);

  // build also the artery null space
  build_artery_block_null_space(solver_, 3);
}

FOUR_C_NAMESPACE_CLOSE
