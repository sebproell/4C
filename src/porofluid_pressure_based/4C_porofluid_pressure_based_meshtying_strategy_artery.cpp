// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_porofluid_pressure_based_meshtying_strategy_artery.hpp"

#include "4C_art_net_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_bio.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_utils_sparse_algebra_print.hpp"
#include "4C_linear_solver_method.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_linear_solver_method_parameters.hpp"
#include "4C_porofluid_pressure_based_algorithm.hpp"
#include "4C_porofluid_pressure_based_elast_scatra_artery_coupling_base.hpp"
#include "4C_porofluid_pressure_based_elast_scatra_utils.hpp"
#include "4C_porofluid_pressure_based_utils.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | constructor                                (public) kremheller 04/18 |
 *----------------------------------------------------------------------*/
PoroPressureBased::MeshtyingArtery::MeshtyingArtery(
    PoroPressureBased::PorofluidAlgorithm* porofluid_algorithm,
    const Teuchos::ParameterList& probparams, const Teuchos::ParameterList& poroparams)
    : porofluid_algorithm_(porofluid_algorithm),
      params_(probparams),
      poroparams_(poroparams),
      vectornormfres_(
          Teuchos::getIntegralValue<PoroPressureBased::VectorNorm>(poroparams_, "VECTORNORM_RESF")),
      vectornorminc_(
          Teuchos::getIntegralValue<PoroPressureBased::VectorNorm>(poroparams_, "VECTORNORM_INC"))
{
  const Teuchos::ParameterList& artdyn = Global::Problem::instance()->arterial_dynamic_params();

  arterydis_ = Global::Problem::instance()->get_dis("artery");

  if (!arterydis_->filled()) arterydis_->fill_complete();

  auto timintscheme =
      Teuchos::getIntegralValue<Inpar::ArtDyn::TimeIntegrationScheme>(artdyn, "DYNAMICTYPE");

  std::shared_ptr<Core::IO::DiscretizationWriter> artery_output = arterydis_->writer();
  artery_output->write_mesh(0, 0.0);

  // build art net time integrator
  artnettimint_ = Arteries::Utils::create_algorithm(timintscheme, arterydis_,
      artdyn.get<int>("LINEAR_SOLVER"), probparams, artdyn, *artery_output);

  // set to false
  artnettimint_->set_solve_scatra(false);

  // initialize
  artnettimint_->init(probparams, artdyn, "artery_scatra");

  // print user info
  if (Core::Communication::my_mpi_rank(porofluid_algorithm->discretization()->get_comm()) == 0)
  {
    std::cout << "\n";
    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    std::cout << "<                                                  >" << std::endl;
    std::cout << "<    Coupling with 1D Artery Network activated     >" << std::endl;
  }

  const bool evaluate_on_lateral_surface =
      poroparams.sublist("ARTERY COUPLING").get<bool>("LATERAL_SURFACE_COUPLING");

  const std::string couplingcondname = std::invoke(
      [&]()
      {
        if (Teuchos::getIntegralValue<
                Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod>(
                Global::Problem::instance()->poro_fluid_multi_phase_dynamic_params().sublist(
                    "ARTERY COUPLING"),
                "ARTERY_COUPLING_METHOD") ==
            Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod::ntp)
        {
          return "ArtPorofluidCouplConNodeToPoint";
        }
        else
        {
          return "ArtPorofluidCouplConNodebased";
        }
      });

  // initialize mesh tying object
  arttoporofluidcoupling_ = PoroPressureBased::create_and_init_artery_coupling_strategy(arterydis_,
      porofluid_algorithm->discretization(), poroparams.sublist("ARTERY COUPLING"),
      couplingcondname, "COUPLEDDOFS_ART", "COUPLEDDOFS_PORO", evaluate_on_lateral_surface);

  // Initialize rhs vector
  rhs_ = std::make_shared<Core::LinAlg::Vector<double>>(*arttoporofluidcoupling_->full_map(), true);

  // Initialize increment vector
  comb_increment_ =
      std::make_shared<Core::LinAlg::Vector<double>>(*arttoporofluidcoupling_->full_map(), true);
  // Initialize phinp vector
  comb_phinp_ =
      std::make_shared<Core::LinAlg::Vector<double>>(*arttoporofluidcoupling_->full_map(), true);

  // initialize Poromultiphase-elasticity-systemmatrix_
  comb_systemmatrix_ =
      std::make_shared<Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>>(
          *arttoporofluidcoupling_->global_extractor(),
          *arttoporofluidcoupling_->global_extractor(), 81, false, true);

  return;
}



/*----------------------------------------------------------------------*
 | prepare time loop                                   kremheller 04/18 |
 *----------------------------------------------------------------------*/
void PoroPressureBased::MeshtyingArtery::prepare_time_loop()
{
  artnettimint_->prepare_time_loop();
  return;
}

/*----------------------------------------------------------------------*
 | setup the variables to do a new time step  (public) kremheller 04/18 |
 *----------------------------------------------------------------------*/
void PoroPressureBased::MeshtyingArtery::prepare_time_step()
{
  artnettimint_->prepare_time_step();
  return;
}

/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                     kremheller 04/18 |
 *----------------------------------------------------------------------*/
void PoroPressureBased::MeshtyingArtery::update()
{
  artnettimint_->time_update();
  return;
}

/*--------------------------------------------------------------------------*
 | initialize the linear solver                            kremheller 07/20 |
 *--------------------------------------------------------------------------*/
void PoroPressureBased::MeshtyingArtery::initialize_linear_solver(
    std::shared_ptr<Core::LinAlg::Solver> solver)
{
  const Teuchos::ParameterList& porofluidparams =
      Global::Problem::instance()->poro_fluid_multi_phase_dynamic_params();
  const int linsolvernumber = porofluidparams.get<int>("LINEAR_SOLVER");
  const Teuchos::ParameterList& solverparams =
      Global::Problem::instance()->solver_params(linsolvernumber);
  const auto solvertype =
      Teuchos::getIntegralValue<Core::LinearSolver::SolverType>(solverparams, "SOLVER");
  // no need to do the rest for direct solvers
  if (solvertype == Core::LinearSolver::SolverType::umfpack or
      solvertype == Core::LinearSolver::SolverType::superlu)
    return;

  if (solvertype != Core::LinearSolver::SolverType::belos)
    FOUR_C_THROW("Iterative solver expected");

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
      FOUR_C_THROW("Block Gauss-Seidel preconditioner expected.");
      break;
  }

  Teuchos::ParameterList& blocksmootherparams1 = solver->params().sublist("Inverse1");
  Core::LinearSolver::Parameters::compute_solver_parameters(
      *porofluid_algorithm_->discretization(), blocksmootherparams1);

  Teuchos::ParameterList& blocksmootherparams2 = solver->params().sublist("Inverse2");
  Core::LinearSolver::Parameters::compute_solver_parameters(*arterydis_, blocksmootherparams2);
}

/*--------------------------------------------------------------------------*
 | solve linear system of equations                        kremheller 04/18 |
 *--------------------------------------------------------------------------*/
void PoroPressureBased::MeshtyingArtery::linear_solve(std::shared_ptr<Core::LinAlg::Solver> solver,
    std::shared_ptr<Core::LinAlg::SparseOperator> sysmat,
    std::shared_ptr<Core::LinAlg::Vector<double>> increment,
    std::shared_ptr<Core::LinAlg::Vector<double>> residual,
    Core::LinAlg::SolverParams& solver_params)
{
  comb_systemmatrix_->complete();

  comb_increment_->put_scalar(0.0);

  // standard solver call
  // system is ready to solve since Dirichlet Boundary conditions have been applied in
  // setup_system_matrix or Evaluate
  solver_params.refactor = true;
  solver->solve(comb_systemmatrix_->epetra_operator(), comb_increment_, rhs_, solver_params);

  return;
}

/*----------------------------------------------------------------------*
 | Calculate problem specific norm                     kremheller 03/18 |
 *----------------------------------------------------------------------*/
void PoroPressureBased::MeshtyingArtery::calculate_norms(std::vector<double>& preresnorm,
    std::vector<double>& incprenorm, std::vector<double>& prenorm,
    const std::shared_ptr<const Core::LinAlg::Vector<double>> increment)
{
  preresnorm.resize(2);
  incprenorm.resize(2);
  prenorm.resize(2);

  prenorm[0] = calculate_vector_norm(vectornorminc_, *porofluid_algorithm_->phinp());
  prenorm[1] = calculate_vector_norm(vectornorminc_, *artnettimint_->pressurenp());

  std::shared_ptr<const Core::LinAlg::Vector<double>> arterypressinc;
  std::shared_ptr<const Core::LinAlg::Vector<double>> porofluidinc;

  arttoporofluidcoupling_->extract_single_field_vectors(
      comb_increment_, porofluidinc, arterypressinc);

  incprenorm[0] = calculate_vector_norm(vectornorminc_, *porofluidinc);
  incprenorm[1] = calculate_vector_norm(vectornorminc_, *arterypressinc);

  std::shared_ptr<const Core::LinAlg::Vector<double>> arterypressrhs;
  std::shared_ptr<const Core::LinAlg::Vector<double>> porofluidrhs;

  arttoporofluidcoupling_->extract_single_field_vectors(rhs_, porofluidrhs, arterypressrhs);

  preresnorm[0] = calculate_vector_norm(vectornormfres_, *porofluidrhs);
  preresnorm[1] = calculate_vector_norm(vectornormfres_, *arterypressrhs);

  return;
}

/*----------------------------------------------------------------------*
 | create result test for this field                   kremheller 04/18 |
 *----------------------------------------------------------------------*/
void PoroPressureBased::MeshtyingArtery::create_field_test()
{
  std::shared_ptr<Core::Utils::ResultTest> arteryresulttest = artnettimint_->create_field_test();
  Global::Problem::instance()->add_field_test(arteryresulttest);
  return;
}

/*----------------------------------------------------------------------*
 |  read restart data                                  kremheller 04/18 |
 -----------------------------------------------------------------------*/
void PoroPressureBased::MeshtyingArtery::read_restart(const int step)
{
  artnettimint_->read_restart(step);

  return;
}

/*----------------------------------------------------------------------*
 | output of solution vector to BINIO                  kremheller 04/18 |
 *----------------------------------------------------------------------*/
void PoroPressureBased::MeshtyingArtery::output()
{
  if (porofluid_algorithm_->step() != 0) artnettimint_->output(false, nullptr);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate matrix and rhs                             kremheller 04/18 |
 *----------------------------------------------------------------------*/
void PoroPressureBased::MeshtyingArtery::evaluate()
{
  arttoporofluidcoupling_->set_solution_vectors(
      porofluid_algorithm_->phinp(), porofluid_algorithm_->phin(), artnettimint_->pressurenp());

  // evaluate the coupling
  arttoporofluidcoupling_->evaluate(comb_systemmatrix_, rhs_);

  // evaluate artery
  artnettimint_->assemble_mat_and_rhs();
  // apply DBC
  artnettimint_->prepare_linear_solve();

  // SetupCoupledArteryPoroFluidSystem();
  arttoporofluidcoupling_->setup_system(comb_systemmatrix_, rhs_,
      porofluid_algorithm_->system_matrix(), artnettimint_->system_matrix(),
      porofluid_algorithm_->rhs(), artnettimint_->rhs(),
      porofluid_algorithm_->get_dbc_map_extractor(), artnettimint_->get_dbc_map_extractor());

  return;
}

/*----------------------------------------------------------------------*
 | extract and update                                  kremheller 04/18 |
 *----------------------------------------------------------------------*/
std::shared_ptr<const Core::LinAlg::Vector<double>>
PoroPressureBased::MeshtyingArtery::extract_and_update_iter(
    const std::shared_ptr<const Core::LinAlg::Vector<double>> inc)
{
  std::shared_ptr<const Core::LinAlg::Vector<double>> arterypressinc;
  std::shared_ptr<const Core::LinAlg::Vector<double>> porofluidinc;

  arttoporofluidcoupling_->extract_single_field_vectors(inc, porofluidinc, arterypressinc);

  artnettimint_->update_iter(arterypressinc);

  return porofluidinc;
}

/*----------------------------------------------------------------------*
 | artery dof row map                                  kremheller 04/18 |
 *----------------------------------------------------------------------*/
std::shared_ptr<const Core::LinAlg::Map> PoroPressureBased::MeshtyingArtery::artery_dof_row_map()
    const
{
  return arttoporofluidcoupling_->artery_dof_row_map();
}

/*-----------------------------------------------------------------------*
 | access to block system matrix of artery poro problem kremheller 04/18 |
 *-----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase>
PoroPressureBased::MeshtyingArtery::artery_porofluid_sysmat() const
{
  return comb_systemmatrix_;
}

/*----------------------------------------------------------------------*
 | return coupled residual                             kremheller 05/18 |
 *----------------------------------------------------------------------*/
std::shared_ptr<const Core::LinAlg::Vector<double>>
PoroPressureBased::MeshtyingArtery::artery_porofluid_rhs() const
{
  return rhs_;
}

/*----------------------------------------------------------------------*
 | extract and update                                  kremheller 04/18 |
 *----------------------------------------------------------------------*/
std::shared_ptr<const Core::LinAlg::Vector<double>>
PoroPressureBased::MeshtyingArtery::combined_increment(
    const std::shared_ptr<const Core::LinAlg::Vector<double>> inc) const
{
  return comb_increment_;
}

/*----------------------------------------------------------------------*
 | check initial fields                                kremheller 06/18 |
 *----------------------------------------------------------------------*/
void PoroPressureBased::MeshtyingArtery::check_initial_fields(
    std::shared_ptr<const Core::LinAlg::Vector<double>> vec_cont) const
{
  arttoporofluidcoupling_->check_initial_fields(vec_cont, artnettimint_->pressurenp());
  return;
}

/*-------------------------------------------------------------------------*
 | set element pairs that are close                       kremheller 03/19 |
 *------------------------------------------------------------------------ */
void PoroPressureBased::MeshtyingArtery::set_nearby_ele_pairs(
    const std::map<int, std::set<int>>* nearbyelepairs)
{
  arttoporofluidcoupling_->set_nearby_ele_pairs(nearbyelepairs);
  return;
}

/*-------------------------------------------------------------------------*
 | setup the strategy                                     kremheller 03/19 |
 *------------------------------------------------------------------------ */
void PoroPressureBased::MeshtyingArtery::setup()
{
  arttoporofluidcoupling_->setup();
  return;
}

/*----------------------------------------------------------------------*
 | apply mesh movement                                 kremheller 06/18 |
 *----------------------------------------------------------------------*/
void PoroPressureBased::MeshtyingArtery::apply_mesh_movement() const
{
  arttoporofluidcoupling_->apply_mesh_movement();
  return;
}

/*----------------------------------------------------------------------*
 | access to blood vessel volume fraction              kremheller 10/19 |
 *----------------------------------------------------------------------*/
std::shared_ptr<const Core::LinAlg::Vector<double>>
PoroPressureBased::MeshtyingArtery::blood_vessel_volume_fraction()
{
  return arttoporofluidcoupling_->blood_vessel_volume_fraction();
}

FOUR_C_NAMESPACE_CLOSE
