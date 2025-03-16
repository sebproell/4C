// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_red_airways_tissue.hpp"

#include "4C_adapter_str_redairway.hpp"
#include "4C_fem_condition.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_red_airways_implicitintegration.hpp"

#include <stdio.h>
#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  Constructor (public)                                 yoshihara 09/12|
 *----------------------------------------------------------------------*/
Airway::RedAirwayTissue::RedAirwayTissue(MPI_Comm comm, const Teuchos::ParameterList& timeparams)
    : Adapter::AlgorithmBase(comm, timeparams)
{
  // Before setting up the structure time integrator, manipulate coupling conditions -> turn them
  // into neumann orthopressure conditions
  std::vector<Core::Conditions::Condition*> surfneumcond;
  std::vector<int> tmp;
  std::shared_ptr<Core::FE::Discretization> structdis =
      Global::Problem::instance()->get_dis("structure");
  if (structdis == nullptr) FOUR_C_THROW("no structure discretization available");

  // First get all Neumann conditions on structure
  structdis->get_condition("SurfaceNeumann", surfneumcond);
  unsigned int numneumcond = surfneumcond.size();
  if (numneumcond == 0) FOUR_C_THROW("no Neumann conditions on structure");

  // Now filter those Neumann conditions that are due to the coupling
  std::vector<Core::Conditions::Condition*> coupcond;
  for (unsigned int i = 0; i < numneumcond; ++i)
  {
    Core::Conditions::Condition* actcond = surfneumcond[i];
    if (actcond->type() == Core::Conditions::RedAirwayTissue) coupcond.push_back(actcond);
  }
  unsigned int numcond = coupcond.size();
  if (numcond == 0) FOUR_C_THROW("no coupling conditions found");

  for (unsigned int i = 0; i < numcond; ++i)
  {
    Core::Conditions::Condition* cond = coupcond[i];
    std::string type = "orthopressure";
    cond->parameters().add("TYPE", type);
    std::vector<int> onoff(6, 0);
    onoff[0] = 1;
    cond->parameters().add("ONOFF", onoff);
    std::vector<double> val(6, 0.0);
    cond->parameters().add("VAL", val);

    int condID = coupcond[i]->parameters().get<int>("coupling_id");
    tmp.push_back(condID);
  }

  std::vector<Core::Conditions::Condition*> nodecond;
  std::shared_ptr<Core::FE::Discretization> redairwaydis =
      Global::Problem::instance()->get_dis("red_airway");
  if (redairwaydis == nullptr) FOUR_C_THROW("no redairway discretization available");

  // First get all redairway prescribed conditions on structure
  redairwaydis->get_condition("RedAirwayPrescribedCond", nodecond);
  unsigned int numnodecond = nodecond.size();
  if (numnodecond == 0)
    FOUR_C_THROW("no redairway prescribed conditions on redairway discretization");

  // Now filter those node conditions that are due to the coupling
  std::vector<Core::Conditions::Condition*> nodecoupcond;
  for (unsigned int i = 0; i < numnodecond; ++i)
  {
    Core::Conditions::Condition* actcond = nodecond[i];
    if (actcond->type() == Core::Conditions::RedAirwayNodeTissue) nodecoupcond.push_back(actcond);
  }
  unsigned int numnodecoupcond = nodecoupcond.size();
  if (numnodecoupcond == 0) FOUR_C_THROW("no coupling conditions found");

  for (unsigned int i = 0; i < numnodecoupcond; ++i)
  {
    Core::Conditions::Condition* cond = nodecoupcond[i];
    std::string bc_data = "flow";
    cond->parameters().add("boundarycond", bc_data);
    std::vector<double> val(1, 0.0);
    cond->parameters().add("VAL", val);
  }


  Epetra_Map redundantmap(
      tmp.size(), tmp.size(), tmp.data(), 0, Core::Communication::as_epetra_comm(comm));
  couppres_ip_ = std::make_shared<Core::LinAlg::Vector<double>>(redundantmap, true);
  couppres_ip_tilde_ = std::make_shared<Core::LinAlg::Vector<double>>(redundantmap, true);
  couppres_im_ = std::make_shared<Core::LinAlg::Vector<double>>(redundantmap, true);
  couppres_im_tilde_ = std::make_shared<Core::LinAlg::Vector<double>>(redundantmap, true);
  couppres_il_ = std::make_shared<Core::LinAlg::Vector<double>>(redundantmap, true);
  omega_np_ = std::make_shared<Core::LinAlg::Vector<double>>(redundantmap, true);
  coupflux_ip_ = std::make_shared<Core::LinAlg::Vector<double>>(redundantmap, true);
  coupflux_im_ = std::make_shared<Core::LinAlg::Vector<double>>(redundantmap, true);
  coupvol_ip_ = std::make_shared<Core::LinAlg::Vector<double>>(redundantmap, true);
  coupvol_im_ = std::make_shared<Core::LinAlg::Vector<double>>(redundantmap, true);

  const Teuchos::ParameterList& sdyn = Global::Problem::instance()->structural_dynamic_params();

  std::shared_ptr<Adapter::StructureBaseAlgorithm> structure =
      std::make_shared<Adapter::StructureBaseAlgorithm>(
          sdyn, const_cast<Teuchos::ParameterList&>(sdyn), structdis);
  structure_ = std::dynamic_pointer_cast<Adapter::StructureRedAirway>(structure->structure_field());
  structure_->setup();

  setup_red_airways();
  const Teuchos::ParameterList& rawdyn =
      Global::Problem::instance()->reduced_d_airway_dynamic_params();

  // Check Time integration parameters
  if (sdyn.get<double>("TIMESTEP") != timeparams.get<double>("TIMESTEP") or
      sdyn.get<int>("NUMSTEP") != timeparams.get<int>("NUMSTEP") or
      sdyn.get<double>("MAXTIME") != timeparams.get<double>("MAXTIME") or
      rawdyn.get<double>("TIMESTEP") != timeparams.get<double>("TIMESTEP") or
      rawdyn.get<int>("NUMSTEP") != timeparams.get<int>("NUMSTEP"))
    FOUR_C_THROW("Parameter(s) for time integrators inconsistent");

  // Check Time integration parameters
  if (sdyn.get<int>("RESTARTEVERY") != rawdyn.get<int>("RESTARTEVERY"))
    FOUR_C_THROW("Parameters for restart inconsistent");

  // Get coupling parameters
  const Teuchos::ParameterList& rawtisdyn =
      Global::Problem::instance()->red_airway_tissue_dynamic_params();
  // Get max iterations
  itermax_ = rawtisdyn.get<int>("MAXITER");

  // Get tolerance for pressure
  tolp_ = rawtisdyn.get<double>("CONVTOL_P");

  // Get tolerance for flux
  tolq_ = rawtisdyn.get<double>("CONVTOL_Q");

  // Dynamic relaxation type
  relaxtype_ =
      Teuchos::getIntegralValue<Inpar::ArteryNetwork::Relaxtype3D0D>(rawtisdyn, "RELAXTYPE");

  // Get normal direction
  // -> if normal == 1.0 : the pressure will be implemented from inside the element to the outside
  // -> if normal ==-1.0 : the pressure will be implemented from outside the element to the inside
  normal_ = rawtisdyn.get<double>("NORMAL");

  // Determine initial volume
  structure_->init_vol();

}  // end of constructor


/*----------------------------------------------------------------------*
 |  Read restart                                              roth 10/13|
 *----------------------------------------------------------------------*/
void Airway::RedAirwayTissue::read_restart(const int step)
{
  structure_->read_restart(step);
  structure_->init_vol();
  redairways_->read_restart(step);

  // Read the coupling variables at restart
  redairways_->get_restart_reader(step)->read_vector(couppres_im_, "couppres_im");
  redairways_->get_restart_reader(step)->read_vector(coupflux_im_, "coupflux_im");
  redairways_->get_restart_reader(step)->read_vector(coupflux_ip_, "coupflux_im");
  redairways_->get_restart_reader(step)->read_vector(coupvol_im_, "coupvol_im");

  // Set timestep and time after restart
  double restartTime = redairways_->time();
  set_time_step(restartTime, step);
}


/*----------------------------------------------------------------------*
 |  Integrate                                            yoshihara 09/12|
 *----------------------------------------------------------------------*/
void Airway::RedAirwayTissue::integrate()
{
  while (not_finished())
  {
    int iter = 0;
    increment_time_and_step();
    redairways_->prepare_time_step();

    do
    {
      do_red_airway_step();
      relax_pressure(iter);
      do_structure_step();
      iter++;
    } while (not_converged(iter) && iter < itermax_);

    if ((iter >= itermax_) && (Core::Communication::my_mpi_rank(couppres_ip_->get_comm()) == 0))
    {
      FOUR_C_THROW("FIELD ITERATION NOT CONVERGED IN {} STEPS AT TIME T={}", itermax_, time());
    }

    update_and_output();
  }
}


/*----------------------------------------------------------------------*
 |  Integrate airways time step and calculate pressures  yoshihara 09/12|
 *----------------------------------------------------------------------*/
void Airway::RedAirwayTissue::do_red_airway_step()
{
  // Scale with -1 (redairway convention: outflow is negative)
  coupflux_ip_->scale(-1.0);

  redairways_->set_airway_flux_from_tissue(*coupflux_ip_);
  redairways_->integrate_step();
  redairways_->extract_pressure(*couppres_ip_tilde_);
  couppres_ip_tilde_->update(0.0, *couppres_ip_tilde_, normal_);
}

/*----------------------------------------------------------------------*
 |  Relax the coupling pressure as described in Kuettler et al.:        |
 |  "Fixed-point-fluid-structure interaction solvers with dynamic       |
 |   relaxation" (2008)                                       roth 10/13|
 *----------------------------------------------------------------------*/
void Airway::RedAirwayTissue::relax_pressure(int iter)
{
  switch (relaxtype_)
  {
    case Inpar::ArteryNetwork::norelaxation:
    {
      printf("No dynamic relaxation \n");
      couppres_ip_->update(1.0, *couppres_ip_tilde_, 0.0);
    }
    break;

    case Inpar::ArteryNetwork::fixedrelaxation:
    {
      // Fixed \omega_i
      omega_ = 0.85;
      printf("Fixed Relaxation Parameter: %f \n", omega_);

      // Dynamic Relaxation formula (35)
      // p^{n+1}_{i+1} = \omega_i \tilde{p}^{n+1}_{i+1} + (1-\omega_i) p^{n+1}_{i}
      // where \tilde{p}^{n+1}_{i+1} = couppres_ip_
      //                 p^{n+1}_{i} = couppres_im_
      couppres_ip_->update(omega_, *couppres_ip_tilde_, (1.0 - omega_), *couppres_im_, 0.0);
    }
    break;

    case Inpar::ArteryNetwork::Aitken:
    {
      // \omega_{i+1} = omega_np_
      // Kuettler (2008), Remark 10: Two previous steps are required to calculate the \omega,
      // thus the relaxation factor is fixed for the first two steps to \omega = 1.0
      if (iter < 2)
      {
        omega_np_->put_scalar(1.0);
      }
      else if (iter >= 2)
      {
        // Relaxation factor formula (41)
        // \omega_{i+1} = (p_{i} - p_{i+1}) / (p_{i} - \tilde{p}_{i+1} - p_{i+1} + \tilde{p}_{i+2})
        //      where p_{i} = couppres_il_
        //          p_{i+1} = couppres_im_
        //  \tilde{p}_{i+1} = couppres_im_tilde_
        //  \tilde{p}_{i+2} = couppres_ip_tilde_
        omega_np_->update(1.0, *couppres_il_, -1.0, *couppres_im_, 0.0);

        Core::LinAlg::Vector<double> denominator(*omega_np_);
        denominator.update(-1.0, *couppres_im_tilde_, +1.0, *couppres_ip_tilde_, 1.0);

        omega_np_->reciprocal_multiply(1.0, denominator, *omega_np_, 0.0);

        // Safety check for \omega_i+1
        for (int i = 0; i < couppres_ip_->get_map().NumMyElements(); ++i)
        {
          if ((*omega_np_)[i] < 0.0)
            (*omega_np_)[i] = 0.0;
          else if ((*omega_np_)[i] > 1.0)
            (*omega_np_)[i] = 1.0;
        }
      }

      // Aitken Relaxation formula (35)
      // p^{n+1}_{i+1} = \omega_i \tilde{p}^{n+1}_{i+1} + (1-\omega_i) p^{n+1}_{i}
      for (int i = 0; i < couppres_ip_->get_map().NumMyElements(); ++i)
      {
        (*couppres_ip_)[i] =
            (*omega_np_)[i] * (*couppres_ip_tilde_)[i] + (1 - (*omega_np_)[i]) * (*couppres_im_)[i];
      }

      // Print relaxation factor \omega_np_
      if (Core::Communication::my_mpi_rank(couppres_ip_->get_comm()) == 0)
      {
        printf("Aitken Relaxation: \n");
        for (int i = 0; i < couppres_ip_->get_map().NumMyElements(); ++i)
          std::cout << "omega_np_[" << i << "]: " << (*omega_np_)[i] << std::endl;
        std::cout << std::endl;
      }
    }
    break;

    // not implemented yet...
    case Inpar::ArteryNetwork::SD:
    {
      FOUR_C_THROW(
          "Currently only two types of RELAXTYPE possible: norelaxation, fixedrelaxation, Aitken");
    }
    break;

    default:
    {
      FOUR_C_THROW(
          "Currently only two types of RELAXTYPE possible: norelaxation, fixedrelaxation, Aitken");
    }
    break;
  }

  /*/Print coupling pressure for validation of dynamic relaxation
  std::cout << "\nSingle couppres_ip_: " << std::endl;
  for (int i=0; i<couppres_ip_->Map().NumMyElements(); ++i)
  {
    printf(" Time: %f couppres_ip_: %6.3e \n",Time(), (*couppres_ip_)[i]);
  }*/
}


/*----------------------------------------------------------------------*
 |  Integrate structure time step and calculate fluxes   yoshihara 09/12|
 *----------------------------------------------------------------------*/
void Airway::RedAirwayTissue::do_structure_step()
{
  structure_->set_pressure(*couppres_ip_);
  structure_->prepare_time_step();
  structure_->solve();
  structure_->calc_flux(*coupflux_ip_, *coupvol_ip_, dt());
}


/*----------------------------------------------------------------------*
 |  Check for convergence between fields                 yoshihara 09/12|
 *----------------------------------------------------------------------*/
/* Note: This has to be done for each field on its own, not as a Norm2
 * over all fields.
 */
bool Airway::RedAirwayTissue::not_converged(int iter)
{
  Core::LinAlg::Vector<double> pres_inc(*couppres_ip_);
  Core::LinAlg::Vector<double> scaled_pres_inc(*couppres_ip_);
  Core::LinAlg::Vector<double> flux_inc(*coupflux_ip_);
  Core::LinAlg::Vector<double> scaled_flux_inc(*coupflux_ip_);

  // Calculate Pressure Norm
  for (int i = 0; i < couppres_ip_->get_map().NumMyElements(); ++i)
  {
    // Calculate pressure increment
    (pres_inc)[i] = abs((*couppres_ip_)[i] - (*couppres_im_)[i]);

    // Calculate scaled pressure increment
    if (abs((*couppres_ip_)[i]) > 1e-05)
      (scaled_pres_inc)[i] = (pres_inc)[i] / abs((*couppres_ip_)[i]);
    else
      (scaled_pres_inc)[i] = (pres_inc)[i];
  }

  // Calculate Flux Norm
  for (int i = 0; i < coupflux_ip_->get_map().NumMyElements(); ++i)
  {
    // Calculate flux increment
    (flux_inc)[i] = abs((*coupflux_ip_)[i] - (*coupflux_im_)[i]);

    // Calculate scaled flux increment
    if (abs((*coupflux_ip_)[i]) > 1e-05)
      (scaled_flux_inc)[i] = (flux_inc)[i] / abs((*coupflux_ip_)[i]);
    else
      (scaled_flux_inc)[i] = (flux_inc)[i];
  }

  // Output
  output_iteration(pres_inc, scaled_pres_inc, flux_inc, scaled_flux_inc, iter);

  // Update values
  couppres_il_->update(1.0, *couppres_im_, 0.0);
  couppres_im_->update(1.0, *couppres_ip_, 0.0);
  couppres_im_tilde_->update(1.0, *couppres_ip_tilde_, 0.0);
  coupflux_im_->update(1.0, *coupflux_ip_, 0.0);
  coupvol_im_->update(1.0, *coupvol_ip_, 0.0);

  double pres_max, flux_max;
  scaled_pres_inc.norm_inf(&pres_max);
  scaled_flux_inc.norm_inf(&flux_max);

  if (pres_max < tolp_ and flux_max < tolq_ and iter > 1) return false;

  return true;
}


/*----------------------------------------------------------------------*
 |  Output of one iteration between fields               yoshihara 09/12|
 *----------------------------------------------------------------------*/
void Airway::RedAirwayTissue::output_iteration(Core::LinAlg::Vector<double>& pres_inc,
    Core::LinAlg::Vector<double>& scaled_pres_inc, Core::LinAlg::Vector<double>& flux_inc,
    Core::LinAlg::Vector<double>& scaled_flux_inc, int iter)
{
  if (Core::Communication::my_mpi_rank(couppres_ip_->get_comm()) == 0)
  {
    printf("\nFIELD ITERATION: %i / %i\n", iter, itermax_);
    printf(
        " Tolerances:                                                                        %4.2e "
        "     %4.2e\n",
        tolp_, tolq_);
    printf(
        " Volume ID      Vol            P             Q            dP            dQ           "
        "dP_scal       dQ_scal\n");
    for (int i = 0; i < couppres_ip_->get_map().NumMyElements(); ++i)
    {
      printf("     %d       %4.3e     %4.3e     %4.3e     %4.3e     %4.3e     %4.3e     %4.3e\n",
          couppres_ip_->get_map().GID(i), (*coupvol_ip_)[i], (*couppres_ip_)[i], (*coupflux_ip_)[i],
          (pres_inc)[i], (flux_inc)[i], (scaled_pres_inc)[i], (scaled_flux_inc)[i]);
    }
    printf(
        "------------------------------------------------------------------------------------------"
        "----------------\n\n");
  }
}


/*----------------------------------------------------------------------*
 |  Update and output                                    yoshihara 09/12|
 *----------------------------------------------------------------------*/
void Airway::RedAirwayTissue::update_and_output()
{
  constexpr bool force_prepare = false;
  structure_->prepare_output(force_prepare);
  structure_->update();
  structure_->output();

  redairways_->time_update();
  redairways_->output();

  // In case of restart write all coupling variables to restart file
  if (redairways_->step() % uprestart_ == 0)
  {
    redairways_->get_output_writer().write_vector("couppres_im", couppres_im_);
    redairways_->get_output_writer().write_vector("coupflux_im", coupflux_im_);
    redairways_->get_output_writer().write_vector("coupvol_im", coupvol_im_);
  }
}


void Airway::RedAirwayTissue::setup_red_airways()
{
  // Access the discretization
  std::shared_ptr<Core::FE::Discretization> actdis = nullptr;
  actdis = Global::Problem::instance()->get_dis("red_airway");

  // Set degrees of freedom in the discretization
  if (!actdis->filled())
  {
    actdis->fill_complete();
  }

  // Context for output and restart
  std::shared_ptr<Core::IO::DiscretizationWriter> output = actdis->writer();
  output->write_mesh(0, 0.0);

  // Set some pointers and variables
  const Teuchos::ParameterList& rawdyn =
      Global::Problem::instance()->reduced_d_airway_dynamic_params();

  // Create a solver
  // Get the solver number
  const int linsolvernumber = rawdyn.get<int>("LINEAR_SOLVER");
  // Check if the present solver has a valid solver number
  if (linsolvernumber == (-1))
    FOUR_C_THROW(
        "no linear solver defined. Please set LINEAR_SOLVER in REDUCED DIMENSIONAL AIRWAYS DYNAMIC "
        "to a valid number!");
  std::unique_ptr<Core::LinAlg::Solver> solver = std::make_unique<Core::LinAlg::Solver>(
      Global::Problem::instance()->solver_params(linsolvernumber), actdis->get_comm(),
      Global::Problem::instance()->solver_params_callback(),
      Teuchos::getIntegralValue<Core::IO::Verbositylevel>(
          Global::Problem::instance()->io_params(), "VERBOSITY"));
  actdis->compute_null_space_if_necessary(solver->params());

  // Set parameters in list required for all schemes
  Teuchos::ParameterList airwaystimeparams;

  // Number of degrees of freedom
  const int ndim = Global::Problem::instance()->n_dim();
  airwaystimeparams.set<int>("number of degrees of freedom", 1 * ndim);

  // Time integration
  // Default time step size
  airwaystimeparams.set<double>("time step size", rawdyn.get<double>("TIMESTEP"));
  // Maximum number of timesteps
  airwaystimeparams.set<int>("max number timesteps", rawdyn.get<int>("NUMSTEP"));

  // Restart and output
  // Restart
  airwaystimeparams.set("write restart every", rawdyn.get<int>("RESTARTEVERY"));
  // Solution output
  airwaystimeparams.set("write solution every", rawdyn.get<int>("RESULTSEVERY"));

  // Solver parameters
  // Solver type
  airwaystimeparams.set(
      "solver type", Teuchos::getIntegralValue<RedAirwaysDyntype>(rawdyn, "SOLVERTYPE"));
  // Tolerance
  airwaystimeparams.set("tolerance", rawdyn.get<double>("TOLERANCE"));
  // Maximum number of iterations
  airwaystimeparams.set("maximum iteration steps", rawdyn.get<int>("MAXITERATIONS"));
  // compute Interdependency
  if (rawdyn.get<bool>("COMPAWACINTER"))
    airwaystimeparams.set("CompAwAcInter", true);
  else
    airwaystimeparams.set("CompAwAcInter", false);
  // Adjust acini volume with pre-stress condition
  if (rawdyn.get<bool>("CALCV0PRESTRESS"))
  {
    airwaystimeparams.set("CalcV0PreStress", true);
    airwaystimeparams.set("transpulmpress", rawdyn.get<double>("TRANSPULMPRESS"));
  }
  else
    airwaystimeparams.set("CalcV0PreStress", false);

  // Get restart timestep
  uprestart_ = rawdyn.get<int>("RESTARTEVERY");

  //------------------------------------------------------------------
  // create all vectors and variables associated with the time
  // integration (call the constructor);
  // the only parameter from the list required here is the number of
  // velocity degrees of freedom
  //------------------------------------------------------------------
  redairways_ = std::make_shared<Airway::RedAirwayImplicitTimeInt>(
      actdis, std::move(solver), airwaystimeparams, *output);

  redairways_->setup_for_coupling();
}

FOUR_C_NAMESPACE_CLOSE
