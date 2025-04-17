// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_porofluid_pressure_based_elast_scatra_partitioned_twoway.hpp"

#include "4C_adapter_art_net.hpp"
#include "4C_adapter_porofluid_pressure_based_wrapper.hpp"
#include "4C_adapter_scatra_base_algorithm.hpp"
#include "4C_adapter_str_wrapper.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_porofluid_pressure_based_elast_base.hpp"
#include "4C_scatra_timint_implicit.hpp"
#include "4C_scatra_timint_meshtying_strategy_artery.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
PoroPressureBased::PoroMultiPhaseScaTraPartitionedTwoWay::PoroMultiPhaseScaTraPartitionedTwoWay(
    MPI_Comm comm, const Teuchos::ParameterList& globaltimeparams)
    : PoroMultiPhaseScaTraPartitioned(comm, globaltimeparams),
      scaincnp_(nullptr),
      structincnp_(nullptr),
      fluidincnp_(nullptr),
      itmax_(-1),
      ittol_(-1),
      artery_coupling_active_(false)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroPressureBased::PoroMultiPhaseScaTraPartitionedTwoWay::init(
    const Teuchos::ParameterList& globaltimeparams, const Teuchos::ParameterList& algoparams,
    const Teuchos::ParameterList& poroparams, const Teuchos::ParameterList& structparams,
    const Teuchos::ParameterList& fluidparams, const Teuchos::ParameterList& scatraparams,
    const std::string& struct_disname, const std::string& fluid_disname,
    const std::string& scatra_disname, bool isale, int nds_disp, int nds_vel, int nds_solidpressure,
    int ndsporofluid_scatra, const std::map<int, std::set<int>>* nearbyelepairs)
{
  // call base class
  PoroPressureBased::PoroMultiPhaseScaTraPartitioned::init(globaltimeparams, algoparams, poroparams,
      structparams, fluidparams, scatraparams, struct_disname, fluid_disname, scatra_disname, isale,
      nds_disp, nds_vel, nds_solidpressure, ndsporofluid_scatra, nearbyelepairs);

  // read input variables
  itmax_ = algoparams.get<int>("ITEMAX");
  ittol_ = algoparams.sublist("PARTITIONED").get<double>("CONVTOL");

  artery_coupling_active_ = algoparams.get<bool>("ARTERY_COUPLING");

  // initialize increment vectors
  scaincnp_ = std::make_shared<Core::LinAlg::Vector<double>>(
      *(scatra_algo()->scatra_field()->discretization()->dof_row_map()));
  structincnp_ =
      std::make_shared<Core::LinAlg::Vector<double>>(*(poro_field()->struct_dof_row_map()));
  fluidincnp_ =
      (std::make_shared<Core::LinAlg::Vector<double>>(*(poro_field()->fluid_dof_row_map())));
  if (artery_coupling_active_)
  {
    arterypressincnp_ = std::make_shared<Core::LinAlg::Vector<double>>(
        *(poro_field()->fluid_field()->artery_dof_row_map()));
    artscaincnp_ =
        std::make_shared<Core::LinAlg::Vector<double>>(*(scatramsht_->art_scatra_dof_row_map()));
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroPressureBased::PoroMultiPhaseScaTraPartitionedTwoWay::setup_system()
{
  poro_field()->setup_system();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroPressureBased::PoroMultiPhaseScaTraPartitionedTwoWay::setup_solver()
{
  poro_field()->setup_solver();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroPressureBased::PoroMultiPhaseScaTraPartitionedTwoWay::do_poro_step()
{
  // Newton-Raphson iteration
  poro_field()->time_step();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroPressureBased::PoroMultiPhaseScaTraPartitionedTwoWay::do_scatra_step()
{
  if (Core::Communication::my_mpi_rank(get_comm()) == 0)
  {
    std::cout << "\n";
    std::cout << "*********************************************************************************"
                 "********************************\n";
    std::cout << "TRANSPORT SOLVER   \n";
    std::cout << "*********************************************************************************"
                 "********************************\n";
  }

  // -------------------------------------------------------------------
  //                  solve nonlinear / linear equation
  // -------------------------------------------------------------------
  scatra_algo()->scatra_field()->solve();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroPressureBased::PoroMultiPhaseScaTraPartitionedTwoWay::print_header_partitioned()
{
  if (Core::Communication::my_mpi_rank(get_comm()) == 0)
  {
    std::cout << "\n";
    std::cout << "********************************************************************************"
              << "***************************************************************\n";
    std::cout << "* PARTITIONED OUTER ITERATION LOOP ----- MULTIPORO  <-------> SCATRA         "
              << "                                                                 *\n";
    std::cout << "* STEP: " << std::setw(5) << std::setprecision(4) << std::scientific << step()
              << "/" << std::setw(5) << std::setprecision(4) << std::scientific << n_step()
              << ", Time: " << std::setw(11) << std::setprecision(4) << std::scientific << time()
              << "/" << std::setw(11) << std::setprecision(4) << std::scientific << max_time()
              << ", Dt: " << std::setw(11) << std::setprecision(4) << std::scientific << dt()
              << "                                                                           *"
              << std::endl;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroPressureBased::PoroMultiPhaseScaTraPartitionedTwoWay::iter_update_states()
{
  // store scalar from first solution for convergence check (like in
  // elch_algorithm: use current values)
  scaincnp_->update(1.0, *scatra_algo()->scatra_field()->phinp(), 0.0);
  structincnp_->update(1.0, *poro_field()->struct_dispnp(), 0.0);
  fluidincnp_->update(1.0, *poro_field()->fluid_phinp(), 0.0);
  if (artery_coupling_active_)
  {
    arterypressincnp_->update(
        1.0, *(poro_field()->fluid_field()->art_net_tim_int()->pressurenp()), 0.0);
    artscaincnp_->update(1.0, *(scatramsht_->art_scatra_field()->phinp()), 0.0);
  }

}  // iter_update_states()

/*----------------------------------------------------------------------*
 | convergence check for both fields (scatra & poro) (copied from tsi)
 *----------------------------------------------------------------------*/
bool PoroPressureBased::PoroMultiPhaseScaTraPartitionedTwoWay::convergence_check(int itnum)
{
  // convergence check based on the scalar increment
  bool stopnonliniter = false;

  //    | scalar increment |_2
  //  -------------------------------- < Tolerance
  //     | scalar+1 |_2

  // variables to save different L2 - Norms
  // define L2-norm of incremental scalar and scalar
  double scaincnorm_L2(0.0);
  double scanorm_L2(0.0);
  double dispincnorm_L2(0.0);
  double dispnorm_L2(0.0);
  double fluidincnorm_L2(0.0);
  double fluidnorm_L2(0.0);
  double artpressincnorm_L2(0.0);
  double artpressnorm_L2(0.0);
  double artscaincnorm_L2(0.0);
  double artscanorm_L2(0.0);

  // build the current scalar increment Inc T^{i+1}
  // \f Delta T^{k+1} = Inc T^{k+1} = T^{k+1} - T^{k}  \f
  scaincnp_->update(1.0, *(scatra_algo()->scatra_field()->phinp()), -1.0);
  structincnp_->update(1.0, *(poro_field()->struct_dispnp()), -1.0);
  fluidincnp_->update(1.0, *(poro_field()->fluid_phinp()), -1.0);
  if (artery_coupling_active_)
  {
    arterypressincnp_->update(
        1.0, *(poro_field()->fluid_field()->art_net_tim_int()->pressurenp()), -1.0);
    artscaincnp_->update(1.0, *(scatramsht_->art_scatra_field()->phinp()), -1.0);
  }

  // build the L2-norm of the scalar increment and the scalar
  scaincnp_->norm_2(&scaincnorm_L2);
  scatra_algo()->scatra_field()->phinp()->norm_2(&scanorm_L2);
  structincnp_->norm_2(&dispincnorm_L2);
  poro_field()->struct_dispnp()->norm_2(&dispnorm_L2);
  fluidincnp_->norm_2(&fluidincnorm_L2);
  poro_field()->fluid_phinp()->norm_2(&fluidnorm_L2);
  if (artery_coupling_active_)
  {
    arterypressincnp_->norm_2(&artpressincnorm_L2);
    poro_field()->fluid_field()->art_net_tim_int()->pressurenp()->norm_2(&artpressnorm_L2);
    artscaincnp_->norm_2(&artscaincnorm_L2);
    poro_field()->fluid_field()->art_net_tim_int()->pressurenp()->norm_2(&artscanorm_L2);
  }

  // care for the case that there is (almost) zero scalar
  if (scanorm_L2 < 1e-6) scanorm_L2 = 1.0;
  if (dispnorm_L2 < 1e-6) dispnorm_L2 = 1.0;
  if (fluidnorm_L2 < 1e-6) fluidnorm_L2 = 1.0;
  if (artpressnorm_L2 < 1e-6) artpressnorm_L2 = 1.0;
  if (artscanorm_L2 < 1e-6) artscanorm_L2 = 1.0;

  // print the incremental based convergence check to the screen
  if (Core::Communication::my_mpi_rank(get_comm()) == 0)
  {
    std::cout << "                                                                                 "
                 "                                                             *\n";
    std::cout << "+--------------------------------------------------------------------------------"
                 "-----------------------------------------+                   *\n";
    std::cout << "| PARTITIONED OUTER ITERATION STEP ----- MULTIPORO  <-------> SCATRA             "
                 "                                         |                   *\n";
    printf(
        "+--------------+---------------------+----------------+----------------+-----"
        "-----------+----------------+----------------+                   *\n");
    printf(
        "|-  step/max  -|-  tol      [norm]  -|-- scalar-inc --|-- disp-inc   --|-- "
        "fluid-inc  --|--  1Dp-inc   --|--  1Ds-inc   --|                   *\n");
    printf(
        "|   %3d/%3d    |  %10.3E[L_2 ]   | %10.3E     | %10.3E     | %10.3E     | "
        "%10.3E     | %10.3E     |",
        itnum, itmax_, ittol_, scaincnorm_L2 / scanorm_L2, dispincnorm_L2 / dispnorm_L2,
        fluidincnorm_L2 / fluidnorm_L2, artpressincnorm_L2 / artpressnorm_L2,
        artscaincnorm_L2 / artscanorm_L2);
    printf("                   *\n");
    printf(
        "+--------------+---------------------+----------------+----------------+-----"
        "-----------+----------------+----------------+                   *\n");
  }

  // converged
  if ((scaincnorm_L2 / scanorm_L2 <= ittol_) and (dispincnorm_L2 / dispnorm_L2 <= ittol_) and
      (fluidincnorm_L2 / fluidnorm_L2 <= ittol_) and
      ((artpressincnorm_L2 / artpressnorm_L2) <= ittol_) and
      ((artscaincnorm_L2 / artscanorm_L2) <= ittol_))
  {
    stopnonliniter = true;
    if (Core::Communication::my_mpi_rank(get_comm()) == 0)
    {
      printf(
          "* MULTIPORO  <-------> SCATRA Outer Iteration loop converged after iteration %3d/%3d !  "
          "                                                      *\n",
          itnum, itmax_);
      printf(
          "****************************************************************************************"
          "*******************************************************\n");
    }
  }

  // break the loop
  // timestep
  if ((itnum == itmax_) and
      ((scaincnorm_L2 / scanorm_L2 > ittol_) or (dispincnorm_L2 / dispnorm_L2 > ittol_) or
          (fluidincnorm_L2 / fluidnorm_L2 > ittol_) or
          ((artpressincnorm_L2 / artpressnorm_L2) > ittol_) or
          ((artscaincnorm_L2 / artscanorm_L2) > ittol_)))
  {
    stopnonliniter = true;
    if ((Core::Communication::my_mpi_rank(get_comm()) == 0))
    {
      printf(
          "* MULTIPORO  <-------> SCATRA Outer Iteration loop not converged in itemax steps        "
          "                                                      *\n");
      printf(
          "****************************************************************************************"
          "*******************************************************\n");
      printf("\n");
      printf("\n");
    }
    handle_divergence();
  }

  return stopnonliniter;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
PoroPressureBased::PoroMultiPhaseScaTraPartitionedTwoWayNested::
    PoroMultiPhaseScaTraPartitionedTwoWayNested(
        MPI_Comm comm, const Teuchos::ParameterList& globaltimeparams)
    : PoroMultiPhaseScaTraPartitionedTwoWay(comm, globaltimeparams)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroPressureBased::PoroMultiPhaseScaTraPartitionedTwoWayNested::init(
    const Teuchos::ParameterList& globaltimeparams, const Teuchos::ParameterList& algoparams,
    const Teuchos::ParameterList& poroparams, const Teuchos::ParameterList& structparams,
    const Teuchos::ParameterList& fluidparams, const Teuchos::ParameterList& scatraparams,
    const std::string& struct_disname, const std::string& fluid_disname,
    const std::string& scatra_disname, bool isale, int nds_disp, int nds_vel, int nds_solidpressure,
    int ndsporofluid_scatra, const std::map<int, std::set<int>>* nearbyelepairs)
{
  // call base class
  PoroPressureBased::PoroMultiPhaseScaTraPartitionedTwoWay::init(globaltimeparams, algoparams,
      poroparams, structparams, fluidparams, scatraparams, struct_disname, fluid_disname,
      scatra_disname, isale, nds_disp, nds_vel, nds_solidpressure, ndsporofluid_scatra,
      nearbyelepairs);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroPressureBased::PoroMultiPhaseScaTraPartitionedTwoWayNested::solve()
{
  int itnum = 0;
  bool stopnonliniter = false;

  print_header_partitioned();

  while (stopnonliniter == false)
  {
    itnum++;

    // update the states to the last solutions obtained
    iter_update_states();

    // set structure-based scalar transport values
    set_scatra_solution();

    // solve structural system
    do_poro_step();

    // set mesh displacement and velocity fields
    set_poro_solution();

    // solve scalar transport equation
    do_scatra_step();

    // check convergence for all fields and stop iteration loop if
    // convergence is achieved overall
    stopnonliniter = convergence_check(itnum);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
PoroPressureBased::PoroMultiPhaseScaTraPartitionedTwoWaySequential::
    PoroMultiPhaseScaTraPartitionedTwoWaySequential(
        MPI_Comm comm, const Teuchos::ParameterList& globaltimeparams)
    : PoroMultiPhaseScaTraPartitionedTwoWay(comm, globaltimeparams)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroPressureBased::PoroMultiPhaseScaTraPartitionedTwoWaySequential::init(
    const Teuchos::ParameterList& globaltimeparams, const Teuchos::ParameterList& algoparams,
    const Teuchos::ParameterList& poroparams, const Teuchos::ParameterList& structparams,
    const Teuchos::ParameterList& fluidparams, const Teuchos::ParameterList& scatraparams,
    const std::string& struct_disname, const std::string& fluid_disname,
    const std::string& scatra_disname, bool isale, int nds_disp, int nds_vel, int nds_solidpressure,
    int ndsporofluid_scatra, const std::map<int, std::set<int>>* nearbyelepairs)
{
  // call base class
  PoroPressureBased::PoroMultiPhaseScaTraPartitionedTwoWay::init(globaltimeparams, algoparams,
      poroparams, structparams, fluidparams, scatraparams, struct_disname, fluid_disname,
      scatra_disname, isale, nds_disp, nds_vel, nds_solidpressure, ndsporofluid_scatra,
      nearbyelepairs);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroPressureBased::PoroMultiPhaseScaTraPartitionedTwoWaySequential::solve()
{
  int itnum = 0;
  bool stopnonliniter = false;

  print_header_partitioned();

  while (stopnonliniter == false)
  {
    itnum++;

    // update the states to the last solutions obtained
    iter_update_states();

    // 1) set scatra and structure solution (on fluid field)
    set_scatra_solution();
    poro_field()->set_struct_solution(
        poro_field()->structure_field()->dispnp(), poro_field()->structure_field()->velnp());

    // 2) solve fluid
    poro_field()->fluid_field()->solve();

    // 3) relaxation
    poro_field()->perform_relaxation(poro_field()->fluid_field()->phinp(), itnum);

    // 4) set relaxed fluid solution on structure field
    poro_field()->set_relaxed_fluid_solution();

    // 5) solve structure
    poro_field()->structure_field()->solve();

    // 6) set mesh displacement and velocity fields on ScaTra
    set_poro_solution();

    // 7) solve scalar transport equation
    do_scatra_step();

    // check convergence for all fields and stop iteration loop if
    // convergence is achieved overall
    stopnonliniter = convergence_check(itnum);
  }
}

FOUR_C_NAMESPACE_CLOSE
