// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_porofluid_pressure_based_elast_scatra_input.hpp"

#include "4C_fem_condition_definition.hpp"
#include "4C_io_input_spec_builders.hpp"
#include "4C_linalg_equilibrate.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN


void PoroMultiPhaseScaTra::set_valid_parameters(std::map<std::string, Core::IO::InputSpec>& list)
{
  using namespace Core::IO::InputSpecBuilders;

  // ----------------------------------------------------------------------
  // (1) general control parameters
  list["POROMULTIPHASESCATRA DYNAMIC"] = all_of({

      // Output type
      parameter<int>(
          "RESTARTEVERY", {.description = "write restart possibility every RESTARTEVERY steps",
                              .default_value = 1}),
      // Time loop control
      parameter<int>(
          "NUMSTEP", {.description = "maximum number of Timesteps", .default_value = 200}),
      parameter<double>(
          "MAXTIME", {.description = "total simulation time", .default_value = 1000.0}),

      parameter<double>("TIMESTEP", {.description = "time step size dt", .default_value = 0.05}),
      parameter<int>(
          "RESULTSEVERY", {.description = "increment for writing solution", .default_value = 1}),
      parameter<int>("ITEMAX",
          {.description = "maximum number of iterations over fields", .default_value = 10}),
      parameter<int>("ITEMIN",
          {.description = "minimal number of iterations over fields", .default_value = 1}),

      // Coupling strategy for poroscatra solvers
      deprecated_selection<SolutionSchemeOverFields>("COUPALGO",
          {
              {"twoway_partitioned_nested", solscheme_twoway_partitioned_nested},
              {"twoway_partitioned_sequential", solscheme_twoway_partitioned_sequential},
              {"twoway_monolithic", solscheme_twoway_monolithic},
          },
          {.description = "Coupling strategies for poroscatra solvers",
              .default_value = solscheme_twoway_partitioned_nested}),

      // coupling with 1D artery network active
      parameter<bool>("ARTERY_COUPLING",
          {.description = "Coupling with 1D blood vessels.", .default_value = false}),

      // no convergence of coupling scheme
      deprecated_selection<DivContAct>("DIVERCONT",
          {
              {"stop", divcont_stop},
              {"continue", divcont_continue},
          },
          {.description =
                  "What to do with time integration when Poromultiphase-Scatra iteration failed",
              .default_value =
                  divcont_stop})});  // ----------------------------------------------------------------------
  // (2) monolithic parameters
  list["POROMULTIPHASESCATRA DYNAMIC/MONOLITHIC"] = all_of({

      deprecated_selection<VectorNorm>("VECTORNORM_RESF",
          {
              {"L1", PoroMultiPhaseScaTra::norm_l1},
              {"L1_Scaled", PoroMultiPhaseScaTra::norm_l1_scaled},
              {"L2", PoroMultiPhaseScaTra::norm_l2},
              {"Rms", PoroMultiPhaseScaTra::norm_rms},
              {"Inf", PoroMultiPhaseScaTra::norm_inf},
          },
          {.description = "type of norm to be applied to residuals",
              .default_value = PoroMultiPhaseScaTra::norm_l2}),

      deprecated_selection<VectorNorm>("VECTORNORM_INC",
          {
              {"L1", PoroMultiPhaseScaTra::norm_l1},
              {"L1_Scaled", PoroMultiPhaseScaTra::norm_l1_scaled},
              {"L2", PoroMultiPhaseScaTra::norm_l2},
              {"Rms", PoroMultiPhaseScaTra::norm_rms},
              {"Inf", PoroMultiPhaseScaTra::norm_inf},
          },
          {.description = "type of norm to be applied to residuals",
              .default_value = PoroMultiPhaseScaTra::norm_l2}),

      // convergence criteria adaptivity --> note ADAPTCONV_BETTER set pretty small
      parameter<bool>("ADAPTCONV",
          {.description =
                  "Switch on adaptive control of linear solver tolerance for nonlinear solution",
              .default_value = false}),
      parameter<double>("ADAPTCONV_BETTER",
          {.description = "The linear solver shall be this much better than the current nonlinear "
                          "residual in the nonlinear convergence limit",
              .default_value = 0.001}),

      // Iterationparameters
      parameter<double>("TOLRES_GLOBAL",
          {.description = "tolerance in the residual norm for the Newton iteration",
              .default_value = 1e-8}),
      parameter<double>("TOLINC_GLOBAL",
          {.description = "tolerance in the increment norm for the Newton iteration",
              .default_value = 1e-8}),

      // number of linear solver used for poroelasticity
      parameter<int>("LINEAR_SOLVER",
          {.description = "number of linear solver used for monolithic poroscatra problems",
              .default_value = -1}),

      // parameters for finite difference check
      deprecated_selection<FdCheck>("FDCHECK",
          {
              {"none", fdcheck_none},
              {"global", fdcheck_global},
          },
          {.description = "flag for finite difference check: none or global",
              .default_value = fdcheck_none}),

      // flag for equilibration of global system of equations
      parameter<Core::LinAlg::EquilibrationMethod>(
          "EQUILIBRATION", {.description = "flag for equilibration of global system of equations",
                               .default_value = Core::LinAlg::EquilibrationMethod::none})});

  // ----------------------------------------------------------------------
  // (3) partitioned parameters
  list["POROMULTIPHASESCATRA DYNAMIC/PARTITIONED"] = all_of({

      // convergence tolerance of outer iteration loop
      parameter<double>(
          "CONVTOL", {.description = "tolerance for convergence check of outer iteration",
                         .default_value = 1e-6})});
}

void PoroMultiPhaseScaTra::set_valid_conditions(
    std::vector<Core::Conditions::ConditionDefinition>& condlist)
{
  using namespace Core::IO::InputSpecBuilders;

  /*--------------------------------------------------------------------*/
  // oxygen partial pressure calculation condition
  {
    // definition of oxygen partial pressure calculation condition
    Core::Conditions::ConditionDefinition oxypartpressline(
        "DESIGN OXYGEN PARTIAL PRESSURE CALCULATION LINE CONDITIONS",
        "PoroMultiphaseScatraOxyPartPressCalcCond",
        "PoroMultiphaseScatra Oxygen Partial Pressure Calculation line condition",
        Core::Conditions::PoroMultiphaseScatraOxyPartPressCalcCond, true,
        Core::Conditions::geometry_type_line);
    Core::Conditions::ConditionDefinition oxypartpresssurf(
        "DESIGN OXYGEN PARTIAL PRESSURE CALCULATION SURF CONDITIONS",
        "PoroMultiphaseScatraOxyPartPressCalcCond",
        "PoroMultiphaseScatra Oxygen Partial Pressure Calculation surface condition",
        Core::Conditions::PoroMultiphaseScatraOxyPartPressCalcCond, true,
        Core::Conditions::geometry_type_surface);
    Core::Conditions::ConditionDefinition oxypartpressvol(
        "DESIGN OXYGEN PARTIAL PRESSURE CALCULATION VOL CONDITIONS",
        "PoroMultiphaseScatraOxyPartPressCalcCond",
        "PoroMultiphaseScatra Oxygen Partial Pressure Calculation volume condition",
        Core::Conditions::PoroMultiphaseScatraOxyPartPressCalcCond, true,
        Core::Conditions::geometry_type_volume);

    const auto make_oxypartpress = [&condlist](Core::Conditions::ConditionDefinition& cond)
    {
      cond.add_component(
          parameter<int>("SCALARID", {.description = "scalar id of oxygen partial pressure"}));
      cond.add_component(parameter<double>("n", {.description = "n"}));
      cond.add_component(parameter<double>("Pb50", {.description = "Pb50"}));
      cond.add_component(parameter<double>("CaO2_max", {.description = "CaO2_max"}));
      cond.add_component(parameter<double>("alpha_bl_eff", {.description = "alpha_bl_eff"}));
      cond.add_component(parameter<double>("rho_oxy", {.description = "rho_oxy"}));
      cond.add_component(parameter<double>("rho_bl", {.description = "rho_bl"}));

      condlist.push_back(cond);
    };

    make_oxypartpress(oxypartpressline);
    make_oxypartpress(oxypartpresssurf);
    make_oxypartpress(oxypartpressvol);
  }
}

FOUR_C_NAMESPACE_CLOSE