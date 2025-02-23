// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_inpar_poromultiphase.hpp"

#include "4C_fem_condition_definition.hpp"
#include "4C_linalg_equilibrate.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN


void Inpar::POROMULTIPHASE::set_valid_parameters(std::map<std::string, Core::IO::InputSpec>& list)
{
  using Teuchos::tuple;

  // ----------------------------------------------------------------------
  // (1) general control parameters
  Core::Utils::SectionSpecs poromultiphasedyn{"POROMULTIPHASE DYNAMIC"};

  // Output type
  Core::Utils::int_parameter(
      "RESTARTEVERY", 1, "write restart possibility every RESTARTEVERY steps", poromultiphasedyn);
  // Time loop control
  Core::Utils::int_parameter("NUMSTEP", 200, "maximum number of Timesteps", poromultiphasedyn);
  Core::Utils::double_parameter("MAXTIME", 1000.0, "total simulation time", poromultiphasedyn);
  Core::Utils::double_parameter("TIMESTEP", -1, "time step size dt", poromultiphasedyn);
  Core::Utils::int_parameter(
      "RESULTSEVERY", 1, "increment for writing solution", poromultiphasedyn);
  Core::Utils::int_parameter(
      "ITEMAX", 10, "maximum number of iterations over fields", poromultiphasedyn);

  // here the computation of the structure can be skipped, this is helpful if only fluid-scatra
  // coupling should be calculated
  Core::Utils::bool_parameter(
      "SOLVE_STRUCTURE", true, "Flag to skip computation of structural field", poromultiphasedyn);


  // Coupling strategy for solvers
  Core::Utils::string_to_integral_parameter<SolutionSchemeOverFields>("COUPALGO",
      "twoway_partitioned", "Coupling strategies for poro multiphase solvers",
      tuple<std::string>("twoway_partitioned", "twoway_monolithic"),
      tuple<SolutionSchemeOverFields>(solscheme_twoway_partitioned, solscheme_twoway_monolithic),
      poromultiphasedyn);

  // coupling with 1D artery network active
  Core::Utils::bool_parameter(
      "ARTERY_COUPLING", false, "Coupling with 1D blood vessels.", poromultiphasedyn);

  poromultiphasedyn.move_into_collection(list);

  // ----------------------------------------------------------------------
  // (2) monolithic parameters
  Core::Utils::SectionSpecs poromultiphasedynmono{poromultiphasedyn, "MONOLITHIC"};

  // convergence tolerances for monolithic coupling
  Core::Utils::double_parameter("TOLRES_GLOBAL", 1e-8,
      "tolerance in the residual norm for the Newton iteration", poromultiphasedynmono);
  Core::Utils::double_parameter("TOLINC_GLOBAL", 1e-8,
      "tolerance in the increment norm for the Newton iteration", poromultiphasedynmono);

  // number of linear solver used for poroelasticity
  Core::Utils::int_parameter("LINEAR_SOLVER", -1,
      "number of linear solver used for poroelasticity problems", poromultiphasedynmono);

  // parameters for finite difference check
  Core::Utils::string_to_integral_parameter<FdCheck>("FDCHECK", "none",
      "flag for finite difference check: none or global",
      tuple<std::string>("none",
          "global"),  // perform finite difference check on time integrator level
      tuple<FdCheck>(fdcheck_none, fdcheck_global), poromultiphasedynmono);

  Core::Utils::string_to_integral_parameter<VectorNorm>("VECTORNORM_RESF", "L2",
      "type of norm to be applied to residuals",
      tuple<std::string>("L1", "L1_Scaled", "L2", "Rms", "Inf"),
      tuple<VectorNorm>(Inpar::POROMULTIPHASE::norm_l1, Inpar::POROMULTIPHASE::norm_l1_scaled,
          Inpar::POROMULTIPHASE::norm_l2, Inpar::POROMULTIPHASE::norm_rms,
          Inpar::POROMULTIPHASE::norm_inf),
      poromultiphasedynmono);

  Core::Utils::string_to_integral_parameter<VectorNorm>("VECTORNORM_INC", "L2",
      "type of norm to be applied to residuals",
      tuple<std::string>("L1", "L1_Scaled", "L2", "Rms", "Inf"),
      tuple<VectorNorm>(Inpar::POROMULTIPHASE::norm_l1, Inpar::POROMULTIPHASE::norm_l1_scaled,
          Inpar::POROMULTIPHASE::norm_l2, Inpar::POROMULTIPHASE::norm_rms,
          Inpar::POROMULTIPHASE::norm_inf),
      poromultiphasedynmono);

  // flag for equilibration of global system of equations
  Core::Utils::string_to_integral_parameter<Core::LinAlg::EquilibrationMethod>("EQUILIBRATION",
      "none", "flag for equilibration of global system of equations",
      tuple<std::string>("none", "rows_full", "rows_maindiag", "columns_full", "columns_maindiag",
          "rowsandcolumns_full", "rowsandcolumns_maindiag"),
      tuple<Core::LinAlg::EquilibrationMethod>(Core::LinAlg::EquilibrationMethod::none,
          Core::LinAlg::EquilibrationMethod::rows_full,
          Core::LinAlg::EquilibrationMethod::rows_maindiag,
          Core::LinAlg::EquilibrationMethod::columns_full,
          Core::LinAlg::EquilibrationMethod::columns_maindiag,
          Core::LinAlg::EquilibrationMethod::rowsandcolumns_full,
          Core::LinAlg::EquilibrationMethod::rowsandcolumns_maindiag),
      poromultiphasedynmono);

  // convergence criteria adaptivity --> note ADAPTCONV_BETTER set pretty small
  Core::Utils::bool_parameter("ADAPTCONV", false,
      "Switch on adaptive control of linear solver tolerance for nonlinear solution",
      poromultiphasedynmono);
  Core::Utils::double_parameter("ADAPTCONV_BETTER", 0.001,
      "The linear solver shall be this much better "
      "than the current nonlinear residual in the nonlinear convergence limit",
      poromultiphasedynmono);

  poromultiphasedynmono.move_into_collection(list);

  // ----------------------------------------------------------------------
  // (3) partitioned parameters
  Core::Utils::SectionSpecs poromultiphasedynpart{poromultiphasedyn, "PARTITIONED"};

  // convergence tolerance of outer iteration loop
  Core::Utils::double_parameter(
      "CONVTOL", 1e-6, "tolerance for convergence check of outer iteration", poromultiphasedynpart);

  // flag for relaxation of partitioned scheme
  Core::Utils::string_to_integral_parameter<RelaxationMethods>("RELAXATION", "none",
      "flag for relaxation of partitioned scheme", tuple<std::string>("none", "Constant", "Aitken"),
      tuple<RelaxationMethods>(relaxation_none, relaxation_constant, relaxation_aitken),
      poromultiphasedynpart);

  // parameters for relaxation of partitioned coupling
  Core::Utils::double_parameter(
      "STARTOMEGA", 1.0, "fixed relaxation parameter", poromultiphasedynpart);
  Core::Utils::double_parameter(
      "MINOMEGA", 0.1, "smallest omega allowed for Aitken relaxation", poromultiphasedynpart);
  Core::Utils::double_parameter(
      "MAXOMEGA", 10.0, "largest omega allowed for Aitken relaxation", poromultiphasedynpart);

  poromultiphasedynpart.move_into_collection(list);
}

FOUR_C_NAMESPACE_CLOSE
