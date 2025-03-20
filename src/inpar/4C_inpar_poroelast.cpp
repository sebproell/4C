// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_inpar_poroelast.hpp"

#include "4C_inpar_fluid.hpp"
#include "4C_linalg_equilibrate.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN


void Inpar::PoroElast::set_valid_parameters(std::map<std::string, Core::IO::InputSpec>& list)
{
  using Teuchos::tuple;
  using namespace Core::IO::InputSpecBuilders;

  Core::Utils::SectionSpecs poroelastdyn{"POROELASTICITY DYNAMIC"};

  // Coupling strategy for (monolithic) porous media solvers
  poroelastdyn.specs.emplace_back(
      deprecated_selection<Inpar::PoroElast::SolutionSchemeOverFields>("COUPALGO",
          {
              {"poro_partitioned", Partitioned},
              {"poro_monolithic", Monolithic},
              {"poro_monolithicstructuresplit", Monolithic_structuresplit},
              {"poro_monolithicfluidsplit", Monolithic_fluidsplit},
              {"poro_monolithicnopenetrationsplit", Monolithic_nopenetrationsplit},
              {"poro_monolithicmeshtying", Monolithic_meshtying},
          },
          {.description = "Coupling strategies for poroelasticity solvers",
              .default_value = Monolithic}));

  // physical type of poro fluid flow (incompressible, varying density, loma, Boussinesq
  // approximation)
  poroelastdyn.specs.emplace_back(deprecated_selection<Inpar::FLUID::PhysicalType>("PHYSICAL_TYPE",
      {
          {"Poro", Inpar::FLUID::poro},
          {"Poro_P1", Inpar::FLUID::poro_p1},
      },
      {.description = "Physical Type of Porofluid", .default_value = Inpar::FLUID::poro}));

  // physical type of poro fluid flow (incompressible, varying density, loma, Boussinesq
  // approximation)
  poroelastdyn.specs.emplace_back(
      deprecated_selection<Inpar::PoroElast::TransientEquationsOfPoroFluid>("TRANSIENT_TERMS",
          {
              {"none", transient_none},
              {"momentum", transient_momentum_only},
              {"continuity", transient_continuity_only},
              {"all", transient_all},
          },
          {.description = "which equation includes transient terms",
              .default_value = transient_all}));

  // Output type
  poroelastdyn.specs.emplace_back(parameter<int>("RESTARTEVERY",
      {.description = "write restart possibility every RESTARTEVERY steps", .default_value = 1}));

  // Time loop control
  poroelastdyn.specs.emplace_back(parameter<int>(
      "NUMSTEP", {.description = "maximum number of Timesteps", .default_value = 200}));
  poroelastdyn.specs.emplace_back(parameter<double>(
      "MAXTIME", {.description = "total simulation time", .default_value = 1000.0}));
  poroelastdyn.specs.emplace_back(
      parameter<double>("TIMESTEP", {.description = "time step size dt", .default_value = 0.05}));
  poroelastdyn.specs.emplace_back(parameter<int>(
      "ITEMAX", {.description = "maximum number of iterations over fields", .default_value = 10}));
  poroelastdyn.specs.emplace_back(parameter<int>(
      "ITEMIN", {.description = "minimal number of iterations over fields", .default_value = 1}));
  poroelastdyn.specs.emplace_back(parameter<int>(
      "RESULTSEVERY", {.description = "increment for writing solution", .default_value = 1}));

  // Iterationparameters
  poroelastdyn.specs.emplace_back(parameter<double>(
      "TOLRES_GLOBAL", {.description = "tolerance in the residual norm for the Newton iteration",
                           .default_value = 1e-8}));
  poroelastdyn.specs.emplace_back(parameter<double>(
      "TOLINC_GLOBAL", {.description = "tolerance in the increment norm for the Newton iteration",
                           .default_value = 1e-8}));
  poroelastdyn.specs.emplace_back(parameter<double>(
      "TOLRES_DISP", {.description = "tolerance in the residual norm for the Newton iteration",
                         .default_value = 1e-8}));
  poroelastdyn.specs.emplace_back(parameter<double>(
      "TOLINC_DISP", {.description = "tolerance in the increment norm for the Newton iteration",
                         .default_value = 1e-8}));
  poroelastdyn.specs.emplace_back(parameter<double>(
      "TOLRES_PORO", {.description = "tolerance in the residual norm for the Newton iteration",
                         .default_value = 1e-8}));
  poroelastdyn.specs.emplace_back(parameter<double>(
      "TOLINC_PORO", {.description = "tolerance in the increment norm for the Newton iteration",
                         .default_value = 1e-8}));
  poroelastdyn.specs.emplace_back(parameter<double>(
      "TOLRES_VEL", {.description = "tolerance in the residual norm for the Newton iteration",
                        .default_value = 1e-8}));
  poroelastdyn.specs.emplace_back(parameter<double>(
      "TOLINC_VEL", {.description = "tolerance in the increment norm for the Newton iteration",
                        .default_value = 1e-8}));
  poroelastdyn.specs.emplace_back(parameter<double>(
      "TOLRES_PRES", {.description = "tolerance in the residual norm for the Newton iteration",
                         .default_value = 1e-8}));
  poroelastdyn.specs.emplace_back(parameter<double>(
      "TOLINC_PRES", {.description = "tolerance in the increment norm for the Newton iteration",
                         .default_value = 1e-8}));
  poroelastdyn.specs.emplace_back(parameter<double>(
      "TOLRES_NCOUP", {.description = "tolerance in the residual norm for the Newton iteration",
                          .default_value = 1e-8}));

  poroelastdyn.specs.emplace_back(deprecated_selection<Inpar::PoroElast::ConvNorm>("NORM_INC",
      {
          {"AbsGlobal", convnorm_abs_global},
          {"AbsSingleFields", convnorm_abs_singlefields},
      },
      {.description = "type of norm for primary variables convergence check",
          .default_value = convnorm_abs_singlefields}));

  poroelastdyn.specs.emplace_back(deprecated_selection<Inpar::PoroElast::ConvNorm>("NORM_RESF",
      {
          {"AbsGlobal", convnorm_abs_global},
          {"AbsSingleFields", convnorm_abs_singlefields},
      },
      {.description = "type of norm for residual convergence check",
          .default_value = convnorm_abs_singlefields}));

  poroelastdyn.specs.emplace_back(
      deprecated_selection<Inpar::PoroElast::BinaryOp>("NORMCOMBI_RESFINC",
          {
              {"And", bop_and},
              {"Or", bop_or},
          },
          {.description = "binary operator to combine primary variables and residual force values",
              .default_value = bop_and}));

  poroelastdyn.specs.emplace_back(
      deprecated_selection<Inpar::PoroElast::VectorNorm>("VECTORNORM_RESF",
          {
              {"L1", norm_l1},
              {"L1_Scaled", norm_l1_scaled},
              {"L2", norm_l2},
              {"Rms", norm_rms},
              {"Inf", norm_inf},
          },
          {.description = "type of norm to be applied to residuals", .default_value = norm_l2}));

  poroelastdyn.specs.emplace_back(
      deprecated_selection<Inpar::PoroElast::VectorNorm>("VECTORNORM_INC",
          {
              {"L1", norm_l1},
              {"L1_Scaled", norm_l1_scaled},
              {"L2", norm_l2},
              {"Rms", norm_rms},
              {"Inf", norm_inf},
          },
          {.description = "type of norm to be applied to residuals", .default_value = norm_l2}));

  poroelastdyn.specs.emplace_back(parameter<bool>("SECONDORDER",
      {.description = "Second order coupling at the interface.", .default_value = true}));

  poroelastdyn.specs.emplace_back(parameter<bool>("CONTIPARTINT",
      {.description = "Partial integration of porosity gradient in continuity equation",
          .default_value = false}));

  poroelastdyn.specs.emplace_back(parameter<bool>("CONTACT_NO_PENETRATION",
      {.description =
              "No-Penetration Condition on active contact surface in case of poro contact problem!",
          .default_value = false}));

  poroelastdyn.specs.emplace_back(
      parameter<bool>("MATCHINGGRID", {.description = "is matching grid", .default_value = true}));

  poroelastdyn.specs.emplace_back(parameter<bool>(
      "CONVECTIVE_TERM", {.description = "convective term ", .default_value = false}));

  // number of linear solver used for poroelasticity
  poroelastdyn.specs.emplace_back(parameter<int>(
      "LINEAR_SOLVER", {.description = "number of linear solver used for poroelasticity problems",
                           .default_value = -1}));

  // flag for equilibration of global system of equations
  poroelastdyn.specs.emplace_back(parameter<Core::LinAlg::EquilibrationMethod>(
      "EQUILIBRATION", {.description = "flag for equilibration of global system of equations",
                           .default_value = Core::LinAlg::EquilibrationMethod::none}));

  poroelastdyn.move_into_collection(list);
}

FOUR_C_NAMESPACE_CLOSE
