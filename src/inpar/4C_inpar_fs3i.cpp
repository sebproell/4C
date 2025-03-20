// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_inpar_fs3i.hpp"

#include "4C_inpar_scatra.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN



void Inpar::FS3I::set_valid_parameters(std::map<std::string, Core::IO::InputSpec>& list)
{
  using Teuchos::tuple;
  using namespace Core::IO::InputSpecBuilders;

  Core::Utils::SectionSpecs fs3idyn{"FS3I DYNAMIC"};

  fs3idyn.specs.emplace_back(
      parameter<double>("TIMESTEP", {.description = "Time increment dt", .default_value = 0.1}));
  fs3idyn.specs.emplace_back(parameter<int>(
      "NUMSTEP", {.description = "Total number of time steps", .default_value = 20}));
  fs3idyn.specs.emplace_back(parameter<double>(
      "MAXTIME", {.description = "Total simulation time", .default_value = 1000.0}));
  fs3idyn.specs.emplace_back(parameter<int>(
      "RESULTSEVERY", {.description = "Increment for writing solution", .default_value = 1}));
  fs3idyn.specs.emplace_back(parameter<int>(
      "RESTARTEVERY", {.description = "Increment for writing restart", .default_value = 1}));
  fs3idyn.specs.emplace_back(deprecated_selection<Inpar::ScaTra::SolverType>("SCATRA_SOLVERTYPE",
      {
          {"linear", Inpar::ScaTra::solvertype_linear_incremental},
          {"nonlinear", Inpar::ScaTra::solvertype_nonlinear},
      },
      {.description = "type of scalar transport solver",
          .default_value = Inpar::ScaTra::solvertype_nonlinear}));
  fs3idyn.specs.emplace_back(parameter<bool>(
      "INF_PERM", {.description = "Flag for infinite permeability", .default_value = true}));
  std::vector<std::string> consthermpress_valid_input = {"No_energy", "No_mass", "Yes"};
  fs3idyn.specs.emplace_back(
      deprecated_selection<std::string>("CONSTHERMPRESS", consthermpress_valid_input,
          {.description = "treatment of thermodynamic pressure in time", .default_value = "Yes"}));

  // number of linear solver used for fs3i problems
  fs3idyn.specs.emplace_back(parameter<int>("COUPLED_LINEAR_SOLVER",
      {.description = "number of linear solver used for fs3i problem", .default_value = -1}));
  fs3idyn.specs.emplace_back(parameter<int>("LINEAR_SOLVER1",
      {.description = "number of linear solver used for fluid problem", .default_value = -1}));
  fs3idyn.specs.emplace_back(parameter<int>("LINEAR_SOLVER2",
      {.description = "number of linear solver used for structural problem", .default_value = -1}));

  fs3idyn.specs.emplace_back(deprecated_selection<Inpar::ScaTra::ConvForm>("STRUCTSCAL_CONVFORM",
      {
          {"convective", Inpar::ScaTra::convform_convective},
          {"conservative", Inpar::ScaTra::convform_conservative},
      },
      {.description = "form of convective term of structure scalar",
          .default_value = Inpar::ScaTra::convform_conservative}));

  fs3idyn.specs.emplace_back(
      deprecated_selection<Inpar::ScaTra::InitialField>("STRUCTSCAL_INITIALFIELD",
          {
              {"zero_field", Inpar::ScaTra::initfield_zero_field},
              {"field_by_function", Inpar::ScaTra::initfield_field_by_function},
          },
          {.description = "Initial Field for structure scalar transport problem",
              .default_value = Inpar::ScaTra::initfield_zero_field}));

  fs3idyn.specs.emplace_back(parameter<int>("STRUCTSCAL_INITFUNCNO",
      {.description = "function number for structure scalar transport initial field",
          .default_value = -1}));

  // Type of coupling strategy between structure and structure-scalar field
  fs3idyn.specs.emplace_back(deprecated_selection<VolumeCoupling>("STRUCTSCAL_FIELDCOUPLING",
      {
          {"volume_matching", coupling_match},
          {"volume_nonmatching", coupling_nonmatch},
      },
      {.description = "Type of coupling strategy between structure and structure-scalar field",
          .default_value = coupling_match}));

  // Type of coupling strategy between fluid and fluid-scalar field
  fs3idyn.specs.emplace_back(deprecated_selection<VolumeCoupling>("FLUIDSCAL_FIELDCOUPLING",
      {
          {"volume_matching", coupling_match},
          {"volume_nonmatching", coupling_nonmatch},
      },
      {.description = "Type of coupling strategy between fluid and fluid-scalar field",
          .default_value = coupling_match}));

  // type of scalar transport
  fs3idyn.specs.emplace_back(deprecated_selection<Inpar::ScaTra::ImplType>("FLUIDSCAL_SCATRATYPE",
      {
          {"Undefined", Inpar::ScaTra::impltype_undefined},
          {"ConvectionDiffusion", Inpar::ScaTra::impltype_std},
          {"Loma", Inpar::ScaTra::impltype_loma},
          {"Advanced_Reaction", Inpar::ScaTra::impltype_advreac},
          {"Chemotaxis", Inpar::ScaTra::impltype_chemo},
          {"Chemo_Reac", Inpar::ScaTra::impltype_chemoreac},
      },
      {.description = "Type of scalar transport problem",
          .default_value = Inpar::ScaTra::impltype_std}));

  // Restart from FSI instead of FS3I
  fs3idyn.specs.emplace_back(parameter<bool>(
      "RESTART_FROM_PART_FSI", {.description = "restart from partitioned fsi problem (e.g. from "
                                               "prestress calculations) instead of fs3i",
                                   .default_value = false}));

  fs3idyn.move_into_collection(list);

  /*----------------------------------------------------------------------*/
  /* parameters for partitioned FS3I */
  /*----------------------------------------------------------------------*/
  Core::Utils::SectionSpecs fs3idynpart{fs3idyn, "PARTITIONED"};

  // Coupling strategy for partitioned FS3I
  fs3idynpart.specs.emplace_back(parameter<SolutionSchemeOverFields>("COUPALGO",
      {.description = "Coupling strategies for FS3I solvers", .default_value = fs3i_IterStagg}));

  // convergence tolerance of outer iteration loop
  fs3idynpart.specs.emplace_back(parameter<double>("CONVTOL",
      {.description = "tolerance for convergence check of outer iteration within partitioned FS3I",
          .default_value = 1e-6}));

  fs3idynpart.specs.emplace_back(parameter<int>(
      "ITEMAX", {.description = "Maximum number of outer iterations", .default_value = 10}));

  fs3idynpart.move_into_collection(list);

  /*----------------------------------------------------------------------  */
  /* parameters for stabilization of the structure-scalar field             */
  /*----------------------------------------------------------------------  */

  /// HACK!
  /// reuse the parameters from scatra
  FOUR_C_ASSERT(list.contains("SCALAR TRANSPORT DYNAMIC/STABILIZATION"),
      "Internal error: Hacky FS3I requires scatra to already have been read.");

  list["FS3I DYNAMIC/STRUCTURE SCALAR STABILIZATION"] =
      list["SCALAR TRANSPORT DYNAMIC/STABILIZATION"];
}

FOUR_C_NAMESPACE_CLOSE
