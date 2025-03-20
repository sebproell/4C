// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_inpar_tsi.hpp"

#include "4C_contact_input.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN



void Inpar::TSI::set_valid_parameters(std::map<std::string, Core::IO::InputSpec>& list)
{
  using Teuchos::tuple;
  using namespace Core::IO::InputSpecBuilders;

  Core::Utils::SectionSpecs tsidyn{"TSI DYNAMIC"};

  // coupling strategy for (partitioned and monolithic) TSI solvers
  tsidyn.specs.emplace_back(deprecated_selection<SolutionSchemeOverFields>("COUPALGO",
      {
          {"tsi_oneway", OneWay},
          {"tsi_sequstagg", SequStagg},
          {"tsi_iterstagg", IterStagg},
          {"tsi_iterstagg_aitken", IterStaggAitken},
          {"tsi_iterstagg_aitkenirons", IterStaggAitkenIrons},
          {"tsi_iterstagg_fixedrelax", IterStaggFixedRel},
          {"tsi_monolithic", Monolithic},
      },
      {.description = "Coupling strategies for TSI solvers", .default_value = Monolithic}));

  tsidyn.specs.emplace_back(
      parameter<bool>("MATCHINGGRID", {.description = "is matching grid", .default_value = true}));

  // output type
  tsidyn.specs.emplace_back(parameter<int>("RESTARTEVERY",
      {.description = "write restart possibility every RESTARTEVERY steps", .default_value = 1}));

  // time loop control
  tsidyn.specs.emplace_back(parameter<int>(
      "NUMSTEP", {.description = "maximum number of Timesteps", .default_value = 200}));
  tsidyn.specs.emplace_back(parameter<double>(
      "MAXTIME", {.description = "total simulation time", .default_value = 1000.0}));
  tsidyn.specs.emplace_back(
      parameter<double>("TIMESTEP", {.description = "time step size dt", .default_value = 0.05}));
  tsidyn.specs.emplace_back(parameter<int>(
      "ITEMAX", {.description = "maximum number of iterations over fields", .default_value = 10}));
  tsidyn.specs.emplace_back(parameter<int>(
      "ITEMIN", {.description = "minimal number of iterations over fields", .default_value = 1}));
  tsidyn.specs.emplace_back(parameter<int>(
      "RESULTSEVERY", {.description = "increment for writing solution", .default_value = 1}));

  tsidyn.specs.emplace_back(deprecated_selection<ConvNorm>("NORM_INC",
      {
          {"Abs", convnorm_abs},
          {"Rel", convnorm_rel},
          {"Mix", convnorm_mix},
      },
      {.description = "type of norm for convergence check of primary variables in TSI",
          .default_value = convnorm_abs}));

  tsidyn.move_into_collection(list);

  /*----------------------------------------------------------------------*/
  /* parameters for monolithic TSI */
  Core::Utils::SectionSpecs tsidynmono{tsidyn, "MONOLITHIC"};

  // convergence tolerance of tsi residual
  tsidynmono.specs.emplace_back(parameter<double>(
      "CONVTOL", {.description = "tolerance for convergence check of TSI", .default_value = 1e-6}));
  // Iterationparameters
  tsidynmono.specs.emplace_back(parameter<double>("TOLINC",
      {.description = "tolerance for convergence check of TSI-increment in monolithic TSI",
          .default_value = 1.0e-6}));

  tsidynmono.specs.emplace_back(deprecated_selection<ConvNorm>("NORM_RESF",
      {
          {"Abs", convnorm_abs},
          {"Rel", convnorm_rel},
          {"Mix", convnorm_mix},
      },
      {.description = "type of norm for residual convergence check",
          .default_value = convnorm_abs}));

  tsidynmono.specs.emplace_back(deprecated_selection<BinaryOp>("NORMCOMBI_RESFINC",
      {
          {"And", bop_and},
          {"Or", bop_or},
          {"Coupl_Or_Single", bop_coupl_or_single},
          {"Coupl_And_Single", bop_coupl_and_single},
          {"And_Single", bop_and_single},
          {"Or_Single", bop_or_single},
      },
      {.description = "binary operator to combine primary variables and residual force values",
          .default_value = bop_coupl_and_single}));

  tsidynmono.specs.emplace_back(deprecated_selection<VectorNorm>("ITERNORM",
      {
          {"L1", norm_l1},
          {"L1_Scaled", norm_l1_scaled},
          {"L2", norm_l2},
          {"Rms", norm_rms},
          {"Inf", norm_inf},
      },
      {.description = "type of norm to be applied to residuals", .default_value = norm_rms}));

  tsidynmono.specs.emplace_back(deprecated_selection<NlnSolTech>("NLNSOL",
      {
          {"fullnewton", soltech_newtonfull},
          {"ptc", soltech_ptc},
      },
      {.description = "Nonlinear solution technique", .default_value = soltech_newtonfull}));

  tsidynmono.specs.emplace_back(
      parameter<double>("PTCDT", {.description = "pseudo time step for pseudo-transient "
                                                 "continuation (PTC) stabilised Newton procedure",
                                     .default_value = 0.1}));

  // number of linear solver used for monolithic TSI
  tsidynmono.specs.emplace_back(parameter<int>(
      "LINEAR_SOLVER", {.description = "number of linear solver used for monolithic TSI problems",
                           .default_value = -1}));

  // convergence criteria adaptivity of monolithic TSI solver
  tsidynmono.specs.emplace_back(parameter<bool>("ADAPTCONV",
      {.description =
              "Switch on adaptive control of linear solver tolerance for nonlinear solution",
          .default_value = false}));
  tsidynmono.specs.emplace_back(parameter<double>("ADAPTCONV_BETTER",
      {.description = "The linear solver shall be this much better than the current nonlinear "
                      "residual in the nonlinear convergence limit",
          .default_value = 0.1}));

  tsidynmono.specs.emplace_back(parameter<bool>("INFNORMSCALING",
      {.description = "Scale blocks of matrix with row infnorm?", .default_value = true}));

  // merge TSI block matrix to enable use of direct solver in monolithic TSI
  // default: "No", i.e. use block matrix
  tsidynmono.specs.emplace_back(parameter<bool>(
      "MERGE_TSI_BLOCK_MATRIX", {.description = "Merge TSI block matrix", .default_value = false}));

  // in case of monolithic TSI nodal values (displacements, temperatures and
  // reaction forces) at fix points of the body can be calculated
  // default: "No", i.e. nothing is calculated
  tsidynmono.specs.emplace_back(parameter<bool>("CALC_NECKING_TSI_VALUES",
      {.description = "Calculate nodal values for evaluation and validation of necking",
          .default_value = false}));

  tsidynmono.specs.emplace_back(deprecated_selection<LineSearch>("TSI_LINE_SEARCH",
      {
          {"none", LS_none},
          {"structure", LS_structure},
          {"thermo", LS_thermo},
          {"and", LS_and},
          {"or", LS_or},
      },
      {.description = "line-search strategy", .default_value = LS_none}));

  tsidynmono.move_into_collection(list);

  /*----------------------------------------------------------------------*/
  /* parameters for partitioned TSI */

  Core::Utils::SectionSpecs tsidynpart{tsidyn, "PARTITIONED"};

  std::vector<std::string> couplvariable_valid_input = {"Displacement", "Temperature"};
  tsidynpart.specs.emplace_back(
      deprecated_selection<std::string>("COUPVARIABLE", couplvariable_valid_input,
          {.description = "Coupling variable", .default_value = "Displacement"}));


  // Solver parameter for relaxation of iterative staggered partitioned TSI
  tsidynpart.specs.emplace_back(parameter<double>("MAXOMEGA",
      {.description = "largest omega allowed for Aitken relaxation (0.0 means no constraint)",
          .default_value = 0.0}));
  tsidynpart.specs.emplace_back(parameter<double>(
      "FIXEDOMEGA", {.description = "fixed relaxation parameter", .default_value = 1.0}));

  // convergence tolerance of outer iteration loop
  tsidynpart.specs.emplace_back(parameter<double>("CONVTOL",
      {.description = "tolerance for convergence check of outer iteraiton within partitioned TSI",
          .default_value = 1e-6}));

  tsidynpart.move_into_collection(list);

  /*----------------------------------------------------------------------*/
  /* parameters for tsi contact */
  Core::Utils::SectionSpecs tsic{"TSI CONTACT"};

  tsic.specs.emplace_back(parameter<double>(
      "HEATTRANSSLAVE", {.description = "Heat transfer parameter for slave side in thermal contact",
                            .default_value = 0.0}));
  tsic.specs.emplace_back(parameter<double>("HEATTRANSMASTER",
      {.description = "Heat transfer parameter for master side in thermal contact",
          .default_value = 0.0}));
  tsic.specs.emplace_back(parameter<double>("TEMP_DAMAGE",
      {.description = "damage temperature at contact interface: friction coefficient zero there",
          .default_value = 1.0e12}));
  tsic.specs.emplace_back(
      parameter<double>("TEMP_REF", {.description = "reference temperature at contact interface: "
                                                    "friction coefficient equals the given value",
                                        .default_value = 0.0}));

  tsic.specs.emplace_back(parameter<double>(
      "NITSCHE_THETA_TSI", {.description = "+1: symmetric, 0: non-symmetric, -1: skew-symmetric",
                               .default_value = 0.0}));

  tsic.specs.emplace_back(parameter<CONTACT::NitscheWeighting>("NITSCHE_WEIGHTING_TSI",
      {.description = "how to weight consistency terms in Nitsche contact formulation",
          .default_value = CONTACT::NitscheWeighting::harmonic}));

  tsic.specs.emplace_back(parameter<bool>("NITSCHE_PENALTY_ADAPTIVE_TSI",
      {.description = "adapt penalty parameter after each converged time step",
          .default_value = true}));

  tsic.specs.emplace_back(parameter<double>("PENALTYPARAM_THERMO",
      {.description = "Penalty parameter for Nitsche solution strategy", .default_value = 0.0}));

  tsic.move_into_collection(list);
}

FOUR_C_NAMESPACE_CLOSE
