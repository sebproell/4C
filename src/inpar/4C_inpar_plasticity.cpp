// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_inpar_plasticity.hpp"

#include "4C_inpar_structure.hpp"
#include "4C_inpar_tsi.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN



void Inpar::Plasticity::set_valid_parameters(std::map<std::string, Core::IO::InputSpec>& list)
{
  using Teuchos::tuple;
  using namespace Core::IO::InputSpecBuilders;

  /*----------------------------------------------------------------------*/
  /* parameters for semi-smooth Newton plasticity algorithm */
  Core::Utils::SectionSpecs iplast{"SEMI-SMOOTH PLASTICITY"};

  iplast.specs.emplace_back(parameter<double>("SEMI_SMOOTH_CPL",
      {.description = "Weighting factor cpl for semi-smooth PDASS", .default_value = 1.0}));
  iplast.specs.emplace_back(parameter<double>("STABILIZATION_S",
      {.description = "Stabilization factor s for semi-smooth PDASS", .default_value = 1.0}));

  // solver convergence test parameters for semi-smooth plasticity formulation
  iplast.specs.emplace_back(deprecated_selection<Inpar::Solid::BinaryOp>(
      "NORMCOMBI_RESFPLASTCONSTR",
      {
          {"And", Inpar::Solid::bop_and},
          {"Or", Inpar::Solid::bop_or},
      },
      {.description = "binary operator to combine plasticity constraints and residual force values",
          .default_value = Inpar::Solid::bop_and}));

  iplast.specs.emplace_back(deprecated_selection<Inpar::Solid::BinaryOp>("NORMCOMBI_DISPPLASTINCR",
      {
          {"And", Inpar::Solid::bop_and},
          {"Or", Inpar::Solid::bop_or},
      },
      {.description = "binary operator to combine displacement increments and plastic flow (Delta "
                      "Lp) increment values",
          .default_value = Inpar::Solid::bop_and}));

  iplast.specs.emplace_back(parameter<double>("TOLPLASTCONSTR",
      {.description = "tolerance in the plastic constraint norm for the newton iteration",
          .default_value = 1.0E-8}));
  iplast.specs.emplace_back(parameter<double>("TOLDELTALP",
      {.description = "tolerance in the plastic flow (Delta Lp) norm for the Newton iteration",
          .default_value = 1.0E-8}));

  iplast.specs.emplace_back(deprecated_selection<Inpar::Solid::BinaryOp>("NORMCOMBI_EASRES",
      {
          {"And", Inpar::Solid::bop_and},
          {"Or", Inpar::Solid::bop_or},
      },
      {.description = "binary operator to combine EAS-residual and residual force values",
          .default_value = Inpar::Solid::bop_and}));

  iplast.specs.emplace_back(deprecated_selection<Inpar::Solid::BinaryOp>("NORMCOMBI_EASINCR",
      {
          {"And", Inpar::Solid::bop_and},
          {"Or", Inpar::Solid::bop_or},
      },
      {.description = "binary operator to combine displacement increments and EAS increment values",
          .default_value = Inpar::Solid::bop_and}));

  iplast.specs.emplace_back(parameter<double>(
      "TOLEASRES", {.description = "tolerance in the EAS residual norm for the newton iteration",
                       .default_value = 1.0E-8}));
  iplast.specs.emplace_back(parameter<double>(
      "TOLEASINCR", {.description = "tolerance in the EAS increment norm for the Newton iteration",
                        .default_value = 1.0E-8}));

  iplast.specs.emplace_back(deprecated_selection<Inpar::TSI::DissipationMode>("DISSIPATION_MODE",
      {
          {"pl_multiplier", Inpar::TSI::pl_multiplier},
          {"pl_flow", Inpar::TSI::pl_flow},
          {"Taylor_Quinney", Inpar::TSI::Taylor_Quinney},
      },
      {.description = "method to calculate the plastic dissipation",
          .default_value = Inpar::TSI::pl_multiplier}));

  iplast.move_into_collection(list);
}

FOUR_C_NAMESPACE_CLOSE
