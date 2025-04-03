// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_inpar_constraint_framework.hpp"

#include "4C_io_input_spec_builders.hpp"

FOUR_C_NAMESPACE_OPEN

/**
 *
 */
void Inpar::CONSTRAINTS::set_valid_parameters(std::map<std::string, Core::IO::InputSpec>& list)
{
  using Teuchos::tuple;
  using namespace Core::IO::InputSpecBuilders;

  list["EMBEDDED MESH COUPLING"] = all_of({

      parameter<EmbeddedMeshCouplingStrategy>(
          "COUPLING_STRATEGY", {.description = "Strategy to couple background and overlapping mesh",
                                   .default_value = EmbeddedMeshCouplingStrategy::none}),


      parameter<SolidToSolidMortarShapefunctions>("MORTAR_SHAPE_FUNCTION",
          {.description = "Shape functions that should be use in case of coupling using the "
                          "Mortar/Lagrange  Multiplier method ",
              .default_value = SolidToSolidMortarShapefunctions::none}),


      parameter<EmbeddedMeshConstraintEnforcement>("CONSTRAINT_ENFORCEMENT",
          {.description = "Apply a constraint enforcement in the embedded mesh coupling strategy",
              .default_value = EmbeddedMeshConstraintEnforcement::none}),

      parameter<double>("CONSTRAINT_ENFORCEMENT_PENALTYPARAM",
          {.description =
                  "Penalty parameter for the constraint enforcement in embedded mesh coupling",
              .default_value = 0.0})});
}

FOUR_C_NAMESPACE_CLOSE