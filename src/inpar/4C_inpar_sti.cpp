// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_inpar_sti.hpp"

#include "4C_fem_condition_definition.hpp"
#include "4C_inpar_scatra.hpp"
#include "4C_linalg_sparseoperator.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*------------------------------------------------------------------------*
 | set valid parameters for scatra-thermo interaction          fang 10/16 |
 *------------------------------------------------------------------------*/
void Inpar::STI::set_valid_parameters(std::map<std::string, Core::IO::InputSpec>& list)
{
  using Teuchos::tuple;
  using namespace Core::IO::InputSpecBuilders;

  Core::Utils::SectionSpecs stidyn{"STI DYNAMIC"};

  // type of scalar transport time integration
  stidyn.specs.emplace_back(deprecated_selection<ScaTraTimIntType>("SCATRATIMINTTYPE",
      {
          {"Standard", ScaTraTimIntType::standard},
          {"Elch", ScaTraTimIntType::elch},
      },
      {.description =
              "scalar transport time integration type is needed to instantiate correct scalar "
              "transport time integration scheme for scatra-thermo interaction problems",
          .default_value = ScaTraTimIntType::standard}));

  // type of coupling between scatra and thermo fields
  stidyn.specs.emplace_back(deprecated_selection<CouplingType>("COUPLINGTYPE",
      {
          {"Undefined", CouplingType::undefined},
          {"Monolithic", CouplingType::monolithic},
          {"OneWay_ScatraToThermo", CouplingType::oneway_scatratothermo},
          {"OneWay_ThermoToScatra", CouplingType::oneway_thermotoscatra},
          {"TwoWay_ScatraToThermo", CouplingType::twoway_scatratothermo},
          {"TwoWay_ScatraToThermo_Aitken", CouplingType::twoway_scatratothermo_aitken},
          {"TwoWay_ScatraToThermo_Aitken_Dofsplit",
              CouplingType::twoway_scatratothermo_aitken_dofsplit},
          {"TwoWay_ThermoToScatra", CouplingType::twoway_thermotoscatra},
          {"TwoWay_ThermoToScatra_Aitken", CouplingType::twoway_thermotoscatra_aitken},
      },
      {.description = "type of coupling between scatra and thermo fields",
          .default_value = CouplingType::undefined}));

  // specification of initial temperature field
  stidyn.specs.emplace_back(deprecated_selection<Inpar::ScaTra::InitialField>("THERMO_INITIALFIELD",
      {
          {"zero_field", Inpar::ScaTra::initfield_zero_field},
          {"field_by_function", Inpar::ScaTra::initfield_field_by_function},
          {"field_by_condition", Inpar::ScaTra::initfield_field_by_condition},
      },
      {.description = "initial temperature field for scatra-thermo interaction problems",
          .default_value = Inpar::ScaTra::initfield_zero_field}));

  // function number for initial temperature field
  stidyn.specs.emplace_back(parameter<int>(
      "THERMO_INITFUNCNO", {.description = "function number for initial temperature field for "
                                           "scatra-thermo interaction problems",
                               .default_value = -1}));

  // ID of linear solver for temperature field
  stidyn.specs.emplace_back(parameter<int>("THERMO_LINEAR_SOLVER",
      {.description = "ID of linear solver for temperature field", .default_value = -1}));

  // flag for double condensation of linear equations associated with temperature field
  stidyn.specs.emplace_back(parameter<bool>("THERMO_CONDENSATION",
      {.description =
              "flag for double condensation of linear equations associated with temperature field",
          .default_value = false}));

  stidyn.move_into_collection(list);

  /*----------------------------------------------------------------------*/
  // valid parameters for monolithic scatra-thermo interaction
  Core::Utils::SectionSpecs stidyn_monolithic{stidyn, "MONOLITHIC"};

  // ID of linear solver for global system of equations
  stidyn_monolithic.specs.emplace_back(parameter<int>("LINEAR_SOLVER",
      {.description = "ID of linear solver for global system of equations", .default_value = -1}));

  // type of global system matrix in global system of equations
  stidyn_monolithic.specs.emplace_back(deprecated_selection<Core::LinAlg::MatrixType>("MATRIXTYPE",
      {
          {"block", Core::LinAlg::MatrixType::block_condition},
          {"sparse", Core::LinAlg::MatrixType::sparse},
      },
      {.description = "type of global system matrix in global system of equations",
          .default_value = Core::LinAlg::MatrixType::block_condition}));

  stidyn_monolithic.move_into_collection(list);

  /*----------------------------------------------------------------------*/
  // valid parameters for partitioned scatra-thermo interaction
  Core::Utils::SectionSpecs stidyn_partitioned{stidyn, "PARTITIONED"};

  // relaxation parameter
  stidyn_partitioned.specs.emplace_back(
      parameter<double>("OMEGA", {.description = "relaxation parameter", .default_value = 1.}));

  // maximum value of Aitken relaxation parameter
  stidyn_partitioned.specs.emplace_back(parameter<double>("OMEGAMAX",
      {.description = "maximum value of Aitken relaxation parameter (0.0 = no constraint)",
          .default_value = 0.}));

  stidyn_partitioned.move_into_collection(list);
}


/*------------------------------------------------------------------------*
 | set valid conditions for scatra-thermo interaction          fang 10/16 |
 *------------------------------------------------------------------------*/
void Inpar::STI::set_valid_conditions(std::vector<Core::Conditions::ConditionDefinition>& condlist)
{
  return;
}

FOUR_C_NAMESPACE_CLOSE
