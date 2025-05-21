// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_structure_new_nln_solver_factory.hpp"

#include "4C_structure_new_nln_solver_nox.hpp"
#include "4C_structure_new_nln_solver_utils.hpp"
#include "4C_structure_new_timint_base.hpp"
#include "4C_utils_enum.hpp"

#include <NOX_Utils.H>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Solid::Nln::SOLVER::Factory::Factory()
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::shared_ptr<Solid::Nln::SOLVER::Generic> Solid::Nln::SOLVER::Factory::build_nln_solver(
    const enum Inpar::Solid::NonlinSolTech& nlnSolType,
    const std::shared_ptr<Solid::TimeInt::BaseDataGlobalState>& gstate,
    const std::shared_ptr<Solid::TimeInt::BaseDataSDyn>& sdyn,
    const std::shared_ptr<Solid::TimeInt::NoxInterface>& noxinterface,
    const std::shared_ptr<Solid::Integrator>& integrator,
    const std::shared_ptr<const Solid::TimeInt::Base>& timint) const
{
  std::shared_ptr<Solid::Nln::SOLVER::Generic> nlnSolver = nullptr;

  const auto it_set = settings.find(nlnSolType);

  Teuchos::ParameterList predefined_nox_settings;

  if (it_set != settings.end())
  {
    predefined_nox_settings = it_set->second(*sdyn);
  }

  // ---------------------------------------------------------------------------
  // STATUS TEST
  // ---------------------------------------------------------------------------
  /* This is only necessary for the special case, that you use no xml-file for
   * the definition of your convergence tests, but you use the input file instead.
   */
  if (not is_xml_status_test_file(sdyn->get_nox_params().sublist("Status Test")))
  {
    std::set<enum NOX::Nln::StatusTest::QuantityType> qtypes;
    Solid::Nln::SOLVER::create_quantity_types(qtypes, *sdyn);

    // remove the unsupported quantity of status test:
    // EAS is removed since it is an element
    // quantity and not a nodal dof
    qtypes.erase(NOX::Nln::StatusTest::quantity_eas);

    Solid::Nln::SOLVER::set_status_test_params(
        predefined_nox_settings.sublist("Status Test"), *sdyn, qtypes);
  }

  return std::make_shared<Solid::Nln::SOLVER::Nox>(
      predefined_nox_settings, gstate, sdyn, noxinterface, integrator, timint);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::shared_ptr<Solid::Nln::SOLVER::Generic> Solid::Nln::SOLVER::build_nln_solver(
    const enum Inpar::Solid::NonlinSolTech& nlnSolType,
    const std::shared_ptr<Solid::TimeInt::BaseDataGlobalState>& gstate,
    const std::shared_ptr<Solid::TimeInt::BaseDataSDyn>& sdyn,
    const std::shared_ptr<Solid::TimeInt::NoxInterface>& noxinterface,
    const std::shared_ptr<Solid::Integrator>& integrator,
    const std::shared_ptr<const Solid::TimeInt::Base>& timint)
{
  Factory factory;
  return factory.build_nln_solver(nlnSolType, gstate, sdyn, noxinterface, integrator, timint);
}

FOUR_C_NAMESPACE_CLOSE
