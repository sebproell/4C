// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_structure_new_nln_solver_generic.hpp"

#include "4C_structure_new_timint_base.hpp"
#include "4C_structure_new_timint_implicit.hpp"
#include "4C_structure_new_timint_noxinterface.hpp"

#include <NOX_Abstract_Group.H>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Solid::Nln::SOLVER::Generic::Generic(
    const std::shared_ptr<Solid::TimeInt::BaseDataGlobalState>& gstate,
    const std::shared_ptr<Solid::TimeInt::BaseDataSDyn>& sdyn,
    const std::shared_ptr<Solid::TimeInt::NoxInterface>& noxinterface,
    const std::shared_ptr<Solid::Integrator>& integr,
    const std::shared_ptr<const Solid::TimeInt::Base>& timint)
    : gstate_ptr_(gstate),
      sdyn_ptr_(sdyn),
      noxinterface_ptr_(noxinterface),
      int_ptr_(integr),
      timint_ptr_(timint),
      group_ptr_(nullptr)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::Abstract::Group>& Solid::Nln::SOLVER::Generic::group_ptr()
{
  return group_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
::NOX::Abstract::Group& Solid::Nln::SOLVER::Generic::get_solution_group()
{
  FOUR_C_ASSERT(group_ptr_, "The group pointer should be initialized beforehand!");

  return *group_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const ::NOX::Abstract::Group& Solid::Nln::SOLVER::Generic::get_solution_group() const
{
  FOUR_C_ASSERT(group_ptr_, "The group pointer should be initialized beforehand!");

  return *group_ptr_;
}

FOUR_C_NAMESPACE_CLOSE
