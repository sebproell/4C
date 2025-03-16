// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_structure_new_nln_solver_factory.hpp"

#include "4C_structure_new_nln_solver_fullnewton.hpp"
#include "4C_structure_new_nln_solver_nox.hpp"
#include "4C_structure_new_nln_solver_ptc.hpp"
#include "4C_structure_new_nln_solver_singlestep.hpp"
#include "4C_structure_new_nln_solver_uzawa.hpp"

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
    const enum Inpar::Solid::NonlinSolTech& nlnSolType) const
{
  std::shared_ptr<Solid::Nln::SOLVER::Generic> nlnSolver = nullptr;

  switch (nlnSolType)
  {
    case Inpar::Solid::soltech_newtonfull:
      nlnSolver = std::make_shared<Solid::Nln::SOLVER::FullNewton>();
      break;
    case Inpar::Solid::soltech_nox_nln:
      nlnSolver = std::make_shared<Solid::Nln::SOLVER::Nox>();
      break;
    case Inpar::Solid::soltech_ptc:
      nlnSolver = std::make_shared<Solid::Nln::SOLVER::PseudoTransient>();
      break;
    case Inpar::Solid::soltech_singlestep:
      nlnSolver = std::make_shared<Solid::Nln::SOLVER::SingleStep>();
      break;
    case Inpar::Solid::soltech_newtonuzawanonlin:
    case Inpar::Solid::soltech_newtonuzawalin:
      //      nlnSolver = Teuchos::rcp(new Solid::Nln::SOLVER::Uzawa());
      //      break;
    default:
      FOUR_C_THROW("Solution technique \"{}\" is not implemented.",
          Inpar::Solid::nonlin_sol_tech_string(nlnSolType).c_str());
      break;
  }

  return nlnSolver;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::shared_ptr<Solid::Nln::SOLVER::Generic> Solid::Nln::SOLVER::build_nln_solver(
    const enum Inpar::Solid::NonlinSolTech& nlnSolType)
{
  Factory factory;
  return factory.build_nln_solver(nlnSolType);
}

FOUR_C_NAMESPACE_CLOSE
