// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_art_net_artery_ele_factory.hpp"

#include "4C_art_net_artery_ele_calc_lin_exp.hpp"
#include "4C_art_net_artery_ele_calc_pres_based.hpp"
#include "4C_art_net_artery_ele_interface.hpp"

FOUR_C_NAMESPACE_OPEN

/*--------------------------------------------------------------------------*
 | (public) kremheller                                                03/18 |
 *--------------------------------------------------------------------------*/
Discret::Elements::ArteryEleInterface* Discret::Elements::ArtNetFactory::provide_impl(
    Core::FE::CellType distype, Inpar::ArtDyn::ImplType problem, const std::string& disname)
{
  switch (distype)
  {
    case Core::FE::CellType::line2:
    {
      return define_problem_type<Core::FE::CellType::line2>(problem, disname);

      break;
    }
      // note by J Kremheller:
      // The current implementation relies on the fact that we only use linear elements on several
      // occasions, for instance, when calculating element volumetric flow and element length
      // I currently do not see any application of higher order elements, since the formulation
      // essentially depends on a linear pressure drop prescribed in each element (Hagen-Poiseuille
      // equation)
      // but if this is ever desired the implementation should be checked carefully
    default:
      FOUR_C_THROW("Only line2 elements available so far");
      break;
  }
  return nullptr;
}


/*--------------------------------------------------------------------------*
 | (public) kremheller                                                03/18 |
 *--------------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::Elements::ArteryEleInterface* Discret::Elements::ArtNetFactory::define_problem_type(
    Inpar::ArtDyn::ImplType problem, const std::string& disname)
{
  switch (problem)
  {
    case Inpar::ArtDyn::ImplType::impltype_lin_exp:
    {
      // 2 dofs per node
      return Discret::Elements::ArteryEleCalcLinExp<distype>::instance(2, disname);
      break;
    }
    case Inpar::ArtDyn::ImplType::impltype_pressure_based:
    {
      // 1 dof per node (only pressure)
      return Discret::Elements::ArteryEleCalcPresBased<distype>::instance(1, disname);
      break;
    }
    default:
    {
      FOUR_C_THROW("Defined problem type {} does not exist!!", problem);
      break;
    }
  }

  return nullptr;
}

FOUR_C_NAMESPACE_CLOSE
