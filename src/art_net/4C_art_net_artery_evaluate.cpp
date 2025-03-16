// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_art_net_artery.hpp"
#include "4C_art_net_artery_ele_action.hpp"
#include "4C_art_net_artery_ele_factory.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_inpar_bio.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_mat_cnst_1d_art.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN


/*---------------------------------------------------------------------*
 //evaluate the element (public)                            ismail 06/09
 *----------------------------------------------------------------------*/
int Discret::Elements::Artery::evaluate(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, Core::Elements::LocationArray& la,
    Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
    Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
    Core::LinAlg::SerialDenseVector& elevec3)
{
  // check for the action parameter
  const auto action = Teuchos::getIntegralValue<Arteries::Action>(params, "action");
  /*
  Here must add the steps for evaluating an element
  */
  std::shared_ptr<Core::Mat::Material> mat = material();

  switch (action)
  {
    case Arteries::calc_sys_matrix_rhs:
    {
      return Discret::Elements::ArtNetFactory::provide_impl(
          shape(), impltype_, discretization.name())
          ->evaluate(
              this, params, discretization, la, elemat1, elemat2, elevec1, elevec2, elevec3, mat);
    }
    break;
    case Arteries::calc_scatra_sys_matrix_rhs:
    {
      return Discret::Elements::ArtNetFactory::provide_impl(
          shape(), impltype_, discretization.name())
          ->scatra_evaluate(this, params, discretization, la[0].lm_, elemat1, elemat2, elevec1,
              elevec2, elevec3, mat);
      break;
    }
    case Arteries::get_initial_artery_state:
    case Arteries::set_term_bc:
    case Arteries::set_scatra_term_bc:
    case Arteries::set_scatra_bc:
    case Arteries::solve_riemann_problem:
    case Arteries::calc_postpro_vals:
    case Arteries::calc_scatra_from_scatra_fb:
    case Arteries::evaluate_wf_wb:
    case Arteries::evaluate_scatra_analytically:
    case Arteries::calc_flow_pressurebased:
    {
      return Discret::Elements::ArtNetFactory::provide_impl(
          shape(), impltype_, discretization.name())
          ->evaluate_service(this, action, params, discretization, la, elemat1, elemat2, elevec1,
              elevec2, elevec3, mat);
    }
    break;
    default:
      FOUR_C_THROW("Unknown type of action {} for Artery", action);
  }  // end of switch(act)


  return 0;
}  // end of Discret::Elements::Artery::Evaluate


int Discret::Elements::Artery::evaluate_neumann(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, Core::Conditions::Condition& condition,
    std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1,
    Core::LinAlg::SerialDenseMatrix* elemat1)
{
  return 0;
}

/*----------------------------------------------------------------------*
 |  do nothing (public)                                     ismail 01/09|
 |                                                                      |
 |  The function is just a dummy.                                       |
 *----------------------------------------------------------------------*/
int Discret::Elements::Artery::evaluate_dirichlet(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, Core::Conditions::Condition& condition,
    std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1)
{
  return 0;
}


// get optimal gaussrule for discretization type
Core::FE::GaussRule1D Discret::Elements::Artery::get_optimal_gaussrule(
    const Core::FE::CellType& distype)
{
  Core::FE::GaussRule1D rule = Core::FE::GaussRule1D::undefined;
  switch (distype)
  {
    case Core::FE::CellType::line2:
      rule = Core::FE::GaussRule1D::line_2point;
      break;
    case Core::FE::CellType::line3:
      rule = Core::FE::GaussRule1D::line_3point;
      break;
    default:
      FOUR_C_THROW("unknown number of nodes for gaussrule initialization");
  }
  return rule;
}


// check, whether higher order derivatives for shape functions (dxdx, dxdy, ...) are necessary
bool Discret::Elements::Artery::is_higher_order_element(const Core::FE::CellType distype) const
{
  bool hoel = true;
  switch (distype)
  {
    case Core::FE::CellType::line3:
      hoel = true;
      break;
    case Core::FE::CellType::line2:
      hoel = false;
      break;
    default:
      FOUR_C_THROW("distype unknown!");
  }
  return hoel;
}

FOUR_C_NAMESPACE_CLOSE
