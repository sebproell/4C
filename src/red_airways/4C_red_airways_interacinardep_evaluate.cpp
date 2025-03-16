// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_mat_list.hpp"
#include "4C_mat_newtonianfluid.hpp"
#include "4C_red_airways_elementbase.hpp"
#include "4C_red_airways_interacinardep_impl.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_SerialDenseSolver.hpp>

FOUR_C_NAMESPACE_OPEN



/*---------------------------------------------------------------------*
 |Evaluate the element (public)                            ismail 09/12|
 *---------------------------------------------------------------------*/
int Discret::Elements::RedInterAcinarDep::evaluate(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
    Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
    Core::LinAlg::SerialDenseVector& elevec3)
{
  Discret::Elements::RedInterAcinarDep::ActionType act = RedInterAcinarDep::none;

  // get the action required
  std::string action = params.get<std::string>("action", "none");
  if (action == "none")
    FOUR_C_THROW("No action supplied");
  else if (action == "calc_sys_matrix_rhs")
    act = RedInterAcinarDep::calc_sys_matrix_rhs;
  else if (action == "calc_sys_matrix_rhs_iad")
    act = RedInterAcinarDep::calc_sys_matrix_rhs_iad;
  else if (action == "get_initial_state")
    act = RedInterAcinarDep::get_initial_state;
  else if (action == "set_bc")
    act = RedInterAcinarDep::set_bc;
  else if (action == "calc_flow_rates")
    act = RedInterAcinarDep::calc_flow_rates;
  else if (action == "calc_elem_volumes")
    act = RedInterAcinarDep::calc_elem_volumes;
  else if (action == "get_coupled_values")
    act = RedInterAcinarDep::get_coupled_values;
  else
  {
    FOUR_C_THROW("Unknown type of action ({}) for inter-acinar linker element", action.c_str());
  }

  /*
    Here one must add the steps for evaluating an element
  */
  std::shared_ptr<Core::Mat::Material> mat = material();

  switch (act)
  {
    case calc_sys_matrix_rhs:
    {
      return Discret::Elements::RedInterAcinarDepImplInterface::impl(this)->evaluate(
          this, params, discretization, lm, elemat1, elemat2, elevec1, elevec2, elevec3, mat);
    }
    break;
    case calc_sys_matrix_rhs_iad:
    {
    }
    break;
    case get_initial_state:
    {
      Discret::Elements::RedInterAcinarDepImplInterface::impl(this)->initial(
          this, params, discretization, lm, elevec3, mat);
    }
    break;
    case set_bc:
    {
      Discret::Elements::RedInterAcinarDepImplInterface::impl(this)->evaluate_terminal_bc(
          this, params, discretization, lm, elevec1, mat);
    }
    break;
    case calc_flow_rates:
    {
    }
    break;
    case calc_elem_volumes:
    {
    }
    break;

    case get_coupled_values:
    {
      Discret::Elements::RedInterAcinarDepImplInterface::impl(this)->get_coupled_values(
          this, params, discretization, lm, mat);
    }
    break;
    default:
      FOUR_C_THROW("Unknown type of action for reduced dimensional acinuss");
      break;
  }  // end of switch(act)

  return 0;
}  // end of Discret::Elements::RedInterAcinarDep::Evaluate


int Discret::Elements::RedInterAcinarDep::evaluate_neumann(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, Core::Conditions::Condition& condition,
    std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1,
    Core::LinAlg::SerialDenseMatrix* elemat1)
{
  return 0;
}

/*----------------------------------------------------------------------*
 |  do nothing (public)                                     ismail 09/12|
 |                                                                      |
 |  The function is just a dummy.                                       |
 *----------------------------------------------------------------------*/
int Discret::Elements::RedInterAcinarDep::evaluate_dirichlet(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, Core::Conditions::Condition& condition,
    std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1)
{
  return 0;
}


/*----------------------------------------------------------------------*
 | Get optimal gaussrule for discretisation type                        |
 |                                                                      |
 *----------------------------------------------------------------------*/
Core::FE::GaussRule1D Discret::Elements::RedInterAcinarDep::get_optimal_gaussrule(
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
      FOUR_C_THROW("Unknown number of nodes for Gaussrule initialization in inter-acinar linker.");
      break;
  }
  return rule;
}


/*----------------------------------------------------------------------*
 | Check, whether higher order derivatives for shape functions          |
 | (dxdx, dxdy, ...) are necessary|                                     |
 *----------------------------------------------------------------------*/
bool Discret::Elements::RedInterAcinarDep::is_higher_order_element(
    const Core::FE::CellType distype) const
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
      break;
  }
  return hoel;
}

FOUR_C_NAMESPACE_CLOSE
