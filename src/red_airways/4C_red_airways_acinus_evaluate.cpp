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
#include "4C_red_airways_acinus_impl.hpp"
#include "4C_red_airways_elementbase.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------*
 |evaluate the element (public)                            ismail 09/12|
 *---------------------------------------------------------------------*/
int Discret::Elements::RedAcinus::evaluate(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
    Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
    Core::LinAlg::SerialDenseVector& elevec3)
{
  Discret::Elements::RedAcinus::ActionType act = RedAcinus::none;

  // get the action required
  std::string action = params.get<std::string>("action", "none");
  if (action == "none")
    FOUR_C_THROW("No action supplied");
  else if (action == "calc_sys_matrix_rhs")
    act = RedAcinus::calc_sys_matrix_rhs;
  else if (action == "calc_sys_matrix_rhs_iad")
    act = RedAcinus::calc_sys_matrix_rhs_iad;
  else if (action == "get_initial_state")
    act = RedAcinus::get_initial_state;
  else if (action == "set_bc")
    act = RedAcinus::set_bc;
  else if (action == "calc_flow_rates")
    act = RedAcinus::calc_flow_rates;
  else if (action == "calc_elem_volumes")
    act = RedAcinus::calc_elem_volumes;
  else if (action == "get_coupled_values")
    act = RedAcinus::get_coupled_values;
  else
  {
    FOUR_C_THROW("Unknown type of action ({}) for reduced dimensional acinus", action.c_str());
  }

  /*
  Here one must add the steps for evaluating an element
  */
  std::shared_ptr<Core::Mat::Material> mat = material();

  switch (act)
  {
    case calc_sys_matrix_rhs:
    {
      return Discret::Elements::RedAcinusImplInterface::impl(this)->evaluate(
          this, params, discretization, lm, elemat1, elemat2, elevec1, elevec2, elevec3, mat);
    }
    break;
    case calc_sys_matrix_rhs_iad:
    {
    }
    break;
    case get_initial_state:
    {
      Discret::Elements::RedAcinusImplInterface::impl(this)->initial(
          this, params, discretization, lm, mat);
    }
    break;
    case set_bc:
    {
      Discret::Elements::RedAcinusImplInterface::impl(this)->evaluate_terminal_bc(
          this, params, discretization, lm, elevec1, mat);
    }
    break;
    case calc_flow_rates:
    {
      Discret::Elements::RedAcinusImplInterface::impl(this)->calc_flow_rates(
          this, params, discretization, lm, mat);
    }
    break;
    case calc_elem_volumes:
    {
      Discret::Elements::RedAcinusImplInterface::impl(this)->calc_elem_volume(
          this, params, discretization, lm, mat);
    }
    break;
    case get_coupled_values:
    {
      Discret::Elements::RedAcinusImplInterface::impl(this)->get_coupled_values(
          this, params, discretization, lm, mat);
    }
    break;
    default:
      FOUR_C_THROW("Unknown type of action for reduced dimensional acinuss");
  }  // end of switch(act)

  return 0;
}  // end of Discret::Elements::RedAcinus::Evaluate


/*----------------------------------------------------------------------*
 |  do nothing (public)                                     ismail 09/12|
 |                                                                      |
 |  The function is just a dummy.                                       |
 *----------------------------------------------------------------------*/
int Discret::Elements::RedAcinus::evaluate_neumann(Teuchos::ParameterList& params,
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
int Discret::Elements::RedAcinus::evaluate_dirichlet(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, Core::Conditions::Condition& condition,
    std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1)
{
  return 0;
}

FOUR_C_NAMESPACE_CLOSE
