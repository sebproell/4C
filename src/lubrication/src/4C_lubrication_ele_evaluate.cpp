// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fem_discretization.hpp"
#include "4C_lubrication_ele.hpp"
#include "4C_lubrication_ele_action.hpp"
#include "4C_lubrication_ele_factory.hpp"
#include "4C_lubrication_ele_interface.hpp"
#include "4C_lubrication_ele_parameter.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                           wirtz 10/15 |
 *----------------------------------------------------------------------*/
int Discret::Elements::Lubrication::evaluate(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, Core::Elements::LocationArray& la,
    Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
    Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
    Core::LinAlg::SerialDenseVector& elevec3)
{
  // check for the action parameter
  const auto action = Teuchos::getIntegralValue<FourC::Lubrication::Action>(params, "action");
  switch (action)
  {
    // all physics-related stuff is included in the implementation class(es) that can
    // be used in principle inside any element (at the moment: only Lubrication element)
    case FourC::Lubrication::calc_mat_and_rhs:
    {
      return Discret::Elements::LubricationFactory::provide_impl(shape(), discretization.name())
          ->evaluate(this, params, discretization, la, elemat1, elemat2, elevec1, elevec2, elevec3);
      break;
    }
    case FourC::Lubrication::calc_lubrication_coupltang:
    {
      return Discret::Elements::LubricationFactory::provide_impl(shape(), discretization.name())
          ->evaluate_ehl_mon(
              this, params, discretization, la, elemat1, elemat2, elevec1, elevec2, elevec3);
      break;
    }
    case FourC::Lubrication::calc_error:
    case FourC::Lubrication::calc_mean_pressures:
    {
      return Discret::Elements::LubricationFactory::provide_impl(shape(), discretization.name())
          ->evaluate_service(
              this, params, discretization, la, elemat1, elemat2, elevec1, elevec2, elevec3);
      break;
    }
    case FourC::Lubrication::set_time_parameter:
    case FourC::Lubrication::set_general_lubrication_parameter:
      // these actions have already been evaluated during element pre-evaluate
      break;
    default:
    {
      FOUR_C_THROW("Unknown type of action '{}' for Lubrication", action);
      break;
    }
  }  // switch(action)

  return 0;
}  // Discret::Elements::Lubrication::Evaluate


/*----------------------------------------------------------------------*
 |  dummy                                                   wirtz 10/15 |
 *----------------------------------------------------------------------*/
int Discret::Elements::Lubrication::evaluate_neumann(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, Core::Conditions::Condition& condition,
    std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1,
    Core::LinAlg::SerialDenseMatrix* elemat1)
{
  //    The function is just a dummy. For Lubrication elements, the integration
  //    integration of volume Neumann conditions (body forces) takes place
  //    in the element. We need it there for potential stabilisation terms! (wirtz)
  return 0;
}

/*---------------------------------------------------------------------*
|  Call the element to set all basic parameter             wirtz 10/15 |
*----------------------------------------------------------------------*/
void Discret::Elements::LubricationType::pre_evaluate(Core::FE::Discretization& dis,
    Teuchos::ParameterList& p, std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix1,
    std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix2,
    std::shared_ptr<Core::LinAlg::Vector<double>> systemvector1,
    std::shared_ptr<Core::LinAlg::Vector<double>> systemvector2,
    std::shared_ptr<Core::LinAlg::Vector<double>> systemvector3)
{
  const auto action = Teuchos::getIntegralValue<FourC::Lubrication::Action>(p, "action");

  switch (action)
  {
    case FourC::Lubrication::set_general_lubrication_parameter:
    {
      Discret::Elements::LubricationEleParameter::instance(dis.name())->set_general_parameters(p);

      break;
    }

    case FourC::Lubrication::set_time_parameter:
    {
      Discret::Elements::LubricationEleParameter::instance(dis.name())->set_time_parameters(p);

      break;
    }
    default:
      // do nothing in all other cases
      break;
  }

  return;
}

FOUR_C_NAMESPACE_CLOSE
