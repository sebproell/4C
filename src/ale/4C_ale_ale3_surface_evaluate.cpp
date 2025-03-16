// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_ale_ale3.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_element_integration_select.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_fem_general_utils_boundary_integration.hpp"
#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_fem_geometry_position_array.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Discret::Elements::Ale3SurfaceImplInterface* Discret::Elements::Ale3SurfaceImplInterface::impl(
    Discret::Elements::Ale3Surface* ele)
{
  switch (ele->shape())
  {
    case Core::FE::CellType::quad4:
    {
      return Discret::Elements::Ale3SurfaceImpl<Core::FE::CellType::quad4>::instance(
          Core::Utils::SingletonAction::create);
    }
    default:
      FOUR_C_THROW("shape {} ({} nodes) not supported", ele->shape(), ele->num_node());
      break;
  }
  return nullptr;
}

template <Core::FE::CellType distype>
Discret::Elements::Ale3SurfaceImpl<distype>* Discret::Elements::Ale3SurfaceImpl<distype>::instance(
    Core::Utils::SingletonAction action)
{
  static auto singleton_owner = Core::Utils::make_singleton_owner(
      []()
      {
        return std::unique_ptr<Discret::Elements::Ale3SurfaceImpl<distype>>(
            new Discret::Elements::Ale3SurfaceImpl<distype>());
      });

  return singleton_owner.instance(action);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
int Discret::Elements::Ale3Surface::evaluate(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
    Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
    Core::LinAlg::SerialDenseVector& elevec3)
{
  const auto act = Teuchos::getIntegralValue<Ale3::ActionType>(params, "action");

  switch (act)
  {
    case Ale3::boundary_calc_ale_node_normal:
    {
      std::shared_ptr<const Core::LinAlg::Vector<double>> dispnp;
      std::vector<double> mydispnp;

      dispnp = discretization.get_state("dispnp");

      if (dispnp != nullptr)
      {
        mydispnp.resize(lm.size());
        mydispnp = Core::FE::extract_values(*dispnp, lm);
      }

      Ale3SurfaceImplInterface::impl(this)->element_node_normal(
          this, params, discretization, lm, elevec1, mydispnp);

      break;
    }
    default:
      FOUR_C_THROW("Unknown type of action '{}' for Ale3Surface", act);
      break;
  }  // end of switch(act)

  return 0;
}  // end of Discret::Elements::Ale3Surface::Evaluate

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
int Discret::Elements::Ale3Surface::evaluate_neumann(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, Core::Conditions::Condition& condition,
    std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1,
    Core::LinAlg::SerialDenseMatrix* elemat1)
{
  return 0;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
template <Core::FE::CellType distype>
inline void Discret::Elements::Ale3SurfaceImpl<distype>::element_node_normal(Ale3Surface* ele,
    Teuchos::ParameterList& params, Core::FE::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseVector& elevec1, std::vector<double>& mydispnp)
{
  Core::FE::element_node_normal<distype>(funct_, deriv_, fac_, unitnormal_, drs_, xsi_, xyze_, ele,
      discretization, elevec1, mydispnp, false, true);
}

FOUR_C_NAMESPACE_CLOSE
