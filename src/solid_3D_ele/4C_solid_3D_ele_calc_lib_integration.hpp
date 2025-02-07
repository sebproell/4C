// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLID_3D_ELE_CALC_LIB_INTEGRATION_HPP
#define FOUR_C_SOLID_3D_ELE_CALC_LIB_INTEGRATION_HPP

#include "4C_config.hpp"

#include "4C_fem_general_cell_type.hpp"
#include "4C_fem_general_element_integration_select.hpp"
#include "4C_fem_general_utils_gausspoints.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

FOUR_C_NAMESPACE_OPEN
namespace Discret::Elements
{
  /*!
   * @brief Compare two Gauss integration rules for equality
   */
  inline bool compare_gauss_integration(const Core::FE::GaussIntegration& integration_a,
      const Core::FE::GaussIntegration& integration_b)
  {
    // currently this simple check is sufficient as we only use the same type of gauss integrations.
    return integration_a.num_points() == integration_b.num_points();
  }

  /*!
   * @brief Get the default Gauss integration rules for different Cell types.
   *
   * @note It follows the rules defined in Discret::Elements::DisTypeToOptGaussRule<celltype>::rule,
   * except for the stiffness matrix of tetrahedral elements.
   *
   */
  /// @{

  template <Core::FE::CellType celltype>
  constexpr auto get_gauss_rule_mass_matrix()
  {
    return Discret::Elements::DisTypeToOptGaussRule<celltype>::rule;
  }

  template <Core::FE::CellType celltype>
  constexpr auto get_gauss_rule_stiffness_matrix()
  {
    return Discret::Elements::DisTypeToOptGaussRule<celltype>::rule;
  }

  template <>
  constexpr auto get_gauss_rule_stiffness_matrix<Core::FE::CellType::tet10>()
  {
    return Core::FE::GaussRule3D::tet_4point;
  }

  template <>
  constexpr auto get_gauss_rule_stiffness_matrix<Core::FE::CellType::tet4>()
  {
    return Core::FE::GaussRule3D::tet_1point;
  }
  /// @}
}  // namespace Discret::Elements
FOUR_C_NAMESPACE_CLOSE

#endif