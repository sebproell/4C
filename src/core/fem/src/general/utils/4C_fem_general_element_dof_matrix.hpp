// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FEM_GENERAL_ELEMENT_DOF_MATRIX_HPP
#define FOUR_C_FEM_GENERAL_ELEMENT_DOF_MATRIX_HPP

#include "4C_config.hpp"

#include "4C_fem_general_cell_type.hpp"
#include "4C_fem_general_cell_type_traits.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_utils_exceptions.hpp"

#include <ranges>
#include <type_traits>



FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  /*!
   * @brief Returns a matrix (num_dof_per_node x num_nodes) from a flattened vector of dof
   * values.
   *
   * @note The flattened vector of dof values needs to be ordered as (x_1, y_1, z_1, x_2, y_2, z_2,
   * ...).
   *
   * @tparam celltype : Celltype of the element
   * @tparam num_dof_per_node : Number of dofs per node
   * @tparam R : Type of the range
   * @param values
   * @return Core::LinAlg::Matrix<num_dof_per_node, Core::FE::num_nodes<celltype>>
   */
  template <Core::FE::CellType celltype, int num_dof_per_node, std::ranges::contiguous_range R>
    requires std::ranges::sized_range<R>
  auto get_element_dof_matrix(R&& values) -> Core::LinAlg::Matrix<num_dof_per_node,
      Core::FE::num_nodes<celltype>, typename std::remove_cvref_t<R>::value_type>
  {
    FOUR_C_ASSERT(std::ranges::size(values) == num_dof_per_node * Core::FE::num_nodes<celltype>,
        "Expecting a size of the span of %d, but got %d",
        num_dof_per_node * Core::FE::num_nodes<celltype>, values.size());

    constexpr bool view = false;
    return Core::LinAlg::Matrix<num_dof_per_node, Core::FE::num_nodes<celltype>,
        typename std::remove_cvref_t<R>::value_type>(std::ranges::data(values), view);
  }

  /*!
   * @copydoc Core::FE::get_element_dof_matrix()
   *
   * @note This function only creates a view of existing data. The parameter @p values must
   outlive
   * the returned matrix object.
   */
  template <Core::FE::CellType celltype, int num_dof_per_node, std::ranges::contiguous_range R>
    requires std::ranges::sized_range<R>
  auto get_element_dof_matrix_view(R&& values) -> Core::LinAlg::Matrix<num_dof_per_node,
      Core::FE::num_nodes<celltype>, typename std::remove_cvref_t<R>::value_type>
  {
    FOUR_C_ASSERT(std::ranges::size(values) == num_dof_per_node * Core::FE::num_nodes<celltype>,
        "Expecting a size of the span of %d, but got %d",
        num_dof_per_node * Core::FE::num_nodes<celltype>, values.size());
    constexpr bool view = true;
    return Core::LinAlg::Matrix<num_dof_per_node, Core::FE::num_nodes<celltype>,
        typename std::remove_cvref_t<R>::value_type>(std::ranges::data(values), view);
  }
}  // namespace Core::FE


FOUR_C_NAMESPACE_CLOSE

#endif
