// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_fem_general_element_dof_matrix.hpp"

#include <vector>

namespace
{
  using namespace FourC;

  template <typename ArrayType>
  ArrayType get_dof_array()
  {
    return ArrayType({1.0, 1.1, 2.0, 2.1, 3.0, 3.1, 4.0, 4.1});
  }

  void expect_dof_matrix_equal(const Core::LinAlg::Matrix<2, 4>& dof_matrix)
  {
    EXPECT_EQ(dof_matrix(0, 0), 1.0);
    EXPECT_EQ(dof_matrix(1, 0), 1.1);

    EXPECT_EQ(dof_matrix(0, 1), 2.0);
    EXPECT_EQ(dof_matrix(1, 1), 2.1);

    EXPECT_EQ(dof_matrix(0, 2), 3.0);
    EXPECT_EQ(dof_matrix(1, 2), 3.1);

    EXPECT_EQ(dof_matrix(0, 3), 4.0);
    EXPECT_EQ(dof_matrix(1, 3), 4.1);
  }

  TEST(ElementDofMatrixTest, TestVector)
  {
    const Core::LinAlg::Matrix<2, 4> dof_matrix =
        Core::FE::get_element_dof_matrix<Core::FE::CellType::quad4, 2>(
            get_dof_array<std::vector<double>>());

    expect_dof_matrix_equal(dof_matrix);
  }

  TEST(ElementDofMatrixTest, TestVectorView)
  {
    const auto vec = get_dof_array<std::vector<double>>();
    const Core::LinAlg::Matrix<2, 4> dof_matrix =
        Core::FE::get_element_dof_matrix_view<Core::FE::CellType::quad4, 2>(vec);

    expect_dof_matrix_equal(dof_matrix);
  }

  TEST(ElementDofMatrixTest, TestArray)
  {
    const Core::LinAlg::Matrix<2, 4> dof_matrix =
        Core::FE::get_element_dof_matrix<Core::FE::CellType::quad4, 2>(
            get_dof_array<std::array<double, 8>>());

    expect_dof_matrix_equal(dof_matrix);
  }

  TEST(ElementDofMatrixTest, TestArrayView)
  {
    const auto arr = get_dof_array<std::array<double, 8>>();
    const Core::LinAlg::Matrix<2, 4> dof_matrix =
        Core::FE::get_element_dof_matrix_view<Core::FE::CellType::quad4, 2>(arr);

    expect_dof_matrix_equal(dof_matrix);
  }

  TEST(ElementDofMatrixTest, ModifyArrayView)
  {
    auto arr = get_dof_array<std::array<double, 8>>();
    Core::LinAlg::Matrix<2, 4> dof_matrix =
        Core::FE::get_element_dof_matrix_view<Core::FE::CellType::quad4, 2>(arr);

    dof_matrix(1, 2) = 12.3;

    EXPECT_EQ(arr[5], 12.3);
  }
}  // namespace
