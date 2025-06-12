// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_config.hpp"

#include "4C_linalg_symmetric_tensor_eigen.hpp"

#include "4C_linalg_symmetric_tensor.hpp"
#include "4C_linalg_tensor_generators.hpp"
#include "4C_unittest_utils_assertions_test.hpp"

FOUR_C_NAMESPACE_OPEN

namespace
{
  TEST(SymmetricTensorEigenTest, eig2x2)
  {
    Core::LinAlg::SymmetricTensor<double, 2, 2> t =
        Core::LinAlg::assume_symmetry(Core::LinAlg::Tensor<double, 2, 2>{
            {{0.9964456203546112, 0.490484665405466}, {0.490484665405466, 0.5611378979071144}}});

    const auto& [eigenvalues, eigenvectors] = Core::LinAlg::eig(t);
    EXPECT_NEAR(eigenvalues[0], 0.242183512545406, 1e-10);
    EXPECT_NEAR(eigenvalues[1], 1.31540000571632, 1e-10);

    Core::LinAlg::Tensor<double, 2, 2> original_matrix =
        eigenvectors * Core::LinAlg::TensorGenerators::diagonal(eigenvalues) *
        Core::LinAlg::transpose(eigenvectors);

    FOUR_C_EXPECT_NEAR(original_matrix, t, 1e-10);
  }

  TEST(SymmetricTensorEigenTest, eig2x2view)
  {
    Core::LinAlg::SymmetricTensor<double, 2, 2> t =
        Core::LinAlg::assume_symmetry(Core::LinAlg::Tensor<double, 2, 2>{
            {{0.9964456203546112, 0.490484665405466}, {0.490484665405466, 0.5611378979071144}}});

    Core::LinAlg::SymmetricTensorView<const double, 2, 2> t_view = t;

    const auto& [eigenvalues, eigenvectors] = Core::LinAlg::eig(t_view);
    EXPECT_NEAR(eigenvalues[0], 0.242183512545406, 1e-10);
    EXPECT_NEAR(eigenvalues[1], 1.31540000571632, 1e-10);

    Core::LinAlg::Tensor<double, 2, 2> original_matrix =
        eigenvectors * Core::LinAlg::TensorGenerators::diagonal(eigenvalues) *
        Core::LinAlg::transpose(eigenvectors);

    FOUR_C_EXPECT_NEAR(original_matrix, t, 1e-10);
  }

  TEST(SymmetricTensorEigenTest, eig3x3)
  {
    Core::LinAlg::SymmetricTensor<double, 3, 3> t =
        Core::LinAlg::assume_symmetry(Core::LinAlg::Tensor<double, 3, 3>{
            {{1.1948716876311152, 0.5399232848096387, 0.6905691418748001},
                {0.5399232848096387, 0.8273468489149256, 0.2538261251935583},
                {0.6905691418748001, 0.2538261251935583, 0.48064209250998646}}});

    const auto& [eigenvalues, eigenvectors] = Core::LinAlg::eig(t);
    EXPECT_NEAR(eigenvalues[0], 0.052879897400611, 1e-10);
    EXPECT_NEAR(eigenvalues[1], 0.516153257193444, 1e-10);
    EXPECT_NEAR(eigenvalues[2], 1.933827474461972, 1e-10);

    Core::LinAlg::Tensor<double, 3, 3> original_matrix =
        eigenvectors * Core::LinAlg::TensorGenerators::diagonal(eigenvalues) *
        Core::LinAlg::transpose(eigenvectors);

    FOUR_C_EXPECT_NEAR(original_matrix, t, 1e-10);
  }

  TEST(SymmetricTensorEigenTest, eig3x3view)
  {
    const Core::LinAlg::SymmetricTensor<double, 3, 3> t =
        Core::LinAlg::assume_symmetry(Core::LinAlg::Tensor<double, 3, 3>{
            {{1.1948716876311152, 0.5399232848096387, 0.6905691418748001},
                {0.5399232848096387, 0.8273468489149256, 0.2538261251935583},
                {0.6905691418748001, 0.2538261251935583, 0.48064209250998646}}});

    Core::LinAlg::SymmetricTensorView<const double, 3, 3> t_view = t;

    const auto& [eigenvalues, eigenvectors] = Core::LinAlg::eig(t_view);
    EXPECT_NEAR(eigenvalues[0], 0.052879897400611, 1e-10);
    EXPECT_NEAR(eigenvalues[1], 0.516153257193444, 1e-10);
    EXPECT_NEAR(eigenvalues[2], 1.933827474461972, 1e-10);

    Core::LinAlg::Tensor<double, 3, 3> original_matrix =
        eigenvectors * Core::LinAlg::TensorGenerators::diagonal(eigenvalues) *
        Core::LinAlg::transpose(eigenvectors);

    FOUR_C_EXPECT_NEAR(original_matrix, t, 1e-10);
  }
}  // namespace

FOUR_C_NAMESPACE_CLOSE