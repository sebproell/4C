// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_config.hpp"

#include "4C_linalg_tensor_svd.hpp"

#include "4C_unittest_utils_assertions_test.hpp"

FOUR_C_NAMESPACE_OPEN

namespace
{
  TEST(TensorSVDTest, svd2x3)
  {
    Core::LinAlg::Tensor<double, 2, 3> t = {
        {{0.771320643266746, 0.0207519493594015, 0.6336482349262754},
            {0.7488038825386119, 0.4985070123025904, 0.22479664553084766}}};

    const auto& [Q, S, VT] = Core::LinAlg::svd(t);
    EXPECT_NEAR(S[0], 1.289138937229215, 1e-10);
    EXPECT_NEAR(S[1], 0.44130158841548733, 1e-10);

    Core::LinAlg::Tensor<double, 2, 3> original_matrix =
        Q * Core::LinAlg::make_rectangular_diagonal_matrix<2, 3>(S) * VT;

    FOUR_C_EXPECT_NEAR(original_matrix, t, 1e-10);
  }

  TEST(TensorSVDTest, svd2x3view)
  {
    Core::LinAlg::Tensor<double, 2, 3> t = {
        {{0.771320643266746, 0.0207519493594015, 0.6336482349262754},
            {0.7488038825386119, 0.4985070123025904, 0.22479664553084766}}};
    Core::LinAlg::TensorView<double, 2, 3> t_view = t;

    const auto& [Q, S, VT] = Core::LinAlg::svd(t_view);
    EXPECT_NEAR(S[0], 1.289138937229215, 1e-10);
    EXPECT_NEAR(S[1], 0.44130158841548733, 1e-10);

    Core::LinAlg::Tensor<double, 2, 3> original_matrix =
        Q * Core::LinAlg::make_rectangular_diagonal_matrix<2, 3>(S) * VT;

    FOUR_C_EXPECT_NEAR(original_matrix, t, 1e-10);
  }

  TEST(TensorSVDTest, svd3x3)
  {
    Core::LinAlg::Tensor<double, 3, 3> t = {
        {{0.771320643266746, 0.0207519493594015, 0.6336482349262754},
            {0.7488038825386119, 0.4985070123025904, 0.22479664553084766},
            {0.19806286475962398, 0.7605307121989587, 0.16911083656253545}}};

    const auto& [Q, S, VT] = Core::LinAlg::svd(t);
    EXPECT_NEAR(S[0], 1.3906212548576882, 1e-10);
    EXPECT_NEAR(S[1], 0.7184380677507587, 1e-10);
    EXPECT_NEAR(S[2], 0.22995629454444302, 1e-10);

    Core::LinAlg::Tensor<double, 3, 3> original_matrix =
        Q * Core::LinAlg::make_rectangular_diagonal_matrix<3, 3>(S) * VT;

    FOUR_C_EXPECT_NEAR(original_matrix, t, 1e-10);
  }

  TEST(TensorSVDTest, svd3x2)
  {
    Core::LinAlg::Tensor<double, 3, 2> t = {{{0.771320643266746, 0.0207519493594015},
        {0.6336482349262754, 0.7488038825386119}, {0.4985070123025904, 0.22479664553084766}}};

    const auto& [Q, S, VT] = Core::LinAlg::svd(t);
    EXPECT_NEAR(S[0], 1.2684609492672299, 1e-10);
    EXPECT_NEAR(S[1], 0.4976274827645476, 1e-10);

    Core::LinAlg::Tensor<double, 3, 2> original_matrix =
        Q * Core::LinAlg::make_rectangular_diagonal_matrix<3, 2>(S) * VT;

    FOUR_C_EXPECT_NEAR(original_matrix, t, 1e-10);
  }
}  // namespace

FOUR_C_NAMESPACE_CLOSE