// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_config.hpp"

#include "4C_linalg_tensor_generators.hpp"

#include "4C_linalg_symmetric_tensor.hpp"
#include "4C_unittest_utils_assertions_test.hpp"

#include <array>

FOUR_C_NAMESPACE_OPEN

namespace
{
  TEST(TensorGeneratorTest, diagonal)
  {
    {
      Core::LinAlg::SymmetricTensor<double, 2, 2> t =
          Core::LinAlg::TensorGenerators::diagonal(std::array{1.0, 2.0});
      EXPECT_EQ(t(0, 0), 1.0);
      EXPECT_EQ(t(1, 1), 2.0);
      EXPECT_EQ(t(0, 1), 0.0);
      EXPECT_EQ(t(1, 0), 0.0);
    }

    {
      Core::LinAlg::SymmetricTensor<double, 2, 2> t =
          Core::LinAlg::TensorGenerators::diagonal<2, 2>(2.0);
      EXPECT_EQ(t(0, 0), 2.0);
      EXPECT_EQ(t(1, 1), 2.0);
      EXPECT_EQ(t(0, 1), 0.0);
      EXPECT_EQ(t(1, 0), 0.0);
    }
  }

  TEST(TensorGeneratorTest, full)
  {
    {
      Core::LinAlg::SymmetricTensor<double, 2, 2> t =
          Core::LinAlg::TensorGenerators::full<2, 2>(1.2);
      EXPECT_EQ(t(0, 0), 1.2);
      EXPECT_EQ(t(1, 1), 1.2);
      EXPECT_EQ(t(0, 1), 1.2);
      EXPECT_EQ(t(1, 0), 1.2);
    }
    {
      Core::LinAlg::Tensor<double, 2> t = Core::LinAlg::TensorGenerators::full<2>(1.2);
      EXPECT_EQ(t(0), 1.2);
      EXPECT_EQ(t(1), 1.2);
    }
  }



  TEST(TensorGeneratorTest, identity)
  {
    {
      Core::LinAlg::SymmetricTensor<double, 2, 2> t =
          Core::LinAlg::TensorGenerators::identity<double, 2, 2>;
      EXPECT_EQ(t(0, 0), 1.0);
      EXPECT_EQ(t(1, 1), 1.0);
      EXPECT_EQ(t(0, 1), 0.0);
      EXPECT_EQ(t(1, 0), 0.0);
    }
    {
      Core::LinAlg::Tensor<double, 2, 2, 2, 2> t =
          Core::LinAlg::TensorGenerators::identity<double, 2, 2, 2, 2>;
      EXPECT_EQ(t(0, 0, 0, 0), 1.0);
      EXPECT_EQ(t(0, 0, 0, 1), 0.0);
      EXPECT_EQ(t(0, 0, 1, 0), 0.0);
      EXPECT_EQ(t(0, 0, 1, 1), 0.0);
      EXPECT_EQ(t(0, 1, 0, 0), 0.0);
      EXPECT_EQ(t(0, 1, 0, 1), 1.0);
      EXPECT_EQ(t(0, 1, 1, 0), 0.0);
      EXPECT_EQ(t(0, 1, 1, 1), 0.0);
      EXPECT_EQ(t(1, 0, 0, 0), 0.0);
      EXPECT_EQ(t(1, 0, 0, 1), 0.0);
      EXPECT_EQ(t(1, 0, 1, 0), 1.0);
      EXPECT_EQ(t(1, 0, 1, 1), 0.0);
      EXPECT_EQ(t(1, 1, 0, 0), 0.0);
      EXPECT_EQ(t(1, 1, 0, 1), 0.0);
      EXPECT_EQ(t(1, 1, 1, 0), 0.0);
      EXPECT_EQ(t(1, 1, 1, 1), 1.0);
    }

    {
      Core::LinAlg::Tensor<double, 3, 3> t = {{
          {1.1, 1.2, 1.3},
          {2.1, 2.2, 2.3},
          {3.1, 3.2, 3.3},
      }};

      Core::LinAlg::Tensor<double, 3, 3> t2 =
          Core::LinAlg::ddot(Core::LinAlg::TensorGenerators::identity<double, 3, 3, 3, 3>, t);

      FOUR_C_EXPECT_NEAR(t, t2, 1e-15);
    }
  }

  TEST(TensorGeneratorTest, symmetricIdentity)
  {
    {
      Core::LinAlg::SymmetricTensor<double, 2, 2, 2, 2> t =
          Core::LinAlg::TensorGenerators::symmetric_identity<double, 2, 2, 2, 2>;
      EXPECT_EQ(t(0, 0, 0, 0), 1.0);
      EXPECT_EQ(t(0, 0, 0, 1), 0.0);
      EXPECT_EQ(t(0, 0, 1, 0), 0.0);
      EXPECT_EQ(t(0, 0, 1, 1), 0.0);
      EXPECT_EQ(t(0, 1, 0, 0), 0.0);
      EXPECT_EQ(t(0, 1, 0, 1), 0.5);
      EXPECT_EQ(t(0, 1, 1, 0), 0.5);
      EXPECT_EQ(t(0, 1, 1, 1), 0.0);
      EXPECT_EQ(t(1, 0, 0, 0), 0.0);
      EXPECT_EQ(t(1, 0, 0, 1), 0.5);
      EXPECT_EQ(t(1, 0, 1, 0), 0.5);
      EXPECT_EQ(t(1, 0, 1, 1), 0.0);
      EXPECT_EQ(t(1, 1, 0, 0), 0.0);
      EXPECT_EQ(t(1, 1, 0, 1), 0.0);
      EXPECT_EQ(t(1, 1, 1, 0), 0.0);
      EXPECT_EQ(t(1, 1, 1, 1), 1.0);
    }

    {
      Core::LinAlg::Tensor<double, 3, 3> t = {{
          {1.1, 1.2, 1.3},
          {2.1, 2.2, 2.3},
          {3.1, 3.2, 3.3},
      }};

      Core::LinAlg::SymmetricTensor<double, 3, 3> t2 = Core::LinAlg::ddot(
          Core::LinAlg::TensorGenerators::symmetric_identity<double, 3, 3, 3, 3>, t);

      FOUR_C_EXPECT_NEAR(0.5 * (t + Core::LinAlg::transpose(t)), t2, 1e-15);
    }
  }
}  // namespace

FOUR_C_NAMESPACE_CLOSE