// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_config.hpp"

#include "4C_linalg_symmetric_tensor.hpp"
#include "4C_linalg_tensor.hpp"

#include <Sacado.hpp>


FOUR_C_NAMESPACE_OPEN

namespace
{
  TEST(TensorFAD, TensorWithFADTypeScale)
  {
    using FadType = Sacado::Fad::DFad<double>;

    const Core::LinAlg::Tensor<FadType, 2> t = {{FadType(2, 0, 1.1), FadType(2, 1, 1.2)}};
    {
      auto scaled_t = 2.3 * t;  // scaled with double
      static_assert(std::is_same_v<typename decltype(scaled_t)::value_type, FadType>);
      EXPECT_NEAR(scaled_t(0).val(), 2.53, 1e-10);
      EXPECT_NEAR(scaled_t(1).val(), 2.76, 1e-10);
      EXPECT_NEAR(scaled_t(0).dx(0), 2.3, 1e-10);
      EXPECT_NEAR(scaled_t(1).dx(0), 0.0, 1e-10);
      EXPECT_NEAR(scaled_t(0).dx(1), 0.0, 1e-10);
      EXPECT_NEAR(scaled_t(1).dx(1), 2.3, 1e-10);
    }

    {
      FadType a(2.3);
      auto scaled_t = 2 * a * t;  // scaled with Fad-operation
      static_assert(std::is_same_v<typename decltype(scaled_t)::value_type, FadType>);
      EXPECT_NEAR(scaled_t(0).val(), 5.06, 1e-10);
      EXPECT_NEAR(scaled_t(1).val(), 5.52, 1e-10);
      EXPECT_NEAR(scaled_t(0).dx(0), 4.6, 1e-10);
      EXPECT_NEAR(scaled_t(1).dx(0), 0.0, 1e-10);
      EXPECT_NEAR(scaled_t(0).dx(1), 0.0, 1e-10);
      EXPECT_NEAR(scaled_t(1).dx(1), 4.6, 1e-10);
    }

    {
      FadType a(2);
      auto scaled_t = 2 * a * t;  // scaled with Fad-operation
      static_assert(std::is_same_v<typename decltype(scaled_t)::value_type, FadType>);
      EXPECT_NEAR(scaled_t(0).val(), 4.4, 1e-10);
      EXPECT_NEAR(scaled_t(1).val(), 4.8, 1e-10);
      EXPECT_NEAR(scaled_t(0).dx(0), 4, 1e-10);
      EXPECT_NEAR(scaled_t(1).dx(0), 0.0, 1e-10);
      EXPECT_NEAR(scaled_t(0).dx(1), 0.0, 1e-10);
      EXPECT_NEAR(scaled_t(1).dx(1), 4, 1e-10);
    }
  }
  TEST(TensorFAD, TensorWithFADTypeAdd)
  {
    using FadType = Sacado::Fad::DFad<double>;

    const Core::LinAlg::Tensor<FadType, 2> t = {{FadType(2, 0, 1.1), FadType(2, 1, 1.2)}};
    const Core::LinAlg::Tensor<FadType, 2> t2 = {{FadType(1.2), FadType(1.3)}};
    {
      auto result = t + t2;
      static_assert(std::is_same_v<typename decltype(result)::value_type, FadType>);
      EXPECT_NEAR(result(0).val(), 2.3, 1e-10);
      EXPECT_NEAR(result(1).val(), 2.5, 1e-10);
      EXPECT_NEAR(result(0).dx(0), 1.0, 1e-10);
      EXPECT_NEAR(result(1).dx(0), 0.0, 1e-10);
      EXPECT_NEAR(result(0).dx(1), 0.0, 1e-10);
      EXPECT_NEAR(result(1).dx(1), 1.0, 1e-10);
    }
  }
  TEST(TensorFAD, TensorWithFADTypeSubtract)
  {
    using FadType = Sacado::Fad::DFad<double>;

    const Core::LinAlg::Tensor<FadType, 2> t = {{FadType(2, 0, 1.1), FadType(2, 1, 1.2)}};
    const Core::LinAlg::Tensor<FadType, 2> t2 = {{FadType(1.2), FadType(1.4)}};
    {
      auto result = t - t2;
      static_assert(std::is_same_v<typename decltype(result)::value_type, FadType>);
      EXPECT_NEAR(result(0).val(), -0.1, 1e-10);
      EXPECT_NEAR(result(1).val(), -0.2, 1e-10);
      EXPECT_NEAR(result(0).dx(0), 1.0, 1e-10);
      EXPECT_NEAR(result(1).dx(0), 0.0, 1e-10);
      EXPECT_NEAR(result(0).dx(1), 0.0, 1e-10);
      EXPECT_NEAR(result(1).dx(1), 1.0, 1e-10);
    }
  }
  TEST(TensorFAD, TensorWithFADTypeDot)
  {
    using FadType = Sacado::Fad::DFad<double>;

    const Core::LinAlg::Tensor<FadType, 2> t = {{FadType(2, 0, 1.1), FadType(2, 1, 1.2)}};
    const Core::LinAlg::Tensor<FadType, 2> t2 = {{FadType(1.2), FadType(1.4)}};
    {
      auto result = t * t2;
      static_assert(std::is_same_v<decltype(result), FadType>);
      EXPECT_NEAR(result.val(), 3.0, 1e-10);
      EXPECT_NEAR(result.dx(0), 1.2, 1e-10);
      EXPECT_NEAR(result.dx(1), 1.4, 1e-10);
    }
  }
  TEST(TensorFAD, TensorWithFADTypeDDot)
  {
    using FadType = Sacado::Fad::DFad<double>;

    const Core::LinAlg::Tensor<FadType, 2, 2> t = {{
        {FadType(4, 0, 1.1), FadType(4, 1, 1.2)},
        {FadType(4, 2, 2.1), FadType(4, 3, 2.2)},
    }};
    const Core::LinAlg::Tensor<FadType, 2, 2> t2 = {{
        {FadType(1.2), FadType(1.4)},
        {FadType(2.2), FadType(2.4)},
    }};
    {
      auto result = Core::LinAlg::ddot(t, t2);
      static_assert(std::is_same_v<decltype(result), FadType>);
      EXPECT_NEAR(result.val(), 12.9, 1e-10);
      EXPECT_NEAR(result.dx(0), 1.2, 1e-10);
      EXPECT_NEAR(result.dx(1), 1.4, 1e-10);
      EXPECT_NEAR(result.dx(2), 2.2, 1e-10);
      EXPECT_NEAR(result.dx(3), 2.4, 1e-10);
    }
  }

  TEST(TensorFAD, SymmetricTensorWithFADTypeScale)
  {
    using FadType = Sacado::Fad::DFad<double>;

    const Core::LinAlg::SymmetricTensor<FadType, 2, 2> t =
        Core::LinAlg::assume_symmetry(Core::LinAlg::Tensor<FadType, 2, 2>{{
            {FadType(3, 0, 1.1), FadType(3, 1, 1.2)},
            {FadType(3, 1, 1.2), FadType(3, 2, 2.2)},
        }});
    {
      auto result = 2.0 * t;
      static_assert(std::is_same_v<typename decltype(result)::value_type, FadType>);
      EXPECT_NEAR(result(0, 0).val(), 2.2, 1e-10);
      EXPECT_NEAR(result(0, 1).val(), 2.4, 1e-10);
      EXPECT_NEAR(result(1, 1).val(), 4.4, 1e-10);

      EXPECT_NEAR(result(0, 0).dx(0), 2, 1e-10);
      EXPECT_NEAR(result(0, 1).dx(0), 0, 1e-10);
      EXPECT_NEAR(result(1, 1).dx(0), 0, 1e-10);

      EXPECT_NEAR(result(0, 0).dx(1), 0, 1e-10);
      EXPECT_NEAR(result(0, 1).dx(1), 2, 1e-10);
      EXPECT_NEAR(result(1, 1).dx(1), 0, 1e-10);

      EXPECT_NEAR(result(0, 0).dx(2), 0, 1e-10);
      EXPECT_NEAR(result(0, 1).dx(2), 0, 1e-10);
      EXPECT_NEAR(result(1, 1).dx(2), 2, 1e-10);
    }
  }
  TEST(TensorFAD, SymmetricTensorWithFADTypeAdd)
  {
    using FadType = Sacado::Fad::DFad<double>;

    const Core::LinAlg::SymmetricTensor<FadType, 2, 2> t =
        Core::LinAlg::assume_symmetry(Core::LinAlg::Tensor<FadType, 2, 2>{{
            {FadType(3, 0, 1.1), FadType(3, 1, 1.2)},
            {FadType(3, 1, 1.2), FadType(3, 2, 2.2)},
        }});
    const Core::LinAlg::SymmetricTensor<FadType, 2, 2> t2 =
        Core::LinAlg::assume_symmetry(Core::LinAlg::Tensor<FadType, 2, 2>{{
            {FadType(1.2), FadType(1.4)},
            {FadType(1.4), FadType(2.4)},
        }});
    {
      auto result = t + t2;
      static_assert(std::is_same_v<typename decltype(result)::value_type, FadType>);
      EXPECT_NEAR(result(0, 0).val(), 2.3, 1e-10);
      EXPECT_NEAR(result(0, 1).val(), 2.6, 1e-10);
      EXPECT_NEAR(result(1, 1).val(), 4.6, 1e-10);

      EXPECT_NEAR(result(0, 0).dx(0), 1, 1e-10);
      EXPECT_NEAR(result(0, 1).dx(0), 0, 1e-10);
      EXPECT_NEAR(result(1, 1).dx(0), 0, 1e-10);

      EXPECT_NEAR(result(0, 0).dx(1), 0, 1e-10);
      EXPECT_NEAR(result(0, 1).dx(1), 1, 1e-10);
      EXPECT_NEAR(result(1, 1).dx(1), 0, 1e-10);

      EXPECT_NEAR(result(0, 0).dx(2), 0, 1e-10);
      EXPECT_NEAR(result(0, 1).dx(2), 0, 1e-10);
      EXPECT_NEAR(result(1, 1).dx(2), 1, 1e-10);
    }
  }
  TEST(TensorFAD, SymmetricTensorWithFADTypeSubtract)
  {
    using FadType = Sacado::Fad::DFad<double>;

    const Core::LinAlg::SymmetricTensor<FadType, 2, 2> t =
        Core::LinAlg::assume_symmetry(Core::LinAlg::Tensor<FadType, 2, 2>{{
            {FadType(3, 0, 1.1), FadType(3, 1, 1.2)},
            {FadType(3, 1, 1.2), FadType(3, 2, 2.2)},
        }});
    const Core::LinAlg::SymmetricTensor<FadType, 2, 2> t2 =
        Core::LinAlg::assume_symmetry(Core::LinAlg::Tensor<FadType, 2, 2>{{
            {FadType(1.2), FadType(1.4)},
            {FadType(1.4), FadType(2.4)},
        }});
    {
      auto result = t - t2;
      static_assert(std::is_same_v<typename decltype(result)::value_type, FadType>);
      EXPECT_NEAR(result(0, 0).val(), -0.1, 1e-10);
      EXPECT_NEAR(result(0, 1).val(), -0.2, 1e-10);
      EXPECT_NEAR(result(1, 1).val(), -0.2, 1e-10);

      EXPECT_NEAR(result(0, 0).dx(0), 1, 1e-10);
      EXPECT_NEAR(result(0, 1).dx(0), 0, 1e-10);
      EXPECT_NEAR(result(1, 1).dx(0), 0, 1e-10);

      EXPECT_NEAR(result(0, 0).dx(1), 0, 1e-10);
      EXPECT_NEAR(result(0, 1).dx(1), 1, 1e-10);
      EXPECT_NEAR(result(1, 1).dx(1), 0, 1e-10);

      EXPECT_NEAR(result(0, 0).dx(2), 0, 1e-10);
      EXPECT_NEAR(result(0, 1).dx(2), 0, 1e-10);
      EXPECT_NEAR(result(1, 1).dx(2), 1, 1e-10);
    }
  }
  TEST(TensorFAD, SymmetricTensorWithFADTypeDot)
  {
    using FadType = Sacado::Fad::DFad<double>;

    const Core::LinAlg::SymmetricTensor<FadType, 2, 2> t =
        Core::LinAlg::assume_symmetry(Core::LinAlg::Tensor<FadType, 2, 2>{{
            {FadType(3, 0, 1.1), FadType(3, 1, 1.2)},
            {FadType(3, 1, 1.2), FadType(3, 2, 2.2)},
        }});
    const Core::LinAlg::SymmetricTensor<FadType, 2, 2> t2 =
        Core::LinAlg::assume_symmetry(Core::LinAlg::Tensor<FadType, 2, 2>{{
            {FadType(1.2), FadType(1.4)},
            {FadType(1.4), FadType(2.4)},
        }});
    {
      auto result = t * t2;
      static_assert(std::is_same_v<typename decltype(result)::value_type, FadType>);
      EXPECT_NEAR(result(0, 0).val(), 3, 1e-10);
      EXPECT_NEAR(result(0, 1).val(), 4.42, 1e-10);
      EXPECT_NEAR(result(1, 0).val(), 4.52, 1e-10);
      EXPECT_NEAR(result(1, 1).val(), 6.96, 1e-10);

      EXPECT_NEAR(result(0, 0).dx(0), 1.2, 1e-10);
      EXPECT_NEAR(result(0, 1).dx(0), 1.4, 1e-10);
      EXPECT_NEAR(result(1, 0).dx(0), 0, 1e-10);
      EXPECT_NEAR(result(1, 1).dx(0), 0, 1e-10);

      EXPECT_NEAR(result(0, 0).dx(1), 1.4, 1e-10);
      EXPECT_NEAR(result(0, 1).dx(1), 2.4, 1e-10);
      EXPECT_NEAR(result(1, 0).dx(1), 1.2, 1e-10);
      EXPECT_NEAR(result(1, 1).dx(1), 1.4, 1e-10);

      EXPECT_NEAR(result(0, 0).dx(2), 0, 1e-10);
      EXPECT_NEAR(result(0, 1).dx(2), 0, 1e-10);
      EXPECT_NEAR(result(1, 0).dx(2), 1.4, 1e-10);
      EXPECT_NEAR(result(1, 1).dx(2), 2.4, 1e-10);
    }
  }
  TEST(TensorFAD, SymmetricTensorWithFADTypeDDot)
  {
    using FadType = Sacado::Fad::DFad<double>;

    const Core::LinAlg::SymmetricTensor<FadType, 2, 2> t =
        Core::LinAlg::assume_symmetry(Core::LinAlg::Tensor<FadType, 2, 2>{{
            {FadType(3, 0, 1.1), FadType(3, 1, 1.2)},
            {FadType(3, 1, 1.2), FadType(3, 2, 2.2)},
        }});
    const Core::LinAlg::SymmetricTensor<FadType, 2, 2> t2 =
        Core::LinAlg::assume_symmetry(Core::LinAlg::Tensor<FadType, 2, 2>{{
            {FadType(1.2), FadType(1.4)},
            {FadType(1.4), FadType(2.4)},
        }});
    {
      auto result = Core::LinAlg::ddot(t, t2);
      static_assert(std::is_same_v<decltype(result), FadType>);
      EXPECT_NEAR(result.val(), 9.96, 1e-10);
      EXPECT_NEAR(result.dx(0), 1.2, 1e-10);
      EXPECT_NEAR(result.dx(1), 2.8, 1e-10);
      EXPECT_NEAR(result.dx(2), 2.4, 1e-10);
    }
  }

}  // namespace

FOUR_C_NAMESPACE_CLOSE