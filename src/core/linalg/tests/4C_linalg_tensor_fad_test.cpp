// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_config.hpp"

#include "4C_linalg_tensor_fad.hpp"

#include "4C_linalg_symmetric_tensor.hpp"
#include "4C_linalg_tensor.hpp"
#include "4C_linalg_tensor_generators.hpp"
#include "4C_unittest_utils_assertions_test.hpp"

#include <Sacado.hpp>

#include <type_traits>


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
      EXPECT_NEAR(result.dx(1), 1.4 * 2, 1e-10);
      EXPECT_NEAR(result.dx(2), 2.4, 1e-10);
    }
  }

  TEST(TensorFAD, TensorDerivativeDerivWrtItself)
  {
    const Core::LinAlg::Tensor<double, 2> t = {{1.1, 1.2}};

    const auto fad_t = Core::LinAlg::make_auto_diff_tensor(t);

    const auto dt_dt =
        Core::LinAlg::extract_derivative<std::remove_cvref_t<decltype(fad_t)>>(fad_t);

    // resulting tensor should be the identity tensor
    constexpr auto id = Core::LinAlg::TensorGenerators::identity<double, 2, 2>;
    FOUR_C_EXPECT_NEAR(dt_dt, id, 1e-10);
  }

  TEST(TensorFAD, TensorDerivativeOperation)
  {
    const Core::LinAlg::Tensor<double, 2> t = {{1.1, 1.2}};

    const auto fad_t = Core::LinAlg::make_auto_diff_tensor(t);

    const auto dt_dt =
        Core::LinAlg::extract_derivative<std::remove_cvref_t<decltype(fad_t)>>(2.5 * fad_t);

    // resulting tensor should be the identity tensor scaled with 2.5
    constexpr auto scaled_id = 2.5 * Core::LinAlg::TensorGenerators::identity<double, 2, 2>;
    FOUR_C_EXPECT_NEAR(dt_dt, scaled_id, 1e-10);
  }

  TEST(TensorFAD, TensorDerivativeDetAndTr)
  {
    const Core::LinAlg::Tensor<double, 2, 2> t = {{
        {1.1, 1.2},
        {1.3, 1.4},
    }};

    const auto fad_t = Core::LinAlg::make_auto_diff_tensor(t);

    const auto ddet_t_dt = Core::LinAlg::extract_derivative<std::remove_cvref_t<decltype(fad_t)>>(
        Core::LinAlg::det(fad_t));

    auto deriv_det = Core::LinAlg::det(t) * Core::LinAlg::transpose(Core::LinAlg::inv(t));

    // compare against the analytical derivative
    FOUR_C_EXPECT_NEAR(ddet_t_dt, deriv_det, 1e-10);

    const auto dtr_t_dt = Core::LinAlg::extract_derivative<std::remove_cvref_t<decltype(fad_t)>>(
        Core::LinAlg::trace(fad_t));


    constexpr auto id = Core::LinAlg::TensorGenerators::identity<double, 2, 2>;
    FOUR_C_EXPECT_NEAR(dtr_t_dt, id, 1e-10);
  }

  TEST(TensorFAD, TensorDerivativeSymmetricDetAndTr)
  {
    const Core::LinAlg::SymmetricTensor<double, 2, 2> t =
        Core::LinAlg::assume_symmetry(Core::LinAlg::Tensor<double, 2, 2>{{
            {1.1, 1.2},
            {1.2, 1.4},
        }});

    const auto fad_t = Core::LinAlg::make_auto_diff_tensor(t);

    const auto ddet_t_dt = Core::LinAlg::extract_derivative<std::remove_cvref_t<decltype(fad_t)>>(
        Core::LinAlg::det(fad_t));

    auto deriv_det_analytical = Core::LinAlg::det(t) * Core::LinAlg::inv(t);

    // compare against the analytical derivative
    FOUR_C_EXPECT_NEAR(ddet_t_dt, deriv_det_analytical, 1e-10);

    const auto dtr_t_dt = Core::LinAlg::extract_derivative<std::remove_cvref_t<decltype(fad_t)>>(
        Core::LinAlg::trace(fad_t));


    constexpr auto id = Core::LinAlg::TensorGenerators::identity<double, 2, 2>;
    FOUR_C_EXPECT_NEAR(dtr_t_dt, id, 1e-10);
  }

  TEST(TensorFAD, TensorDerivativeTensorWrtScalar)
  {
    const Core::LinAlg::Tensor<double, 2, 2> t = {{{1.1, 1.2}, {1.3, 1.4}}};

    Sacado::Fad::DFad<double> a = Sacado::Fad::DFad<double>(1, 0, 2.1);

    const auto dat_dt =
        Core::LinAlg::extract_derivative<std::remove_cvref_t<decltype(a)>>(2 * a * t);

    FOUR_C_EXPECT_NEAR(dat_dt, 2 * t, 1e-10);
  }


  TEST(TensorFAD, TensorDerivativeSymmetricTensorWrtScalar)
  {
    const Core::LinAlg::SymmetricTensor<double, 2, 2> t =
        Core::LinAlg::assume_symmetry(Core::LinAlg::Tensor<double, 2, 2>{{
            {1.1, 1.2},
            {1.2, 1.4},
        }});

    Sacado::Fad::DFad<double> a = Sacado::Fad::DFad<double>(1, 0, 2.1);

    const auto dat_dt =
        Core::LinAlg::extract_derivative<std::remove_cvref_t<decltype(a)>>(2 * a * t);

    FOUR_C_EXPECT_NEAR(dat_dt, 2 * t, 1e-10);
  }


  TEST(TensorFAD, TensorDerivativeSymmetricTensorWrtSymmetricTensor)
  {
    const Core::LinAlg::SymmetricTensor<double, 2, 2> t =
        Core::LinAlg::assume_symmetry(Core::LinAlg::Tensor<double, 2, 2>{{
            {1.1, 1.2},
            {1.2, 1.4},
        }});

    const auto fad_t = Core::LinAlg::make_auto_diff_tensor(t);

    const auto dt_dt =
        Core::LinAlg::extract_derivative<std::remove_cvref_t<decltype(fad_t)>>(fad_t);

    auto symmetric_identity =
        Core::LinAlg::TensorGenerators::symmetric_identity<double, 2, 2, 2, 2>;
    FOUR_C_EXPECT_NEAR(dt_dt, symmetric_identity, 1e-10);



    const auto dscaled_t_dt =
        Core::LinAlg::extract_derivative<std::remove_cvref_t<decltype(fad_t)>>(2.5 * fad_t);
    FOUR_C_EXPECT_NEAR(dscaled_t_dt, 2.5 * symmetric_identity, 1e-10);
  }
}  // namespace

FOUR_C_NAMESPACE_CLOSE