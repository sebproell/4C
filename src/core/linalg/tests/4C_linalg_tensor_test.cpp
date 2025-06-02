// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_config.hpp"

#include "4C_linalg_tensor.hpp"

#include "4C_linalg_fixedsizematrix.hpp"

FOUR_C_NAMESPACE_OPEN

namespace
{


  TEST(TensorTest, TestConstructability)
  {
    // Tensor should be default constructible
    static_assert(std::is_default_constructible_v<Core::LinAlg::Tensor<double, 3>>);
  }

  TEST(TensorTest, DefaultConstruct1Tensor)
  {
    Core::LinAlg::Tensor<double, 3> t{};
    EXPECT_DOUBLE_EQ(t(0), 0.0);
    EXPECT_DOUBLE_EQ(t(0), 0.0);
    EXPECT_DOUBLE_EQ(t(0), 0.0);
  }

  TEST(TensorTest, DefaultConstruct2Tensor)
  {
    Core::LinAlg::Tensor<double, 2, 2> t{};
    EXPECT_DOUBLE_EQ(t(0, 0), 0.0);
    EXPECT_DOUBLE_EQ(t(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(t(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(t(1, 1), 0.0);
  }

  TEST(TensorTest, InitializationAndIndexing1Tensor)
  {
    Core::LinAlg::Tensor<double, 3> t = {{0.0, 1.0, 2.0}};

    EXPECT_DOUBLE_EQ(t(0), 0.0);
    EXPECT_DOUBLE_EQ(t(1), 1.0);
    EXPECT_DOUBLE_EQ(t(2), 2.0);
  }

  TEST(TensorTest, InitializationAndIndexing2Tensor)
  {
    Core::LinAlg::Tensor<double, 3, 2> t = {{{0.0, 0.1}, {1.0, 1.1}, {2.0, 2.1}}};

    EXPECT_DOUBLE_EQ(t(0, 0), 0.0);
    EXPECT_DOUBLE_EQ(t(0, 1), 0.1);
    EXPECT_DOUBLE_EQ(t(1, 0), 1.0);
    EXPECT_DOUBLE_EQ(t(1, 1), 1.1);
    EXPECT_DOUBLE_EQ(t(2, 0), 2.0);
    EXPECT_DOUBLE_EQ(t(2, 1), 2.1);
  }

  TEST(TensorTest, InitializationAndIndexing3Tensor)
  {
    Core::LinAlg::Tensor<double, 3, 2, 2> t = {{
        {{0.00, 0.01}, {0.10, 0.11}},
        {{1.00, 1.01}, {1.10, 1.11}},
        {{2.00, 2.01}, {2.10, 2.11}},
    }};

    EXPECT_DOUBLE_EQ(t(0, 0, 0), 0.00);
    EXPECT_DOUBLE_EQ(t(0, 0, 1), 0.01);

    EXPECT_DOUBLE_EQ(t(0, 1, 0), 0.10);
    EXPECT_DOUBLE_EQ(t(0, 1, 1), 0.11);

    EXPECT_DOUBLE_EQ(t(1, 0, 0), 1.00);
    EXPECT_DOUBLE_EQ(t(1, 0, 1), 1.01);

    EXPECT_DOUBLE_EQ(t(1, 1, 0), 1.10);
    EXPECT_DOUBLE_EQ(t(1, 1, 1), 1.11);

    EXPECT_DOUBLE_EQ(t(2, 0, 0), 2.00);
    EXPECT_DOUBLE_EQ(t(2, 0, 1), 2.01);

    EXPECT_DOUBLE_EQ(t(2, 1, 0), 2.10);
    EXPECT_DOUBLE_EQ(t(2, 1, 1), 2.11);
  }

  TEST(TensorTest, InitializationAndIndexing4Tensor)
  {
    Core::LinAlg::Tensor<double, 3, 2, 2, 2> t = {{
        {{{0.000, 0.001}, {0.010, 0.011}}, {{0.100, 0.101}, {0.110, 0.111}}},
        {{{1.000, 1.001}, {1.010, 1.011}}, {{1.100, 1.101}, {1.110, 1.111}}},
        {{{2.000, 2.001}, {2.010, 2.011}}, {{2.100, 2.101}, {2.110, 2.111}}},
    }};

    EXPECT_DOUBLE_EQ(t(0, 0, 0, 0), 0.000);
    EXPECT_DOUBLE_EQ(t(0, 0, 0, 1), 0.001);
    EXPECT_DOUBLE_EQ(t(0, 0, 1, 0), 0.010);
    EXPECT_DOUBLE_EQ(t(0, 0, 1, 1), 0.011);
    EXPECT_DOUBLE_EQ(t(0, 1, 0, 0), 0.100);
    EXPECT_DOUBLE_EQ(t(0, 1, 0, 1), 0.101);
    EXPECT_DOUBLE_EQ(t(0, 1, 1, 0), 0.110);
    EXPECT_DOUBLE_EQ(t(0, 1, 1, 1), 0.111);

    EXPECT_DOUBLE_EQ(t(1, 0, 0, 0), 1.000);
    EXPECT_DOUBLE_EQ(t(1, 0, 0, 1), 1.001);
    EXPECT_DOUBLE_EQ(t(1, 0, 1, 0), 1.010);
    EXPECT_DOUBLE_EQ(t(1, 0, 1, 1), 1.011);
    EXPECT_DOUBLE_EQ(t(1, 1, 0, 0), 1.100);
    EXPECT_DOUBLE_EQ(t(1, 1, 0, 1), 1.101);
    EXPECT_DOUBLE_EQ(t(1, 1, 1, 0), 1.110);
    EXPECT_DOUBLE_EQ(t(1, 1, 1, 1), 1.111);

    EXPECT_DOUBLE_EQ(t(2, 0, 0, 0), 2.000);
    EXPECT_DOUBLE_EQ(t(2, 0, 0, 1), 2.001);
    EXPECT_DOUBLE_EQ(t(2, 0, 1, 0), 2.010);
    EXPECT_DOUBLE_EQ(t(2, 0, 1, 1), 2.011);
    EXPECT_DOUBLE_EQ(t(2, 1, 0, 0), 2.100);
    EXPECT_DOUBLE_EQ(t(2, 1, 0, 1), 2.101);
    EXPECT_DOUBLE_EQ(t(2, 1, 1, 0), 2.110);
    EXPECT_DOUBLE_EQ(t(2, 1, 1, 1), 2.111);
  }

  TEST(TensorTest, AlignmentSameAsMatrix)
  {
    Core::LinAlg::Tensor<double, 2, 2> t = {{{0.0, 0.1}, {1.0, 1.1}}};
    Core::LinAlg::Matrix<2, 2> m(t.data());

    EXPECT_DOUBLE_EQ(m(0, 0), 0.0);
    EXPECT_DOUBLE_EQ(m(0, 1), 0.1);
    EXPECT_DOUBLE_EQ(m(1, 0), 1.0);
    EXPECT_DOUBLE_EQ(m(1, 1), 1.1);
  }

  TEST(TensorTest, PrintTensor)
  {
    {
      Core::LinAlg::Tensor<double, 2> t1 = {{0.0, 0.1}};
      std::stringstream ss;
      ss << t1;
      EXPECT_EQ(ss.str(), "Tensor<double, 2>[0, 0.1]");
    }

    {
      Core::LinAlg::Tensor<double, 2, 3> t2 = {{{0.0, 0.1, 0.2}, {1.0, 1.1, 1.2}}};
      std::stringstream ss;
      ss << t2;
      EXPECT_EQ(ss.str(), R"(Tensor<double, 2, 3>[
 [0, 0.1, 0.2],
 [1, 1.1, 1.2]
])");
    }

    {
      Core::LinAlg::Tensor<double, 2, 3, 2> t3 = {
          {{{0.0, 0.01}, {0.1, 0.11}, {0.2, 0.21}}, {{1.0, 1.01}, {1.1, 1.11}, {1.2, 1.21}}}};
      std::stringstream ss;
      ss << t3;
      EXPECT_EQ(ss.str(), R"(Tensor<double, 2, 3, 2>[
 [
  [0, 0.01],
  [0.1, 0.11],
  [0.2, 0.21]
 ],
 [
  [1, 1.01],
  [1.1, 1.11],
  [1.2, 1.21]
 ]
])");
    }
  }


  TEST(TensorViewTest, TestConstructability)
  {
    // TensorView should not be default constructable
    static_assert(!std::is_default_constructible_v<Core::LinAlg::TensorView<double, 3>>);

    // TensorView should be constructable from Tensor&
    static_assert(std::is_constructible_v<Core::LinAlg::TensorView<double, 3>,
        Core::LinAlg::Tensor<double, 3>&>);

    // const TensorView should be constructable from const Tensor&
    static_assert(std::is_constructible_v<Core::LinAlg::TensorView<const double, 3>,
        const Core::LinAlg::Tensor<double, 3>&>);

    // TensorView should not be constructable from Tensor&&
    static_assert(!std::is_constructible_v<Core::LinAlg::TensorView<double, 3>,
        Core::LinAlg::Tensor<double, 3>&&>);

    // TensorView should NOT be constructable from const Tensor&
    static_assert(!std::is_constructible_v<Core::LinAlg::TensorView<double, 3>,
        const Core::LinAlg::Tensor<double, 3>&>);

    // TensorView should be constructible from another TensorView
    static_assert(std::is_constructible_v<Core::LinAlg::TensorView<double, 3>,
        Core::LinAlg::TensorView<double, 3>>);

    // TensorView should NOT be constructible from another const TensorView
    static_assert(!std::is_constructible_v<Core::LinAlg::TensorView<double, 3>,
        Core::LinAlg::TensorView<const double, 3>>);

    // const TensorView should be constructible from another const TensorView
    static_assert(std::is_constructible_v<Core::LinAlg::TensorView<const double, 3>,
        Core::LinAlg::TensorView<const double, 3>>);
  }

  TEST(TensorViewTest, MutableTensorView)
  {
    Core::LinAlg::Tensor<double, 3> t{};

    Core::LinAlg::TensorView<double, 3> t_view(t);
    EXPECT_DOUBLE_EQ(t(0), 0.0);
    EXPECT_DOUBLE_EQ(t(0), 0.0);
    EXPECT_DOUBLE_EQ(t(0), 0.0);

    t_view(1) = 2.0;

    EXPECT_DOUBLE_EQ(t(1), 2.0);

    // this is the actual ultimate test whether they work on the same data
    EXPECT_EQ(t_view.data(), t.data());
  }

  TEST(TensorViewTest, ConstantTensorView)
  {
    const Core::LinAlg::Tensor<double, 3> t{};
    Core::LinAlg::TensorView<const double, 3> t_view(t);

    // this is the actual ultimate test whether they work on the same data
    EXPECT_EQ(t_view.data(), t.data());
  }


  TEST(TensorViewTest, PrintTensor)
  {
    {
      Core::LinAlg::Tensor<double, 2> t1 = {{0.0, 0.1}};
      Core::LinAlg::TensorView<double, 2> t1_view(t1);
      std::stringstream ss;
      ss << t1_view;
      EXPECT_EQ(ss.str(), "TensorView<double, 2>[0, 0.1]");
    }

    {
      Core::LinAlg::Tensor<double, 2, 3> t2 = {{{0.0, 0.1, 0.2}, {1.0, 1.1, 1.2}}};
      Core::LinAlg::TensorView<double, 2, 3> t2_view(t2);
      std::stringstream ss;
      ss << t2_view;

      EXPECT_EQ(ss.str(), R"(TensorView<double, 2, 3>[
 [0, 0.1, 0.2],
 [1, 1.1, 1.2]
])");
    }



    {
      Core::LinAlg::Tensor<double, 2, 3, 2> t3 = {
          {{{0.0, 0.01}, {0.1, 0.11}, {0.2, 0.21}}, {{1.0, 1.01}, {1.1, 1.11}, {1.2, 1.21}}}};
      Core::LinAlg::TensorView<double, 2, 3, 2> t3_view(t3);
      std::stringstream ss;
      ss << t3_view;
      EXPECT_EQ(ss.str(), R"(TensorView<double, 2, 3, 2>[
 [
  [0, 0.01],
  [0.1, 0.11],
  [0.2, 0.21]
 ],
 [
  [1, 1.01],
  [1.1, 1.11],
  [1.2, 1.21]
 ]
])");
    }
  }

}  // namespace

FOUR_C_NAMESPACE_CLOSE