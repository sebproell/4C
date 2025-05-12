// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_config.hpp"

#include "4C_linalg_tensor.hpp"

#include <type_traits>

FOUR_C_NAMESPACE_OPEN

namespace
{
  namespace math = Core::LinAlg;

  enum class StorageType
  {
    owning,
    view
  };

  template <StorageType storage_type, typename T, std::size_t... n>
  class TensorHolder;

  template <typename T, std::size_t... n>
  class TensorHolder<StorageType::owning, T, n...>
  {
   public:
    template <typename... U>
    TensorHolder(U&&... u) : tens_(std::forward<U>(u)...)
    {
    }

    math::Tensor<T, n...>& get_tensor() { return tens_; };

   private:
    math::Tensor<T, n...> tens_;
  };

  template <typename T, std::size_t... n>
  class TensorHolder<StorageType::view, T, n...>
  {
   public:
    template <typename... U>
    TensorHolder(U&&... u) : tens_(std::forward<U>(u)...), tens_view_(tens_)
    {
    }

    math::TensorView<T, n...>& get_tensor() { return tens_view_; };


   private:
    math::Tensor<T, n...> tens_;
    math::TensorView<T, n...> tens_view_;
  };

  template <typename T>
  class TensorOperationsTest : public testing::Test
  {
   public:
    static constexpr StorageType STORAGE_TYPE = T();
  };

  using MyTypes = ::testing::Types<std::integral_constant<StorageType, StorageType::owning>,
      std::integral_constant<StorageType, StorageType::view>>;
  TYPED_TEST_SUITE(TensorOperationsTest, MyTypes);


  TYPED_TEST(TensorOperationsTest, Determinant)
  {
    TensorHolder<TestFixture::STORAGE_TYPE, double, 3, 3> t{};
    EXPECT_EQ(math::det(t.get_tensor()), 0.0);

    TensorHolder<TestFixture::STORAGE_TYPE, double, 2, 2> t2{
        math::Tensor<double, 2, 2>{{{1.0, 2.0}, {3.0, 4.0}}}};
    EXPECT_EQ(math::det(t2.get_tensor()), -2.0);
  }


  TYPED_TEST(TensorOperationsTest, Trace)
  {
    TensorHolder<TestFixture::STORAGE_TYPE, double, 3, 3> t{};
    EXPECT_EQ(math::trace(t.get_tensor()), 0.0);

    TensorHolder<TestFixture::STORAGE_TYPE, double, 2, 2> t2{
        math::Tensor<double, 2, 2>{{{1.0, 2.0}, {3.0, 4.0}}}};
    EXPECT_EQ(math::trace(t2.get_tensor()), 5.0);
  }

  TYPED_TEST(TensorOperationsTest, Inverse)
  {
    TensorHolder<TestFixture::STORAGE_TYPE, double, 2, 2> t{
        math::Tensor<double, 2, 2>{{{1.0, 2.0}, {3.0, 4.0}}}};

    math::Tensor<double, 2, 2> t_inv = math::inv(t.get_tensor());

    EXPECT_EQ(t_inv(0, 0), -2.0);
    EXPECT_EQ(t_inv(1, 1), -0.5);
    EXPECT_EQ(t_inv(0, 1), 1.0);
    EXPECT_EQ(t_inv(1, 0), 1.5);
  }

  TYPED_TEST(TensorOperationsTest, Transpose)
  {
    TensorHolder<TestFixture::STORAGE_TYPE, double, 3, 2> t{math::Tensor<double, 3, 2>{{
        {1.0, 2.0},
        {3.0, 4.0},
        {5.0, 6.0},
    }}};

    math::Tensor<double, 2, 3> t_t = math::transpose(t.get_tensor());

    EXPECT_EQ(t_t.at(0, 0), 1.0);
    EXPECT_EQ(t_t.at(0, 1), 3.0);
    EXPECT_EQ(t_t.at(0, 2), 5.0);
    EXPECT_EQ(t_t.at(1, 0), 2.0);
    EXPECT_EQ(t_t.at(1, 1), 4.0);
    EXPECT_EQ(t_t.at(1, 2), 6.0);
  }

  TYPED_TEST(TensorOperationsTest, Vec_dot_Vec)
  {
    TensorHolder<TestFixture::STORAGE_TYPE, double, 2> a = math::Tensor<double, 2>{{2.0, 3.0}};
    TensorHolder<TestFixture::STORAGE_TYPE, double, 2> b = math::Tensor<double, 2>{{1.0, 2.0}};

    double axb = math::dot(a.get_tensor(), b.get_tensor());

    EXPECT_EQ(axb, 8.0);


    double axb2 = a.get_tensor() * b.get_tensor();

    EXPECT_EQ(axb2, 8.0);
  }

  TYPED_TEST(TensorOperationsTest, Mat_dot_Vec)
  {
    TensorHolder<TestFixture::STORAGE_TYPE, double, 3, 2> A = math::Tensor<double, 3, 2>{{
        {1.0, 2.0},
        {3.0, 4.0},
        {5.0, 6.0},
    }};
    TensorHolder<TestFixture::STORAGE_TYPE, double, 2> b = math::Tensor<double, 2>{{1.0, 2.0}};

    math::Tensor<double, 3> Axb = math::dot(A.get_tensor(), b.get_tensor());

    EXPECT_EQ(Axb.at(0), 5.0);
    EXPECT_EQ(Axb.at(1), 11.0);
    EXPECT_EQ(Axb.at(2), 17.0);
  }

  TYPED_TEST(TensorOperationsTest, Mat_dot_Mat)
  {
    TensorHolder<TestFixture::STORAGE_TYPE, double, 3, 2> A = math::Tensor<double, 3, 2>{{
        {1.0, 2.0},
        {3.0, 4.0},
        {5.0, 6.0},
    }};
    TensorHolder<TestFixture::STORAGE_TYPE, double, 2, 3> B = math::Tensor<double, 2, 3>{{
        {1.0, 2.0, 3.0},
        {3.0, 4.0, 5.0},
    }};

    math::Tensor<double, 3, 3> AxB = math::dot(A.get_tensor(), B.get_tensor());

    EXPECT_EQ(AxB.at(0, 0), 7.0);
    EXPECT_EQ(AxB.at(0, 1), 10.0);
    EXPECT_EQ(AxB.at(0, 2), 13.0);

    EXPECT_EQ(AxB.at(1, 0), 15.0);
    EXPECT_EQ(AxB.at(1, 1), 22.0);
    EXPECT_EQ(AxB.at(1, 2), 29.0);

    EXPECT_EQ(AxB.at(2, 0), 23.0);
    EXPECT_EQ(AxB.at(2, 1), 34.0);
    EXPECT_EQ(AxB.at(2, 2), 45.0);
  }



  TYPED_TEST(TensorOperationsTest, Vec_dot_Mat)
  {
    TensorHolder<TestFixture::STORAGE_TYPE, double, 3> a = math::Tensor<double, 3>{{1.0, 2.0, 3.0}};
    TensorHolder<TestFixture::STORAGE_TYPE, double, 3, 2> B = math::Tensor<double, 3, 2>{{
        {1.0, 2.0},
        {3.0, 4.0},
        {5.0, 6.0},
    }};

    math::Tensor<double, 2> a_dot_B = math::dot(a.get_tensor(), B.get_tensor());

    EXPECT_EQ(a_dot_B.at(0), 22.0);
    EXPECT_EQ(a_dot_B.at(1), 28.0);
  }

  TYPED_TEST(TensorOperationsTest, Mat_ddot_Mat)
  {
    TensorHolder<TestFixture::STORAGE_TYPE, double, 3, 2> A = math::Tensor<double, 3, 2>{{
        {1.0, 2.0},
        {3.0, 4.0},
        {5.0, 6.0},
    }};
    TensorHolder<TestFixture::STORAGE_TYPE, double, 3, 2> B = math::Tensor<double, 3, 2>{{
        {2.0, 3.0},
        {4.0, 5.0},
        {6.0, 7.0},
    }};

    double A_ddot_B = math::ddot(A.get_tensor(), B.get_tensor());

    EXPECT_EQ(A_ddot_B, 112.0);
  }

  TYPED_TEST(TensorOperationsTest, Mat_dyadic_Mat)
  {
    TensorHolder<TestFixture::STORAGE_TYPE, double, 3, 2> A = math::Tensor<double, 3, 2>{{
        {1.0, 2.0},
        {3.0, 4.0},
        {5.0, 6.0},
    }};
    TensorHolder<TestFixture::STORAGE_TYPE, double, 2, 3> B = math::Tensor<double, 2, 3>{{
        {1.0, 2.0, 3.0},
        {3.0, 4.0, 5.0},
    }};

    math::Tensor<double, 3, 2, 2, 3> AxB = math::dyadic(A.get_tensor(), B.get_tensor());

    for (std::size_t i = 0; i < A.get_tensor().template extent<0>(); ++i)
    {
      for (std::size_t j = 0; j < A.get_tensor().template extent<1>(); ++j)
      {
        for (std::size_t k = 0; k < B.get_tensor().template extent<0>(); ++k)
        {
          for (std::size_t l = 0; l < B.get_tensor().template extent<1>(); ++l)
          {
            EXPECT_EQ(AxB.at(i, j, k, l), A.get_tensor()(i, j) * B.get_tensor()(k, l));
          }
        }
      }
    }
  }

  TYPED_TEST(TensorOperationsTest, Mat_dyadic_Vec)
  {
    TensorHolder<TestFixture::STORAGE_TYPE, double, 3, 2> A = math::Tensor<double, 3, 2>{{
        {1.0, 2.0},
        {3.0, 4.0},
        {5.0, 6.0},
    }};
    TensorHolder<TestFixture::STORAGE_TYPE, double, 2> b = math::Tensor<double, 2>{
        {1.0, 2.0},
    };

    math::Tensor<double, 3, 2, 2> Axb = math::dyadic(A.get_tensor(), b.get_tensor());

    for (std::size_t i = 0; i < A.get_tensor().template extent<0>(); ++i)
    {
      for (std::size_t j = 0; j < A.get_tensor().template extent<1>(); ++j)
      {
        for (std::size_t k = 0; k < b.get_tensor().template extent<0>(); ++k)
        {
          EXPECT_EQ(Axb.at(i, j, k), A.get_tensor()(i, j) * b.get_tensor()(k));
        }
      }
    }
  }

  TYPED_TEST(TensorOperationsTest, Vec_dyadic_Vec)
  {
    TensorHolder<TestFixture::STORAGE_TYPE, double, 3> a = math::Tensor<double, 3>{
        {1.0, 2.0, 3.0},
    };
    TensorHolder<TestFixture::STORAGE_TYPE, double, 2> b = math::Tensor<double, 2>{
        {1.0, 2.0},
    };

    math::Tensor<double, 3, 2> axb = math::dyadic(a.get_tensor(), b.get_tensor());

    for (std::size_t i = 0; i < a.get_tensor().template extent<0>(); ++i)
    {
      for (std::size_t j = 0; j < b.get_tensor().template extent<0>(); ++j)
      {
        EXPECT_EQ(axb.at(i, j), a.get_tensor()(i) * b.get_tensor()(j));
      }
    }
  }

  TYPED_TEST(TensorOperationsTest, Addition)
  {
    TensorHolder<TestFixture::STORAGE_TYPE, double, 3> a = math::Tensor<double, 3>{
        {1.0, 2.0, 3.0},
    };
    TensorHolder<TestFixture::STORAGE_TYPE, double, 3> b = math::Tensor<double, 3>{
        {2.0, 3.0, 4.0},
    };

    math::Tensor<double, 3> a_plus_b = math::add(a.get_tensor(), b.get_tensor());

    EXPECT_EQ(a_plus_b(0), 3.0);
    EXPECT_EQ(a_plus_b(1), 5.0);
    EXPECT_EQ(a_plus_b(2), 7.0);

    math::Tensor<double, 3> a_plus_b2 = a.get_tensor() + b.get_tensor();

    EXPECT_EQ(a_plus_b(0), 3.0);
    EXPECT_EQ(a_plus_b(1), 5.0);
    EXPECT_EQ(a_plus_b(2), 7.0);
  }

  TYPED_TEST(TensorOperationsTest, Subtract)
  {
    TensorHolder<TestFixture::STORAGE_TYPE, double, 3> a = math::Tensor<double, 3>{
        {1.0, 2.0, 3.0},
    };
    TensorHolder<TestFixture::STORAGE_TYPE, double, 3> b = math::Tensor<double, 3>{
        {2.0, 4.0, 6.0},
    };

    math::Tensor<double, 3> a_negative_b = math::subtract(a.get_tensor(), b.get_tensor());

    EXPECT_EQ(a_negative_b(0), -1.0);
    EXPECT_EQ(a_negative_b(1), -2.0);
    EXPECT_EQ(a_negative_b(2), -3.0);

    math::Tensor<double, 3> a_negative_b2 = a.get_tensor() - b.get_tensor();

    EXPECT_EQ(a_negative_b2(0), -1.0);
    EXPECT_EQ(a_negative_b2(1), -2.0);
    EXPECT_EQ(a_negative_b2(2), -3.0);

    a.get_tensor() -= b.get_tensor();

    EXPECT_EQ(a.get_tensor()(0), -1.0);
    EXPECT_EQ(a.get_tensor()(1), -2.0);
    EXPECT_EQ(a.get_tensor()(2), -3.0);
  }

  TYPED_TEST(TensorOperationsTest, Division)
  {
    TensorHolder<TestFixture::STORAGE_TYPE, double, 3> a = math::Tensor<double, 3>{
        {1.0, 2.0, 3.0},
    };

    math::Tensor<double, 3> a_divide = a.get_tensor() / 3.0;

    EXPECT_EQ(a_divide(0), 1.0 / 3.0);
    EXPECT_EQ(a_divide(1), 2.0 / 3.0);
    EXPECT_EQ(a_divide(2), 1.0);

    a.get_tensor() /= 3.0;

    EXPECT_EQ(a.get_tensor()(0), 1.0 / 3.0);
    EXPECT_EQ(a.get_tensor()(1), 2.0 / 3.0);
    EXPECT_EQ(a.get_tensor()(2), 1.0);
  }

  TEST(TensorOperationsTest, Equality)
  {
    Core::LinAlg::Tensor<double, 3> a = {
        {1.0, 2.0, 3.0},
    };
    Core::LinAlg::Tensor<double, 3> a1 = {
        {1.0, 2.0, 3.0},
    };
    Core::LinAlg::Tensor<double, 3> b = {
        {1.0, 2.0, 3.1},
    };

    Core::LinAlg::TensorView<double, 3> a_view = a;
    Core::LinAlg::TensorView<double, 3> b_view = b;

    EXPECT_TRUE(a == a1);
    EXPECT_TRUE(a1 == a);
    EXPECT_TRUE(a_view == a1);
    EXPECT_TRUE(a1 == a_view);
    EXPECT_FALSE(a_view == b);
    EXPECT_FALSE(b == a_view);
    EXPECT_FALSE(a == b_view);
    EXPECT_FALSE(b_view == a);

    EXPECT_FALSE(a != a1);
    EXPECT_FALSE(a1 != a);
    EXPECT_FALSE(a_view != a1);
    EXPECT_FALSE(a1 != a_view);
    EXPECT_TRUE(a_view != b);
    EXPECT_TRUE(b != a_view);
    EXPECT_TRUE(a != b_view);
    EXPECT_TRUE(b_view != a);
  }

  TEST(TensorOperationsTest, ReorderAxis)
  {
    Core::LinAlg::Tensor<double, 3, 2> A = {{
        {1.0, 2.0},
        {4.0, 5.0},
        {6.0, 7.0},
    }};

    Core::LinAlg::Tensor<double, 2, 3> B = Core::LinAlg::reorder_axis<1, 0>(A);

    EXPECT_DOUBLE_EQ(B(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(B(0, 1), 4.0);
    EXPECT_DOUBLE_EQ(B(0, 2), 6.0);
    EXPECT_DOUBLE_EQ(B(1, 0), 2.0);
    EXPECT_DOUBLE_EQ(B(1, 1), 5.0);
    EXPECT_DOUBLE_EQ(B(1, 2), 7.0);


    Core::LinAlg::Tensor<double, 4, 3, 2> A2 = {{
        {{0.00, 0.01}, {0.10, 0.11}},
        {{1.00, 1.01}, {1.10, 1.11}},
        {{2.00, 2.01}, {2.10, 2.11}},
        {{3.00, 3.01}, {3.10, 3.11}},
    }};

    Core::LinAlg::Tensor<double, 3, 4, 2> B2 = Core::LinAlg::reorder_axis<1, 0, 2>(A2);

    for (int i = 0; i < 4; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        for (int k = 0; k < 2; ++k)
        {
          EXPECT_DOUBLE_EQ(B2.at(j, i, k), A2.at(i, j, k));
        }
      }
    }

    Core::LinAlg::Tensor<double, 2, 4, 3> B3 = Core::LinAlg::reorder_axis<2, 0, 1>(A2);

    for (int i = 0; i < 4; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        for (int k = 0; k < 2; ++k)
        {
          EXPECT_DOUBLE_EQ(B3.at(k, i, j), A2.at(i, j, k));
        }
      }
    }
  }

}  // namespace

FOUR_C_NAMESPACE_CLOSE