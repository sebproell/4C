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

  TYPED_TEST(TensorOperationsTest, VectorDotVector)
  {
    TensorHolder<TestFixture::STORAGE_TYPE, double, 2> a = math::Tensor<double, 2>{{2.0, 3.0}};
    TensorHolder<TestFixture::STORAGE_TYPE, double, 2> b = math::Tensor<double, 2>{{1.0, 2.0}};

    double axb = math::dot(a.get_tensor(), b.get_tensor());

    EXPECT_EQ(axb, 8.0);


    double axb2 = a.get_tensor() * b.get_tensor();

    EXPECT_EQ(axb2, 8.0);
  }

  TYPED_TEST(TensorOperationsTest, MatrixDotVector)
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

  TYPED_TEST(TensorOperationsTest, MatrixDotMatrix)
  {
    TensorHolder<TestFixture::STORAGE_TYPE, double, 3, 2> A =
        math::Tensor<double, 3, 2>{{{0.9770807259226818, 0.4407738249006665},
            {0.3182728054789512, 0.5197969858753801}, {0.5781364298824675, 0.8539337505004864}}};
    TensorHolder<TestFixture::STORAGE_TYPE, double, 2, 4> B = math::Tensor<double, 2, 4>{
        {{0.06809727353795003, 0.46453080777933253, 0.7819491186191484, 0.7186028103822503},
            {0.5860219800531759, 0.037094413234407875, 0.350656391283133, 0.563190684492745}}};

    math::Tensor<double, 3, 4> AxB = math::dot(A.get_tensor(), B.get_tensor());

    EXPECT_EQ(AxB.at(0, 0), 0.3248396830857161);
    EXPECT_EQ(AxB.at(0, 1), 0.47023434528225583);
    EXPECT_EQ(AxB.at(0, 2), 0.91858757126673);
    EXPECT_EQ(AxB.at(0, 3), 0.95037266777066);

    EXPECT_EQ(AxB.at(1, 0), 0.3262859691827539);
    EXPECT_EQ(AxB.at(1, 1), 0.1671290876153926);
    EXPECT_EQ(AxB.at(1, 2), 0.43114327499162);
    EXPECT_EQ(AxB.at(1, 3), 0.5214565527578386);

    EXPECT_EQ(AxB.at(2, 0), 0.5397934619104899);
    EXPECT_EQ(AxB.at(2, 1), 0.3002383541958349);
    EXPECT_EQ(AxB.at(2, 2), 0.7515105991335884);
    EXPECT_EQ(AxB.at(2, 3), 0.8963779967537278);
  }



  TYPED_TEST(TensorOperationsTest, VectorDotMatrix)
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

  TYPED_TEST(TensorOperationsTest, MatrixDdotMatrix)
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

  TYPED_TEST(TensorOperationsTest, ThreeTensorDdotMatrix)
  {
    TensorHolder<TestFixture::STORAGE_TYPE, double, 4, 3, 2> A = math::Tensor<double, 4, 3, 2>{
        {
            {{0.771320643266746, 0.0207519493594015}, {0.6336482349262754, 0.7488038825386119},
                {0.4985070123025904, 0.22479664553084766}},
            {{0.19806286475962398, 0.7605307121989587}, {0.16911083656253545, 0.08833981417401027},
                {0.6853598183677972, 0.9533933461949365}},
            {{0.003948266327914451, 0.5121922633857766}, {0.8126209616521135, 0.6125260668293881},
                {0.7217553174317995, 0.29187606817063316}},
            {{0.9177741225129434, 0.7145757833976906}, {0.5425443680112613, 0.14217004760152696},
                {0.3733407600514692, 0.6741336150663453}},
        },
    };
    TensorHolder<TestFixture::STORAGE_TYPE, double, 3, 2> B =
        math::Tensor<double, 3, 2>{{{0.4418331744229961, 0.4340139933332937},
            {0.6177669784693172, 0.5131382425543909}, {0.6503971819314672, 0.6010389534045444}}};

    Core::LinAlg::Tensor<double, 4> A_ddot_B = math::ddot(A.get_tensor(), B.get_tensor());

    EXPECT_NEAR(A_ddot_B(0), 1.58482764506344, 1e-10);
    EXPECT_NEAR(A_ddot_B(1), 1.5861759767039691, 1e-10);
    EXPECT_NEAR(A_ddot_B(2), 1.6852205412426562, 1e-10);
    EXPECT_NEAR(A_ddot_B(3), 1.7717581672187814, 1e-10);
  }

  TYPED_TEST(TensorOperationsTest, MatrixDdotThreeTensor)
  {
    TensorHolder<TestFixture::STORAGE_TYPE, double, 3, 2> A = math::Tensor<double, 3, 2>{{
        {0.8052231968327465, 0.5216471523936341},
        {0.9086488808086682, 0.3192360889885453},
        {0.09045934927090737, 0.30070005663620336},
    }};
    TensorHolder<TestFixture::STORAGE_TYPE, double, 3, 2, 4> B = math::Tensor<double, 3, 2, 4>{{
        {{0.11398436186354977, 0.8286813263076767, 0.04689631938924976, 0.6262871483113925},
            {0.5475861559192435, 0.8192869956700687, 0.1989475396788123, 0.8568503024577332}},
        {{0.3516526394320879, 0.7546476915298572, 0.2959617068796787, 0.8839364795611863},
            {0.3255116378322488, 0.16501589771914849, 0.3925292439465873, 0.0934603745586503}},
        {{0.8211056578369285, 0.15115201964256386, 0.3841144486921996, 0.9442607122388011},
            {0.9876254749018722, 0.4563045470947841, 0.8261228438427398, 0.25137413420705934}},
    }};

    Core::LinAlg::Tensor<double, 4> A_ddot_B = math::ddot(A.get_tensor(), B.get_tensor());

    EXPECT_NEAR(A_ddot_B(0), 1.172129170338301, 1e-10);
    EXPECT_NEAR(A_ddot_B(1), 1.9839248816243478, 1e-10);
    EXPECT_NEAR(A_ddot_B(2), 0.8189391251432785, 1e-10);
    EXPECT_NEAR(A_ddot_B(3), 1.9453037032761435, 1e-10);
  }

  TYPED_TEST(TensorOperationsTest, ThreeTensDdotThreeTensor)
  {
    TensorHolder<TestFixture::STORAGE_TYPE, double, 4, 3, 2> A = math::Tensor<double, 4, 3, 2>{{
        {{0.5973716482308843, 0.9028317603316274}, {0.5345579488018151, 0.5902013629854229},
            {0.03928176722538734, 0.3571817586345363}},
        {{0.07961309015596418, 0.30545991834281827}, {0.330719311982132, 0.7738302962105958},
            {0.039959208689977266, 0.42949217843163834}},
        {{0.3149268718426883, 0.6364911430675446}, {0.34634715008003303, 0.04309735620499444},
            {0.879915174517916, 0.763240587143681}},
        {{0.8780966427248583, 0.41750914383926696}, {0.6055775643937568, 0.5134666274082884},
            {0.5978366479629736, 0.2622156611319503}},
    }};
    TensorHolder<TestFixture::STORAGE_TYPE, double, 3, 2, 2> B = math::Tensor<double, 3, 2, 2>{{
        {{0.30087130894070724, 0.025399782050106068}, {0.30306256065103476, 0.24207587540352737}},
        {{0.5575781886626442, 0.5655070198881675}, {0.47513224741505056, 0.2927979762895091}},
        {{0.06425106069482445, 0.9788191457576426}, {0.33970784363786366, 0.4950486308824543}},
    }};

    auto A_ddot_B = math::ddot(A.get_tensor(), B.get_tensor());

    EXPECT_NEAR(A_ddot_B(0, 0), 1.1556893879240193, 1e-10);
    EXPECT_NEAR(A_ddot_B(0, 1), 0.92410502208987, 1e-10);
    EXPECT_NEAR(A_ddot_B(1, 0), 0.8170696456976152, 1e-10);
    EXPECT_NEAR(A_ddot_B(1, 1), 0.7412990229546594, 1e-10);
    EXPECT_NEAR(A_ddot_B(2, 0), 0.8170559534228106, 1e-10);
    EXPECT_NEAR(A_ddot_B(2, 1), 1.6096778150801772, 1e-10);
    EXPECT_NEAR(A_ddot_B(3, 0), 1.0998352261678221, 1e-10);
    EXPECT_NEAR(A_ddot_B(3, 1), 1.3311561690778437, 1e-10);
  }

  TYPED_TEST(TensorOperationsTest, FourTensorDdotMatrix)
  {
    TensorHolder<TestFixture::STORAGE_TYPE, double, 2, 4, 3, 2> A = math::Tensor<double, 2, 4, 3,
        2>{{
        {
            {{0.771320643266746, 0.0207519493594015}, {0.6336482349262754, 0.7488038825386119},
                {0.4985070123025904, 0.22479664553084766}},
            {{0.19806286475962398, 0.7605307121989587}, {0.16911083656253545, 0.08833981417401027},
                {0.6853598183677972, 0.9533933461949365}},
            {{0.003948266327914451, 0.5121922633857766}, {0.8126209616521135, 0.6125260668293881},
                {0.7217553174317995, 0.29187606817063316}},
            {{0.9177741225129434, 0.7145757833976906}, {0.5425443680112613, 0.14217004760152696},
                {0.3733407600514692, 0.6741336150663453}},
        },
        {
            {{0.4418331744229961, 0.4340139933332937}, {0.6177669784693172, 0.5131382425543909},
                {0.6503971819314672, 0.6010389534045444}},
            {{0.8052231968327465, 0.5216471523936341}, {0.9086488808086682, 0.3192360889885453},
                {0.09045934927090737, 0.30070005663620336}},
            {{0.11398436186354977, 0.8286813263076767}, {0.04689631938924976, 0.6262871483113925},
                {0.5475861559192435, 0.8192869956700687}},
            {{0.1989475396788123, 0.8568503024577332}, {0.3516526394320879, 0.7546476915298572},
                {0.2959617068796787, 0.8839364795611863}},
        },
    }};
    TensorHolder<TestFixture::STORAGE_TYPE, double, 3, 2> B = math::Tensor<double, 3, 2>{
        {
            {0.3255116378322488, 0.16501589771914849},
            {0.3925292439465873, 0.0934603745586503},
            {0.8211056578369285, 0.15115201964256386},
        },
    };

    Core::LinAlg::Tensor<double, 2, 4> A_ddot_B = math::ddot(A.get_tensor(), B.get_tensor());

    EXPECT_NEAR(A_ddot_B(0, 0), 1.0165125966071782, 1e-10);
    EXPECT_NEAR(A_ddot_B(0, 1), 0.971468800965393, 1e-10);
    EXPECT_NEAR(A_ddot_B(0, 2), 1.0987845120181288, 1e-10);
    EXPECT_NEAR(A_ddot_B(0, 3), 1.0513631864533974, 1e-10);

    EXPECT_NEAR(A_ddot_B(1, 0), 1.1307838039469182, 1e-10);
    EXPECT_NEAR(A_ddot_B(1, 1), 0.8544248817704513, 1e-10);
    EXPECT_NEAR(A_ddot_B(1, 2), 0.8242530123983, 1e-10);
    EXPECT_NEAR(A_ddot_B(1, 3), 0.7913418780962238, 1e-10);
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

  TYPED_TEST(TensorOperationsTest, MatrixDyadicVector)
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

  TYPED_TEST(TensorOperationsTest, VectorDyadicVector)
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