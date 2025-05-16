// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LINALG_TENSOR_HPP
#define FOUR_C_LINALG_TENSOR_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_tensor_internals.hpp"
#include "4C_linalg_tensor_meta_utils.hpp"

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <cstring>
#include <span>
#include <tuple>
#include <type_traits>
#include <utility>


FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  template <typename Tensor>
  concept TensorConcept = requires(const Tensor& t) {
    { t.data() } -> std::convertible_to<const typename Tensor::value_type*>;
  };

  template <typename Tensor>
  concept Rank1TensorConcept = TensorConcept<Tensor> && Tensor::rank() == 1;

  template <typename Tensor>
  concept Rank2TensorConcept = TensorConcept<Tensor> && Tensor::rank() == 2;

  template <typename Tensor>
  concept SquareTensor =
      Rank2TensorConcept<Tensor> && Tensor::template extent<0>() == Tensor::template extent<1>();

  /*!
   * @brief An owning, dense tensor of arbitrary rank
   *
   * @copydetails TensorInternal
   */
  template <typename T, std::size_t... n>
  using Tensor = TensorInternal<T, TensorStorageType::owning, n...>;  // owns the memory

  /*!
   * @brief A dense view to a Tensor of arbitrary rank
   *
   * @copydetails TensorInternal
   */
  template <typename T, std::size_t... n>
  using TensorView = TensorInternal<T, TensorStorageType::view, n...>;  // view onto's other
                                                                        // memory

  // Implementation of tensor operations

  namespace Internal
  {
    template <typename T, typename S1>
    struct SameShapeTensorResult;

    template <typename T, std::size_t... s>
    struct SameShapeTensorResult<T, std::integer_sequence<std::size_t, s...>>
    {
      using type = TensorInternal<T, TensorStorageType::owning, s...>;
    };


    template <typename T, typename S1, typename S2>
    struct DyadicProductTensorResult;

    template <typename T, std::size_t... s1, std::size_t... s2>
    struct DyadicProductTensorResult<T, std::integer_sequence<std::size_t, s1...>,
        std::integer_sequence<std::size_t, s2...>>
    {
      using type = TensorInternal<T, TensorStorageType::owning, s1..., s2...>;
    };

    template <typename T, typename Shape, std::size_t... new_order>
    struct ReorderAxisTensorHelper;

    template <typename T, std::size_t... shape, std::size_t... new_order>
    struct ReorderAxisTensorHelper<T, std::integer_sequence<std::size_t, shape...>, new_order...>
    {
      static consteval std::size_t get_axis(std::size_t t)
      {
        constexpr std::array<std::size_t, sizeof...(shape)> tensor_shape = {shape...};
        return tensor_shape[t];
      }

      static consteval std::tuple<decltype(shape)...> get_new_index(
          std::tuple<decltype(shape)...> old_index)
      {
        auto make_array_from_tuple = [](auto&&... args) constexpr
        { return std::array<std::size_t, sizeof...(args)>{args...}; };

        std::array<std::size_t, sizeof...(shape)> old_index_array =
            std::apply(make_array_from_tuple, old_index);

        std::array<std::size_t, sizeof...(shape)> new_order_array = {new_order...};
        std::array<std::size_t, sizeof...(shape)> new_index_array;
        for (std::size_t i = 0; i < sizeof...(shape); ++i)
        {
          auto where_is_it = std::ranges::find(new_order_array, i);

          new_index_array[i] = old_index_array[std::distance(new_order_array.begin(), where_is_it)];
        }

        return std::tuple_cat(new_index_array);
      }

      static consteval std::array<std::size_t, (shape * ...)> get_index_mapping()
      {
        std::array<std::size_t, (shape * ...)> index_mapping = {0};

        for (std::size_t i = 0; i < (shape * ...); ++i)
        {
          auto old_index = get_md_index<OrderType::column_major, TensorBoundCheck::no_check,
              get_axis(new_order)...>(i);

          auto new_index = get_new_index(old_index);

          index_mapping[i] = std::apply(
              get_flat_index<OrderType::column_major, TensorBoundCheck::no_check, shape...>,
              new_index);
        }

        return index_mapping;
      }

      using result_type = TensorInternal<T, TensorStorageType::owning, get_axis(new_order)...>;
    };
  }  // namespace Internal

  /*!
   * @brief Computes and returns the determinant of a square tensor of rank 2
   *
   * @param A
   * @return auto
   */
  constexpr auto det(const SquareTensor auto& A);

  /*!
   * @brief Computes and returns the trace of a square tensor of rank 2
   *
   * @param A
   * @return auto
   */
  constexpr auto trace(const SquareTensor auto& A);

  /*!
   * @brief Computes and returns the inverse of a square tensor of rank 2
   *
   * @param A
   * @return auto
   */
  auto inv(const SquareTensor auto& A);

  /*!
   * @brief Computes and returns the transpose of a rank 2 tensor
   *
   * @param A
   * @return auto
   */
  constexpr auto transpose(const Rank2TensorConcept auto& A);

  /*!
   * @brief Computes the dot product of two rank-1 tensors
   *
   * This function performs the inner product between two vectors a * b. Both sizes need to match.
   *
   * @tparam TensorLeft The type of the rank-1 tensor (vector).
   * @tparam TensorRight The type of the rank-1 tensor (vector).
   * @param a The rank-1 tensor (vector).
   * @param b The rank-1 tensor (vector).
   * @return scalar resulting from the dot product.
   */
  template <typename TensorLeft, typename TensorRight>
    requires(Rank1TensorConcept<TensorLeft> && Rank1TensorConcept<TensorRight> &&
             TensorLeft::template extent<0>() == TensorRight::template extent<0>())
  constexpr auto dot(const TensorLeft& a, const TensorRight& b);

  /*!
   * @brief Computes the dot product of a rank-2 tensor and a rank-1 tensor.
   *
   * This function performs the matrix-vector multiplication A * b, where A is a rank-2 tensor
   * (matrix) and b is a rank-1 tensor (vector). The number of columns in A must match the size of
   * b.
   *
   * @tparam TensorLeft The type of the rank-2 tensor (matrix).
   * @tparam TensorRight The type of the rank-1 tensor (vector).
   * @param A The rank-2 tensor (matrix).
   * @param b The rank-1 tensor (vector).
   * @return A rank-1 tensor (vector) resulting from the dot product.
   */
  template <typename TensorLeft, typename TensorRight>
    requires(Rank2TensorConcept<TensorLeft> && Rank1TensorConcept<TensorRight> &&
             TensorLeft::template extent<1>() == TensorRight::template extent<0>())
  constexpr auto dot(const TensorLeft& A, const TensorRight& b);

  /*!
   * @brief Computes the dot product of a rank-1 tensor and a rank-2 tensor.
   *
   * This function performs the vector-matrix multiplication a * B, where a is a rank-1 tensor
   * (vector) and B is a rank-2 tensor (matrix). The size of a must match the number of rows in B.
   *
   * @tparam TensorLeft The type of the rank-1 tensor (vector).
   * @tparam TensorRight The type of the rank-2 tensor (matrix).
   * @param a The rank-1 tensor (vector).
   * @param B The rank-2 tensor (matrix).
   * @return A rank-1 tensor (vector) resulting from the dot product.
   */
  template <typename TensorLeft, typename TensorRight>
    requires(Rank1TensorConcept<TensorLeft> && Rank2TensorConcept<TensorRight> &&
             TensorLeft::template extent<0>() == TensorRight::template extent<0>())
  constexpr auto dot(const TensorLeft& a, const TensorRight& B);

  /*!
   * @brief Computes the dot product of two rank-2 tensors.
   *
   * This function performs the matrix-matrix multiplication A * B, where A and B are rank-2 tensors
   * (matrices). The number of columns in A must match the number of rows in B.
   *
   * @tparam TensorLeft The type of the first rank-2 tensor (matrix).
   * @tparam TensorRight The type of the second rank-2 tensor (matrix).
   * @param A The first rank-2 tensor (matrix).
   * @param B The second rank-2 tensor (matrix).
   * @return A rank-2 tensor (matrix) resulting from the dot product.
   */
  template <typename TensorLeft, typename TensorRight>
    requires(Rank2TensorConcept<TensorLeft> && Rank2TensorConcept<TensorRight> &&
             TensorLeft::template extent<1>() == TensorRight::template extent<0>())
  constexpr auto dot(const TensorLeft& A, const TensorRight& B);

  /*!
   * @brief Computes the double dot product of two rank-2 tensors.
   *
   * This function computes the double dot product of two rank-2 tensors (matrices) A and B, which
   * is defined as the sum of the element-wise products of the corresponding elements in A and B.
   *
   * @tparam TensorLeft The type of the first rank-2 tensor (matrix).
   * @tparam TensorRight The type of the second rank-2 tensor (matrix).
   * @param A The first rank-2 tensor (matrix).
   * @param B The second rank-2 tensor (matrix).
   * @return A scalar value representing the double dot product.
   */
  template <typename TensorLeft, typename TensorRight>
    requires(Rank2TensorConcept<TensorLeft> && Rank2TensorConcept<TensorRight> &&
             TensorLeft::template extent<0>() == TensorRight::template extent<0>() &&
             TensorLeft::template extent<1>() == TensorRight::template extent<1>())
  constexpr auto ddot(const TensorLeft& A, const TensorRight& B);

  /*!
   * @brief Adds two tensors element-wise and returns the resulting tensor.
   *
   * This function performs an element-wise addition of two tensors with the same shape.
   * The resulting tensor will have the same shape as the input tensors.
   *
   * @tparam TensorLeft The type of the first tensor.
   * @tparam TensorRight The type of the second tensor.
   * @param A The first tensor.
   * @param B The second tensor.
   * @return A tensor containing the element-wise sum of the input tensors.
   *
   * @note The input tensors must have the same shape, as enforced by the `requires` clause.
   */
  template <typename TensorLeft, typename TensorRight>
    requires(std::is_same_v<typename TensorLeft::shape_type, typename TensorRight::shape_type>)
  constexpr auto add(const TensorLeft& A, const TensorRight& B);

  /*!
   * @brief Subtracts two tensors element-wise and returns the resulting tensor.
   *
   * This function performs an element-wise subtraction of two tensors with the same shape.
   * The resulting tensor will have the same shape as the input tensors.
   *
   * @tparam TensorLeft The type of the first tensor.
   * @tparam TensorRight The type of the second tensor.
   * @param A The first tensor.
   * @param B The second tensor.
   * @return A tensor containing the element-wise subtraction of the input tensors.
   *
   * @note The input tensors must have the same shape, as enforced by the `requires` clause.
   */
  template <typename TensorLeft, typename TensorRight>
    requires(std::is_same_v<typename TensorLeft::shape_type, typename TensorRight::shape_type>)
  constexpr auto subtract(const TensorLeft& A, const TensorRight& B);

  /*!
   * @brief Return a tensor scaled by a scalar
   *
   * @tparam Tensor
   * @tparam Scalar
   */
  template <typename Tensor, typename Scalar>
    requires(std::is_arithmetic_v<Scalar>)
  constexpr auto scale(const Tensor& tensor, const Scalar& b);

  /*!
   * @brief Returns the dyadic product of two tensors with arbitrary rank
   *
   * @tparam TensorLeft
   * @tparam TensorRight
   */
  template <typename TensorLeft, typename TensorRight>
    requires(TensorConcept<TensorLeft> && TensorConcept<TensorRight>)
  constexpr auto dyadic(const TensorLeft& A, const TensorRight& B);

  /*!
   * @brief Reorders the axes of a tensor (generalization of transpose)
   *
   * Given the tensor A of shape (2, 3, 4), @p reorder_axis<0, 2, 1>(A) will return a tensor of the
   * shape (2, 4, 3), i.e., the the order given in the template brackets is the new order of the
   * axes.
   *
   * @tparam new_order
   * @return A tensor with the axes reordered according to @p new_order
   */
  template <std::size_t... new_order>
  constexpr auto reorder_axis(auto tensor)
    requires(TensorConcept<decltype(tensor)> && sizeof...(new_order) == decltype(tensor)::rank());


  /*!
   * @brief Add another tensor onto this tensor
   *
   * @tparam OtherTensor
   */
  template <typename OtherTensor, typename Number, TensorStorageType storage_type, std::size_t... n>
    requires(std::is_same_v<typename OtherTensor::shape_type, typename OtherTensor::shape_type>)
  constexpr TensorInternal<Number, storage_type, n...>& operator+=(
      TensorInternal<Number, storage_type, n...>& tensor, const OtherTensor& B);

  /*!
   * @brief Subtract another tensor from this tensor
   *
   * @tparam OtherTensor
   */
  template <typename OtherTensor, typename Number, TensorStorageType storage_type, std::size_t... n>
    requires(std::is_same_v<typename OtherTensor::shape_type, typename OtherTensor::shape_type>)
  constexpr TensorInternal<Number, storage_type, n...>& operator-=(
      TensorInternal<Number, storage_type, n...>& tensor, const OtherTensor& B);

  /*!
   * @brief Scale the tensor with a scalar value
   *
   * @tparam OtherTensor
   */
  template <typename Scalar, typename Number, TensorStorageType storage_type, std::size_t... n>
    requires(std::is_arithmetic_v<Scalar>)
  constexpr TensorInternal<Number, storage_type, n...>& operator*=(
      TensorInternal<Number, storage_type, n...>& tensor, const Scalar b);

  /*!
   * @brief Scales the tensor with the inverse of the scalar value b
   *
   * @tparam Scalar
   * @tparam Number
   * @tparam storage_type
   * @tparam n
   */
  template <typename Scalar, typename Number, TensorStorageType storage_type, std::size_t... n>
    requires(std::is_arithmetic_v<Scalar>)
  constexpr TensorInternal<Number, storage_type, n...>& operator/=(
      TensorInternal<Number, storage_type, n...>& tensor, const Scalar b);

  template <typename TensorLeft, typename TensorRight>
    requires(std::is_same_v<typename TensorLeft::shape_type, typename TensorRight::shape_type>)
  constexpr auto operator+(const TensorLeft& A, const TensorRight& B)
  {
    return add(A, B);
  }
  template <typename TensorLeft, typename TensorRight>
    requires(std::is_same_v<typename TensorLeft::shape_type, typename TensorRight::shape_type>)
  constexpr auto operator-(const TensorLeft& A, const TensorRight& B)
  {
    return subtract(A, B);
  }

  template <typename Scalar, typename Tensor>
    requires(std::is_arithmetic_v<Scalar>)
  constexpr auto operator+(const Scalar b, const Tensor& tensor)
  {
    return add_scalar(tensor, b);
  }

  template <typename Scalar, typename Tensor>
    requires(std::is_arithmetic_v<Scalar>)
  constexpr auto operator*(const Tensor& tensor, const Scalar b)
  {
    return scale(tensor, b);
  }

  template <typename Scalar, typename Tensor>
    requires(std::is_arithmetic_v<Scalar>)
  constexpr auto operator*(const Scalar b, const Tensor& tensor)
  {
    return scale(tensor, b);
  }

  template <typename Scalar, typename Tensor>
    requires(std::is_arithmetic_v<Scalar>)
  constexpr auto operator/(const Tensor& tensor, const Scalar b)
  {
    return scale(tensor, Scalar(1) / b);
  }

  template <typename TensorLeft, typename TensorRight>
    requires(TensorConcept<TensorLeft> && TensorConcept<TensorRight>)
  constexpr auto operator*(const TensorLeft& A, const TensorRight& B)
  {
    return dot(A, B);
  }

  template <typename TensorLeft, typename TensorRight>
    requires(std::is_same_v<typename TensorLeft::shape_type, typename TensorRight::shape_type>)
  constexpr bool operator==(const TensorLeft& A, const TensorRight& B)
  {
    return std::ranges::equal(A.container(), B.container());
  }

  template <typename TensorLeft, typename TensorRight>
    requires(std::is_same_v<typename TensorLeft::shape_type, typename TensorRight::shape_type>)
  constexpr bool operator!=(const TensorLeft& A, const TensorRight& B)
  {
    return !std::ranges::equal(A.container(), B.container());
  }


  // Actual implementations
  constexpr auto det(const SquareTensor auto& A)
  {
    using Tensor = std::remove_cvref_t<decltype(A)>;
    using value_type = Tensor::value_type;
    constexpr std::size_t n = Tensor::template extent<0>();

    return DenseFunctions::determinant<value_type, n, n>(A.data());
  }

  constexpr auto trace(const SquareTensor auto& A)
  {
    using Tensor = std::remove_cvref_t<decltype(A)>;
    using value_type = Tensor::value_type;
    constexpr std::size_t n = Tensor::template extent<0>();

    return DenseFunctions::trace<value_type, n, n>(A.data());
  }

  auto inv(const SquareTensor auto& A)
  {
    using Tensor = std::remove_cvref_t<decltype(A)>;
    using value_type = Tensor::value_type;
    constexpr std::size_t n = Tensor::template extent<0>();

    Core::LinAlg::Tensor<value_type, n, n> dest;
    DenseFunctions::invert<value_type, n, n>(dest.data(), A.data());
    return dest;
  }

  constexpr auto transpose(const Rank2TensorConcept auto& A) { return reorder_axis<1, 0>(A); }

  template <typename TensorLeft, typename TensorRight>
    requires(Rank1TensorConcept<TensorLeft> && Rank1TensorConcept<TensorRight> &&
             TensorLeft::template extent<0>() == TensorRight::template extent<0>())
  constexpr auto dot(const TensorLeft& a, const TensorRight& b)
  {
    using value_type = decltype(std::declval<typename TensorLeft::value_type>() *
                                std::declval<typename TensorRight::value_type>());
    constexpr std::size_t m = TensorLeft::template extent<0>();

    return DenseFunctions::dot<value_type, m, 1>(a.data(), b.data());
  }

  template <typename TensorLeft, typename TensorRight>
    requires(Rank2TensorConcept<TensorLeft> && Rank1TensorConcept<TensorRight> &&
             TensorLeft::template extent<1>() == TensorRight::template extent<0>())
  constexpr auto dot(const TensorLeft& A, const TensorRight& b)
  {
    using value_type = decltype(std::declval<typename TensorLeft::value_type>() *
                                std::declval<typename TensorRight::value_type>());
    constexpr std::size_t m = TensorLeft::template extent<0>();
    constexpr std::size_t n = TensorLeft::template extent<1>();

    Tensor<value_type, m> dest;
    DenseFunctions::multiply<value_type, m, n, 1>(dest.data(), A.data(), b.data());
    return dest;
  }

  template <typename TensorLeft, typename TensorRight>
    requires(Rank1TensorConcept<TensorLeft> && Rank2TensorConcept<TensorRight> &&
             TensorLeft::template extent<0>() == TensorRight::template extent<0>())
  constexpr auto dot(const TensorLeft& a, const TensorRight& B)
  {
    using value_type = decltype(std::declval<typename TensorLeft::value_type>() *
                                std::declval<typename TensorRight::value_type>());
    constexpr std::size_t m = TensorRight::template extent<0>();
    constexpr std::size_t n = TensorRight::template extent<1>();

    Tensor<value_type, n> dest;
    DenseFunctions::multiply<value_type, 1, m, n>(dest.data(), a.data(), B.data());
    return dest;
  }

  template <typename TensorLeft, typename TensorRight>
    requires(Rank2TensorConcept<TensorLeft> && Rank2TensorConcept<TensorRight> &&
             TensorLeft::template extent<1>() == TensorRight::template extent<0>())
  constexpr auto dot(const TensorLeft& A, const TensorRight& B)
  {
    using value_type = decltype(std::declval<typename TensorLeft::value_type>() *
                                std::declval<typename TensorRight::value_type>());
    constexpr std::size_t m = TensorLeft::template extent<0>();
    constexpr std::size_t k = TensorLeft::template extent<1>();
    constexpr std::size_t n = TensorRight::template extent<1>();

    Core::LinAlg::Tensor<value_type, n, m> dest;
    DenseFunctions::multiply<value_type, m, k, n>(dest.data(), A.data(), B.data());
    return dest;
  }

  template <typename TensorLeft, typename TensorRight>
    requires(Rank2TensorConcept<TensorLeft> && Rank2TensorConcept<TensorRight> &&
             TensorLeft::template extent<0>() == TensorRight::template extent<0>() &&
             TensorLeft::template extent<1>() == TensorRight::template extent<1>())
  constexpr auto ddot(const TensorLeft& A, const TensorRight& B)
  {
    using value_type = decltype(std::declval<typename TensorLeft::value_type>() *
                                std::declval<typename TensorRight::value_type>());
    constexpr std::size_t m = TensorLeft::template extent<0>();
    constexpr std::size_t n = TensorLeft::template extent<1>();

    return DenseFunctions::dot<value_type, m, n>(A.data(), B.data());
  }

  template <typename TensorLeft, typename TensorRight>
    requires(std::is_same_v<typename TensorLeft::shape_type, typename TensorRight::shape_type>)
  constexpr auto add(const TensorLeft& A, const TensorRight& B)
  {
    using result_value_type = decltype(std::declval<typename TensorLeft::value_type>() +
                                       std::declval<typename TensorRight::value_type>());

    typename Internal::SameShapeTensorResult<result_value_type,
        typename TensorLeft::shape_type>::type tens_out{};
    DenseFunctions::update<result_value_type, TensorLeft::size(), 1>(
        tens_out.data(), A.data(), B.data());

    return tens_out;
  }

  template <typename TensorLeft, typename TensorRight>
    requires(std::is_same_v<typename TensorLeft::shape_type, typename TensorRight::shape_type>)
  constexpr auto subtract(const TensorLeft& A, const TensorRight& B)
  {
    using result_value_type = decltype(std::declval<typename TensorLeft::value_type>() +
                                       std::declval<typename TensorRight::value_type>());

    typename Internal::SameShapeTensorResult<result_value_type,
        typename TensorLeft::shape_type>::type tens_out{};
    DenseFunctions::update<result_value_type, TensorLeft::size(), 1>(
        tens_out.data(), 1.0, A.data(), -1.0, B.data());

    return tens_out;
  }

  template <typename Tensor, typename Scalar>
    requires(std::is_arithmetic_v<Scalar>)
  constexpr auto scale(const Tensor& tensor, const Scalar& b)
  {
    using result_value_type =
        decltype(std::declval<Scalar>() + std::declval<typename Tensor::value_type>());

    typename Internal::SameShapeTensorResult<result_value_type, typename Tensor::shape_type>::type
        tens_out = tensor;

    tens_out *= b;

    return tens_out;
  }

  template <typename TensorLeft, typename TensorRight>
    requires(TensorConcept<TensorLeft> && TensorConcept<TensorRight>)
  constexpr auto dyadic(const TensorLeft& A, const TensorRight& B)
  {
    using value_type = decltype(std::declval<typename TensorLeft::value_type>() *
                                std::declval<typename TensorRight::value_type>());

    typename Internal::DyadicProductTensorResult<value_type, typename TensorLeft::shape_type,
        typename TensorRight::shape_type>::type dest{};

    DenseFunctions::multiply<value_type, TensorLeft::size(), 1, TensorRight::size()>(
        dest.data(), A.data(), B.data());

    return dest;
  }

  template <std::size_t... new_order>
  constexpr auto reorder_axis(auto tensor)
    requires(TensorConcept<decltype(tensor)> && sizeof...(new_order) == decltype(tensor)::rank())
  {
    constexpr std::array new_order_array = {new_order...};
    static_assert(std::ranges::max(new_order_array) < tensor.rank(),
        "Invalid tensor axis reordering. An axis index is out of bounds.");
    static_assert(
        [new_order_array]() constexpr
        {
          std::array sorted_indices = new_order_array;
          std::ranges::sort(sorted_indices);
          return std::ranges::adjacent_find(sorted_indices) == std::end(sorted_indices);
        }(),
        "All indices must be unique during reordering!");

    using Tensor = decltype(tensor);
    typename Internal::ReorderAxisTensorHelper<typename Tensor::value_type,
        typename Tensor::shape_type, new_order...>::result_type result;

    constexpr std::array<std::size_t, result.size()> index_reorder =
        Internal::ReorderAxisTensorHelper<typename Tensor::value_type, typename Tensor::shape_type,
            new_order...>::get_index_mapping();

    std::transform(index_reorder.begin(), index_reorder.end(), result.data(),
        [&tensor](const auto& i) { return tensor.container()[i]; });

    return result;
  }

  template <typename OtherTensor, typename Number, TensorStorageType storage_type, std::size_t... n>
    requires(std::is_same_v<typename OtherTensor::shape_type, typename OtherTensor::shape_type>)
  constexpr TensorInternal<Number, storage_type, n...>& operator+=(
      TensorInternal<Number, storage_type, n...>& tensor, const OtherTensor& B)
  {
    DenseFunctions::update<typename TensorInternal<Number, storage_type, n...>::value_type,
        OtherTensor::size(), 1>(1.0, tensor.data(), 1.0, B.data());
    return tensor;
  }

  template <typename OtherTensor, typename Number, TensorStorageType storage_type, std::size_t... n>
    requires(std::is_same_v<typename OtherTensor::shape_type, typename OtherTensor::shape_type>)
  constexpr TensorInternal<Number, storage_type, n...>& operator-=(
      TensorInternal<Number, storage_type, n...>& tensor, const OtherTensor& B)
  {
    DenseFunctions::update<typename TensorInternal<Number, storage_type, n...>::value_type,
        OtherTensor::size(), 1>(1.0, tensor.data(), -1.0, B.data());
    return tensor;
  }

  template <typename Scalar, typename Number, TensorStorageType storage_type, std::size_t... n>
    requires(std::is_arithmetic_v<Scalar>)
  constexpr TensorInternal<Number, storage_type, n...>& operator*=(
      TensorInternal<Number, storage_type, n...>& tensor, const Scalar b)
  {
    for (auto& value : tensor.container()) value *= b;
    return tensor;
  }

  template <typename Scalar, typename Number, TensorStorageType storage_type, std::size_t... n>
    requires(std::is_arithmetic_v<Scalar>)
  constexpr TensorInternal<Number, storage_type, n...>& operator/=(
      TensorInternal<Number, storage_type, n...>& tensor, const Scalar b)
  {
    tensor *= Scalar(1.0) / b;
    return tensor;
  }
}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE

#endif