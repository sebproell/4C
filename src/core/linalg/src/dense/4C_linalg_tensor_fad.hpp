// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LINALG_TENSOR_FAD_HPP
#define FOUR_C_LINALG_TENSOR_FAD_HPP

#include "4C_config.hpp"

#include "4C_linalg_symmetric_tensor.hpp"
#include "4C_linalg_tensor.hpp"
#include "4C_linalg_tensor_internals.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_fad_meta.hpp"

#include <Sacado.hpp>

#include <cstring>
#include <type_traits>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  namespace Internal
  {
    template <typename Tensor>
    struct TensorCastType;


    template <typename ScalarType, TensorStorageType storage_type, typename Compression,
        std::size_t... n>
    struct TensorCastType<TensorInternal<ScalarType, storage_type, Compression, n...>>
    {
      template <typename NewScalarType>
      using type = TensorInternal<NewScalarType, TensorStorageType::owning, Compression, n...>;
    };
  }  // namespace Internal

  /*!
   * @brief Create a tensor of FAD-type values from a tensor of arithmetic values for
   * differentiating an expression w.r.t. this tensor.
   */
  auto make_auto_diff_tensor(const auto& tensor)
    requires(is_tensor<decltype(tensor)> &&
             std::is_arithmetic_v<typename std::remove_cvref_t<decltype(tensor)>::value_type>)
  {
    using OldValueType =
        std::remove_cvref_t<typename std::remove_cvref_t<decltype(tensor)>::value_type>;
    typename Internal::TensorCastType<typename std::remove_cvref_t<decltype(tensor)>>::
        template type<Sacado::Fad::DFad<OldValueType>>
            result;

    constexpr std::size_t compressed_size = std::remove_cvref_t<decltype(tensor)>::compressed_size;

    for (std::size_t i = 0; i < compressed_size; ++i)
    {
      result.container()[i] =
          Sacado::Fad::DFad<OldValueType>(compressed_size, i, tensor.container()[i]);
    }
    return result;
  }

  /*!
   * @brief Computes the resulting tensor derivative of the tensorial expression @p tensor w.r.t.
   * the tensor previously selected via @p Core::LinAlg::make_auto_diff_tensor(...).
   *
   * @note It is necessary to pass the original tensor type (the type of the tensor passed to @p
   * Core::LinAlg::make_auto_diff_tensor(...)) as a template parameter to this function so that it
   * can determine the correct shape of the derivative tensor.
   *
   * @tparam OriginalTensorType
   */
  template <typename OriginalTensorType>
  auto extract_derivative(const auto& tensor)
    requires(is_tensor<decltype(tensor)> && is_tensor<OriginalTensorType> &&
             FADUtils::SacadoFadType<typename std::remove_cvref_t<decltype(tensor)>::value_type> &&
             !is_compressed_tensor<OriginalTensorType> && !is_compressed_tensor<decltype(tensor)>)
  {
    using ResultingScalarType = decltype(std::declval<
        std::remove_cvref_t<typename std::remove_cvref_t<decltype(tensor)>::value_type>>()
            .dx(0));

    using ValueTensorShape = typename std::remove_cvref_t<decltype(tensor)>::shape_type;
    using OriginalTensorShape = typename OriginalTensorType::shape_type;
    using ResultingTensorType = typename Internal::DyadicProductTensorResult<ResultingScalarType,
        ValueTensorShape, OriginalTensorShape>::type;

    ResultingTensorType result;

    // Create a view on the arbitrary-rank result tensor and reinterpret it as a rank-2 tensor
    auto result_view = make_tensor_view<std::integer_sequence<std::size_t,
        std::remove_cvref_t<decltype(tensor)>::size(), OriginalTensorType::size()>>(result.data());

    // Ensure that the tensor FAD types is long enough that it can hold the derivative w.r.t. the
    // original tensor
    FOUR_C_ASSERT(std::ranges::all_of(tensor.container(),
                      [](const auto& fad_type)
                      {
                        return static_cast<std::size_t>(fad_type.length()) >=
                               OriginalTensorType::compressed_size;
                      }),
        "The size of the tensor FAD type does not match the expected size! Expecting at least {}.",
        0, OriginalTensorType::compressed_size);

    for (std::size_t i = 0; i < std::remove_cvref_t<decltype(tensor)>::compressed_size; ++i)
    {
      for (std::size_t j = 0; j < OriginalTensorType::compressed_size; ++j)
      {
        result_view(i, j) =
            tensor.container()[i].dx(j);  // This is the derivative of the i-th component w.r.t. j
      }
    }

    return result;
  }

  /*!
   * @brief Computes the resulting tensor derivative of the symmetric tensorial expression @p tensor
   * w.r.t. the symmetric tensor previously selected via @p
   * Core::LinAlg::make_auto_diff_tensor(...).
   *
   * @note It is necessary to pass the original tensor type (the type of the tensor passed to @p
   * Core::LinAlg::make_auto_diff_tensor(...)) as a template parameter to this function so that it
   * can determine the correct shape of the derivative tensor.
   *
   * @tparam OriginalTensorType
   */
  template <typename OriginalTensorType>
  auto extract_derivative(const auto& tensor)
    requires(is_tensor<decltype(tensor)> && is_tensor<OriginalTensorType> &&
             FADUtils::SacadoFadType<typename std::remove_cvref_t<decltype(tensor)>::value_type> &&
             is_symmetric_tensor<OriginalTensorType> && is_symmetric_tensor<decltype(tensor)> &&
             OriginalTensorType::rank() == 2 && std::remove_cvref_t<decltype(tensor)>::rank() == 2)
  {
    using ResultingScalarType = decltype(std::declval<
        std::remove_cvref_t<typename std::remove_cvref_t<decltype(tensor)>::value_type>>()
            .dx(0));

    using ValueTensorShape = typename std::remove_cvref_t<decltype(tensor)>::shape_type;
    using OriginalTensorShape = typename OriginalTensorType::shape_type;
    using ResultingTensorType = decltype(assume_symmetry(
        std::declval<typename Internal::DyadicProductTensorResult<ResultingScalarType,
            ValueTensorShape, OriginalTensorShape>::type>()));
    ResultingTensorType result;

    // Create a view on the arbitrary-rank result tensor and reinterpret it as a rank-2 tensor
    auto result_view = make_tensor_view<
        std::integer_sequence<std::size_t, std::remove_cvref_t<decltype(tensor)>::compressed_size,
            OriginalTensorType::compressed_size>>(result.data());

    // Ensure that the tensor FAD types is long enough that it can hold the derivative w.r.t. the
    // original tensor
    FOUR_C_ASSERT(std::ranges::all_of(tensor.container(),
                      [](const auto& fad_type)
                      {
                        return static_cast<std::size_t>(fad_type.length()) >=
                               OriginalTensorType::compressed_size;
                      }),
        "The size of the tensor FAD type does not match the expected size! Expecting at least {}.",
        OriginalTensorType::compressed_size);

    for (std::size_t j = 0; j < OriginalTensorType::compressed_size; ++j)
    {
      // Note: It is necessary to unscale the off-diagonal elements of symmtric tensors since the
      // derivative is computed only with the upper triangle, but we return the the full tensor as
      // derivative.
      const double unscale_factor = (j < OriginalTensorType::template extent<0>()) ? 1.0 : 0.5;
      for (std::size_t i = 0; i < std::remove_cvref_t<decltype(tensor)>::compressed_size; ++i)
      {
        result_view(i, j) =
            unscale_factor *
            tensor.container()[i].dx(j);  // This is the derivative of the i-th component w.r.t. j
      }
    }

    return result;
  }

  /*!
   * @brief Derivative of a scalar expression w.r.t. a tensor previously selected via @p
   * Core::LinAlg::make_auto_diff_tensor(...).
   *
   * @note It is necessary to pass the original tensor type (the type of the tensor passed to @p
   * Core::LinAlg::make_auto_diff_tensor(...)) as a template parameter to this function so it can
   * determine the correct shape of the derivative tensor.
   *
   * @tparam OriginalTensorType
   */
  template <typename OriginalTensorType>
  auto extract_derivative(const auto& scalar)
    requires(is_scalar<decltype(scalar)> && is_tensor<OriginalTensorType> &&
             FADUtils::SacadoFadType<std::remove_cvref_t<decltype(scalar)>> &&
             (!is_compressed_tensor<OriginalTensorType> || OriginalTensorType::rank() == 2))
  {
    using ResultingScalarType =
        decltype(std::declval<std::remove_cvref_t<decltype(scalar)>>().dx(0));

    using ResultingTensorType =
        Internal::TensorCastType<OriginalTensorType>::template type<ResultingScalarType>;

    ResultingTensorType result;

    // Create a view on the arbitrary-rank result tensor and reinterpret it as a rank-1 tensor
    auto result_view =
        make_tensor_view<std::integer_sequence<std::size_t, OriginalTensorType::compressed_size>>(
            result.data());

    // Ensure that the tensor FAD types is long enough that it can hold the derivative w.r.t. the
    // original tensor
    FOUR_C_ASSERT(static_cast<std::size_t>(scalar.length()) >= OriginalTensorType::compressed_size,
        "The size of the tensor FAD type does not match the expected size! Expecting at least {}, "
        "but only got {}.",
        OriginalTensorType::compressed_size, scalar.length());

    for (std::size_t j = 0; j < OriginalTensorType::compressed_size; ++j)
    {
      result_view(j) = scalar.dx(j);
    }

    if constexpr (is_symmetric_tensor<OriginalTensorType>)
    {
      // we need to unscale the off-diagonal elements since the result is a symmetric tensor and the
      // off-diagonal elements are multiplied by 2
      std::for_each(result_view.container().begin() + OriginalTensorType::template extent<0>(),
          result_view.container().end(), [](auto& value) { value *= 0.5; });
    }

    return result;
  }

  /*!
   * @brief Derivative of a tensorial expression @p tensor w.r.t. a scalar
   *
   * @tparam ScalarType
   */
  template <typename ScalarType>
  auto extract_derivative(const auto& tensor)
    requires(is_scalar<ScalarType> && is_tensor<decltype(tensor)> &&
             FADUtils::SacadoFadType<typename std::remove_cvref_t<decltype(tensor)>::value_type>)
  {
    using ResultingScalarType = decltype(std::declval<
        std::remove_cvref_t<typename std::remove_cvref_t<decltype(tensor)>::value_type>>()
            .dx(0));

    using ResultingTensorType = Internal::TensorCastType<
        std::remove_cvref_t<decltype(tensor)>>::template type<ResultingScalarType>;

    ResultingTensorType result;

    // Create a view on the arbitrary-rank result tensor and reinterpret it as a rank-1 tensor
    auto result_view = make_tensor_view<
        std::integer_sequence<std::size_t, std::remove_cvref_t<decltype(tensor)>::compressed_size>>(
        result.data());

    // Ensure that the tensor FAD types is long enough that it can hold the derivative w.r.t. the
    // original tensor
    FOUR_C_ASSERT(std::ranges::all_of(tensor.container(), [](const auto& fad_type)
                      { return static_cast<std::size_t>(fad_type.length()) >= 1; }),
        "The size of the tensor FAD type does not match the expected size! Expecting at least {}.",
        1);

    for (std::size_t j = 0; j < std::remove_cvref_t<decltype(tensor)>::compressed_size; ++j)
    {
      result_view(j) = tensor.container()[j].dx(0);
    }

    return result;
  }
}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE

#endif