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


  namespace Internal
  {

    enum class StorageType
    {
      owning,
      view
    };

    /*!
     * @brief General tensor class for dense tensors of arbitrary rank
     *
     * This class is designed for small tensors typically used in physics simulations, such as
     * strains, stresses, and linearizations at the Gauss point level. The tensor dimensions are
     * usually small (e.g., 2x2 or 3x3), with a rank of 1 or 2, and rarely rank 4. This
     * results in typically 9 or fewer entries (rarely 81 for 4-tensors in 3D). The data is
     * stored on the stack avoiding dynamic memory allocation and ensuring efficient memory access
     * via typically high locality.
     *
     * Key characteristics:
     * - **Stack Allocation**: The tensor's data (if it is owning) is stored on the stack, avoiding
     * the overhead of dynamic memory allocation.
     * - **Temporary Objects**: Operations on tensors may create temporary objects, but these do not
     *   involve dynamic memory allocation, ensuring efficient computation.
     * - **Use Case**: This class is ideal for small-scale tensor operations at the Gauss point
     * level in finite element simulations. It allows to write physical equations close to their
     * pendant on paper.
     *
     * Comparison with @p LinAlg::Matrix:
     * - **Tensor**: Designed for small, fixed-size data stored on the stack. Suitable for
     * operations at the Gauss point level (e.g., strains, stresses).
     * - **Matrix**: Designed for larger, dynamically allocated data. Suitable for element-level
     *   calculations (e.g., stiffness and mass matrices). For example, a stiffness matrix for a
     *   hex27 element may have dimensions 81x81 with 6561 entries, requiring dynamic memory
     * allocation.
     *
     * Shared Backend:
     * Both @p LinAlg::Tensor and @p LinAlg::Matrix share the same backend for matrix and vector
     * operations, ensuring consistency and avoiding duplication of functionality.
     *
     * @tparam Number Type of the tensor elements (e.g., `double`).
     * @tparam storage_type Specifies whether the tensor owns its data (`StorageType::owning`) or is
     * a view (`StorageType::view`).
     * @tparam n Dimensions of the tensor, specified as a variadic template parameter.
     */
    template <typename Number, StorageType storage_type, std::size_t... n>
    class Tensor
    {
     public:
      using value_type = Number;
      using shape_type = std::integer_sequence<std::size_t, n...>;

      static constexpr std::size_t rank_ = sizeof...(n);
      static constexpr std::size_t size_ = (n * ...);

      static_assert(storage_type == StorageType::view or
                        std::is_same_v<std::remove_cv_t<value_type>, value_type>,
          "Owning Tensor must have a non-const, non-volatile value_type");
      static_assert(storage_type == StorageType::owning or
                        std::is_same_v<std::remove_volatile_t<value_type>, value_type>,
          "TensorView must have a non-volatile value_type");
      static_assert(sizeof...(n) > 0, "Tensor must have at least rank one");
      static_assert(
          []() constexpr
          {
            std::array shape = {n...};
            return std::ranges::all_of(shape, [](auto dim) { return dim > 0; });
          }(),
          "The extents of all axes must larger than zero");

      // Container is either a std::array for an owning tensor or a std::span for a view
      using ContainerType = std::conditional_t<storage_type == StorageType::owning,
          std::array<Number, size_>, std::span<Number, size_>>;

     private:
      ContainerType data_;

     public:
      /*!
       * @brief Default constructor
       *
       * @note Only an owning tensor has the default constructor
       */
      Tensor()
        requires(std::is_default_constructible_v<ContainerType>)
      = default;

      Tensor(const Tensor&) = default;

      /*!
       * @brief Move constructor
       *
       * @note Copy operation of a view is cheap, so we don't need to have a move constructor
       */
      Tensor(Tensor&&) noexcept
        requires(storage_type == StorageType::owning)
      = default;
      Tensor& operator=(const Tensor&) noexcept = default;

      /*!
       * @brief Move assignment operator
       *
       * @note Copy operation of a view is cheap, so we don't need to have a move assignment
       * operator
       */
      Tensor& operator=(Tensor&&) noexcept
        requires(storage_type == StorageType::owning)
      = default;
      ~Tensor() = default;

      /*!
       * @brief Construct an owning tensor from given data
       *
       * You can assign the tensor via
       *
       * @code
       * Tensor<double, 3, 4> A = {{
       *     {1.0, 2.0, 3.0, 4.0},
       *     {5.0, 6.0, 7.0, 8.0},
       *     {9.0, 10.0, 11.0, 12.0}
       * }};
       * @endcode
       *
       */
      Tensor(const TensorInitializerList<Number, n...>::type& lst)
        requires(std::is_default_constructible_v<ContainerType>);

      /*!
       * @brief Create a view on an owning tensor
       */
      Tensor(Tensor<Number, StorageType::owning, n...>& view_on)
        requires(storage_type == StorageType::view)
          : data_(view_on.container())
      {
      }

      /*!
       * @brief Create a constant view on a const owning tensor
       */
      Tensor(const Tensor<std::remove_cv_t<Number>, StorageType::owning, n...>& view_on)
        requires(storage_type == StorageType::view && std::is_const_v<Number>)
          : data_(view_on.container())
      {
      }

      /*!
       * @brief Create a copy of a view
       */
      Tensor(const Tensor<std::remove_cv_t<Number>, StorageType::view, n...>& other)
        requires(storage_type == StorageType::owning)
      {
        std::ranges::copy(other.container(), data_.begin());
      }

      /*!
       * @brief Access to the underlying raw data of the tensor
       *
       * @return Number*
       */
      [[nodiscard]] Number* data() { return data_.data(); }

      /*!
       * @brief Access to the underlying raw data of the tensor for readonly access
       *
       * @return Number*
       */
      [[nodiscard]] const Number* data() const { return data_.data(); }

      /*!
       * @brief Access to the underlying container (std::array or std::span, depending on whether it
       * is owning or a view)
       *
       * @return Number*
       */
      [[nodiscard]] ContainerType& container() { return data_; }

      /*!
       * @brief Access to the underlying readonly container (std::array or std::span, depending on
       * whether it is owning or a view)
       *
       * @return Number*
       */
      [[nodiscard]] const ContainerType& container() const { return data_; }

      /*!
       * @brief Indexing operator to access individual values of the tensor (without bound checks)
       *
       * @param i
       * @return Number&
       */
      [[nodiscard]] Number& operator()(decltype(n)... i);

      /*!
       * @brief Indexing operator to access individual values of the tensor in readonly mode
       * (without bound checks)
       *
       * @param i
       * @return Number&
       */
      [[nodiscard]] const Number& operator()(decltype(n)... i) const;

      /*!
       * @brief Indexing operator to access individual values of the tensor (with bound checks)
       *
       * @param i
       * @return Number&
       */
      [[nodiscard]] Number& at(decltype(n)... i);

      /*!
       * @brief Indexing operator to access individual values of the tensor in readonly mode
       * (with bound checks)
       *
       * @param i
       * @return Number&
       */
      [[nodiscard]] const Number& at(decltype(n)... i) const;

      [[nodiscard]] static constexpr std::size_t rank() { return rank_; }
      [[nodiscard]] static constexpr std::size_t size() { return size_; }

      /*!
       * @brief Returns the dimension of the i-th axis of the tensor
       *
       * @tparam i
       */
      template <std::size_t i>
        requires(i < rank())
      [[nodiscard]] static constexpr auto extent()
      {
        return std::get<i>(std::make_tuple(n...));
      }

      /*!
       * @brief Returns a std::tuple of the shape of the tensor
       *
       * @return constexpr auto
       */
      [[nodiscard]] static constexpr auto shape() { return std::make_tuple(n...); }

      /*!
       * @brief Fills the tensor with the given value
       *
       * @param value
       */
      constexpr void fill(const Number& value) { std::ranges::fill(data_, value); }
    };

    // actual implementations

    template <typename Number, Internal::StorageType storage_type, std::size_t... n>
    Tensor<Number, storage_type, n...>::Tensor(const TensorInitializerList<Number, n...>::type& lst)
      requires(std::is_default_constructible_v<ContainerType>)
    {
      constexpr std::array<std::size_t, size()> index_mapping = order_type_mapping<n...>();
      for (std::size_t i = 0; i < size(); ++i)
      {
        data_[index_mapping[i]] = *(get_view_to_first_element<Number, n...>(lst) + i);
      }
    }

    template <typename Number, Internal::StorageType storage_type, std::size_t... n>
    Number& Tensor<Number, storage_type, n...>::operator()(decltype(n)... i)
    {
      return data_[get_flat_index<OrderType::column_major, TensorBoundCheck::no_check, n...>(i...)];
    }

    template <typename Number, Internal::StorageType storage_type, std::size_t... n>
    const Number& Tensor<Number, storage_type, n...>::operator()(decltype(n)... i) const
    {
      return data_[get_flat_index<OrderType::column_major, TensorBoundCheck::no_check, n...>(i...)];
    }

    template <typename Number, Internal::StorageType storage_type, std::size_t... n>
    Number& Tensor<Number, storage_type, n...>::at(decltype(n)... i)
    {
      return data_[get_flat_index<OrderType::column_major, TensorBoundCheck::check, n...>(i...)];
    }

    template <typename Number, Internal::StorageType storage_type, std::size_t... n>
    const Number& Tensor<Number, storage_type, n...>::at(decltype(n)... i) const
    {
      return data_[get_flat_index<OrderType::column_major, TensorBoundCheck::check, n...>(i...)];
    }
  }  // namespace Internal

  /*!
   * @brief An owning, dense tensor of arbitrary rank
   *
   * @copydetails Internal::Tensor
   */
  template <typename T, std::size_t... n>
  using Tensor = Internal::Tensor<T, Internal::StorageType::owning, n...>;  // owns the memory

  /*!
   * @brief A dense view to a Tensor of arbitrary rank
   *
   * @copydetails Internal::Tensor
   */
  template <typename T, std::size_t... n>
  using TensorView = Internal::Tensor<T, Internal::StorageType::view, n...>;  // view onto's other
                                                                              // memory

  // Implementation of tensor operations

  namespace Internal
  {
    template <typename T, typename S1>
    struct SameShapeTensorResult;

    template <typename T, std::size_t... s>
    struct SameShapeTensorResult<T, std::integer_sequence<std::size_t, s...>>
    {
      using type = Tensor<T, StorageType::owning, s...>;
    };


    template <typename T, typename S1, typename S2>
    struct DyadicProductTensorResult;

    template <typename T, std::size_t... s1, std::size_t... s2>
    struct DyadicProductTensorResult<T, std::integer_sequence<std::size_t, s1...>,
        std::integer_sequence<std::size_t, s2...>>
    {
      using type = Tensor<T, StorageType::owning, s1..., s2...>;
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

      using result_type = Tensor<T, StorageType::owning, get_axis(new_order)...>;
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

  namespace Internal
  {
    // The underlying class of Tensor is in the Internal namespace, so we need to define these
    // operators also there.

    /*!
     * @brief Add another tensor onto this tensor
     *
     * @tparam OtherTensor
     */
    template <typename OtherTensor, typename Number, StorageType storage_type, std::size_t... n>
      requires(std::is_same_v<typename OtherTensor::shape_type, typename OtherTensor::shape_type>)
    constexpr Tensor<Number, storage_type, n...>& operator+=(
        Tensor<Number, storage_type, n...>& tensor, const OtherTensor& B);

    /*!
     * @brief Subtract another tensor from this tensor
     *
     * @tparam OtherTensor
     */
    template <typename OtherTensor, typename Number, StorageType storage_type, std::size_t... n>
      requires(std::is_same_v<typename OtherTensor::shape_type, typename OtherTensor::shape_type>)
    constexpr Tensor<Number, storage_type, n...>& operator-=(
        Tensor<Number, storage_type, n...>& tensor, const OtherTensor& B);

    /*!
     * @brief Scale the tensor with a scalar value
     *
     * @tparam OtherTensor
     */
    template <typename Scalar, typename Number, StorageType storage_type, std::size_t... n>
      requires(std::is_arithmetic_v<Scalar>)
    constexpr Tensor<Number, storage_type, n...>& operator*=(
        Tensor<Number, storage_type, n...>& tensor, const Scalar b);

    /*!
     * @brief Scales the tensor with the inverse of the scalar value b
     *
     * @tparam Scalar
     * @tparam Number
     * @tparam storage_type
     * @tparam n
     */
    template <typename Scalar, typename Number, StorageType storage_type, std::size_t... n>
      requires(std::is_arithmetic_v<Scalar>)
    constexpr Tensor<Number, storage_type, n...>& operator/=(
        Tensor<Number, storage_type, n...>& tensor, const Scalar b);

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
  }  // namespace Internal


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

  namespace Internal
  {
    template <typename OtherTensor, typename Number, StorageType storage_type, std::size_t... n>
      requires(std::is_same_v<typename OtherTensor::shape_type, typename OtherTensor::shape_type>)
    constexpr Tensor<Number, storage_type, n...>& operator+=(
        Tensor<Number, storage_type, n...>& tensor, const OtherTensor& B)
    {
      DenseFunctions::update<typename Tensor<Number, storage_type, n...>::value_type,
          OtherTensor::size(), 1>(1.0, tensor.data(), 1.0, B.data());
      return tensor;
    }

    template <typename OtherTensor, typename Number, StorageType storage_type, std::size_t... n>
      requires(std::is_same_v<typename OtherTensor::shape_type, typename OtherTensor::shape_type>)
    constexpr Tensor<Number, storage_type, n...>& operator-=(
        Tensor<Number, storage_type, n...>& tensor, const OtherTensor& B)
    {
      DenseFunctions::update<typename Tensor<Number, storage_type, n...>::value_type,
          OtherTensor::size(), 1>(1.0, tensor.data(), -1.0, B.data());
      return tensor;
    }

    template <typename Scalar, typename Number, StorageType storage_type, std::size_t... n>
      requires(std::is_arithmetic_v<Scalar>)
    constexpr Tensor<Number, storage_type, n...>& operator*=(
        Tensor<Number, storage_type, n...>& tensor, const Scalar b)
    {
      for (auto& value : tensor.container()) value *= b;
      return tensor;
    }

    template <typename Scalar, typename Number, StorageType storage_type, std::size_t... n>
      requires(std::is_arithmetic_v<Scalar>)
    constexpr Tensor<Number, storage_type, n...>& operator/=(
        Tensor<Number, storage_type, n...>& tensor, const Scalar b)
    {
      tensor *= Scalar(1.0) / b;
      return tensor;
    }
  }  // namespace Internal
}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE

#endif