// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LINALG_TENSOR_SVD_HPP
#define FOUR_C_LINALG_TENSOR_SVD_HPP

#include "4C_config.hpp"

#include "4C_linalg_tensor.hpp"
#include "4C_linalg_tensor_internals.hpp"

#include <Teuchos_LAPACK.hpp>

#include <type_traits>


FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  /*!
   * @brief A struct holding the decomposition of a 2-tensor into its singular value decomposition
   * (SVD) such that tensor = Q * make_rectangular_diagonal_matrix(S) * VT.
   */
  template <typename Scalar, std::size_t rows, std::size_t cols>
  struct SVDDecomposition
  {
    //! Orthogonal matrix Q
    Core::LinAlg::Tensor<Scalar, rows, rows> Q{};

    //!< Array of singular values
    std::array<Scalar, std::min(rows, cols)> S{};

    //!< Orthogonal matrix VT
    Core::LinAlg::Tensor<Scalar, cols, cols> VT{};
  };

  /*!
   * @brief Computes the singular value decomposition of a 2-tensor
   *
   * @tparam Tensor
   * @return @p SVDDecomposition holding @p Q, @p s and @p VT, where @p Q and @p VT are orthogonal
   * 2-tensors and @p s is a std::array containing the singular valurs, such that t = Q *
   * make_rectangular_diagonal_matrix(s) * VT.
   */
  template <typename Tensor>
    requires(!is_compressed_tensor<Tensor> && Tensor::rank() == 2)
  auto svd(const Tensor& t)
  {
    using ValueType = std::remove_cvref_t<typename Tensor::value_type>;
    constexpr std::size_t rows = Tensor::template extent<0>();
    constexpr std::size_t cols = Tensor::template extent<1>();

    SVDDecomposition<ValueType, rows, cols> svd_composition{};

    Core::LinAlg::Tensor<ValueType, rows, cols> t_copy = t;

    constexpr int lwork =
        std::max(3 * std::min(rows, cols) + std::max(rows, cols), 5 * std::min(rows, cols));
    std::array<ValueType, lwork> work;
    ValueType rwork;
    int info;
    Teuchos::LAPACK<int, ValueType> lapack;
    lapack.GESVD('A', 'A', rows, cols, t_copy.data(), rows, svd_composition.S.data(),
        svd_composition.Q.data(), rows, svd_composition.VT.data(), cols, work.data(), lwork, &rwork,
        &info);
    FOUR_C_ASSERT_ALWAYS(info == 0, "Singular value decomposition failed with error code {}", info);

    return svd_composition;
  }

  /**
   * @brief Creates a rectangular diagonal matrix from a given diagonal array.
   *
   * Constructs a matrix of dimensions `rows` x `cols` where the diagonal elements
   * are initialized from the provided `diag` array, and all other elements are default-initialized.
   * The length of the diagonal is the minimum of `rows` and `cols`.
   *
   * @tparam rows Number of rows in the resulting matrix.
   * @tparam cols Number of columns in the resulting matrix.
   * @tparam ValueType Type of the matrix elements.
   * @param diag Array containing the diagonal elements. Its size must be `std::min(rows, cols)`.
   * @return Core::LinAlg::Tensor<ValueType, rows, cols> The resulting rectangular diagonal matrix.
   */
  template <std::size_t rows, std::size_t cols, typename ValueType>
  auto make_rectangular_diagonal_matrix(const std::array<ValueType, std::min(rows, cols)>& diag)
  {
    Core::LinAlg::Tensor<ValueType, rows, cols> result{};
    for (std::size_t i = 0; i < std::min(rows, cols); ++i)
    {
      result(i, i) = diag[i];
    }
    return result;
  }
}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE

#endif