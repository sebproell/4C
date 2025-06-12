// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LINALG_SYMMETRIC_TENSOR_EIGEN_HPP
#define FOUR_C_LINALG_SYMMETRIC_TENSOR_EIGEN_HPP

#include "4C_config.hpp"

#include "4C_linalg_symmetric_tensor.hpp"

#include <Teuchos_LAPACK.hpp>


FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  /*!
   * @brief Computes the eigenvalue decomposition of a symmetric 2-tensor
   *
   * @tparam Tensor
   * @return A tuple of a std::array containing the eigenvalues and a 2-tensor containing the
   * eigenvectors.
   */
  template <typename Tensor>
    requires(is_symmetric_tensor<Tensor> && Tensor::rank() == 2)
  auto eig(const Tensor& t)
  {
    // Lapack uses this tensor to compute the eigenvectors
    auto eigenvectors = get_full(t);

    constexpr std::size_t size = Tensor::template extent<0>();


    std::array<double, size> eigenvalues;
    // ----- perform eigendecomposition ----- //
    const int lwork = 2 * size * size + 6 * size + 1;
    std::array<double, lwork> work;
    int info;
    Teuchos::LAPACK<int, double> lapack;
    lapack.SYEV(
        'V', 'U', size, eigenvectors.data(), size, eigenvalues.data(), work.data(), lwork, &info);
    FOUR_C_ASSERT_ALWAYS(info == 0, "Eigenvalue decomposition failed with error code {}", info);

    return std::make_tuple(eigenvalues, eigenvectors);
  }
}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE

#endif