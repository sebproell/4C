// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LINALG_UTILS_DENSEMATRIX_DETERMINANT_HPP
#define FOUR_C_LINALG_UTILS_DENSEMATRIX_DETERMINANT_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_serialdensematrix.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  /**
   *\brief Explicit determinant of a nonsymmetric 2x2 matrix.
   *
   * @param A (in) Matrix A.
   * @return determinant of matrix A.
   */
  template <typename T>
  T determinant(const Core::LinAlg::Matrix<2, 2, T>& A)
  {
    T b00 = A(0, 0);
    T b01 = A(0, 1);
    T b10 = A(1, 0);
    T b11 = A(1, 1);
    T det = b00 * b11 - b01 * b10;
    return det;
  }

  /**
   *\brief Explicit determinant of a nonsymmetric 3x3 matrix.
   *
   * @param A (in) Matrix A.
   * @return determinant of matrix A.
   */
  template <typename T>
  T determinant(const Core::LinAlg::Matrix<3, 3, T>& A)
  {
    T b00 = A(0, 0);
    T b01 = A(0, 1);
    T b02 = A(0, 2);
    T b10 = A(1, 0);
    T b11 = A(1, 1);
    T b12 = A(1, 2);
    T b20 = A(2, 0);
    T b21 = A(2, 1);
    T b22 = A(2, 2);
    T a = b11 * b22 - b21 * b12;
    T b = -b10 * b22 + b20 * b12;
    T c = b10 * b21 - b20 * b11;
    T det = b00 * a + b01 * b + b02 * c;
    return det;
  }

  /**
   *\brief Explicit determinant of a nonsymmetric 4x4 matrix.
   *
   * @param A (in) Matrix A.
   * @return determinant of matrix A.
   */
  template <typename T>
  T determinant(const Core::LinAlg::Matrix<4, 4, T>& A)
  {
    T a00 = A(0, 0);
    T a01 = A(0, 1);
    T a02 = A(0, 2);
    T a03 = A(0, 3);

    T a10 = A(1, 0);
    T a11 = A(1, 1);
    T a12 = A(1, 2);
    T a13 = A(1, 3);

    T a20 = A(2, 0);
    T a21 = A(2, 1);
    T a22 = A(2, 2);
    T a23 = A(2, 3);

    T a30 = A(3, 0);
    T a31 = A(3, 1);
    T a32 = A(3, 2);
    T a33 = A(3, 3);

    T det = a00 * a11 * a33 * a22 - a00 * a11 * a32 * a23 - a00 * a31 * a13 * a22 +
            a00 * a32 * a21 * a13 + a00 * a31 * a12 * a23 - a00 * a33 * a21 * a12 +
            a11 * a30 * a02 * a23 - a11 * a33 * a20 * a02 + a11 * a32 * a20 * a03 -
            a11 * a30 * a03 * a22 - a31 * a10 * a02 * a23 - a33 * a22 * a10 * a01 +
            a33 * a21 * a10 * a02 + a33 * a20 * a01 * a12 + a30 * a03 * a21 * a12 +
            a31 * a13 * a20 * a02 + a31 * a10 * a03 * a22 + a30 * a01 * a13 * a22 +
            a32 * a23 * a10 * a01 - a32 * a21 * a10 * a03 - a32 * a20 * a01 * a13 -
            a30 * a02 * a21 * a13 - a31 * a12 * a20 * a03 - a30 * a01 * a12 * a23;
    return det;
  }

  /*!
  \brief determinant of a nonsymmetric matrix using LU factorization

  \return the determinant of the matrix A
  */
  double determinant_lu(const Core::LinAlg::SerialDenseMatrix& A);

}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE

#endif
