// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_linalg_utils_densematrix_determinant.hpp"

#include "4C_utils_exceptions.hpp"

#include <Teuchos_LAPACK.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double Core::LinAlg::determinant_lu(const Core::LinAlg::SerialDenseMatrix& A)
{
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (A.numRows() != A.numCols()) FOUR_C_THROW("Matrix is not square");
#endif
  Core::LinAlg::SerialDenseMatrix tmp(A);
  const int n = tmp.numCols();
  const int m = tmp.numRows();
  std::vector<int> ipiv(n);
  int info;

  Teuchos::LAPACK<int, double> lapack;
  lapack.GETRF(m, n, tmp.values(), tmp.stride(), ipiv.data(), &info);

  if (info < 0)
    FOUR_C_THROW("Lapack's dgetrf returned {}", info);
  else if (info > 0)
    return 0.0;
  double d = tmp(0, 0);
  for (int i = 1; i < n; ++i) d *= tmp(i, i);
  // swapping rows of A changes the sign of the determinant, so we have to
  // undo lapack's permutation w.r.t. the determinant
  // note the fortran indexing convention in ipiv
  for (int i = 0; i < n; ++i)
    if (ipiv[i] != i + 1) d *= -1.0;
  return d;
}

FOUR_C_NAMESPACE_CLOSE
