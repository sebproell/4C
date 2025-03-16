// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fem_general_utils_bspline.hpp"

#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

//--------------------------------------------------
// Constructor
//--------------------------------------------------
Core::FE::Nurbs::BsplinePolynomial::BsplinePolynomial(
    const int degree, const Core::LinAlg::SerialDenseVector local_knotvector)
    : myknotvector_(local_knotvector),
      bspline_(degree + 1),
      degree_(degree),
      degree_plus_one_(degree + 1)
{
  return;
}

//--------------------------------------------------
// Copy constructor
//--------------------------------------------------
Core::FE::Nurbs::BsplinePolynomial::BsplinePolynomial(const BsplinePolynomial& old)
    : degree_(old.degree_)
{
  myknotvector_ = old.myknotvector_;
  return;
}


FOUR_C_NAMESPACE_CLOSE
