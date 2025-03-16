// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FEM_GENERAL_UTILS_GDER2_HPP
#define FOUR_C_FEM_GENERAL_UTILS_GDER2_HPP

#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_linalg_fixedsizematrix.hpp"

FOUR_C_NAMESPACE_OPEN


namespace Core::FE
{
  /*---------------------------------------------------------------------*/
  /*!
   * \brief calculate second global derivatives w.r.t x/y/z at point r,s,t
   */
  /*                                                          gammi 07/07
  |
  | To compute the second derivative, the first derivative is derived once more.
  | This is an application of the product rule and the chain rule.
  | The trick to compute it is to realize that the transformation dr/dx must be
  | derived w.r.t. dx (but dr/dx = (dx/dr)^-1), which gives the term -(dx/dr)^{-2}
  | plus some inner derivatives.
  |
  | This formula is valid in all dimension if one interprets x and r as tensors
  | and interprets multiplication with (dx/dr)^-1 from the right is as a matrix
  | multiplication with (dx/dr)^-T.
  |
  |             +-            -+
  |  d^2N    d  |  /dx\-1   dN |   d   /dx\-1   dN     /dx\    d  dN
  |  ----  = -- | | -- |  * -- | = -- | -- |  * --  + | -- | * -- --
  |  dx^2    dx |  \dr/     dr |   dx  \dr/     dr     \dr/    dx dr
  |             +-            -+
  |
  |             /dx\-2 dN d  dx    /dx\-1   d^2N    /dx\-1
  |        = - | -- | *--*-- -- + | -- |  * ---- * | -- |
  |             \dr/   dr dx dr    \dr/     dr^2    \dr/

  |             /dx\-1 dN d^2x  /dx\-1     /dx\-1 d^2N  /dx\-1
  |        = - | -- | *--*----*| -- |   + | -- | *----*| -- |
  |             \dr/   dr dr^2  \dr/       \dr/   dr^2  \dr/

  |                    +-                -+
  |           /dx\-1   |   dN d^2x   d^2N |    /dx\-1
  |        = | -- |  * | - --*---- + ---- | * | -- |
  |           \dr/     |   dx dr^2   dr^2 |    \dr/
  |                    +-                -+
  *----------------------------------------------------------------------*/
  // for enriched elements (e.g. xwall), num_node may be larger than the number of element nodes
  // for all other elements, num_node==numnode
  template <Core::FE::CellType distype, int num_node, int prob_dim>
  void gder2(const Core::LinAlg::Matrix<prob_dim, prob_dim>& xjm,
      const Core::LinAlg::Matrix<prob_dim, num_node>& derxy,
      const Core::LinAlg::Matrix<Core::FE::DisTypeToNumDeriv2<distype>::numderiv2, num_node>&
          deriv2,
      const Core::LinAlg::Matrix<prob_dim, num_node>& xyze,
      Core::LinAlg::Matrix<Core::FE::DisTypeToNumDeriv2<distype>::numderiv2, num_node>& derxy2)
  {
    // some numbers already known during compilation
    const int numnode = Core::FE::num_nodes<distype>;
    FOUR_C_ASSERT(numnode <= num_node, "Expect at least numNodePerElement matrix columns");
    const int nsd = prob_dim;
    const int numderiv2 = Core::FE::DisTypeToNumDeriv2<distype>::numderiv2;

    // compute d^2x/dr^2
    double xder2[numderiv2 * nsd];
    Core::LinAlg::DenseFunctions::multiply_nt<double, numderiv2, numnode, nsd>(
        xder2, deriv2.data(), xyze.data());

    // compute -(dN/dx)*(d^2x/dr^2)
    Core::LinAlg::DenseFunctions::multiply<double, numderiv2, nsd, numnode>(
        derxy2.data(), -1.0, xder2, derxy.data());

    // compute -(dN/dx)*(d^2x/dr^2) + (d^2N/dr^2)
    derxy2.update(1.0, deriv2, 1.0);

    // finally multiply by (dx/dr)^-1 from the left and (dx/dr)^-T from the right
    // write out the products by hand because derxy2 only stores the symmetric part
    // of the second derivative tensor (which means we cannot use matrix-matrix products)
    double xjiData[9];
    double* xji[3];
    for (int i = 0; i < nsd; ++i) xji[i] = &xjiData[nsd * i];
    Core::LinAlg::DenseFunctions::invert<double, nsd, nsd>(xjiData, xjm.data());
    for (int node = 0; node < numnode; ++node)
    {
      double tmp[3][3];
      switch (nsd)
      {
        case 3:
          for (int j = 0; j < nsd; ++j)
          {
            tmp[0][j] = derxy2(0, node) * xji[0][j] + derxy2(3, node) * xji[1][j] +
                        derxy2(4, node) * xji[2][j];
            tmp[1][j] = derxy2(3, node) * xji[0][j] + derxy2(1, node) * xji[1][j] +
                        derxy2(5, node) * xji[2][j];
            tmp[2][j] = derxy2(4, node) * xji[0][j] + derxy2(5, node) * xji[1][j] +
                        derxy2(2, node) * xji[2][j];
          }
          derxy2(0, node) = xji[0][0] * tmp[0][0] + xji[1][0] * tmp[1][0] + xji[2][0] * tmp[2][0];
          derxy2(1, node) = xji[0][1] * tmp[0][1] + xji[1][1] * tmp[1][1] + xji[2][1] * tmp[2][1];
          derxy2(2, node) = xji[0][2] * tmp[0][2] + xji[1][2] * tmp[1][2] + xji[2][2] * tmp[2][2];
          derxy2(3, node) = xji[0][0] * tmp[0][1] + xji[1][0] * tmp[1][1] + xji[2][0] * tmp[2][1];
          derxy2(4, node) = xji[0][0] * tmp[0][2] + xji[1][0] * tmp[1][2] + xji[2][0] * tmp[2][2];
          derxy2(5, node) = xji[0][1] * tmp[0][2] + xji[1][1] * tmp[1][2] + xji[2][1] * tmp[2][2];
          break;
        case 2:
          for (int j = 0; j < nsd; ++j)
          {
            tmp[0][j] = derxy2(0, node) * xji[0][j] + derxy2(2, node) * xji[1][j];
            tmp[1][j] = derxy2(2, node) * xji[0][j] + derxy2(1, node) * xji[1][j];
          }
          derxy2(0, node) = xji[0][0] * tmp[0][0] + xji[1][0] * tmp[1][0];
          derxy2(1, node) = xji[0][1] * tmp[0][1] + xji[1][1] * tmp[1][1];
          derxy2(2, node) = xji[0][0] * tmp[0][1] + xji[1][0] * tmp[1][1];
          break;
        case 1:
          derxy2(0, node) *= xji[0][0] * xji[0][0];
          break;
        default:
          FOUR_C_THROW("Illegal number of space dimensions: {}", nsd);
          break;
      }
    }

    return;

  }  // Core::FE::gder2

  template <Core::FE::CellType distype, int num_node>
  void gder2(const Core::LinAlg::Matrix<Core::FE::dim<distype>, Core::FE::dim<distype>>& xjm,
      const Core::LinAlg::Matrix<Core::FE::dim<distype>, num_node>& derxy,
      const Core::LinAlg::Matrix<Core::FE::DisTypeToNumDeriv2<distype>::numderiv2, num_node>&
          deriv2,
      const Core::LinAlg::Matrix<Core::FE::dim<distype>, num_node>& xyze,
      Core::LinAlg::Matrix<Core::FE::DisTypeToNumDeriv2<distype>::numderiv2, num_node>& derxy2)
  {
    gder2<distype, num_node, Core::FE::dim<distype>>(xjm, derxy, deriv2, xyze, derxy2);
  }


}  // namespace Core::FE

FOUR_C_NAMESPACE_CLOSE

#endif
