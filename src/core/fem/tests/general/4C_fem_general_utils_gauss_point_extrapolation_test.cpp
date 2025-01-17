// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_fem_general_utils_gauss_point_extrapolation.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_fem_general_element_integration_select.hpp"
#include "4C_fem_general_utils_gausspoints.hpp"
#include "4C_fem_general_utils_integration.hpp"

#include <vector>

namespace
{
  using namespace FourC;

  TEST(ElementServiceTest, TestGaussPointProjectionMatrixHex8)
  {
    constexpr Core::FE::CellType distype = Core::FE::CellType::hex8;
    constexpr int nsd = 3;

    Core::FE::IntPointsAndWeights<nsd> intpoints(
        Discret::Elements::DisTypeToOptGaussRule<distype>::rule);

    // format as Discret::Utils::GaussIntegration
    std::shared_ptr<Core::FE::CollectedGaussPoints> gp =
        std::make_shared<Core::FE::CollectedGaussPoints>();

    std::array<double, nsd> xi{};
    for (int i = 0; i < intpoints.ip().nquad; ++i)
    {
      for (int d = 0; d < nsd; ++d) xi[d] = intpoints.ip().qxg[i][d];
      gp->append(xi[0], xi[1], xi[2], intpoints.ip().qwgt[i]);
    }

    // save default integration rule
    Core::FE::GaussIntegration integration(gp);
    Core::LinAlg::SerialDenseMatrix m =
        Core::FE::evaluate_gauss_points_to_nodes_extrapolation_matrix<distype>(integration);
  }

}  // namespace
