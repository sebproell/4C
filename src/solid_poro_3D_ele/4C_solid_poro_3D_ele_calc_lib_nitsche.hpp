// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLID_PORO_3D_ELE_CALC_LIB_NITSCHE_HPP
#define FOUR_C_SOLID_PORO_3D_ELE_CALC_LIB_NITSCHE_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_solid_3D_ele_calc_lib_nitsche.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret::Elements
{

  template <int dim>
  struct SolidPoroCauchyNDirLinearizations
  {
    /// all pure solid linearizations
    CauchyNDirLinearizations<dim> solid{};

    /// first derivative w.r.t. pressure
    Core::LinAlg::SerialDenseMatrix* d_cauchyndir_dp = nullptr;
  };
}  // namespace Discret::Elements

FOUR_C_NAMESPACE_CLOSE
#endif