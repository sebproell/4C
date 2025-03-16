// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_beam3_triad_interpolation.hpp"

#include "4C_beam3_triad_interpolation_local_rotation_vectors.hpp"
#include "4C_utils_exceptions.hpp"

#include <Sacado.hpp>

#include <memory>

FOUR_C_NAMESPACE_OPEN

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
LargeRotations::TriadInterpolation<T>::TriadInterpolation()
{
  // empty constructor
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
std::shared_ptr<LargeRotations::TriadInterpolation<T>>
LargeRotations::TriadInterpolation<T>::create(unsigned int numnodes)
{
  // so far, the only implemented variant is the one based on local rotation vectors

  switch (numnodes)
  {
    case 2:
    {
      return std::make_shared<LargeRotations::TriadInterpolationLocalRotationVectors<2, T>>();
    }
    case 3:
    {
      return std::make_shared<LargeRotations::TriadInterpolationLocalRotationVectors<3, T>>();
    }
    case 4:
    {
      return std::make_shared<LargeRotations::TriadInterpolationLocalRotationVectors<4, T>>();
    }
    case 5:
    {
      return std::make_shared<LargeRotations::TriadInterpolationLocalRotationVectors<5, T>>();
    }
    default:
    {
      FOUR_C_THROW("{} is no valid number of nodes used for triad interpolation! choose 2,3,4 or 5",
          numnodes);
      break;
    }
  }

  return nullptr;
}

// explicit template instantiations
template class LargeRotations::TriadInterpolation<double>;
template class LargeRotations::TriadInterpolation<Sacado::Fad::DFad<double>>;

FOUR_C_NAMESPACE_CLOSE
