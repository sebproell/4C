// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_UTILS_STD23_UNREACHABLE_HPP
#define FOUR_C_UTILS_STD23_UNREACHABLE_HPP

#include "4C_config.hpp"

FOUR_C_NAMESPACE_OPEN

namespace std23
{
  /**
   * Marks a point in the code as unreachable.
   *
   * See https://en.cppreference.com/w/cpp/utility/unreachable for more information.
   */
  [[noreturn]] inline void unreachable()
  {
#if defined(__GNUC__) || defined(__clang__)
    __builtin_unreachable();
#endif
  }
}  // namespace std23

FOUR_C_NAMESPACE_CLOSE

#endif