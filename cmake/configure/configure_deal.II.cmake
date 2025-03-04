# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

# At this point we deliberately turned on the option to use deal.II, so fail if it is not available.
find_package(deal.II REQUIRED)

# Assume we always get full Trilinos from deal.II
if(NOT DEAL_II_WITH_TRILINOS)
  message(
    FATAL_ERROR
      "deal.II was not compiled with Trilinos support. Please recompile deal.II with Trilinos support."
    )
endif()

target_link_libraries(four_c_all_enabled_external_dependencies INTERFACE dealii::dealii)

# TODO installation
