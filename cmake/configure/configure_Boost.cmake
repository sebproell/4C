# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

# We are using the following Boost components:
#  - graph (a compiled library, so we need to list the component)
#  - stacktrace (a header-only library, so we don't need to list the component)
find_package(
  Boost
  COMPONENTS graph
  REQUIRED
  )

# post-process found targets
if(Boost_FOUND)
  message(STATUS "Boost component libraries: ${Boost_LIBRARIES}")

  target_compile_definitions(
    Boost::graph
    INTERFACE "-DBOOST_MAJOR_VERSION=${Boost_MAJOR_VERSION}"
              "-DBOOST_MINOR_VERSION=${Boost_MINOR_VERSION}"
    )

  target_link_libraries(four_c_all_enabled_external_dependencies INTERFACE Boost::graph)

  configure_file(
    ${PROJECT_SOURCE_DIR}/cmake/templates/Boost.cmake.in
    ${PROJECT_BINARY_DIR}/cmake/templates/Boost.cmake
    @ONLY
    )
endif()
