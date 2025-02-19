# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

message(STATUS "Fetch content for ryml")
fetchcontent_declare(
  ryml
  GIT_REPOSITORY https://github.com/biojppm/rapidyaml.git
  GIT_TAG e65999dc4d65368e1b3eb770e28d775885b0f525 # version 0.8.0
  )
set(RYML_INSTALL
    ON
    CACHE BOOL "Turn on ryml install" FORCE
    )
fetchcontent_makeavailable(ryml)
set(FOUR_C_RYML_ROOT "${CMAKE_INSTALL_PREFIX}/lib/cmake/ryml")

four_c_add_external_dependency(four_c_all_enabled_external_dependencies ryml::ryml)

configure_file(
  ${CMAKE_SOURCE_DIR}/cmake/templates/ryml.cmake.in
  ${CMAKE_BINARY_DIR}/cmake/templates/ryml.cmake
  @ONLY
  )
