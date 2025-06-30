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
  GIT_TAG 47ec2fa184209687c20fd5bc05621e1cb1200311 # version 0.9.0
  )
set(RYML_INSTALL
    ON
    CACHE BOOL "Turn on ryml install" FORCE
    )
fetchcontent_makeavailable(ryml)
set(FOUR_C_RYML_ROOT "${CMAKE_INSTALL_PREFIX}")

four_c_add_external_dependency(four_c_all_enabled_external_dependencies ryml::ryml)

configure_file(
  ${PROJECT_SOURCE_DIR}/cmake/templates/ryml.cmake.in
  ${PROJECT_BINARY_DIR}/cmake/templates/ryml.cmake
  @ONLY
  )
