# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

set(OBJS_FOUR_C_MAIN
    # cmake-format: sortable
    4C_global_full_cal_control.cpp
    4C_global_full_control.cpp
    4C_global_full_init_control.cpp
    4C_global_full_inp_control.cpp
    4C_global_full_main.cpp
    )

add_executable(${FOUR_C_EXECUTABLE_NAME} ${OBJS_FOUR_C_MAIN})
set_target_properties(
  ${FOUR_C_EXECUTABLE_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  )
four_c_set_up_executable(${FOUR_C_EXECUTABLE_NAME})

if(FOUR_C_ENABLE_METADATA_GENERATION)
  if(Python3_VERSION VERSION_GREATER_EQUAL "3.10")
    add_custom_command(
      TARGET ${FOUR_C_EXECUTABLE_NAME}
      POST_BUILD
      COMMAND
        ${CMAKE_COMMAND} -E rm -f ${PROJECT_BINARY_DIR}/4C_metadata.yaml
        ${PROJECT_BINARY_DIR}/4C_schema.json ${PROJECT_BINARY_DIR}/4C_schema_partial.json
      COMMAND
        ${FOUR_C_ENABLE_ADDRESS_SANITIZER_TEST_OPTIONS} "$<TARGET_FILE:${FOUR_C_EXECUTABLE_NAME}>"
        -p > ${PROJECT_BINARY_DIR}/4C_metadata.yaml
      COMMAND
        ${FOUR_C_PYTHON_VENV_BUILD}/bin/python ${PROJECT_SOURCE_DIR}/utilities/create_json_schema.py
        ${PROJECT_BINARY_DIR}/4C_metadata.yaml ${PROJECT_BINARY_DIR}/4C_schema.json
      COMMENT "Creating metadata and JSON schema files"
      )
  else()
    message(
      WARNING
        "Your python version ${Python3_VERSION} is too old to generate metadata and JSON schema files"
      )
    add_custom_command(
      TARGET ${FOUR_C_EXECUTABLE_NAME}
      POST_BUILD
      COMMAND
        ${CMAKE_COMMAND} -E echo
        "Skipping metadata and JSON schema generation due to outdated python version"
      )
  endif()
endif()
