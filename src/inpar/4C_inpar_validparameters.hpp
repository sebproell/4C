// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_INPAR_VALIDPARAMETERS_HPP
#define FOUR_C_INPAR_VALIDPARAMETERS_HPP

#include "4C_config.hpp"

#include "4C_io_input_spec.hpp"

#include <iostream>
#include <map>
#include <memory>
#include <string>

FOUR_C_NAMESPACE_OPEN

namespace Input
{
  /**
   * The valid parameters per section.
   */
  std::map<std::string, Core::IO::InputSpec> valid_parameters();

}  // namespace Input


/*! print help message */
void print_help_message();

/*! print flag sections of dat file with default flags */
void print_default_dat_header();


FOUR_C_NAMESPACE_CLOSE

#endif
