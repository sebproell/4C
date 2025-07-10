// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fem_general_element_definition.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_io_input_parameter_container.hpp"
#include "4C_io_input_spec.hpp"

FOUR_C_NAMESPACE_OPEN


Core::Elements::ElementDefinition::ElementDefinition()
{
  Core::Communication::ParObjectFactory::instance().setup_element_definition(definitions);
}

const Core::IO::InputSpec& Core::Elements::ElementDefinition::get(
    const std::string& element_name, Core::FE::CellType cell_type) const
{
  auto it = definitions.find(element_name);
  if (it == definitions.end()) FOUR_C_THROW("No element '{}' found.", element_name);
  auto it2 = it->second.find(cell_type);
  if (it2 == it->second.end())
    FOUR_C_THROW("Element '{}' does not seem to know cell type '{}'.", element_name, cell_type);
  return it2->second;
}

FOUR_C_NAMESPACE_CLOSE
