// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FEM_GENERAL_ELEMENT_DEFINITION_HPP
#define FOUR_C_FEM_GENERAL_ELEMENT_DEFINITION_HPP

#include "4C_config.hpp"

#include "4C_fem_general_cell_type.hpp"
#include "4C_io_input_spec.hpp"

#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN


namespace Core::Elements
{
  /**
   * Collection of valid element input.
   *
   * This glass gathers all valid element definitions from the global state. This is possible
   * since the ElementType subclass register themselves in the global Factory.
   */
  struct ElementDefinition
  {
    //! Gather all valid element definitions from global state
    ElementDefinition();

    /**
     * Convenience access with check for existence of element definition.
     */
    const Core::IO::InputSpec& get(
        const std::string& element_name, const std::string& cell_type) const;

    //! Map from physics to cell type to InputSpec.
    std::map<std::string, std::map<std::string, Core::IO::InputSpec>> definitions;
  };

}  // namespace Core::Elements


FOUR_C_NAMESPACE_CLOSE

#endif
