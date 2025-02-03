// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_GLOBAL_LEGACY_MODULE_HPP
#define FOUR_C_GLOBAL_LEGACY_MODULE_HPP

#include "4C_config.hpp"

#include "4C_module_registry_callbacks.hpp"

FOUR_C_NAMESPACE_OPEN

/**
 * Expose callbacks for the global legacy module.
 *
 * @note This function exists for historic reasons. Its internals need to be split up to remove
 * forced dependencies on the main apps.
 */
[[nodiscard]] ModuleCallbacks global_legacy_module_callbacks();

/**
 * Write metadata for input files to @p out.
 *
 * @note This function exists for historic reasons. Its internals need to be split up to remove
 * forced dependencies on the main apps.
 */
void write_input_metadata(std::ostream& out);

FOUR_C_NAMESPACE_CLOSE

#endif
