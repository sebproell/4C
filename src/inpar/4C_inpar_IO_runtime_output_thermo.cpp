// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_inpar_IO_runtime_output_thermo.hpp"

#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Inpar
{
  namespace IORuntimeOutput
  {
    namespace Thermo
    {
      void set_valid_parameters(Teuchos::ParameterList& list)
      {
        using Teuchos::setStringToIntegralParameter;
        using Teuchos::tuple;

        // related sublist
        Teuchos::ParameterList& sublist_IO = list.sublist("IO", false, "");
        Teuchos::ParameterList& sublist_IO_output =
            sublist_IO.sublist("RUNTIME VTK OUTPUT", false, "");
        Teuchos::ParameterList& sublist_IO_output_thermo =
            sublist_IO_output.sublist("THERMO", false, "");

        // whether to write output for thermo
        Core::Utils::bool_parameter(
            "OUTPUT_THERMO", "No", "write thermo output", &sublist_IO_output_thermo);

        // whether to write temperature state
        Core::Utils::bool_parameter(
            "TEMPERATURE", "No", "write temperature output", &sublist_IO_output_thermo);

        // whether to write heatflux state
        Core::Utils::bool_parameter(
            "HEATFLUX", "No", "write heatflux output", &sublist_IO_output_thermo);

        // whether to write temperature gradient state
        Core::Utils::bool_parameter(
            "TEMPGRAD", "No", "write temperature gradient output", &sublist_IO_output_thermo);

        // whether to write energy state
        Core::Utils::bool_parameter(
            "ENERGY", "No", "write energy output", &sublist_IO_output_thermo);

        // whether to write element owner
        Core::Utils::bool_parameter(
            "ELEMENT_OWNER", "No", "write element owner", &sublist_IO_output_thermo);

        // whether to write element GIDs
        Core::Utils::bool_parameter(
            "ELEMENT_GID", "No", "write 4C internal element GIDs", &sublist_IO_output_thermo);

        // whether to write node GIDs
        Core::Utils::bool_parameter(
            "NODE_GID", "No", "write 4C internal node GIDs", &sublist_IO_output_thermo);
      }
    }  // namespace Thermo
  }  // namespace IORuntimeOutput
}  // namespace Inpar

FOUR_C_NAMESPACE_CLOSE
