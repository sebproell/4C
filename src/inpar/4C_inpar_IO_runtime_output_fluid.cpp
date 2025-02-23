// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_inpar_IO_runtime_output_fluid.hpp"

#include "4C_utils_parameter_list.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Inpar
{
  namespace IORuntimeOutput
  {
    namespace FLUID
    {
      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      void set_valid_parameters(std::map<std::string, Core::IO::InputSpec>& list)
      {
        using Teuchos::tuple;

        // related sublist
        Core::Utils::SectionSpecs sublist_IO{"IO"};
        Core::Utils::SectionSpecs sublist_IO_output{sublist_IO, "RUNTIME VTK OUTPUT"};
        Core::Utils::SectionSpecs sublist_IO_output_fluid{sublist_IO_output, "FLUID"};

        // whether to write output for fluid
        Core::Utils::bool_parameter(
            "OUTPUT_FLUID", false, "write fluid output", sublist_IO_output_fluid);

        // whether to write velocity state
        Core::Utils::bool_parameter(
            "VELOCITY", false, "write velocity output", sublist_IO_output_fluid);

        // whether to write pressure state
        Core::Utils::bool_parameter(
            "PRESSURE", false, "write pressure output", sublist_IO_output_fluid);

        // whether to write acceleration state
        Core::Utils::bool_parameter(
            "ACCELERATION", false, "write acceleration output", sublist_IO_output_fluid);

        // whether to write displacement state
        Core::Utils::bool_parameter(
            "DISPLACEMENT", false, "write displacement output", sublist_IO_output_fluid);

        // whether to write displacement state
        Core::Utils::bool_parameter(
            "GRIDVELOCITY", false, "write grid velocity output", sublist_IO_output_fluid);

        // whether to write element owner
        Core::Utils::bool_parameter(
            "ELEMENT_OWNER", false, "write element owner", sublist_IO_output_fluid);

        // whether to write element GIDs
        Core::Utils::bool_parameter(
            "ELEMENT_GID", false, "write 4C internal element GIDs", sublist_IO_output_fluid);

        // whether to write node GIDs
        Core::Utils::bool_parameter(
            "NODE_GID", false, "write 4C internal node GIDs", sublist_IO_output_fluid);

        sublist_IO_output_fluid.move_into_collection(list);
      }
    }  // namespace FLUID
  }  // namespace IORuntimeOutput
}  // namespace Inpar

FOUR_C_NAMESPACE_CLOSE
