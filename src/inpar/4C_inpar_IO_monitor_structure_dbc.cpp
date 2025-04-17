// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_inpar_IO_monitor_structure_dbc.hpp"

#include "4C_io_input_spec_builders.hpp"


FOUR_C_NAMESPACE_OPEN

namespace Inpar
{
  namespace IOMonitorStructureDBC
  {
    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    void set_valid_parameters(std::map<std::string, Core::IO::InputSpec>& list)
    {
      using namespace Core::IO::InputSpecBuilders;

      // related sublist
      list["IO/MONITOR STRUCTURE DBC"] = group("IO/MONITOR STRUCTURE DBC",
          {

              // output interval regarding steps: write output every INTERVAL_STEPS steps
              parameter<int>("INTERVAL_STEPS",
                  {.description = "write reaction force output every INTERVAL_STEPS steps",
                      .default_value = -1}),

              // precision for file
              parameter<int>("PRECISION_FILE",
                  {.description = "precision for written file", .default_value = 16}),

              // precision for screen
              parameter<int>("PRECISION_SCREEN",
                  {.description = "precision for written screen output", .default_value = 5}),

              // type of written output file
              deprecated_selection<Inpar::IOMonitorStructureDBC::FileType>("FILE_TYPE",
                  {
                      {"csv", Inpar::IOMonitorStructureDBC::csv},
                      {"CSV", Inpar::IOMonitorStructureDBC::csv},
                      {"Csv", Inpar::IOMonitorStructureDBC::csv},
                      {"data", Inpar::IOMonitorStructureDBC::data},
                      {"Data", Inpar::IOMonitorStructureDBC::data},
                      {"DATA", Inpar::IOMonitorStructureDBC::data},
                  },
                  {.description = "type of written output file",
                      .default_value = Inpar::IOMonitorStructureDBC::csv}),

              // whether to write output in every iteration of the nonlinear solver
              parameter<bool>("WRITE_HEADER",
                  {.description =
                          "write information about monitored boundary condition to output file",
                      .default_value = false})},
          {.defaultable = true});
    }

    std::string to_string(FileType type)
    {
      switch (type)
      {
        case FileType::csv:
          return "csv";
        case FileType::data:
          return "data";
        default:
          FOUR_C_THROW("Unknown file type");
      }
    }
  }  // namespace IOMonitorStructureDBC
}  // namespace Inpar

FOUR_C_NAMESPACE_CLOSE