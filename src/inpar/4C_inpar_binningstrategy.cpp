// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_inpar_binningstrategy.hpp"

#include "4C_binstrategy.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN



void Inpar::BINSTRATEGY::set_valid_parameters(std::map<std::string, Core::IO::InputSpec>& list)
{
  using Teuchos::tuple;
  using namespace Core::IO::InputSpecBuilders;

  Core::Utils::SectionSpecs binningstrategy{"BINNING STRATEGY"};


  Core::Utils::double_parameter("BIN_SIZE_LOWER_BOUND", -1.0,
      "Lower bound for bin size. Exact bin size is computed via (Domain edge "
      "length)/BIN_SIZE_LOWER_BOUND. This also determines the number of bins in each spatial "
      "direction",
      binningstrategy);

  Core::Utils::string_parameter("BIN_PER_DIR", "-1 -1 -1",
      "Number of bins per direction (x, y, z) in particle simulations. Either Define this value or "
      "BIN_SIZE_LOWER_BOUND",
      binningstrategy);

  Core::Utils::string_parameter("PERIODICONOFF", "0 0 0",
      "Turn on/off periodic boundary conditions in each spatial direction", binningstrategy);

  Core::Utils::string_parameter("DOMAINBOUNDINGBOX", "1e12 1e12 1e12 1e12 1e12 1e12",
      "Bounding box for computational domain using binning strategy. Specify diagonal corner "
      "points",
      binningstrategy);

  Core::Utils::string_to_integral_parameter<Core::Binstrategy::WriteBins>("WRITEBINS", "none",
      "Write none, row or column bins for visualization",
      tuple<std::string>("none", "rows", "cols"),
      tuple<Core::Binstrategy::WriteBins>(Core::Binstrategy::WriteBins::none,
          Core::Binstrategy::WriteBins::rows, Core::Binstrategy::WriteBins::cols),
      binningstrategy);

  binningstrategy.move_into_collection(list);
}

FOUR_C_NAMESPACE_CLOSE
