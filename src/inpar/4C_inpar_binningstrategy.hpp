// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_INPAR_BINNINGSTRATEGY_HPP
#define FOUR_C_INPAR_BINNINGSTRATEGY_HPP

#include "4C_config.hpp"

#include "4C_io_input_spec.hpp"

#include <map>

FOUR_C_NAMESPACE_OPEN

namespace Inpar
{
  namespace BINSTRATEGY
  {

    /// set the binning strategy parameters
    void set_valid_parameters(std::map<std::string, Core::IO::InputSpec>& list);

  }  // namespace BINSTRATEGY

}  // namespace Inpar

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
