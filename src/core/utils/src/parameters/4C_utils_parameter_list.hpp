// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_UTILS_PARAMETER_LIST_HPP
#define FOUR_C_UTILS_PARAMETER_LIST_HPP

#include "4C_config.hpp"

#include "4C_io_input_spec_builders.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core
{
  namespace Utils
  {
    //! add entry as item of enum class @p value to @p list with name @p parameter_name
    template <class EnumType>
    void add_enum_class_to_parameter_list(
        const std::string& parameter_name, const EnumType value, Teuchos::ParameterList& list)
    {
      const std::string docu = "";
      const std::string value_name = "val";
      Teuchos::setStringToIntegralParameter<EnumType>(parameter_name, value_name, docu,
          Teuchos::tuple<std::string>(value_name), Teuchos::tuple<EnumType>(value), &list);
    }

  }  // namespace Utils
}  // namespace Core


FOUR_C_NAMESPACE_CLOSE

#endif
