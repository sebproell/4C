// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_IO_INPUT_PARAMETER_CONTAINER_TEMPLATES_HPP
#define FOUR_C_IO_INPUT_PARAMETER_CONTAINER_TEMPLATES_HPP

#include "4C_config.hpp"

#include "4C_io_input_parameter_container.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

// The add() function requires expensive-to-include Teuchos headers, but it is rarely needed, since
// most user code only reads from the InputParameterContainer. Therefore, we put the implementation
// into this special templates file.

template <typename T>
void Core::IO::InputParameterContainer::add(const std::string& name, const T& data)
{
  entries_[name] = {.data = std::any{data}};
  ensure_type_action_registered<T>();
}

template <typename T>
void Core::IO::InputParameterContainer::ensure_type_action_registered()
{
  auto& type_actions = get_type_actions();
  if (!type_actions.contains(typeid(T)))
  {
    type_actions[typeid(T)] = {
        .print = Internal::InputParameterContainerImplementation::PrintHelper<T>{},
        .write_to_pl = [](Teuchos::ParameterList& pl, const std::string& name, const std::any& data)
        { pl.set<T>(name, std::any_cast<T>(data)); },
    };
  }
}

// --- Declare that these templates are instantiated for some types --- //

extern template void Core::IO::InputParameterContainer::add(
    const std::string& name, const int& data);
extern template void Core::IO::InputParameterContainer::add(
    const std::string& name, const double& data);
extern template void Core::IO::InputParameterContainer::add(
    const std::string& name, const std::string& data);
extern template void Core::IO::InputParameterContainer::add(
    const std::string& name, const bool& data);

extern template void Core::IO::InputParameterContainer::add(
    const std::string& name, const std::vector<int>& data);
extern template void Core::IO::InputParameterContainer::add(
    const std::string& name, const std::vector<double>& data);
extern template void Core::IO::InputParameterContainer::add(
    const std::string& name, const std::vector<std::string>& data);
extern template void Core::IO::InputParameterContainer::add(
    const std::string& name, const std::vector<bool>& data);

FOUR_C_NAMESPACE_CLOSE

#endif
