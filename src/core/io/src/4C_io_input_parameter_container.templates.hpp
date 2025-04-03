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
#include "4C_utils_enum.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

// The add() function requires expensive-to-include Teuchos headers, but it is rarely needed, since
// most user code only reads from the InputParameterContainer. Therefore, we put the implementation
// into this special templates file.
namespace Core::IO::Internal::InputParameterContainerImplementation
{
  // Default printer if not printable.
  template <typename T>
  struct PrintHelper
  {
    void operator()(std::ostream& os, const std::any& data) { os << "<not printable> "; }
  };

  template <typename T>
  concept StreamInsertable = requires(std::ostream& os, const T& t) { os << t; };

  // Specialization for stream insert.
  template <StreamInsertable T>
  struct PrintHelper<T>
  {
    void operator()(std::ostream& os, const std::any& data) { os << std::any_cast<T>(data) << " "; }
  };

  template <StreamInsertable T>
  struct PrintHelper<std::optional<T>>
  {
    void operator()(std::ostream& os, const std::any& data)
    {
      auto val = std::any_cast<std::optional<T>>(data);
      if (val.has_value())
        PrintHelper<T>{}(os, *val);
      else
        os << "none ";
    }
  };

  // Specialization for vectors.
  template <typename T>
  struct PrintHelper<std::vector<T>>
  {
    void operator()(std::ostream& os, const std::any& data)
    {
      FOUR_C_ASSERT(typeid(std::vector<T>) == data.type(), "Implementation error.");
      const auto& vec = std::any_cast<std::vector<T>>(data);
      for (const auto& v : vec)
      {
        PrintHelper<T>{}(os, v);
      }
    }
  };

  // Specialization for maps.
  template <typename Key, typename Value>
  struct PrintHelper<std::map<Key, Value>>
  {
    void operator()(std::ostream& os, const std::any& data)
    {
      FOUR_C_ASSERT(typeid(std::map<Key, Value>) == data.type(), "Implementation error.");
      const auto& map = std::any_cast<std::map<Key, Value>>(data);
      for (const auto& [key, value] : map)
      {
        os << key << " : ";
        PrintHelper<Value>{}(os, value);
      }
    }
  };
}  // namespace Core::IO::Internal::InputParameterContainerImplementation



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
