// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_IO_YAML_EMITTER_HPP
#define FOUR_C_IO_YAML_EMITTER_HPP

#include "4C_config.hpp"

#include <ryml.hpp>
#include <ryml_std.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::IO
{
  /**
   * This class wraps the yaml implementation in our own namespace. Most code does not need the
   * implementation details of the yaml library, so we hide it here and forward declare the class
   * in relevant headers.
   */
  class YamlEmitter
  {
   public:
    /**
     * Construct a new YamlEmitter object as a wrapper around the given node.
     */
    YamlEmitter(ryml::NodeRef node) : node(node) {}

    /**
     * The node to which the emitting code is supposed to write.
     */
    ryml::NodeRef node;
  };

  template <typename T>
  concept YamlSupportedType = requires(ryml::NodeRef node, const T& value) {
    { node << value };
  };

  template <YamlSupportedType T>
  void emit_value_as_yaml(ryml::NodeRef node, const T& value)
  {
    node << value;
  }

  inline void emit_value_as_yaml(ryml::NodeRef node, const std::string& value)
  {
    node << ryml::to_csubstr(value);
  }

  inline void emit_value_as_yaml(ryml::NodeRef node, const bool& value)
  {
    node << (value ? "true" : "false");
  }


  template <typename T, typename U>
  void emit_value_as_yaml(ryml::NodeRef node, const std::pair<T, U>& value)
  {
    node |= ryml::SEQ | ryml::FLOW_SL;
    auto first = node.append_child();
    emit_value_as_yaml(first, value.first);
    auto second = node.append_child();
    emit_value_as_yaml(second, value.second);
  }

  template <typename T>
  void emit_value_as_yaml(ryml::NodeRef node, const std::vector<T>& value)
  {
    node |= ryml::SEQ | ryml::FLOW_SL;
    for (const auto& v : value)
    {
      auto child = node.append_child();
      emit_value_as_yaml(child, v);
    }
  }
}  // namespace Core::IO

FOUR_C_NAMESPACE_CLOSE

#endif
