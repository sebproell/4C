// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_IO_YAML_HPP
#define FOUR_C_IO_YAML_HPP

#include "4C_config.hpp"

#include "4C_utils_exceptions.hpp"

#include <magic_enum/magic_enum.hpp>
#include <ryml.hpp>      // IWYU pragma: export
#include <ryml_std.hpp>  // IWYU pragma: export

#include <filesystem>
#include <format>
#include <optional>

FOUR_C_NAMESPACE_OPEN

namespace Core::IO
{
  namespace Internal
  {
    /**
     * This class wraps the yaml implementation in our own namespace. Most code does not need the
     * implementation details of the yaml library, so we hide it here and forward declare the class
     * in relevant headers. This has the additional benefit that we can attach more data that is
     * useful when traversing yaml trees, e.g. an associated file path. Since ryml treats const-ness
     * with different node types, we model this behavior with a template parameter selecting between
     * a const and a non-const node.
     *
     * @note This class uses CRTP to allow wrap() to return the concrete derived type.
     */
    template <bool is_const, typename Derived>
    class YamlNodeRefImpl
    {
     public:
      using NodeRef = std::conditional_t<is_const, ryml::ConstNodeRef, ryml::NodeRef>;

     private:
      YamlNodeRefImpl(NodeRef node, std::filesystem::path associated_file)
          : node(node), associated_file(std::move(associated_file))
      {
      }


     public:
      /**
       * When traversing the yaml tree, we frequently need to process child nodes of the wrapped
       * node. This function wraps the @p next node while keeping any other information stored in
       * this wrapper object.
       */
      [[nodiscard]] Derived wrap(NodeRef next) const { return Derived(next, associated_file); }

      /**
       * The wrapped ryml node.
       */
      NodeRef node;

      /**
       * The file associated with this node. May be used for error messages and to correctly resolve
       * relative paths.
       */
      std::filesystem::path associated_file;

      /**
       * Make the derived class a friend to allow it to call the constructor.
       */
      friend Derived;
    };
  }  // namespace Internal

  /**
   * Refers to a non-const yaml node.
   */
  class YamlNodeRef : public Internal::YamlNodeRefImpl<false, YamlNodeRef>
  {
   public:
    YamlNodeRef(NodeRef node, std::filesystem::path associated_file)
        : Internal::YamlNodeRefImpl<false, YamlNodeRef>(node, std::move(associated_file))
    {
    }
  };

  /**
   * Refers to a const yaml node.
   */
  class ConstYamlNodeRef : public Internal::YamlNodeRefImpl<true, ConstYamlNodeRef>
  {
   public:
    ConstYamlNodeRef(ryml::ConstNodeRef node, std::filesystem::path associated_file)
        : Internal::YamlNodeRefImpl<true, ConstYamlNodeRef>(node, std::move(associated_file))
    {
    }
  };


  /**
   * The types that we support for yaml (de)serialization.
   */
  template <typename T>
  concept YamlSupportedType =
      std::same_as<T, int> || std::same_as<T, double> || std::same_as<T, bool> ||
      std::same_as<T, std::string> || std::same_as<T, std::filesystem::path> || std::is_enum_v<T>;

  /**
   * An exception thrown when an error occurs during yaml processing..
   */
  class YamlException : public Core::Exception
  {
   public:
    explicit YamlException(const std::string& message);
  };


  /**
   * Initialize a ryml tree which throws YamlExceptions on parse errors.
   */
  [[nodiscard]] ryml::Tree init_yaml_tree_with_exceptions();

  void emit_value_as_yaml(ryml::NodeRef node, const int& value);

  void emit_value_as_yaml(ryml::NodeRef node, const double& value);

  void emit_value_as_yaml(ryml::NodeRef node, const std::string& value);
  void emit_value_as_yaml(ryml::NodeRef node, const std::string_view& value);

  void emit_value_as_yaml(ryml::NodeRef node, const bool& value);

  void emit_value_as_yaml(ryml::NodeRef node, const std::filesystem::path& value);

  template <typename T>
    requires(std::is_enum_v<T>)
  void emit_value_as_yaml(ryml::NodeRef node, const T& value);

  template <typename T>
  void emit_value_as_yaml(ryml::NodeRef node, const std::optional<T>& value);

  template <typename T>
  void emit_value_as_yaml(ryml::NodeRef node, const std::map<std::string, T>& value);

  template <typename T>
  void emit_value_as_yaml(ryml::NodeRef node, const std::vector<T>& value);

  template <YamlSupportedType T>
    requires(!std::is_enum_v<T>)
  void read_value_from_yaml(ConstYamlNodeRef node, T& value)
  {
    FOUR_C_ASSERT_ALWAYS(node.node.has_val(), "Expected a value node.");
    node.node >> value;
  }

  void read_value_from_yaml(ConstYamlNodeRef node, double& value);

  /**
   * Reading a bool currently supports the following options:
   *
   * - truthy values: "true", "yes", "on", "1"
   * - falsy values: "false", "no", "off", "0"
   *
   * All of them are case-insensitive.
   */
  void read_value_from_yaml(ConstYamlNodeRef node, bool& value);

  template <typename T>
    requires(std::is_enum_v<T>)
  void read_value_from_yaml(ConstYamlNodeRef node, T& value)
  {
    FOUR_C_ASSERT_ALWAYS(node.node.has_val(), "Expected a value node.");
    auto substr = node.node.val();
    auto val = magic_enum::enum_cast<T>(std::string_view(substr.data(), substr.size()));
    if (val)
    {
      value = *val;
    }
    else
    {
      FOUR_C_THROW("Could not parse value '{}' as an enum constant of type '{}'.",
          std::string_view(substr.data(), substr.size()), magic_enum::enum_type_name<T>());
    }
  }

  template <typename T>
  void read_value_from_yaml(ConstYamlNodeRef node, std::optional<T>& value)
  {
    FOUR_C_ASSERT_ALWAYS(node.node.has_val(), "Expected a value node.");
    if (node.node.val_is_null())
    {
      value.reset();
    }
    else
    {
      value = T{};
      read_value_from_yaml(node, value.value());
    }
  }

  inline void read_value_from_yaml(ConstYamlNodeRef node, std::filesystem::path& value)
  {
    FOUR_C_ASSERT_ALWAYS(node.node.has_val(), "Expected a value node.");
    std::string path;
    read_value_from_yaml(node, path);

    value = std::filesystem::path(path);
    if (value.is_relative())
    {
      value = node.associated_file.parent_path() / value;
    }
  }

  template <typename U>
  void read_value_from_yaml(ConstYamlNodeRef node, std::map<std::string, U>& value)
  {
    FOUR_C_ASSERT_ALWAYS(node.node.is_map(), "Expected a map node for a map.");
    value.clear();
    for (auto child : node.node.children())
    {
      FOUR_C_ASSERT_ALWAYS(child.has_key(), "Expected a key in the map node.");
      std::string key(child.key().data(), child.key().size());
      read_value_from_yaml(node.wrap(child), value[key]);
    }
  }

  template <typename T>
  void read_value_from_yaml(ConstYamlNodeRef node, std::vector<T>& value)
  {
    FOUR_C_ASSERT_ALWAYS(node.node.is_seq(), "Expected a sequence node for a vector.");
    value.clear();
    for (auto child : node.node.children())
    {
      T v;
      read_value_from_yaml(node.wrap(child), v);
      value.push_back(v);
    }
  }
}  // namespace Core::IO

template <typename T>
  requires(std::is_enum_v<T>)
void Core::IO::emit_value_as_yaml(ryml::NodeRef node, const T& value)
{
  std::string_view str = magic_enum::enum_name(value);
  node << ryml::csubstr(str.data(), str.size());
}

template <typename T>
void Core::IO::emit_value_as_yaml(ryml::NodeRef node, const std::optional<T>& value)
{
  if (value.has_value())
  {
    emit_value_as_yaml(node, value.value());
  }
  else
  {
    node = "null";
  }
}

template <typename T>
void Core::IO::emit_value_as_yaml(ryml::NodeRef node, const std::map<std::string, T>& value)
{
  node |= ryml::MAP;
  for (const auto& [key, v] : value)
  {
    auto child = node.append_child();
    child << ryml::key(key);
    emit_value_as_yaml(child, v);
  }
}

template <typename T>
void Core::IO::emit_value_as_yaml(ryml::NodeRef node, const std::vector<T>& value)
{
  node |= ryml::SEQ | ryml::FLOW_SL;
  for (const auto& v : value)
  {
    auto child = node.append_child();
    emit_value_as_yaml(child, v);
  }
}

FOUR_C_NAMESPACE_CLOSE

#endif
