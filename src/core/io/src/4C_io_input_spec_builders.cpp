// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_io_input_spec_builders.hpp"

#include <set>
#include <unordered_map>
#include <utility>

FOUR_C_NAMESPACE_OPEN

namespace
{
  void parse_in_arbitrary_order(Core::IO::ValueParser& parser,
      const std::vector<Core::IO::InputSpec>& line, Core::IO::InputParameterContainer& container)
  {
    std::unordered_map<std::string, const Core::IO::InputSpec*> name_to_entry_map;
    std::set<const Core::IO::InputSpec*> unnamed_entries;
    for (const auto& entry : line)
    {
      const auto& name = entry.impl().name();
      if (name.empty())
      {
        unnamed_entries.insert(&entry);
      }
      else
      {
        name_to_entry_map[name] = &entry;
      }
    }

    const auto are_entries_left = [&]()
    { return !name_to_entry_map.empty() || !unnamed_entries.empty(); };


    // Parse as long as there are tokens and we expect more entries.
    while (!parser.at_end() && are_entries_left())
    {
      // Only peek at the next value, do not consume it yet.
      const std::string next(parser.peek());

      if (next == "//")
      {
        parser.consume_comment(next);
        break;
      }

      // The typical case: peeking reveals the key of the next value.
      const auto it = name_to_entry_map.find(next);
      if (it != name_to_entry_map.end())
      {
        // Will consume the name as well as the value.
        it->second->impl().parse(parser, container);

        // Drop the entry from the map: we do not want to parse the same value twice. This also
        // allows to check if all required values have been parsed.
        name_to_entry_map.erase(it);
        continue;
      }

      // We might find a parseable unnamed component.
      if (!unnamed_entries.empty())
      {
        bool parsed = false;
        Core::IO::ValueParser::BacktrackScope backtrack_scope(parser);
        for (const auto& entry : unnamed_entries)
        {
          try
          {
            entry->impl().parse(parser, container);
          }
          catch (const Core::Exception&)
          {
            // Try the next component.
            parser.backtrack();
          }

          // Drop the entry from the set: we do not want to parse the same value twice.
          unnamed_entries.erase(entry);
          parsed = true;
          break;
        }

        if (parsed) continue;
      }

      // If we reach this point, we could not parse the next value. What remains must therefore
      // belong to other specs.
      break;
    }

    // Consume a potential trailing comment
    {
      const auto next = parser.peek();
      if (next == "//") parser.consume_comment("//");
    }

    for (const auto& entry : unnamed_entries)
    {
      // Unnamed entries contain a useful description, which indicates what is missing.
      FOUR_C_ASSERT_ALWAYS(!entry->impl().required(), "Required '%s' not found in input line",
          entry->impl().description().c_str());
    }

    // Check if all required values have been parsed, i.e., any remaining component must be optional
    for (const auto& [name, entry] : name_to_entry_map)
    {
      if (entry->impl().required())
      {
        FOUR_C_THROW("Required value '%s' not found in input line", name.c_str());
      }
      else if (entry->impl().has_default_value())
      {
        entry->impl().set_default_value(container);
      }
    }
  }

  void assert_unique_or_empty_names(const std::vector<Core::IO::InputSpec>& specs)
  {
    std::set<std::string> names;
    for (const auto& spec : specs)
    {
      if (!spec.impl().name().empty())
        FOUR_C_ASSERT_ALWAYS(names.insert(spec.impl().name()).second,
            "Duplicate component name '%s' found in input line.", spec.impl().name().c_str());
    }
  }

  bool all_have_default_values(const std::vector<Core::IO::InputSpec>& specs)
  {
    return std::all_of(specs.begin(), specs.end(),
        [](const Core::IO::InputSpec& spec) { return spec.impl().has_default_value(); });
  }

  [[nodiscard]] const std::string& describe(const Core::IO::InputSpec& spec)
  {
    return (spec.impl().name().empty()) ? spec.impl().description() : spec.impl().name();
  }

  [[nodiscard]] std::string describe(const std::vector<Core::IO::InputSpec>& specs)
  {
    if (specs.empty()) return "{}";

    std::string description = "{";
    for (const auto& spec : specs)
    {
      // Unnamed InputSpecs are created internally and will have a description that is useful for
      // error messages.
      description += describe(spec);
      description += ", ";
    }
    description.pop_back();
    description.pop_back();
    description += "}";

    return description;
  }


}  // namespace

void Core::IO::InputSpecBuilders::Internal::GroupSpec::parse(
    Core::IO::ValueParser& parser, Core::IO::InputParameterContainer& container) const
{
  // Support the special case of unnamed groups.
  if (!name.empty()) parser.consume(name);

  // Parse into a separate container to avoid side effects if parsing fails.
  Core::IO::InputParameterContainer subcontainer;
  parse_in_arbitrary_order(parser, specs, subcontainer);

  if (name.empty())
    container.merge(subcontainer);
  else
    container.group(name) = subcontainer;
}


void Core::IO::InputSpecBuilders::Internal::GroupSpec::set_default_value(
    Core::IO::InputParameterContainer& container) const
{
  auto& subcontainer = (name.empty()) ? container : container.group(name);
  for (const auto& spec : specs)
  {
    if (spec.impl().has_default_value()) spec.impl().set_default_value(subcontainer);
  }
}


void Core::IO::InputSpecBuilders::Internal::GroupSpec::print(
    std::ostream& stream, std::size_t indent) const
{
  if (!name.empty())
  {
    stream << "// " << std::string(indent, ' ') << name << ":\n";
    indent += 2;
  }
  // Only print an all_of if it is not present at the top-level.
  else if (indent > 0)
  {
    stream << "// " << std::string(indent, ' ') << "<all_of>:\n";
    indent += 2;
  }

  for (const auto& spec : specs)
  {
    spec.impl().print(stream, indent);
  }
}

void Core::IO::InputSpecBuilders::Internal::GroupSpec::emit_metadata(ryml::NodeRef node) const
{
  node |= ryml::MAP;
  if (!name.empty())
  {
    node["name"] << name;
  }

  node["type"] << (name.empty() ? "all_of" : "group");
  if (!data.description.empty())
  {
    node["description"] << data.description;
  }
  emit_value_as_yaml(node["required"], data.required);
  node["specs"] |= ryml::SEQ;
  {
    for (const auto& spec : specs)
    {
      auto child = node["specs"].append_child();
      child |= ryml::MAP;
      spec.impl().emit_metadata(child);
    };
  }
}

void Core::IO::InputSpecBuilders::Internal::OneOfSpec::parse(
    Core::IO::ValueParser& parser, Core::IO::InputParameterContainer& container) const
{
  ValueParser::BacktrackScope backtrack_scope(parser);

  // Try to parse a component and backtrack if it fails.
  const auto try_parse = [&](const Core::IO::InputSpec& spec) -> bool
  {
    try
    {
      spec.impl().parse(parser, container);
      return true;
    }
    catch (const Core::Exception&)
    {
      parser.backtrack();
      return false;
    }
  };

  auto component = specs.begin();
  for (; component != specs.end(); ++component)
  {
    const auto success = try_parse(*component);

    // Now check that no other component can be parsed.
    if (success)
    {
      ValueParser::BacktrackScope backtrack_scope_after_success(parser);
      for (auto other = specs.begin(); other != specs.end(); ++other)
      {
        if (other == component) continue;
        if (try_parse(*other))
        {
          // Backtrack to the original position.
          parser.backtrack();
          FOUR_C_THROW(
              "Ambiguous input: both '%s' and '%s' could be parsed, but only one of them is "
              "expected.",
              describe(*component).c_str(), describe(*other).c_str());
        }
      }

      // If we reach this point, we have successfully parsed the input.
      std::size_t index = std::distance(specs.begin(), component);
      if (on_parse_callback) on_parse_callback(parser, container, index);

      return;
    }
  }

  // Convert remainder into a null-terminated string for the error message.
  std::string remainder(parser.get_unparsed_remainder());
  FOUR_C_THROW("While parsing '%s'.\nNone of the specs fit the input. Expected %s",
      remainder.c_str(), data.description.c_str());
}


void Core::IO::InputSpecBuilders::Internal::OneOfSpec::set_default_value(
    Core::IO::InputParameterContainer& container) const
{
  FOUR_C_THROW("Implementation error: OneOfSpec cannot have a default value.");
}


void Core::IO::InputSpecBuilders::Internal::OneOfSpec::print(
    std::ostream& stream, std::size_t indent) const
{
  stream << "// " << std::string(indent, ' ') << "<one_of>:\n";
  for (const auto& spec : specs)
  {
    spec.impl().print(stream, indent + 2);
  }
}

void Core::IO::InputSpecBuilders::Internal::OneOfSpec::emit_metadata(ryml::NodeRef node) const
{
  node |= ryml::MAP;

  node["type"] << "one_of";
  if (!data.description.empty())
  {
    node["description"] << data.description;
  }
  emit_value_as_yaml(node["required"], data.required);
  node["specs"] |= ryml::SEQ;
  for (const auto& spec : specs)
  {
    auto child = node["specs"].append_child();
    spec.impl().emit_metadata(child);
  };
}


Core::IO::InputSpec Core::IO::InputSpecBuilders::tag(std::string name, ScalarData<bool> data)
{
  Internal::sanitize_required_default(data);
  return user_defined<bool>(
      name, data,
      // If we encounter the tag, we set the value to true.
      [name](ValueParser& parser, InputParameterContainer& container)
      {
        parser.consume(name);
        container.add(name, true);
      },
      // Print the tag as a comment.
      [name, data](std::ostream& out, std::size_t indent)
      {
        out << "// " << std::string(indent, ' ') << name;
        out << " <tag>";
        if (!data.required.value()) out << " (optional)";
        if (!data.description.empty()) out << " " << std::quoted(data.description);
        out << "\n";
      });
}

namespace
{
  // An all_of is only used for grouping purposes and does not have a name. If it is used inside
  // another group, it can be flattened by stealing its specs and adding them to the parent group.
  std::vector<Core::IO::InputSpec> flatten_all_ofs(std::vector<Core::IO::InputSpec> specs)
  {
    std::vector<Core::IO::InputSpec> flattened_specs;
    for (auto&& spec : specs)
    {
      if (auto* group = dynamic_cast<Core::IO::Internal::InputSpecTypeErasedImplementation<
              Core::IO::InputSpecBuilders::Internal::GroupSpec>*>(&spec.impl());
          // must be an all_of to steal the specs
          group && group->wrapped.name.empty())
      {
        for (auto&& sub_spec : group->wrapped.specs)
        {
          flattened_specs.emplace_back(std::move(sub_spec));
        }
      }
      else
      {
        flattened_specs.emplace_back(std::move(spec));
      }
    }

    return flattened_specs;
  }
}  // namespace


Core::IO::InputSpec Core::IO::InputSpecBuilders::group(
    std::string name, std::vector<InputSpec> specs, Core::IO::InputSpecBuilders::GroupData data)
{
  auto flattened_specs = flatten_all_ofs(std::move(specs));

  assert_unique_or_empty_names(flattened_specs);

  IO::Internal::InputSpecTypeErasedBase::CommonData common_data{
      .name = name,
      .description = data.description,
      .required = data.required,
      .has_default_value = all_have_default_values(flattened_specs),
  };

  return IO::Internal::make_spec(
      Internal::GroupSpec{
          .name = name, .data = std::move(data), .specs = std::move(flattened_specs)},
      common_data);
}


Core::IO::InputSpec Core::IO::InputSpecBuilders::all_of(std::vector<InputSpec> specs)
{
  auto flattened_specs = flatten_all_ofs(std::move(specs));

  assert_unique_or_empty_names(flattened_specs);

  const bool any_required =
      std::ranges::any_of(flattened_specs, [](const auto& spec) { return spec.impl().required(); });

  // Generate a description of the form "group {a, b, c}".
  std::string description = "all_of " + describe(flattened_specs);

  IO::Internal::InputSpecTypeErasedBase::CommonData common_data{
      .name = "",
      .description = description,
      .required = any_required,
      .has_default_value = all_have_default_values(flattened_specs),
  };

  return IO::Internal::make_spec(Internal::GroupSpec{.name = "",
                                     .data = {.description = description, .required = any_required},
                                     .specs = std::move(flattened_specs)},
      common_data);
}


Core::IO::InputSpec Core::IO::InputSpecBuilders::one_of(std::vector<InputSpec> specs,
    std::function<void(ValueParser& parser, InputParameterContainer& container, std::size_t index)>
        on_parse_callback)
{
  FOUR_C_ASSERT_ALWAYS(!specs.empty(), "`one_of` must contain at least one entry.");

  assert_unique_or_empty_names(specs);

  // Assert that all specs are required.
  const bool all_required =
      std::ranges::all_of(specs, [](const auto& spec) { return spec.impl().required(); });
  if (!all_required)
  {
    std::string non_required;
    for (const auto& spec : specs)
    {
      if (!spec.impl().required()) non_required += "    " + describe(spec) + "\n";
    }

    FOUR_C_THROW(
        "All specs in a 'one_of' must be required to avoid confusion. The following InputSpecs "
        "are not required:\n%s",
        non_required.c_str());
  }

  std::string description = "one_of " + describe(specs);
  IO::Internal::InputSpecTypeErasedBase::CommonData common_data{
      .name = "",
      .description = description,
      .required = true,
      .has_default_value = false,
  };

  GroupData group_data{
      .description = description,
      .required = true,
  };

  return IO::Internal::make_spec(Internal::OneOfSpec{.data = std::move(group_data),
                                     .specs = std::move(specs),
                                     .on_parse_callback = std::move(on_parse_callback)},
      std::move(common_data));
}


FOUR_C_NAMESPACE_CLOSE