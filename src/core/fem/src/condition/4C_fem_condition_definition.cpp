// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fem_condition_definition.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_io_input_file.hpp"
#include "4C_io_input_file_utils.hpp"
#include "4C_io_input_spec_builders.hpp"
#include "4C_utils_exceptions.hpp"

#include <algorithm>
#include <iterator>
#include <utility>

FOUR_C_NAMESPACE_OPEN



/* -----------------------------------------------------------------------------------------------*
 | Class ConditionDefinition                                                                      |
 * -----------------------------------------------------------------------------------------------*/

Core::Conditions::ConditionDefinition::ConditionDefinition(std::string sectionname,
    std::string conditionname, std::string description, Core::Conditions::ConditionType condtype,
    bool buildgeometry, Core::Conditions::GeometryType gtype)
    : sectionname_(std::move(sectionname)),
      conditionname_(std::move(conditionname)),
      description_(std::move(description)),
      condtype_(condtype),
      buildgeometry_(buildgeometry),
      gtype_(gtype)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::Conditions::ConditionDefinition::add_component(Core::IO::InputSpec&& spec)
{
  specs_.emplace_back(std::move(spec));
}


void Core::Conditions::ConditionDefinition::add_component(const Core::IO::InputSpec& spec)
{
  specs_.emplace_back(spec);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::Conditions::ConditionDefinition::read(Core::IO::InputFile& input,
    std::multimap<int, std::shared_ptr<Core::Conditions::Condition>>& cmap) const
{
  using namespace Core::IO::InputSpecBuilders;
  auto condition_spec = all_of({
      parameter<int>("E", {.description = "ID of the condition. This ID refers to the respective "
                                          "topological entity of the condition."}),
      all_of(specs_),
  });

  for (const auto& fragment : input.in_section(section_name()))
  {
    std::optional<Core::IO::InputParameterContainer> data;
    try
    {
      data = fragment.match(condition_spec);
    }
    catch (const Core::Exception& e)
    {
      FOUR_C_THROW("Failed to match condition specification in section '%s'. The error was:\n%s.",
          section_name().c_str(), e.what());
    }
    if (!data)
    {
      FOUR_C_THROW(
          "Failed to match condition specification in section '%s'.", sectionname_.c_str());
    }

    // Read a one-based condition number but convert it to zero-based for internal use.
    const int dobjid = data->get<int>("E") - 1;

    std::shared_ptr<Core::Conditions::Condition> condition =
        std::make_shared<Core::Conditions::Condition>(dobjid, condtype_, buildgeometry_, gtype_);
    condition->parameters() = *std::move(data);

    //------------------------------- put condition in map of conditions
    cmap.insert(std::pair<int, std::shared_ptr<Core::Conditions::Condition>>(dobjid, condition));
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::ostream& Core::Conditions::ConditionDefinition::print(std::ostream& stream)
{
  Core::IO::print_section_header(stream, sectionname_);

  using namespace Core::IO::InputSpecBuilders;
  auto condition_spec = all_of({
      parameter<int>("E"),
      all_of(specs_),
  });

  condition_spec.print_as_dat(stream);

  stream << "\n";
  return stream;
}

FOUR_C_NAMESPACE_CLOSE
