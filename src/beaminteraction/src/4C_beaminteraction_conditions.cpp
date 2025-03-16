// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_beaminteraction_conditions.hpp"

#include "4C_beaminteraction_beam_to_beam_contact_condition.hpp"
#include "4C_beaminteraction_beam_to_beam_point_coupling_pair.hpp"
#include "4C_beaminteraction_beam_to_beam_point_coupling_pair_condition.hpp"
#include "4C_beaminteraction_beam_to_solid_conditions.hpp"
#include "4C_beaminteraction_beam_to_solid_surface_contact_params.hpp"
#include "4C_beaminteraction_beam_to_solid_surface_meshtying_params.hpp"
#include "4C_beaminteraction_beam_to_solid_volume_meshtying_params.hpp"
#include "4C_beaminteraction_contact_pair.hpp"
#include "4C_beaminteraction_contact_params.hpp"
#include "4C_fem_condition.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_geometry_pair_element.hpp"
#include "4C_geometry_pair_evaluation_data_base.hpp"
#include "4C_inpar_beam_to_solid.hpp"

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
BeamInteraction::BeamInteractionConditionBase::BeamInteractionConditionBase(
    const std::shared_ptr<const Core::Conditions::Condition>& condition_line)
    : condition_line_(condition_line), line_ids_()
{
}

/**
 *
 */
void BeamInteraction::BeamInteractionConditionBase::build_id_sets(
    const std::shared_ptr<const Core::FE::Discretization>& discretization)
{
  // Set the IDs of the line elements.
  std::vector<int> line_ids;
  condition_to_element_ids(*condition_line_, line_ids);
  line_ids_ = std::set<int>(line_ids.begin(), line_ids.end());
}

/**
 *
 */
void BeamInteraction::BeamInteractionConditionBase::setup(
    const std::shared_ptr<const Core::FE::Discretization>& discret)
{
}

/**
 *
 */
void BeamInteraction::BeamInteractionConditionBase::clear() {}

/**
 *
 */
BeamInteraction::BeamInteractionConditions::BeamInteractionConditions() {}

/**
 *
 */
void BeamInteraction::BeamInteractionConditions::set_beam_interaction_conditions(
    const Core::FE::Discretization& discret, const BeamContactParams& params_ptr)
{
  condition_map_.clear();

  // Get all available interaction types.
  std::vector<Inpar::BeamInteraction::BeamInteractionConditions> interaction_types;
  Inpar::BeamInteraction::beam_interaction_conditions_get_all(interaction_types);

  // Loop over interaction types.
  for (const auto& interaction_type : interaction_types)
  {
    if (interaction_type == Inpar::BeamInteraction::BeamInteractionConditions::beam_to_beam_contact)
    {
      // Add all beam-to-beam conditions.
      std::vector<std::shared_ptr<BeamInteractionConditionBase>>& interaction_vector =
          condition_map_[interaction_type];

      // Get the names for the conditions of this type.
      std::string condition_name = "BeamToBeamContact";

      // Get the line conditions from the discretization.
      std::vector<std::shared_ptr<Core::Conditions::Condition>> condition_lines;
      discret.get_condition(condition_name, condition_lines);

      // Match the coupling IDs from the input line.
      std::map<int, std::pair<std::shared_ptr<const Core::Conditions::Condition>,
                        std::shared_ptr<const Core::Conditions::Condition>>>
          coupling_id_map;
      for (const auto& condition : condition_lines)
      {
        const int coupling_id = condition->parameters().get<int>("COUPLING_ID");
        auto& condition_1 = coupling_id_map[coupling_id].first;
        auto& condition_2 = coupling_id_map[coupling_id].second;
        if (condition_1 == nullptr)
          condition_1 = condition;
        else if (condition_2 == nullptr)
          condition_2 = condition;
        else
          FOUR_C_THROW("There can not be three different beam-to-beam coupling conditions.");
      }

      for (const auto& map_item : coupling_id_map)
      {
        const auto& condition_1 = map_item.second.first;
        const auto& condition_2 = map_item.second.second;
        if (condition_1 != nullptr && condition_2 != nullptr)
        {
          // We found the matching conditions, now create the beam-to-beam condition objects.
          interaction_vector.push_back(
              std::make_shared<BeamInteraction::BeamToBeamContactCondition>(
                  condition_1, condition_2));
        }
        else
          FOUR_C_THROW("Could not find both conditions ({}) for the COUPLING_ID {}",
              condition_name.c_str(), map_item.first);
      }

      // Check that all conditions were added, i.e. that there are no double definitions of
      // COUPLING_ID.
      if (2 * interaction_vector.size() != condition_lines.size())
        FOUR_C_THROW("There are multiple definitions of the same COUPLING_ID for {}",
            condition_name.c_str());
    }
    else if (interaction_type == Inpar::BeamInteraction::BeamInteractionConditions::
                                     beam_to_solid_volume_meshtying or
             interaction_type == Inpar::BeamInteraction::BeamInteractionConditions::
                                     beam_to_solid_surface_meshtying or
             interaction_type ==
                 Inpar::BeamInteraction::BeamInteractionConditions::beam_to_solid_surface_contact)
    {
      // Add all beam-to-solid conditions.
      std::vector<std::shared_ptr<BeamInteractionConditionBase>>& interaction_vector =
          condition_map_[interaction_type];

      // Get the names for the conditions of this type.
      std::array<std::string, 2> condition_names;
      Inpar::BeamToSolid::beam_to_solid_interaction_get_string(interaction_type, condition_names);

      // Get the conditions from the discretization.
      std::vector<std::shared_ptr<Core::Conditions::Condition>> condition_line;
      std::vector<std::shared_ptr<Core::Conditions::Condition>> condition_other;
      discret.get_condition(condition_names[0], condition_line);
      discret.get_condition(condition_names[1], condition_other);

      // There has to be an equal number of sections for lines and surfaces / volumes.
      if (condition_line.size() != condition_other.size())
        FOUR_C_THROW("There are {} {} sections and {} {} sections. The numbers have to match!",
            condition_line.size(), condition_names[0].c_str(), condition_other.size(),
            condition_names[1].c_str());

      // Match the coupling IDs from the input line.
      std::map<int, std::pair<std::shared_ptr<const Core::Conditions::Condition>,
                        std::shared_ptr<const Core::Conditions::Condition>>>
          coupling_id_map;
      for (const auto& condition : condition_line)
        coupling_id_map[condition->parameters().get<int>("COUPLING_ID")].first = condition;
      for (const auto& condition : condition_other)
        coupling_id_map[condition->parameters().get<int>("COUPLING_ID")].second = condition;
      for (const auto& map_item : coupling_id_map)
      {
        if (map_item.second.first != nullptr && map_item.second.second != nullptr)
        {
          // We found the matching conditions, now create the beam-to-solid condition objects.
          std::shared_ptr<BeamInteractionConditionBase> new_condition;
          if (interaction_type ==
              Inpar::BeamInteraction::BeamInteractionConditions::beam_to_solid_volume_meshtying)
            new_condition = std::make_shared<BeamInteraction::BeamToSolidConditionVolumeMeshtying>(
                map_item.second.first, map_item.second.second,
                params_ptr.beam_to_solid_volume_meshtying_params());
          else if (interaction_type == Inpar::BeamInteraction::BeamInteractionConditions::
                                           beam_to_solid_surface_meshtying)
            new_condition = std::make_shared<BeamInteraction::BeamToSolidConditionSurface>(
                map_item.second.first, map_item.second.second,
                params_ptr.beam_to_solid_surface_meshtying_params(), true);
          else if (interaction_type ==
                   Inpar::BeamInteraction::BeamInteractionConditions::beam_to_solid_surface_contact)
            new_condition = std::make_shared<BeamInteraction::BeamToSolidConditionSurface>(
                map_item.second.first, map_item.second.second,
                params_ptr.beam_to_solid_surface_contact_params(), false);
          else
            FOUR_C_THROW("Got unexpected interaction type.");
          interaction_vector.push_back(new_condition);
        }
        else
          FOUR_C_THROW("Could not find both conditions ({}, {}) for the COUPLING_ID {}",
              condition_names[0].c_str(), condition_names[1].c_str(), map_item.first);
      }

      // Check that all conditions were added, i.e. that there are no double definitions of
      // COUPLING_ID.
      if (interaction_vector.size() != condition_line.size())
        FOUR_C_THROW("There are multiple definitions of the same COUPLING_ID for {} and {}",
            condition_names[0].c_str(), condition_names[1].c_str());
    }
    else if (interaction_type ==
             Inpar::BeamInteraction::BeamInteractionConditions::beam_to_beam_point_coupling)
    {
      std::vector<std::shared_ptr<BeamInteractionConditionBase>>& interaction_vector =
          condition_map_[interaction_type];

      // Get the conditions from the discretization.
      std::vector<std::shared_ptr<Core::Conditions::Condition>> coupling_conditions;
      discret.get_condition("PenaltyPointCouplingCondition", coupling_conditions);
      for (const auto& condition : coupling_conditions)
      {
        // We found the matching conditions, now create the beam-to-beam coupling condition object
        std::shared_ptr<BeamInteractionConditionBase> new_condition;

        new_condition = std::make_shared<BeamInteraction::BeamToBeamPointCouplingCondition>(
            condition, condition->parameters().get<double>("POSITIONAL_PENALTY_PARAMETER"),
            condition->parameters().get<double>("ROTATIONAL_PENALTY_PARAMETER"));

        interaction_vector.push_back(new_condition);
      }
    }
    else
      FOUR_C_THROW("Got unexpected interaction type.");
  }
}

/**
 *
 */
void BeamInteraction::BeamInteractionConditions::build_id_sets(
    std::shared_ptr<Core::FE::Discretization> discretization)
{
  for (auto const& map_pair : condition_map_)
    for (auto const& condition : map_pair.second) condition->build_id_sets(discretization);
}

/**
 *
 */
void BeamInteraction::BeamInteractionConditions::set_state(
    const std::shared_ptr<const Core::FE::Discretization>& discret,
    const std::shared_ptr<const Solid::ModelEvaluator::BeamInteractionDataState>&
        beaminteraction_data_state)
{
  for (auto const& map_pair : condition_map_)
    for (auto const& condition : map_pair.second)
      condition->set_state(discret, beaminteraction_data_state);
}

/**
 *
 */
void BeamInteraction::BeamInteractionConditions::setup(
    const std::shared_ptr<const Core::FE::Discretization>& discret)
{
  for (auto const& map_pair : condition_map_)
    for (auto const& condition : map_pair.second) condition->setup(discret);
}

/**
 *
 */
void BeamInteraction::BeamInteractionConditions::clear()
{
  for (auto const& map_pair : condition_map_)
    for (auto const& condition : map_pair.second) condition->clear();
}

/**
 *
 */
std::shared_ptr<BeamInteraction::BeamContactPair>
BeamInteraction::BeamInteractionConditions::create_contact_pair(
    const std::vector<Core::Elements::Element const*>& ele_ptrs)
{
  std::shared_ptr<BeamInteraction::BeamContactPair> new_pair;
  for (auto& map_pair : condition_map_)
  {
    for (auto& condition : map_pair.second)
    {
      new_pair = condition->create_contact_pair(ele_ptrs);
      if (new_pair != nullptr) return new_pair;
    }
  }

  // Default return value, i.e. the pair was not found in any of the conditions.
  return nullptr;
}

/**
 *
 */
void BeamInteraction::BeamInteractionConditions::create_indirect_assembly_managers(
    const std::shared_ptr<const Core::FE::Discretization>& discret,
    std::vector<std::shared_ptr<SUBMODELEVALUATOR::BeamContactAssemblyManager>>& assembly_managers)
{
  std::shared_ptr<BeamInteraction::SUBMODELEVALUATOR::BeamContactAssemblyManager>
      condition_assembly_manager = nullptr;
  for (auto& map_pair : condition_map_)
  {
    for (auto& condition : map_pair.second)
    {
      condition_assembly_manager = condition->create_indirect_assembly_manager(discret);
      if (not(condition_assembly_manager == nullptr))
        assembly_managers.push_back(condition_assembly_manager);
    }
  }
}

/**
 *
 */
void BeamInteraction::condition_to_element_ids(
    const Core::Conditions::Condition& condition, std::vector<int>& element_ids)
{
  // Loop over the elements in the condition and get the "real" element by comparing the node IDs.
  element_ids.clear();
  element_ids.reserve(condition.geometry().size());
  for (const auto& item : condition.geometry())
  {
    int n_nodes = item.second->num_node();

    // Create the node sets and store the node IDs from the condition element in it.
    std::set<int> nodes_condition;
    std::set<int> nodes_element;
    for (int i = 0; i < n_nodes; i++) nodes_condition.insert(item.second->nodes()[i]->id());

    // Loop over all elements connected to a node and check if the nodal IDs are the same. Use the
    // last node, since if there are nodes connected to fewer elements, those are usually at the
    // end of the list.
    int local_node_id = n_nodes - 1;
    Core::Elements::Element** elements = item.second->nodes()[local_node_id]->elements();
    for (int i_element = 0; i_element < item.second->nodes()[local_node_id]->num_element();
        i_element++)
    {
      if (elements[i_element]->num_node() != n_nodes) continue;

      // Fill up the node ID map.
      nodes_element.clear();
      for (int i_nodes = 0; i_nodes < n_nodes; i_nodes++)
        nodes_element.insert(elements[i_element]->nodes()[i_nodes]->id());

      // Check if the maps are equal.
      if (std::equal(nodes_condition.begin(), nodes_condition.end(), nodes_element.begin()))
        element_ids.push_back(elements[i_element]->id());
    }
  }

  // Check if all elements were found.
  if (condition.geometry().size() != element_ids.size())
    FOUR_C_THROW("Could not find the IDs of all elements!");
}

FOUR_C_NAMESPACE_CLOSE
