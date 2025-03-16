// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_beaminteraction_beam_to_solid_visualization_output_writer_base.hpp"

#include "4C_beaminteraction_beam_to_solid_visualization_output_writer_visualization.hpp"
#include "4C_utils_exceptions.hpp"

#include <utility>

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
BeamInteraction::BeamToSolidVisualizationOutputWriterBase::BeamToSolidVisualizationOutputWriterBase(
    const std::string& base_output_name, Core::IO::VisualizationParameters visualization_params)
    : base_output_name_(base_output_name), visualization_params_(std::move(visualization_params))
{
}

/**
 *
 */
std::shared_ptr<BeamInteraction::BeamToSolidOutputWriterVisualization>
BeamInteraction::BeamToSolidVisualizationOutputWriterBase::add_visualization_writer(
    const std::string& writer_name, const std::string& writer_name_key)
{
  const auto& it = visualization_writers_.find(writer_name_key);
  if (it != visualization_writers_.end())
  {
    FOUR_C_THROW(
        "The output writer key '{}' you want to add already exists.", writer_name_key.c_str());
  }
  else
  {
    auto new_writer = std::make_shared<BeamInteraction::BeamToSolidOutputWriterVisualization>(

        base_output_name_ + "-" + writer_name, visualization_params_);
    visualization_writers_[writer_name_key] = new_writer;
    return new_writer;
  }
}

/**
 *
 */
std::shared_ptr<BeamInteraction::BeamToSolidOutputWriterVisualization>
BeamInteraction::BeamToSolidVisualizationOutputWriterBase::add_visualization_writer(
    const std::string& writer_name)
{
  return add_visualization_writer(writer_name, writer_name);
}

/**
 *
 */
std::shared_ptr<BeamInteraction::BeamToSolidOutputWriterVisualization>
BeamInteraction::BeamToSolidVisualizationOutputWriterBase::get_visualization_writer(
    const std::string& writer_name)
{
  const auto& it = visualization_writers_.find(writer_name);
  if (it != visualization_writers_.end())
    return it->second;
  else
    return nullptr;
}

/**
 *
 */
void BeamInteraction::BeamToSolidVisualizationOutputWriterBase::write(
    const unsigned int timestep_number, const double time)
{
  for (auto& it : visualization_writers_) it.second->write(timestep_number, time);
}

FOUR_C_NAMESPACE_CLOSE
