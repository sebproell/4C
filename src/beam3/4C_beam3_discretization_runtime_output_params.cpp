// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_beam3_discretization_runtime_output_params.hpp"

#include "4C_utils_exceptions.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
Discret::Elements::BeamRuntimeOutputParams::BeamRuntimeOutputParams()
    : isinit_(false),
      issetup_(false),
      output_displacement_state_(false),
      use_absolute_positions_visualizationpoint_coordinates_(true),
      write_internal_energy_element_(false),
      write_kinetic_energy_element_(false),
      write_triads_visualizationpoints_(false),
      write_material_crosssection_strains_gausspoints_(false),
      write_material_crosssection_strains_continuous_(false),
      write_material_crosssection_stresses_gausspoints_(false),
      write_spatial_crosssection_stresses_gausspoints_(false),
      write_filament_condition_(false),
      write_orientation_parameter_(false),
      write_rve_crosssection_forces_(false),
      write_ref_length_(false),
      write_element_gid_(false),
      write_element_ghosting_(false),
      n_subsegments_(0)
{
  // empty constructor
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void Discret::Elements::BeamRuntimeOutputParams::init(
    const Teuchos::ParameterList& IO_vtk_structure_beams_paramslist)
{
  // We have to call setup() after init()
  issetup_ = false;

  // initialize the parameter values

  output_displacement_state_ = IO_vtk_structure_beams_paramslist.get<bool>("DISPLACEMENT");

  use_absolute_positions_visualizationpoint_coordinates_ =
      IO_vtk_structure_beams_paramslist.get<bool>("USE_ABSOLUTE_POSITIONS");

  write_internal_energy_element_ =
      IO_vtk_structure_beams_paramslist.get<bool>("INTERNAL_ENERGY_ELEMENT");

  write_kinetic_energy_element_ =
      IO_vtk_structure_beams_paramslist.get<bool>("KINETIC_ENERGY_ELEMENT");

  write_triads_visualizationpoints_ =
      IO_vtk_structure_beams_paramslist.get<bool>("TRIAD_VISUALIZATIONPOINT");

  write_material_crosssection_strains_gausspoints_ =
      IO_vtk_structure_beams_paramslist.get<bool>("STRAINS_GAUSSPOINT");

  write_material_crosssection_strains_continuous_ =
      IO_vtk_structure_beams_paramslist.get<bool>("STRAINS_CONTINUOUS");

  write_material_crosssection_stresses_gausspoints_ =
      IO_vtk_structure_beams_paramslist.get<bool>("MATERIAL_FORCES_GAUSSPOINT");

  write_material_crosssection_strains_continuous_ =
      IO_vtk_structure_beams_paramslist.get<bool>("MATERIAL_FORCES_CONTINUOUS");

  write_spatial_crosssection_stresses_gausspoints_ =
      IO_vtk_structure_beams_paramslist.get<bool>("SPATIAL_FORCES_GAUSSPOINT");

  write_orientation_parameter_ =
      IO_vtk_structure_beams_paramslist.get<bool>("ORIENTATION_PARAMETER");

  write_rve_crosssection_forces_ =
      IO_vtk_structure_beams_paramslist.get<bool>("RVE_CROSSSECTION_FORCES");

  write_ref_length_ = IO_vtk_structure_beams_paramslist.get<bool>("REF_LENGTH");

  write_element_gid_ = IO_vtk_structure_beams_paramslist.get<bool>("ELEMENT_GID");

  write_element_ghosting_ = IO_vtk_structure_beams_paramslist.get<bool>("ELEMENT_GHOSTING");

  n_subsegments_ = IO_vtk_structure_beams_paramslist.get<int>("NUMBER_SUBSEGMENTS");
  if (n_subsegments_ < 1)
    FOUR_C_THROW("The number of subsegments has to be at least 1. Got {}", n_subsegments_);

  isinit_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void Discret::Elements::BeamRuntimeOutputParams::setup()
{
  if (not is_init()) FOUR_C_THROW("init() has not been called, yet!");

  // Nothing to do here at the moment

  issetup_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void Discret::Elements::BeamRuntimeOutputParams::check_init_setup() const
{
  if (not is_init() or not is_setup()) FOUR_C_THROW("Call init() and setup() first!");
}

FOUR_C_NAMESPACE_CLOSE
