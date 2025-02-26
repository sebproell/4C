// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_global_data.hpp"
#include "4C_mat_maxwell_0d_acinus.hpp"
#include "4C_red_airways_elementbase.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
| Read in the RED_AIRWAY elements                                       |
*-----------------------------------------------------------------------*/
bool Discret::Elements::RedAirway::read_element(const std::string& eletype,
    const std::string& distype, const Core::IO::InputParameterContainer& container)
{
  const int ndim = Global::Problem::instance()->n_dim();
  if (ndim != 3)
    FOUR_C_THROW("Problem defined as %dd, but found Reduced dimensional AIRWAY element.", ndim);

  // Read number of material model
  int material_id = container.get<int>("MAT");
  set_material(0, Mat::factory(material_id));

  // Read the element type, the element specific variables and store them to airwayParams_
  elem_type_ = container.get<std::string>("TYPE");
  resistance_ = container.get<std::string>("Resistance");
  elemsolving_type_ = container.get<std::string>("ElemSolvingType");

  double velPow = container.get<double>("PowerOfVelocityProfile");
  double Ew = container.get<double>("WallElasticity");
  double nu = container.get<double>("PoissonsRatio");
  double Ts = container.get<double>("ViscousTs");
  double Phis = container.get<double>("ViscousPhaseShift");
  double tw = container.get<double>("WallThickness");
  double A = container.get<double>("Area");
  int generation = container.get<int>("Generation");

  if (container.get<std::optional<double>>("AirwayColl").has_value())
  {
    airway_params_.airway_coll = *container.get<std::optional<double>>("AirwayColl");
    airway_params_.s_close = *container.get<std::optional<double>>("S_Close");
    airway_params_.s_open = *container.get<std::optional<double>>("S_Open");
    airway_params_.p_crit_open = *container.get<std::optional<double>>("Pcrit_Open");
    airway_params_.p_crit_close = *container.get<std::optional<double>>("Pcrit_Close");
    airway_params_.open_init = *container.get<std::optional<double>>("Open_Init");
  }

  // Correct the velocity profile power
  // this is because the 2.0 is the minimum energy consumtive laminar profile
  if (velPow < 2.0) velPow = 2.0;
  airway_params_.power_velocity_profile = velPow;
  airway_params_.wall_elasticity = Ew;
  airway_params_.poisson_ratio = nu;
  airway_params_.wall_thickness = tw;
  airway_params_.area = A;
  airway_params_.viscous_Ts = Ts;
  airway_params_.viscous_phase_shift = Phis;
  airway_params_.generation = generation;
  airway_params_.branch_length = container.get<std::optional<double>>("BranchLength").value_or(-1);

  return true;
}


/*----------------------------------------------------------------------*
| Read in the RED_ACINUS elements                                       |
*-----------------------------------------------------------------------*/
bool Discret::Elements::RedAcinus::read_element(const std::string& eletype,
    const std::string& distype, const Core::IO::InputParameterContainer& container)
{
  const int ndim = Global::Problem::instance()->n_dim();
  if (ndim != 3)
    FOUR_C_THROW("Problem defined as %dd, but found Reduced dimensional ACINUS element.", ndim);

  // Read number of material model
  int material_id = container.get<int>("MAT");
  set_material(0, Mat::factory(material_id));

  // Read the element type, the element specific variables and store them to acinusParams_
  elem_type_ = container.get<std::string>("TYPE");

  acinus_params_.volume_relaxed = container.get<double>("AcinusVolume");
  acinus_params_.alveolar_duct_volume = container.get<double>("AlveolarDuctVolume");
  acinus_params_.volume_init = acinus_params_.volume_relaxed;
  acinus_params_.generation = -1;

  // Setup material, calls overloaded function setup(linedef) for each Maxwell_0d_acinus material
  std::shared_ptr<Core::Mat::Material> mat = material();
  std::shared_ptr<Mat::Maxwell0dAcinus> acinus_mat =
      std::dynamic_pointer_cast<Mat::Maxwell0dAcinus>(material());
  acinus_mat->setup(container);

  return true;
}


/*----------------------------------------------------------------------*
| Read in the RED_ACINAR_INTER_DEP elements                             |
*-----------------------------------------------------------------------*/
bool Discret::Elements::RedInterAcinarDep::read_element(const std::string& eletype,
    const std::string& distype, const Core::IO::InputParameterContainer& container)
{
  const int ndim = Global::Problem::instance()->n_dim();
  if (ndim != 3)
    FOUR_C_THROW(
        "Problem defined as %dd, but found Reduced dimensional INTER ACINAR DEPENDENCE element.",
        ndim);

  // set generation
  const int generation = -2;
  generation_ = generation;

  // Read number of material model
  int material_id = container.get<int>("MAT");
  set_material(0, Mat::factory(material_id));


  return true;
}


FOUR_C_NAMESPACE_CLOSE
