// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_material_base.hpp"
#include "4C_shell7p_ele_calc_interface.hpp"
#include "4C_shell7p_ele_factory.hpp"
#include "4C_shell7p_ele_neumann_evaluator.hpp"
#include "4C_shell7p_ele_scatra.hpp"
#include "4C_shell7p_ele_scatra_preevaluator.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_function_of_time.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace
{
  inline std::vector<char>& get_mutable_stress_data(
      const Discret::Elements::Shell7pScatra& ele, const Teuchos::ParameterList& params)
  {
    if (ele.is_params_interface())
    {
      return *ele.str_params_interface().stress_data_ptr();
    }
    else
    {
      return *params.get<std::shared_ptr<std::vector<char>>>("stress");
    }
  }

  inline std::vector<char>& get_mutable_strain_data(
      const Discret::Elements::Shell7pScatra& ele, const Teuchos::ParameterList& params)
  {
    if (ele.is_params_interface())
    {
      return *ele.str_params_interface().strain_data_ptr();
    }
    else
    {
      return *params.get<std::shared_ptr<std::vector<char>>>("strain");
    }
  }

  inline Inpar::Solid::StressType get_io_stress_type(
      const Discret::Elements::Shell7pScatra& ele, const Teuchos::ParameterList& params)
  {
    if (ele.is_params_interface())
    {
      return ele.str_params_interface().get_stress_output_type();
    }
    else
    {
      return Teuchos::getIntegralValue<Inpar::Solid::StressType>(params, "iostress");
    }
  }

  inline Inpar::Solid::StrainType get_io_strain_type(
      const Discret::Elements::Shell7pScatra& ele, const Teuchos::ParameterList& params)
  {
    if (ele.is_params_interface())
    {
      return ele.str_params_interface().get_strain_output_type();
    }
    else
    {
      return Teuchos::getIntegralValue<Inpar::Solid::StrainType>(params, "iostrain");
    }
  }
}  // namespace

int Discret::Elements::Shell7pScatra::evaluate(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, Core::Elements::LocationArray& la,
    Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
    Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
    Core::LinAlg::SerialDenseVector& elevec3)

{
  // get params interface pointer
  set_params_interface_ptr(params);

  const Core::Elements::ActionType action = std::invoke(
      [&]()
      {
        if (is_params_interface())
          return str_params_interface().get_action_type();
        else
          return Core::Elements::string_to_action_type(params.get<std::string>("action", "none"));
      });

  // in some cases we need to write/change some data before evaluating
  Discret::Elements::Shell::pre_evaluate_scatra_by_element(*this, params, discretization, la);

  // what should the element do
  switch (action)
  {
    // nonlinear stiffness and internal force vector
    case Core::Elements::struct_calc_nlnstiff:
    {
      shell_interface_->evaluate_nonlinear_force_stiffness_mass(*this, *solid_material(),
          discretization, nodal_directors_, la[0].lm_, params, &elevec1, &elemat1, nullptr);
    }
    break;
    case Core::Elements::struct_calc_linstiff:
    {
      shell_interface_->evaluate_nonlinear_force_stiffness_mass(*this, *solid_material(),
          discretization, nodal_directors_, la[0].lm_, params, &elevec1, &elemat1, &elemat2);
    }
    break;
    case Core::Elements::struct_calc_internalforce:
    {
      shell_interface_->evaluate_nonlinear_force_stiffness_mass(*this, *solid_material(),
          discretization, nodal_directors_, la[0].lm_, params, &elevec1, nullptr, nullptr);
    }
    break;
    case Core::Elements::struct_calc_linstiffmass:
    {
      FOUR_C_THROW("Case not yet implemented: struct_calc_linstiffmass");
    }
    case Core::Elements::struct_calc_nlnstiffmass:   // do mass, stiffness and internal forces
    case Core::Elements::struct_calc_nlnstifflmass:  // do lump mass, stiffness and internal forces
    {
      shell_interface_->evaluate_nonlinear_force_stiffness_mass(*this, *solid_material(),
          discretization, nodal_directors_, la[0].lm_, params, &elevec1, &elemat1, &elemat2);
      if (action == Core::Elements::struct_calc_nlnstifflmass)
        Solid::Utils::Shell::lump_mass_matrix(elemat2);
    }
    break;
    case Core::Elements::struct_calc_recover:
    {
      shell_interface_->recover(*this, discretization, la[0].lm_, params, str_params_interface());
    }
    break;
    case Core::Elements::struct_calc_stress:
    {
      shell_interface_->calculate_stresses_strains(*this, *solid_material(),
          ShellStressIO{get_io_stress_type(*this, params), get_mutable_stress_data(*this, params)},
          ShellStrainIO{get_io_strain_type(*this, params), get_mutable_strain_data(*this, params)},
          discretization, nodal_directors_, la[0].lm_, params);
    }
    break;
    case Core::Elements::struct_calc_energy:
    {
      double int_energy = shell_interface_->calculate_internal_energy(
          *this, *solid_material(), discretization, nodal_directors_, la[0].lm_, params);

      if (is_params_interface())
      {
        // new structural time integration
        str_params_interface().add_contribution_to_energy_type(int_energy, Solid::internal_energy);
      }
      else
      {
        // old structural time integration
        // check length of elevec1
        if (elevec1.length() < 1) FOUR_C_THROW("The given result vector is too short.");
        elevec1(0) = int_energy;
      }
    }
    break;
    case Core::Elements::struct_calc_update_istep:
    {
      shell_interface_->update(
          *this, *solid_material(), discretization, nodal_directors_, la[0].lm_, params);
    }
    break;
    case Core::Elements::struct_calc_reset_istep:
    {
      // Reset of history (if needed)
      shell_interface_->reset_to_last_converged(*this, *solid_material());
    }
    break;
    case Core::Elements::struct_calc_predict:
    case Core::Elements::struct_create_backup:
    case Core::Elements::struct_recover_from_backup:
    {
      // do nothing for now
    }
    break;
    default:
      FOUR_C_THROW("The element action {} is not yet implemented for the Shell element yet",
          action_type_to_string(action).c_str());
  }
  return 0;
}

// Integrate a Surface Neumann boundary condition
int Discret::Elements::Shell7pScatra::evaluate_neumann(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, Core::Conditions::Condition& condition,
    std::vector<int>& la, Core::LinAlg::SerialDenseVector& elevec1,
    Core::LinAlg::SerialDenseMatrix* elemat1)
{
  set_params_interface_ptr(params);

  const double time = std::invoke(
      [&]()
      {
        if (is_params_interface())
          return str_params_interface().get_total_time();
        else
          return params.get("total time", -1.0);
      });

  Discret::Elements::Shell::evaluate_neumann_by_element(
      *this, discretization, condition, la, elevec1, elemat1, time);
  return 0;
}

FOUR_C_NAMESPACE_CLOSE
