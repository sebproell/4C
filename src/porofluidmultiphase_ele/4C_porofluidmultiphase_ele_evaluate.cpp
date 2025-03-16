// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fem_discretization.hpp"
#include "4C_porofluidmultiphase_ele.hpp"
#include "4C_porofluidmultiphase_ele_action.hpp"
#include "4C_porofluidmultiphase_ele_factory.hpp"
#include "4C_porofluidmultiphase_ele_interface.hpp"
#include "4C_porofluidmultiphase_ele_parameter.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                           vuong 08/16 |
 *----------------------------------------------------------------------*/
int Discret::Elements::PoroFluidMultiPhase::evaluate(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, Core::Elements::LocationArray& la,
    Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
    Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
    Core::LinAlg::SerialDenseVector& elevec3)
{
  // we assume here, that numdofpernode is equal for every node within
  // the element and does not change during the computations
  const int numdofpernode = num_dof_per_node(*(nodes()[0]));

  // check for the action parameter
  const auto action = Teuchos::getIntegralValue<POROFLUIDMULTIPHASE::Action>(params, "action");
  switch (action)
  {
    // all physics-related stuff is included in the implementation class(es) that can
    // be used in principle inside any element (at the moment: only PoroFluidMultiPhase element)
    case POROFLUIDMULTIPHASE::calc_mat_and_rhs:
    case POROFLUIDMULTIPHASE::calc_fluid_struct_coupl_mat:
    case POROFLUIDMULTIPHASE::calc_fluid_scatra_coupl_mat:
    case POROFLUIDMULTIPHASE::calc_error:
    case POROFLUIDMULTIPHASE::calc_pres_and_sat:
    case POROFLUIDMULTIPHASE::calc_solidpressure:
    case POROFLUIDMULTIPHASE::calc_porosity:
    case POROFLUIDMULTIPHASE::recon_flux_at_nodes:
    case POROFLUIDMULTIPHASE::calc_phase_velocities:
    case POROFLUIDMULTIPHASE::calc_initial_time_deriv:
    case POROFLUIDMULTIPHASE::calc_valid_dofs:
    case POROFLUIDMULTIPHASE::calc_domain_integrals:
    {
      std::vector<Core::LinAlg::SerialDenseMatrix*> elemat(2);
      elemat[0] = &elemat1;
      elemat[1] = &elemat2;

      std::vector<Core::LinAlg::SerialDenseVector*> elevec(3);
      elevec[0] = &elevec1;
      elevec[1] = &elevec2;
      elevec[2] = &elevec3;

      return Discret::Elements::PoroFluidMultiPhaseFactory::provide_impl(
          shape(), numdofpernode, discretization.name())
          ->evaluate(this, params, discretization, la, elemat, elevec);
      break;
    }
    case POROFLUIDMULTIPHASE::set_timestep_parameter:
    case POROFLUIDMULTIPHASE::set_general_parameter:
      // these actions have already been evaluated during element pre-evaluate
      break;
    default:
    {
      FOUR_C_THROW("Unknown type of action '{}' for PoroFluidMultiPhase", action);
      break;
    }
  }  // switch(action)

  return 0;
}  // Discret::Elements::PoroFluidMultiPhase::Evaluate


/*----------------------------------------------------------------------*
 |  dummy                                                   vuong 08/16 |
 *----------------------------------------------------------------------*/
int Discret::Elements::PoroFluidMultiPhase::evaluate_neumann(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, Core::Conditions::Condition& condition,
    std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1,
    Core::LinAlg::SerialDenseMatrix* elemat1)
{
  FOUR_C_THROW("evaluate_neumann for PoroFluidMultiPhase  not yet implemented!");
  //    The function is just a dummy. For PoroFluidMultiPhase elements, the integration
  //    integration of volume Neumann conditions (body forces) takes place
  //    in the element. We need it there for potential stabilisation terms!
  return 0;
}

/*---------------------------------------------------------------------*
|  Call the element to set all basic parameter             vuong 08/16 |
*----------------------------------------------------------------------*/
void Discret::Elements::PoroFluidMultiPhaseType::pre_evaluate(Core::FE::Discretization& dis,
    Teuchos::ParameterList& p, std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix1,
    std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix2,
    std::shared_ptr<Core::LinAlg::Vector<double>> systemvector1,
    std::shared_ptr<Core::LinAlg::Vector<double>> systemvector2,
    std::shared_ptr<Core::LinAlg::Vector<double>> systemvector3)
{
  const auto action = Teuchos::getIntegralValue<POROFLUIDMULTIPHASE::Action>(p, "action");

  switch (action)
  {
    case POROFLUIDMULTIPHASE::set_general_parameter:
    {
      Discret::Elements::PoroFluidMultiPhaseEleParameter::instance(dis.name())
          ->set_general_parameters(p);

      break;
    }

    case POROFLUIDMULTIPHASE::set_timestep_parameter:
    {
      Discret::Elements::PoroFluidMultiPhaseEleParameter::instance(dis.name())
          ->set_time_step_parameters(p);

      break;
    }
    default:
      // do nothing in all other cases
      break;
  }

  return;
}

FOUR_C_NAMESPACE_CLOSE
