// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_scatra_ele_calc_artery.hpp"

#include "4C_fem_general_extract_values.hpp"
#include "4C_mat_scatra.hpp"
#include "4C_scatra_ele_parameter_std.hpp"
#include "4C_scatra_ele_parameter_timint.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
Discret::Elements::ScaTraEleCalcArtery<distype, probdim>::ScaTraEleCalcArtery(
    const int numdofpernode, const int numscal, const std::string& disname)
    : my::ScaTraEleCalc(numdofpernode, numscal, disname)
{
  // safety check
  // Note: if higher order 1D elements should be used, the approach with adding the length from
  // below
  //        will not work anymore
  if (nen_ != 2)
  {
    FOUR_C_THROW(
        "Only line2 elements supported so far, you have {} nodes, if called with 2D or 3D element, "
        "think again",
        nen_);
  }
  // replace internal variable manager by internal variable manager for artery
  my::scatravarmanager_ =
      std::make_shared<ScaTraEleInternalVariableManagerArtery<nsd_, nen_>>(my::numscal_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
Discret::Elements::ScaTraEleCalcArtery<distype, probdim>*
Discret::Elements::ScaTraEleCalcArtery<distype, probdim>::instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = Core::Utils::make_singleton_map<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleCalcArtery<distype, probdim>>(
            new ScaTraEleCalcArtery<distype, probdim>(numdofpernode, numscal, disname));
      });

  return singleton_map[disname].instance(
      Core::Utils::SingletonAction::create, numdofpernode, numscal, disname);
}

/*----------------------------------------------------------------------*
 | setup element evaluation                            kremheller 03/18 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
int Discret::Elements::ScaTraEleCalcArtery<distype, probdim>::setup_calc(
    Core::Elements::Element* ele, Core::FE::Discretization& discretization)
{
  // base class
  my::setup_calc(ele, discretization);

  // set the artery material in the variable manager
  var_manager()->set_artery_material(ele);

  return 0;
}

/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)              kremheller 04/18 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::Elements::ScaTraEleCalcArtery<distype, probdim>::materials(
    const std::shared_ptr<const Core::Mat::Material> material,  //!< pointer to current material
    const int k,                                                //!< id of current scalar
    double& densn,                                              //!< density at t_(n)
    double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,  //!< density at t_(n+alpha_M)
    double& visc,    //!< fluid viscosity
    const int iquad  //!< id of current gauss point

)
{
  switch (material->material_type())
  {
    case Core::Materials::m_scatra:
    {
      const std::shared_ptr<const Mat::ScatraMat>& actmat =
          std::dynamic_pointer_cast<const Mat::ScatraMat>(material);

      densn = 1.0;
      densam = 1.0;
      densnp = 1.0;

      my::diffmanager_->set_isotropic_diff(actmat->diffusivity(), k);

      break;
    }
    default:
    {
      FOUR_C_THROW("Material type {} is not supported!", material->material_type());
      break;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | set internal variables                              kremheller 03/18 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::Elements::ScaTraEleCalcArtery<distype,
    probdim>::set_internal_variables_for_mat_and_rhs()
{
  var_manager()->set_internal_variables_artery(my::funct_, my::derxy_, my::deriv_, my::xjm_,
      my::ephinp_, my::ephin_, my::ehist_, earterypressurenp_);

  return;
}

/*----------------------------------------------------------------------*
 | extract element based or nodal values               kremheller 03/18 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::Elements::ScaTraEleCalcArtery<distype, probdim>::extract_element_and_node_values(
    Core::Elements::Element* ele, Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, Core::Elements::LocationArray& la)
{
  //---------------------------------------------------------------------------------------------
  //                                 SCATRA
  //---------------------------------------------------------------------------------------------

  // extract local values from the global vectors
  std::shared_ptr<const Core::LinAlg::Vector<double>> hist = discretization.get_state("hist");
  std::shared_ptr<const Core::LinAlg::Vector<double>> phinp = discretization.get_state("phinp");
  if (hist == nullptr || phinp == nullptr)
    FOUR_C_THROW("Cannot get state vector 'hist' and/or 'phinp'");

  // values of scatra field are always in first dofset
  const std::vector<int>& lm = la[0].lm_;
  Core::FE::extract_my_values<Core::LinAlg::Matrix<nen_, 1>>(*hist, my::ehist_, lm);
  Core::FE::extract_my_values<Core::LinAlg::Matrix<nen_, 1>>(*phinp, my::ephinp_, lm);

  if (my::scatraparatimint_->is_gen_alpha() and not my::scatraparatimint_->is_incremental())
  {
    // extract additional local values from global vector
    std::shared_ptr<const Core::LinAlg::Vector<double>> phin = discretization.get_state("phin");
    if (phin == nullptr) FOUR_C_THROW("Cannot get state vector 'phin'");
    Core::FE::extract_my_values<Core::LinAlg::Matrix<nen_, 1>>(*phin, my::ephin_, lm);
  }

  //---------------------------------------------------------------------------------------------
  //                                 CURRENT LENGTH
  //---------------------------------------------------------------------------------------------
  // extract element and node values of the artery
  if (discretization.has_state(1, "curr_seg_lengths"))
  {
    std::shared_ptr<const Core::LinAlg::Vector<double>> curr_seg_lengths =
        discretization.get_state(1, "curr_seg_lengths");
    std::vector<double> seglengths = Core::FE::extract_values(*curr_seg_lengths, la[1].lm_);

    const double curr_ele_length = std::accumulate(seglengths.begin(), seglengths.end(), 0.0);

    static Core::LinAlg::Matrix<probdim, 1> arteryrefpos0;
    for (unsigned int d = 0; d < probdim; ++d) arteryrefpos0(d) = my::xyze_(d, 0);
    static Core::LinAlg::Matrix<probdim, 1> arteryrefpos1;
    for (unsigned int d = 0; d < probdim; ++d) arteryrefpos1(d) = my::xyze_(d, 1);

    static Core::LinAlg::Matrix<probdim, 1> dist0;
    dist0.update(-1.0, arteryrefpos0, 1.0, arteryrefpos1, 0.0);
    const double arteryreflength = dist0.norm2();

    // this is a hack
    // will not work for anything else but line2 elements
    // change in length is simply added to displacement of second node
    for (unsigned int d = 0; d < probdim; ++d)
    {
      my::edispnp_(d, 0) = 0.0;
      my::edispnp_(d, 1) = (curr_ele_length / arteryreflength - 1.0) * dist0(d);
    }

    my::update_node_coordinates();
  }

  int ndsscatra_artery = 1;
  if (discretization.num_dof_sets() == 3) ndsscatra_artery = 2;

  //---------------------------------------------------------------------------------------------
  //                                 ARTERY
  //---------------------------------------------------------------------------------------------
  // extract element and node values of the artery
  if (discretization.has_state(ndsscatra_artery, "one_d_artery_pressure"))
  {
    std::shared_ptr<const Core::LinAlg::Vector<double>> arterypn =
        discretization.get_state(ndsscatra_artery, "one_d_artery_pressure");
    // values of scatra field are always in first dofset
    const std::vector<int>& lm_artery = la[ndsscatra_artery].lm_;
    Core::FE::extract_my_values<Core::LinAlg::Matrix<nen_, 1>>(
        *arterypn, earterypressurenp_, lm_artery);
  }
  else
    FOUR_C_THROW("Something went wrong here, scatra-dis does not have artery primary variable");

  // ---------------------------------------------------------------------
  // call routine for calculation of body force in element nodes
  // (time n+alpha_F for generalized-alpha scheme, at time n+1 otherwise)
  // ---------------------------------------------------------------------
  my::body_force(ele);
  //--------------------------------------------------------------------------------
  // further node-based source terms not given via Neumann volume condition
  // i.e., set special body force for homogeneous isotropic turbulence
  //--------------------------------------------------------------------------------
  my::other_node_based_source_terms(lm, discretization, params);
}

/*-----------------------------------------------------------------------------*
 |  calculation of convective element matrix                                    |
 |  in convective form (OD fluid)                              kremheller 05/18 |
 *-----------------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::Elements::ScaTraEleCalcArtery<distype, probdim>::calc_mat_conv_od_fluid(
    Core::LinAlg::SerialDenseMatrix& emat, const int k, const int ndofpernodefluid,
    const double timefacfac, const double densnp, const Core::LinAlg::Matrix<nsd_, 1>& gradphi)
{
  const double prefac = timefacfac * var_manager()->diam() * var_manager()->diam() / 32.0 /
                        var_manager()->visc() * (-1.0);
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    const int fvi = vi * my::numdofpernode_ + k;
    const double v = prefac * my::funct_(vi);


    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      // get correct factor
      double laplawf(0.0);
      for (unsigned j = 0; j < nsd_; j++) laplawf += my::derxy_(j, ui) * gradphi(j);
      const unsigned fui = ui;
      emat(fvi, fui) += v * laplawf;
    }
  }
  return;
}

// template classes

// 1D elements
template class Discret::Elements::ScaTraEleCalcArtery<Core::FE::CellType::line2, 1>;
template class Discret::Elements::ScaTraEleCalcArtery<Core::FE::CellType::line2, 2>;
template class Discret::Elements::ScaTraEleCalcArtery<Core::FE::CellType::line2, 3>;
template class Discret::Elements::ScaTraEleCalcArtery<Core::FE::CellType::line3, 1>;
// template class Discret::Elements::ScaTraEleCalcStd<Core::FE::CellType::line3,2>;
// template class Discret::Elements::ScaTraEleCalcStd<Core::FE::CellType::line3,3>;

// 2D elements
template class Discret::Elements::ScaTraEleCalcArtery<Core::FE::CellType::tri3, 2>;
template class Discret::Elements::ScaTraEleCalcArtery<Core::FE::CellType::tri3, 3>;
template class Discret::Elements::ScaTraEleCalcArtery<Core::FE::CellType::tri6, 2>;
template class Discret::Elements::ScaTraEleCalcArtery<Core::FE::CellType::quad4, 2>;
template class Discret::Elements::ScaTraEleCalcArtery<Core::FE::CellType::quad4, 3>;
// template class Discret::Elements::ScaTraEleCalcStd<Core::FE::CellType::quad8>;
template class Discret::Elements::ScaTraEleCalcArtery<Core::FE::CellType::quad9, 2>;
// template class Discret::Elements::ScaTraEleCalcStd<Core::FE::CellType::quad9,3>;
template class Discret::Elements::ScaTraEleCalcArtery<Core::FE::CellType::nurbs9, 2>;
// template class Discret::Elements::ScaTraEleCalcStd<Core::FE::CellType::nurbs9,3>;

// 3D elements
template class Discret::Elements::ScaTraEleCalcArtery<Core::FE::CellType::hex8, 3>;
// template class Discret::Elements::ScaTraEleCalcStd<Core::FE::CellType::hex20>;
template class Discret::Elements::ScaTraEleCalcArtery<Core::FE::CellType::hex27, 3>;
template class Discret::Elements::ScaTraEleCalcArtery<Core::FE::CellType::tet4, 3>;
template class Discret::Elements::ScaTraEleCalcArtery<Core::FE::CellType::tet10, 3>;
// template class Discret::Elements::ScaTraEleCalcStd<Core::FE::CellType::wedge6>;
template class Discret::Elements::ScaTraEleCalcArtery<Core::FE::CellType::pyramid5, 3>;
// template class Discret::Elements::ScaTraEleCalcStd<Core::FE::CellType::nurbs27>;s

FOUR_C_NAMESPACE_CLOSE
