/*----------------------------------------------------------------------*/
/*! \file

\brief Internal implementation of poro Fluid element (standard poro fluid)

\level 2


*/
/*----------------------------------------------------------------------*/

#include "fluid_ele_poro.H"
#include "lib_element_integration_select.H"
#include "fluid_ele_action.H"
#include "fluid_ele_parameter_poro.H"

#include "fluid_rotsym_periodicbc.H"

#include "discretization_fem_general_utils_gder2.H"
#include "discretization_fem_general_utils_nurbs_shapefunctions.H"

#include "geometry_position_array.H"

#include "mat_fluidporo.H"
#include "mat_structporo.H"
#include "mat_stvenantkirchhoff.H"
#include "mat_elasthyper.H"

#include "lib_globalproblem.H"

#include "nurbs_discret_nurbs_utils.H"

#include "fluid_ele_calc_poro.H"

#define STAB

template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleCalcPoro<distype>* DRT::ELEMENTS::FluidEleCalcPoro<distype>::Instance(
    ::UTILS::SingletonAction action)
{
  static auto singleton_owner = ::UTILS::MakeSingletonOwner(
      []()
      {
        return std::unique_ptr<DRT::ELEMENTS::FluidEleCalcPoro<distype>>(
            new DRT::ELEMENTS::FluidEleCalcPoro<distype>());
      });

  return singleton_owner.Instance(action);
}


template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleCalcPoro<distype>::FluidEleCalcPoro()
    : DRT::ELEMENTS::FluidEleCalc<distype>::FluidEleCalc(),
      N_XYZ_(true),
      N_XYZ2_(true),
      N_XYZ2full_(true),
      xyze0_(true),
      xyzeold_(true),
      hist_con_(true),
      porosity_(0.0),
      grad_porosity_(true),
      gridvel_int_(true),
      gridvel_n_int_(true),
      convvel_(true),
      gridvel_div_(0.0),
      J_(0.0),
      press_(0.0),
      press_dot_(0.0),
      refgrad_press_(true),
      mat_reac_tensor_(true),
      reac_tensor_(true),
      reac_tensor_linOD_vel_(true),
      reac_tensor_linOD_grid_vel_(true),
      reac_tensor_vel_(true),
      reac_tensor_gridvel_(true),
      reac_tensor_convvel_(true),
      dtau_dphi_(true),
      tau_struct_(0.0),
      mixres_(true),
      struct_mat_(Teuchos::null),
      const_permeability_(true),
      kintype_(INPAR::STR::kinem_vague)
{
  anisotropic_permeability_directions_.resize(nsd_, std::vector<double>(nsd_, 0.0));
  anisotropic_permeability_nodal_coeffs_.resize(nsd_, std::vector<double>(nen_, 0.0));

  // change pointer to parameter list in base class to poro parameters
  Base::fldpara_ = DRT::ELEMENTS::FluidEleParameterPoro::Instance();
  // this is just for convenience. The same pointer as above to circumvent casts when accessing poro
  // specific paramters
  porofldpara_ = DRT::ELEMENTS::FluidEleParameterPoro::Instance();
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::PreEvaluate(
    Teuchos::ParameterList& params, DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization)
{
  // do nothing
}

template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcPoro<distype>::EvaluateService(DRT::ELEMENTS::Fluid* ele,
    Teuchos::ParameterList& params, Teuchos::RCP<MAT::Material>& mat,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3)
{
  // get the action required
  const auto act = DRT::INPUT::get<FLD::Action>(params, "action");

  switch (act)
  {
    case FLD::calc_volume:
    {
      return ComputeVolume(params, ele, discretization, lm, elevec1);
      break;
    }
    case FLD::calc_fluid_error:
    {
      return ComputeError(ele, params, mat, discretization, lm, elevec1);
      break;
    }
    default:
      dserror("unknown action for EvaluateService() in poro fluid element");
      break;
  }
  return -1;
}

template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcPoro<distype>::Evaluate(DRT::ELEMENTS::Fluid* ele,
    DRT::Discretization& discretization, const std::vector<int>& lm, Teuchos::ParameterList& params,
    Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elemat1_epetra,
    Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
    Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra,
    bool offdiag)
{
  Teuchos::RCP<const MAT::FluidPoro> actmat = Teuchos::rcp_static_cast<const MAT::FluidPoro>(mat);
  const_permeability_ = (actmat->PermeabilityFunction() == MAT::PAR::constant);

  auto* poroele = dynamic_cast<DRT::ELEMENTS::FluidPoro*>(ele);

  if (poroele)
  {
    kintype_ = poroele->KinematicType();
    anisotropic_permeability_directions_ = poroele->GetAnisotropicPermeabilityDirections();
    anisotropic_permeability_nodal_coeffs_ = poroele->GetAnisotropicPermeabilityNodalCoeffs();
  }

  if (not offdiag)
  {  // evaluate diagonal block (pure fluid block)
    return Evaluate(ele, discretization, lm, params, mat, elemat1_epetra, elemat2_epetra,
        elevec1_epetra, elevec2_epetra, elevec3_epetra, Base::intpoints_);
  }
  else
  {  // evaluate off diagonal block (coupling block)
    return EvaluateOD(ele, discretization, lm, params, mat, elemat1_epetra, elemat2_epetra,
        elevec1_epetra, elevec2_epetra, elevec3_epetra, Base::intpoints_);
  }
}

template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcPoro<distype>::Evaluate(DRT::ELEMENTS::Fluid* ele,
    DRT::Discretization& discretization, const std::vector<int>& lm, Teuchos::ParameterList& params,
    Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elemat1_epetra,
    Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
    Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra,
    const DRT::UTILS::GaussIntegration& intpoints)
{
  //----------------------------------------------------------------
  // Now do the nurbs specific stuff (for isogeometric elements)
  //----------------------------------------------------------------
  if (Base::isNurbs_)
  {
    // access knots and weights for this element
    bool zero_size =
        DRT::NURBS::GetMyNurbsKnotsAndWeights(discretization, ele, Base::myknots_, Base::weights_);

    // if we have a zero sized element due to a interpolated point -> exit here
    if (zero_size) return (0);
  }  // Nurbs specific stuff

  // set element id
  Base::eid_ = ele->Id();
  // get structure material
  GetStructMaterial(ele);

  // rotationally symmetric periodic bc's: do setup for current element
  // (only required to be set up for routines "ExtractValuesFromGlobalVector")
  Base::rotsymmpbc_->Setup(ele);

  // ---------------------------------------------------------------------
  // call routine for calculation of body force in element nodes,
  // with pressure gradient prescribed as body force included for turbulent
  // channel flow and with scatra body force included for variable-density flow
  // (evaluation at time n+alpha_F for generalized-alpha scheme,
  //  and at time n+1 otherwise)
  // ---------------------------------------------------------------------
  static LINALG::Matrix<nsd_, nen_> ebofoaf(true);
  ebofoaf.Clear();
  static LINALG::Matrix<nsd_, nen_> eprescpgaf(true);
  eprescpgaf.Clear();
  static LINALG::Matrix<nen_, 1> escabofoaf(true);
  escabofoaf.Clear();
  Base::BodyForce(ele, ebofoaf, eprescpgaf, escabofoaf);

  // ---------------------------------------------------------------------
  // get all general state vectors: velocity/pressure, acceleration
  // and history
  // velocity/pressure values are at time n+alpha_F/n+alpha_M for
  // generalized-alpha scheme and at time n+1/n for all other schemes
  // acceleration values are at time n+alpha_M for
  // generalized-alpha scheme and at time n+1 for all other schemes
  // ---------------------------------------------------------------------
  // fill the local element vector/matrix with the global values
  // af_genalpha: velocity/pressure at time n+alpha_F
  // np_genalpha: velocity at time n+alpha_F, pressure at time n+1
  // ost:         velocity/pressure at time n+1
  static LINALG::Matrix<nsd_, nen_> evelaf(true);
  evelaf.Clear();
  static LINALG::Matrix<nen_, 1> epreaf(true);
  epreaf.Clear();
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, &evelaf, &epreaf, "velaf");

  // np_genalpha: additional vector for velocity at time n+1
  static LINALG::Matrix<nsd_, nen_> evelnp(true);
  evelnp.Clear();
  static LINALG::Matrix<nen_, 1> eprenp(true);
  eprenp.Clear();
  if (Base::fldparatimint_->IsGenalphaNP())
    Base::ExtractValuesFromGlobalVector(
        discretization, lm, *Base::rotsymmpbc_, &evelnp, &eprenp, "velnp");

  static LINALG::Matrix<nsd_, nen_> emhist(true);
  emhist.Clear();
  static LINALG::Matrix<nen_, 1> echist(true);
  echist.Clear();
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, &emhist, &echist, "hist");

  static LINALG::Matrix<nsd_, nen_> eaccam(true);
  static LINALG::Matrix<nen_, 1> epressam_timederiv(true);
  eaccam.Clear();
  epressam_timederiv.Clear();

  if (Base::fldparatimint_->IsGenalpha())
    Base::ExtractValuesFromGlobalVector(
        discretization, lm, *Base::rotsymmpbc_, &eaccam, &epressam_timederiv, "accam");

  static LINALG::Matrix<nen_, 1> epressn_timederiv(true);
  epressn_timederiv.Clear();
  if (Base::fldparatimint_->IsGenalpha())
    Base::ExtractValuesFromGlobalVector(
        discretization, lm, *Base::rotsymmpbc_, nullptr, &epressn_timederiv, "accn");

  static LINALG::Matrix<nen_, 1> epren(true);
  epren.Clear();
  static LINALG::Matrix<nsd_, nen_> eveln(true);
  eveln.Clear();
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, &eveln, &epren, "veln");

  static LINALG::Matrix<nen_, 1> epressnp_timederiv(true);
  epressnp_timederiv.Clear();
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, nullptr, &epressnp_timederiv, "accnp");

  static LINALG::Matrix<nen_, 1> escaaf(true);
  escaaf.Clear();
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, nullptr, &escaaf, "scaaf");

  // ---------------------------------------------------------------------
  // get additional state vectors for ALE case: grid displacement and vel.
  // ---------------------------------------------------------------------
  static LINALG::Matrix<nsd_, nen_> edispnp(true);
  edispnp.Clear();
  static LINALG::Matrix<nsd_, nen_> egridv(true);
  egridv.Clear();
  static LINALG::Matrix<nsd_, nen_> egridvn(true);
  egridvn.Clear();
  static LINALG::Matrix<nsd_, nen_> edispn(true);
  edispn.Clear();

  LINALG::Matrix<nen_, 1> eporositynp(true);

  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, &edispnp, nullptr, "dispnp");
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, &egridv, nullptr, "gridv");
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, &egridvn, nullptr, "gridvn");
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, &edispn, nullptr, "dispn");

  // get node coordinates and number of elements per node
  GEO::fillInitialPositionArray<distype, nsd_, LINALG::Matrix<nsd_, nen_>>(ele, Base::xyze_);

  // construct views
  LINALG::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_> elemat1(elemat1_epetra, true);
  // LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_>
  // elemat2(elemat2_epetra,true);
  LINALG::Matrix<(nsd_ + 1) * nen_, 1> elevec1(elevec1_epetra, true);
  // elevec2 and elevec3 are currently not in use

  PreEvaluate(params, ele, discretization);

  // call inner evaluate (does not know about DRT element or discretization object)
  int result = Evaluate(params, ebofoaf, elemat1, elevec1, evelaf, epreaf, evelnp, eveln, eprenp,
      epren, emhist, echist, epressnp_timederiv, epressam_timederiv, epressn_timederiv, eaccam,
      edispnp, edispn, egridv, egridvn, escaaf, nullptr, nullptr, nullptr, mat, ele->IsAle(),
      intpoints);

  return result;
}

template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcPoro<distype>::EvaluateOD(DRT::ELEMENTS::Fluid* ele,
    DRT::Discretization& discretization, const std::vector<int>& lm, Teuchos::ParameterList& params,
    Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elemat1_epetra,
    Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
    Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra,
    const DRT::UTILS::GaussIntegration& intpoints)
{
  //----------------------------------------------------------------
  // Now do the nurbs specific stuff (for isogeometric elements)
  //----------------------------------------------------------------
  if (Base::isNurbs_)
  {
    // access knots and weights for this element
    bool zero_size =
        DRT::NURBS::GetMyNurbsKnotsAndWeights(discretization, ele, Base::myknots_, Base::weights_);

    // if we have a zero sized element due to a interpolated point -> exit here
    if (zero_size) return (0);
  }  // Nurbs specific stuff

  // set element id
  Base::eid_ = ele->Id();

  // get structure material
  GetStructMaterial(ele);

  // rotationally symmetric periodic bc's: do setup for current element
  // (only required to be set up for routines "ExtractValuesFromGlobalVector")
  Base::rotsymmpbc_->Setup(ele);

  // construct views
  LINALG::Matrix<(nsd_ + 1) * nen_, nsd_ * nen_> elemat1(elemat1_epetra, true);
  //  LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_>
  //  elemat2(elemat2_epetra,true);
  LINALG::Matrix<(nsd_ + 1) * nen_, 1> elevec1(elevec1_epetra, true);
  // elevec2 and elevec3 are currently not in use

  // ---------------------------------------------------------------------
  // call routine for calculation of body force in element nodes,
  // with pressure gradient prescribed as body force included for turbulent
  // channel flow and with scatra body force included for variable-density flow
  // (evaluation at time n+alpha_F for generalized-alpha scheme,
  //  and at time n+1 otherwise)
  // ---------------------------------------------------------------------
  static LINALG::Matrix<nsd_, nen_> ebofoaf(true);
  static LINALG::Matrix<nsd_, nen_> eprescpgaf(true);
  static LINALG::Matrix<nen_, 1> escabofoaf(true);
  ebofoaf.Clear();
  eprescpgaf.Clear();
  escabofoaf.Clear();
  Base::BodyForce(ele, ebofoaf, eprescpgaf, escabofoaf);

  // ---------------------------------------------------------------------
  // get all general state vectors: velocity/pressure, acceleration
  // and history
  // velocity/pressure values are at time n+alpha_F/n+alpha_M for
  // generalized-alpha scheme and at time n+1/n for all other schemes
  // acceleration values are at time n+alpha_M for
  // generalized-alpha scheme and at time n+1 for all other schemes
  // ---------------------------------------------------------------------
  // fill the local element vector/matrix with the global values
  // af_genalpha: velocity/pressure at time n+alpha_F
  // np_genalpha: velocity at time n+alpha_F, pressure at time n+1
  // ost:         velocity/pressure at time n+1
  static LINALG::Matrix<nsd_, nen_> evelaf(true);
  static LINALG::Matrix<nen_, 1> epreaf(true);
  evelaf.Clear();
  epreaf.Clear();
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, &evelaf, &epreaf, "velaf");

  // np_genalpha: additional vector for velocity at time n+1
  static LINALG::Matrix<nsd_, nen_> evelnp(true);
  static LINALG::Matrix<nen_, 1> eprenp(true);
  evelnp.Clear();
  eprenp.Clear();
  if (Base::fldparatimint_->IsGenalphaNP())
    Base::ExtractValuesFromGlobalVector(
        discretization, lm, *Base::rotsymmpbc_, &evelnp, &eprenp, "velnp");

  static LINALG::Matrix<nsd_, nen_> eveln(true);
  static LINALG::Matrix<nen_, 1> epren(true);
  eveln.Clear();
  epren.Clear();
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, &eveln, &epren, "veln");

  static LINALG::Matrix<nsd_, nen_> eaccam(true);
  static LINALG::Matrix<nen_, 1> epressam_timederiv(true);
  eaccam.Clear();
  epressam_timederiv.Clear();

  if (Base::fldparatimint_->IsGenalpha())
    Base::ExtractValuesFromGlobalVector(
        discretization, lm, *Base::rotsymmpbc_, &eaccam, &epressam_timederiv, "accam");

  static LINALG::Matrix<nen_, 1> epressn_timederiv(true);
  epressn_timederiv.Clear();
  if (Base::fldparatimint_->IsGenalpha())
    Base::ExtractValuesFromGlobalVector(
        discretization, lm, *Base::rotsymmpbc_, nullptr, &epressn_timederiv, "accn");

  static LINALG::Matrix<nen_, 1> epressnp_timederiv(true);
  epressnp_timederiv.Clear();
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, nullptr, &epressnp_timederiv, "accnp");

  static LINALG::Matrix<nen_, 1> escaaf(true);
  epressnp_timederiv.Clear();
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, nullptr, &escaaf, "scaaf");

  static LINALG::Matrix<nsd_, nen_> emhist(true);
  static LINALG::Matrix<nen_, 1> echist(true);
  emhist.Clear();
  echist.Clear();
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, &emhist, &echist, "hist");

  // ---------------------------------------------------------------------
  // get additional state vectors for ALE case: grid displacement and vel.
  // ---------------------------------------------------------------------
  static LINALG::Matrix<nsd_, nen_> edispnp(true);
  edispnp.Clear();
  static LINALG::Matrix<nsd_, nen_> edispn(true);
  edispn.Clear();
  static LINALG::Matrix<nsd_, nen_> egridv(true);
  egridv.Clear();
  static LINALG::Matrix<nsd_, nen_> egridvn(true);
  egridvn.Clear();

  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, &edispnp, nullptr, "dispnp");
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, &edispn, nullptr, "dispn");
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, &egridv, nullptr, "gridv");
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, &egridvn, nullptr, "gridvn");

  // get node coordinates and number of elements per node
  GEO::fillInitialPositionArray<distype, nsd_, LINALG::Matrix<nsd_, nen_>>(ele, Base::xyze_);

  PreEvaluate(params, ele, discretization);

  // call inner evaluate (does not know about DRT element or discretization object)
  return EvaluateOD(params, ebofoaf, elemat1, elevec1, evelaf, epreaf, evelnp, eveln, eprenp, epren,
      epressnp_timederiv, epressam_timederiv, epressn_timederiv, eaccam, edispnp, edispn, egridv,
      egridvn, escaaf, emhist, echist, nullptr, mat, ele->IsAle(), intpoints);
}

template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcPoro<distype>::Evaluate(Teuchos::ParameterList& params,
    const LINALG::Matrix<nsd_, nen_>& ebofoaf,
    LINALG::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& elemat1,
    LINALG::Matrix<(nsd_ + 1) * nen_, 1>& elevec1, const LINALG::Matrix<nsd_, nen_>& evelaf,
    const LINALG::Matrix<nen_, 1>& epreaf, const LINALG::Matrix<nsd_, nen_>& evelnp,
    const LINALG::Matrix<nsd_, nen_>& eveln, const LINALG::Matrix<nen_, 1>& eprenp,
    const LINALG::Matrix<nen_, 1>& epren, const LINALG::Matrix<nsd_, nen_>& emhist,
    const LINALG::Matrix<nen_, 1>& echist, const LINALG::Matrix<nen_, 1>& epressnp_timederiv,
    const LINALG::Matrix<nen_, 1>& epressam_timederiv,
    const LINALG::Matrix<nen_, 1>& epressn_timederiv, const LINALG::Matrix<nsd_, nen_>& eaccam,
    const LINALG::Matrix<nsd_, nen_>& edispnp, const LINALG::Matrix<nsd_, nen_>& edispn,
    const LINALG::Matrix<nsd_, nen_>& egridv, const LINALG::Matrix<nsd_, nen_>& egridvn,
    const LINALG::Matrix<nen_, 1>& escaaf, const LINALG::Matrix<nen_, 1>* eporositynp,
    const LINALG::Matrix<nen_, 1>* eporositydot, const LINALG::Matrix<nen_, 1>* eporositydotn,
    Teuchos::RCP<MAT::Material> mat, bool isale, const DRT::UTILS::GaussIntegration& intpoints)
{
  // flag for higher order elements
  Base::is_higher_order_ele_ = IsHigherOrder<distype>::ishigherorder;
  // overrule higher_order_ele if input-parameter is set
  // this might be interesting for fast (but slightly
  // less accurate) computations
  if (Base::fldpara_->IsInconsistent()) Base::is_higher_order_ele_ = false;

  // ---------------------------------------------------------------------
  // call routine for calculating element matrix and right hand side
  // ---------------------------------------------------------------------
  Sysmat(params, ebofoaf, evelaf, evelnp, eveln, epreaf, eprenp, epren, eaccam, emhist, echist,
      epressnp_timederiv, epressam_timederiv, epressn_timederiv, edispnp, edispn, egridv, egridvn,
      escaaf, eporositynp, eporositydot, eporositydotn, elemat1, elevec1, mat, isale, intpoints);

  return 0;
}

template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcPoro<distype>::EvaluateOD(Teuchos::ParameterList& params,
    const LINALG::Matrix<nsd_, nen_>& ebofoaf,
    LINALG::Matrix<(nsd_ + 1) * nen_, nsd_ * nen_>& elemat1,
    LINALG::Matrix<(nsd_ + 1) * nen_, 1>& elevec1, const LINALG::Matrix<nsd_, nen_>& evelaf,
    const LINALG::Matrix<nen_, 1>& epreaf, const LINALG::Matrix<nsd_, nen_>& evelnp,
    const LINALG::Matrix<nsd_, nen_>& eveln, const LINALG::Matrix<nen_, 1>& eprenp,
    const LINALG::Matrix<nen_, 1>& epren, const LINALG::Matrix<nen_, 1>& epressnp_timederiv,
    const LINALG::Matrix<nen_, 1>& epressam_timederiv,
    const LINALG::Matrix<nen_, 1>& epressn_timederiv, const LINALG::Matrix<nsd_, nen_>& eaccam,
    const LINALG::Matrix<nsd_, nen_>& edispnp, const LINALG::Matrix<nsd_, nen_>& edispn,
    const LINALG::Matrix<nsd_, nen_>& egridv, const LINALG::Matrix<nsd_, nen_>& egridvn,
    const LINALG::Matrix<nen_, 1>& escaaf, const LINALG::Matrix<nsd_, nen_>& emhist,
    const LINALG::Matrix<nen_, 1>& echist, const LINALG::Matrix<nen_, 1>* eporositynp,
    Teuchos::RCP<MAT::Material> mat, bool isale, const DRT::UTILS::GaussIntegration& intpoints)
{
  // flag for higher order elements
  Base::is_higher_order_ele_ = IsHigherOrder<distype>::ishigherorder;
  // overrule higher_order_ele if input-parameter is set
  // this might be interesting for fast (but slightly
  // less accurate) computations
  if (Base::fldpara_->IsInconsistent()) Base::is_higher_order_ele_ = false;

  // stationary formulation does not support ALE formulation
  // if (isale and Base::fldparatimint_->IsStationary())
  //  dserror("No ALE support within stationary fluid solver.");

  // ---------------------------------------------------------------------
  // call routine for calculating element matrix and right hand side
  // ---------------------------------------------------------------------
  SysmatOD(params, ebofoaf, evelaf, evelnp, eveln, epreaf, eprenp, epren, eaccam,
      epressnp_timederiv, epressam_timederiv, epressn_timederiv, edispnp, edispn, egridv, egridvn,
      escaaf, emhist, echist, eporositynp, elemat1, elevec1, mat, isale, intpoints);

  return 0;
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::Sysmat(Teuchos::ParameterList& params,
    const LINALG::Matrix<nsd_, nen_>& ebofoaf, const LINALG::Matrix<nsd_, nen_>& evelaf,
    const LINALG::Matrix<nsd_, nen_>& evelnp, const LINALG::Matrix<nsd_, nen_>& eveln,
    const LINALG::Matrix<nen_, 1>& epreaf, const LINALG::Matrix<nen_, 1>& eprenp,
    const LINALG::Matrix<nen_, 1>& epren, const LINALG::Matrix<nsd_, nen_>& eaccam,
    const LINALG::Matrix<nsd_, nen_>& emhist, const LINALG::Matrix<nen_, 1>& echist,
    const LINALG::Matrix<nen_, 1>& epressnp_timederiv,
    const LINALG::Matrix<nen_, 1>& epressam_timederiv,
    const LINALG::Matrix<nen_, 1>& epressn_timederiv, const LINALG::Matrix<nsd_, nen_>& edispnp,
    const LINALG::Matrix<nsd_, nen_>& edispn, const LINALG::Matrix<nsd_, nen_>& egridv,
    const LINALG::Matrix<nsd_, nen_>& egridvn, const LINALG::Matrix<nen_, 1>& escaaf,
    const LINALG::Matrix<nen_, 1>* eporositynp, const LINALG::Matrix<nen_, 1>* eporositydot,
    const LINALG::Matrix<nen_, 1>* eporositydotn,
    LINALG::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& estif,
    LINALG::Matrix<(nsd_ + 1) * nen_, 1>& eforce, Teuchos::RCP<const MAT::Material> material,
    bool isale, const DRT::UTILS::GaussIntegration& intpoints)
{
  //------------------------------------------------------------------------
  //  preliminary definitions and evaluations
  //------------------------------------------------------------------------
  // definition of matrices (static to avoid unnecessary reallocation of memory)
  static LINALG::Matrix<nen_ * nsd_, nen_ * nsd_> estif_u(true);
  static LINALG::Matrix<nen_ * nsd_, nen_> estif_p_v(true);
  static LINALG::Matrix<nen_, nen_ * nsd_> estif_q_u(true);
  static LINALG::Matrix<nen_, nen_> ppmat(true);

  estif_u.Clear();
  estif_p_v.Clear();
  estif_q_u.Clear();
  ppmat.Clear();

  // definition of vectors (static to avoid unnecessary reallocation of memory)
  LINALG::Matrix<nen_, 1> preforce(true);
  LINALG::Matrix<nsd_, nen_> velforce(true);
  preforce.Clear();
  velforce.Clear();

  // material coordinates xyze0
  xyze0_ = Base::xyze_;
  xyzeold_ = Base::xyze_;

  {
    Base::xyze_ += edispnp;
    xyzeold_ += edispn;
  }

  //------------------------------------------------------------------------
  // potential evaluation of material parameters, subgrid viscosity
  // and/or stabilization parameters at element center
  //------------------------------------------------------------------------
  // evaluate shape functions and derivatives at element center
  Base::EvalShapeFuncAndDerivsAtEleCenter();

  //------------------------------------------------------------------------
  //  start loop over integration points
  //------------------------------------------------------------------------
  GaussPointLoop(params, ebofoaf, evelaf, evelnp, eveln, epreaf, eprenp, epren, epressnp_timederiv,
      epressam_timederiv, epressn_timederiv, eaccam, edispnp, edispn, egridv, egridvn, escaaf,
      emhist, echist, eporositynp, eporositydot, eporositydotn, estif_u, estif_p_v, estif_q_u,
      ppmat, preforce, velforce, material, intpoints);

  //------------------------------------------------------------------------
  //  end loop over integration points
  //------------------------------------------------------------------------

  //------------------------------------------------------------------------
  //  add contributions to element matrix and right-hand-side vector
  //------------------------------------------------------------------------
  // add pressure part to right-hand-side vector
  for (int vi = 0; vi < nen_; ++vi)
  {
    eforce(Base::numdofpernode_ * vi + nsd_) += preforce(vi);
  }

  // add velocity part to right-hand-side vector
  for (int vi = 0; vi < nen_; ++vi)
  {
    for (int idim = 0; idim < nsd_; ++idim)
    {
      eforce(Base::numdofpernode_ * vi + idim) += velforce(idim, vi);
    }
  }

  // add pressure-pressure part to matrix
  for (int ui = 0; ui < nen_; ++ui)
  {
    const int fuipp = Base::numdofpernode_ * ui + nsd_;

    for (int vi = 0; vi < nen_; ++vi)
    {
      const int numdof_vi_p_nsd = Base::numdofpernode_ * vi + nsd_;

      estif(numdof_vi_p_nsd, fuipp) += ppmat(vi, ui);
    }
  }

  // add velocity-velocity part to matrix
  for (int ui = 0; ui < nen_; ++ui)
  {
    const int numdof_ui = Base::numdofpernode_ * ui;
    const int nsd_ui = nsd_ * ui;

    for (int jdim = 0; jdim < nsd_; ++jdim)
    {
      const int numdof_ui_jdim = numdof_ui + jdim;
      const int nsd_ui_jdim = nsd_ui + jdim;

      for (int vi = 0; vi < nen_; ++vi)
      {
        const int numdof_vi = Base::numdofpernode_ * vi;
        const int nsd_vi = nsd_ * vi;

        for (int idim = 0; idim < nsd_; ++idim)
        {
          estif(numdof_vi + idim, numdof_ui_jdim) += estif_u(nsd_vi + idim, nsd_ui_jdim);
        }
      }
    }
  }

  // add velocity-pressure part to matrix
  for (int ui = 0; ui < nen_; ++ui)
  {
    const int numdof_ui_nsd = Base::numdofpernode_ * ui + nsd_;

    for (int vi = 0; vi < nen_; ++vi)
    {
      const int nsd_vi = nsd_ * vi;
      const int numdof_vi = Base::numdofpernode_ * vi;

      for (int idim = 0; idim < nsd_; ++idim)
      {
        estif(numdof_vi + idim, numdof_ui_nsd) += estif_p_v(nsd_vi + idim, ui);
      }
    }
  }

  // add pressure-velocity part to matrix
  for (int ui = 0; ui < nen_; ++ui)
  {
    const int numdof_ui = Base::numdofpernode_ * ui;
    const int nsd_ui = nsd_ * ui;

    for (int jdim = 0; jdim < nsd_; ++jdim)
    {
      const int numdof_ui_jdim = numdof_ui + jdim;
      const int nsd_ui_jdim = nsd_ui + jdim;

      for (int vi = 0; vi < nen_; ++vi)
        estif(Base::numdofpernode_ * vi + nsd_, numdof_ui_jdim) += estif_q_u(vi, nsd_ui_jdim);
    }
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::SysmatOD(Teuchos::ParameterList& params,
    const LINALG::Matrix<nsd_, nen_>& ebofoaf, const LINALG::Matrix<nsd_, nen_>& evelaf,
    const LINALG::Matrix<nsd_, nen_>& evelnp, const LINALG::Matrix<nsd_, nen_>& eveln,
    const LINALG::Matrix<nen_, 1>& epreaf, const LINALG::Matrix<nen_, 1>& eprenp,
    const LINALG::Matrix<nen_, 1>& epren, const LINALG::Matrix<nsd_, nen_>& eaccam,
    const LINALG::Matrix<nen_, 1>& epressnp_timederiv,
    const LINALG::Matrix<nen_, 1>& epressam_timederiv,
    const LINALG::Matrix<nen_, 1>& epressn_timederiv, const LINALG::Matrix<nsd_, nen_>& edispnp,
    const LINALG::Matrix<nsd_, nen_>& edispn, const LINALG::Matrix<nsd_, nen_>& egridv,
    const LINALG::Matrix<nsd_, nen_>& egridvn, const LINALG::Matrix<nen_, 1>& escaaf,
    const LINALG::Matrix<nsd_, nen_>& emhist, const LINALG::Matrix<nen_, 1>& echist,
    const LINALG::Matrix<nen_, 1>* eporositynp,
    LINALG::Matrix<(nsd_ + 1) * nen_, nsd_ * nen_>& ecoupl,
    LINALG::Matrix<(nsd_ + 1) * nen_, 1>& eforce, Teuchos::RCP<const MAT::Material> material,
    bool isale, const DRT::UTILS::GaussIntegration& intpoints)
{
  //------------------------------------------------------------------------
  //  preliminary definitions and evaluations
  //------------------------------------------------------------------------
  // definition of matrices
  static LINALG::Matrix<nen_ * nsd_, nen_ * nsd_> ecoupl_u(
      true);  // coupling matrix for momentum equation
  static LINALG::Matrix<nen_, nen_ * nsd_> ecoupl_p(
      true);  // coupling matrix for continuity equation
  // LINALG::Matrix<(nsd_ + 1) * nen_, nen_ * nsd_> emesh(true); //
  // linearisation of mesh motion

  ecoupl_u.Clear();
  ecoupl_p.Clear();

  // material coordinates xyze0
  xyze0_ = Base::xyze_;
  xyzeold_ = Base::xyze_;

  // add displacement when fluid nodes move in the ALE case
  // if (isale)
  // if(kintype_!=INPAR::STR::kinem_linear)
  {
    Base::xyze_ += edispnp;
    xyzeold_ += edispn;
  }

  //------------------------------------------------------------------------
  // potential evaluation of material parameters, subgrid viscosity
  // and/or stabilization parameters at element center
  //------------------------------------------------------------------------
  // evaluate shape functions and derivatives at element center
  Base::EvalShapeFuncAndDerivsAtEleCenter();

  //------------------------------------------------------------------------
  //  start loop over integration points
  //------------------------------------------------------------------------
  GaussPointLoopOD(params, ebofoaf, evelaf, evelnp, eveln, epreaf, eprenp, epren,
      epressnp_timederiv, epressam_timederiv, epressn_timederiv, eaccam, edispnp, edispn, egridv,
      egridvn, escaaf, emhist, echist, eporositynp, eforce, ecoupl_u, ecoupl_p, material,
      intpoints);
  //------------------------------------------------------------------------
  //  end loop over integration points
  //------------------------------------------------------------------------

  //------------------------------------------------------------------------
  //  add contributions to element matrix
  //------------------------------------------------------------------------

  // add fluid velocity-structure displacement part to matrix
  for (int ui = 0; ui < nen_; ++ui)
  {
    const int nsd_ui = nsd_ * ui;

    for (int jdim = 0; jdim < nsd_; ++jdim)
    {
      const int nsd_ui_jdim = nsd_ui + jdim;

      for (int vi = 0; vi < nen_; ++vi)
      {
        const int numdof_vi = Base::numdofpernode_ * vi;
        const int nsd_vi = nsd_ * vi;

        for (int idim = 0; idim < nsd_; ++idim)
        {
          ecoupl(numdof_vi + idim, nsd_ui_jdim) += ecoupl_u(nsd_vi + idim, nsd_ui_jdim);
        }
      }
    }
  }

  // add fluid pressure-structure displacement part to matrix
  for (int ui = 0; ui < nen_; ++ui)
  {
    const int nsd_ui = nsd_ * ui;

    for (int jdim = 0; jdim < nsd_; ++jdim)
    {
      const int nsd_ui_jdim = nsd_ui + jdim;

      for (int vi = 0; vi < nen_; ++vi)
      {
        ecoupl(Base::numdofpernode_ * vi + nsd_, nsd_ui_jdim) += ecoupl_p(vi, nsd_ui_jdim);
      }
    }
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::EvaluatePressureEquation(
    Teuchos::ParameterList& params, const double& timefacfacpre, const double& rhsfac,
    const double& dphi_dp, const double& dphi_dJ, const double& dphi_dJdp, const double& dphi_dpp,
    const LINALG::Matrix<nen_, 1>* eporositydot, const LINALG::Matrix<nen_, 1>* eporositydotn,
    const LINALG::Matrix<nen_, 1>& echist, const LINALG::Matrix<nsd_, nen_>& dgradphi_dp,
    LINALG::Matrix<nen_, nen_ * nsd_>& estif_q_u, LINALG::Matrix<nen_, nen_>& ppmat,
    LINALG::Matrix<nen_, 1>& preforce)
{
  // first evaluate terms without porosity time derivative
  EvaluatePressureEquationNonTransient(params, timefacfacpre, rhsfac, dphi_dp, dphi_dJ, dphi_dJdp,
      dphi_dpp, dgradphi_dp, estif_q_u, ppmat, preforce);

  // now the porosity time derivative (different for standard poro and other poro elements)
  if (!porofldpara_->IsStationaryConti())
  {
    // inertia terms on the right hand side for instationary fluids
    if (Base::fldparatimint_->IsGenalpha())
    {
      for (int vi = 0; vi < nen_; ++vi)
        preforce(vi) -= timefacfacpre * (press_dot_ * dphi_dp) * Base::funct_(vi);
    }
    else  // one step theta
    {
      for (int vi = 0; vi < nen_; ++vi)
        preforce(vi) -= Base::fac_ * (press_ * dphi_dp) * Base::funct_(vi);

      const double rhsfac_rhscon = rhsfac * dphi_dp * Base::rhscon_;
      for (int vi = 0; vi < nen_; ++vi)
      {
        /* additional rhs term of continuity equation */
        preforce(vi) += rhsfac_rhscon * Base::funct_(vi);
      }
    }

    for (int vi = 0; vi < nen_; ++vi)
      preforce(vi) -= rhsfac * Base::funct_(vi) * dphi_dJ * J_ * gridvel_div_;

    // additional left hand side term as history values are multiplied by dphi_dp^(n+1)
    for (int vi = 0; vi < nen_; ++vi)
    {
      const double factor = timefacfacpre * Base::funct_(vi) * Base::rhscon_ * dphi_dpp;
      for (int ui = 0; ui < nen_; ++ui)
      {
        ppmat(vi, ui) -= factor * Base::funct_(ui);
      }
    }

    // in case of reactive porous medium : additional rhs term
    double refporositydot = struct_mat_->RefPorosityTimeDeriv();
    for (int vi = 0; vi < nen_; ++vi)
    {
      preforce(vi) -= rhsfac * refporositydot * Base::funct_(vi);
    }
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::EvaluatePressureEquationNonTransient(
    Teuchos::ParameterList& params, const double& timefacfacpre, const double& rhsfac,
    const double& dphi_dp, const double& dphi_dJ, const double& dphi_dJdp, const double& dphi_dpp,
    const LINALG::Matrix<nsd_, nen_>& dgradphi_dp, LINALG::Matrix<nen_, nen_ * nsd_>& estif_q_u,
    LINALG::Matrix<nen_, nen_>& ppmat, LINALG::Matrix<nen_, 1>& preforce)
{
  double vel_grad_porosity = 0.0;
  for (int idim = 0; idim < nsd_; ++idim)
    vel_grad_porosity += grad_porosity_(idim) * Base::velint_(idim);

  double grad_porosity_gridvelint = 0.0;
  for (int j = 0; j < nsd_; j++) grad_porosity_gridvelint += grad_porosity_(j) * gridvel_int_(j);

  if (!porofldpara_->PoroContiPartInt())
  {
    for (int idim = 0; idim < nsd_; ++idim)
    {
      const double grad_porosity_idim = grad_porosity_(idim);
      for (int ui = 0; ui < nen_; ++ui)
      {
        const int fui = nsd_ * ui;
        const double funct_ui = Base::funct_(ui);
        const double derxy_idim_ui = Base::derxy_(idim, ui);

        for (int vi = 0; vi < nen_; ++vi)
        {
          /* continuity term */
          /*
               /                      \
              |                        |
              | phi * nabla o Du  , q  |
              |                        |
               \                      /
          */
          /* porosity gradient term */
          /*
               /                   \
              |                     |
              | grad(phi)* Du  , q  |
              |                     |
               \                   /
          */
          estif_q_u(vi, fui + idim) += timefacfacpre * Base::funct_(vi) *
                                       (porosity_ * derxy_idim_ui + grad_porosity_idim * funct_ui);
        }
      }
    }

    // auxiliary variables
    static LINALG::Matrix<1, nen_> dgradphi_dp_gridvel;
    static LINALG::Matrix<1, nen_> dgradphi_dp_velint;
    dgradphi_dp_gridvel.MultiplyTN(gridvel_int_, dgradphi_dp);
    dgradphi_dp_velint.MultiplyTN(Base::velint_, dgradphi_dp);

    // pressure terms on left-hand side
    /* poroelasticity term */
    /*
         /                            \
        |                   n+1        |
        | d(grad(phi))/dp* u    Dp, q  |
        |                   (i)        |
         \                            /

         /                            \
        |                  n+1        |
     +  | d(phi)/dp * div u    Dp, q  |
        |                  (i)        |
         \                            /
    */

    for (int vi = 0; vi < nen_; ++vi)
    {
      const double v = timefacfacpre * Base::funct_(vi);

      for (int ui = 0; ui < nen_; ++ui)
      {
        ppmat(vi, ui) += v * (dphi_dp * Base::vdiv_ * Base::funct_(ui) + dgradphi_dp_velint(ui));
      }
    }

    // right-hand side
    const double rhsfac_vdiv = rhsfac * Base::vdiv_;
    for (int vi = 0; vi < nen_; ++vi)
    {
      // velocity term on right-hand side
      preforce(vi) -= rhsfac_vdiv * porosity_ * Base::funct_(vi) +
                      rhsfac * vel_grad_porosity * Base::funct_(vi);
    }

    // transient porosity terms
    /*
         /                             \      /                                             \
        |                   n+1         |    |                        /   n+1  \             |
      - | d(grad(phi))/dp* vs    Dp, q  |  + | d(phi)/(dJdp) * J *div| vs       |  * Dp , q  |
        |                   (i)         |    |                        \  (i)   /             |
         \                             /      \                                             /

         /                    \     /                                \
        |                      |   |                    n+1           |
      + | d(phi)/dp *  Dp , q  | + | d^2(phi)/(dp)^2 * p   *  Dp , q  |
        |                      |   |                    (i)           |
         \                    /     \                                /

    */

    if (!Base::fldparatimint_->IsStationary())
    {
      for (int vi = 0; vi < nen_; ++vi)
      {
        const double v = timefacfacpre * Base::funct_(vi);
        const double w = Base::fac_ * Base::funct_(vi);
        for (int ui = 0; ui < nen_; ++ui)
        {
          const double funct_ui = Base::funct_(ui);
          ppmat(vi, ui) += -v * dgradphi_dp_gridvel(ui) +
                           v * (dphi_dJdp * J_ * gridvel_div_) * funct_ui + w * funct_ui * dphi_dp +
                           w * dphi_dpp * funct_ui * press_;
        }
      }

      // coupling term on right hand side
      for (int vi = 0; vi < nen_; ++vi)
      {
        preforce(vi) -= rhsfac * Base::funct_(vi) * (-grad_porosity_gridvelint);
      }
    }
  }
  else  // Base::fldpara_->PoroContiPartInt() == true
  {
    for (int idim = 0; idim < nsd_; ++idim)
    {
      for (int ui = 0; ui < nen_; ++ui)
      {
        const int fui = nsd_ * ui;
        const double val = -1.0 * timefacfacpre * porosity_ * Base::funct_(ui);

        for (int vi = 0; vi < nen_; ++vi)
        {
          /* porosity convective term */
          /*
               /                   \
              |                     |
            - | phi * Du       , q  |
              |                     |
               \                   /
          */
          estif_q_u(vi, fui + idim) += val * Base::derxy_(idim, vi);
        }
      }
    }

    LINALG::Matrix<1, nen_> deriv_vel;
    deriv_vel.MultiplyTN(Base::velint_, Base::derxy_);
    // stationary right-hand side
    for (int vi = 0; vi < nen_; ++vi)
    {
      // velocity term on right-hand side
      preforce(vi) -= -1.0 * rhsfac * porosity_ * deriv_vel(vi);
    }

    // pressure terms on left-hand side
    /*
         /                                   \
        |                   n+1               |
        | -d(phi)/dp      * u    Dp, grad(q)  |
        |                   (i)               |
         \                                   /
    */

    for (int vi = 0; vi < nen_; ++vi)
    {
      const double factor = -timefacfacpre * dphi_dp * deriv_vel(vi);
      for (int ui = 0; ui < nen_; ++ui)
      {
        ppmat(vi, ui) += factor * Base::funct_(ui);
      }
    }

    if (!Base::fldparatimint_->IsStationary())
    {
      // transient porosity terms
      /*
          /                                             \
         |                        /   n+1  \             |
         | d(phi)/(dJdp) * J *div| vs       |  * Dp , q  |
         |                        \  (i)   /             |
          \                                             /

           /                    \       /                                \
          |                      |   |                    n+1           |
        + | d(phi)/dp *  Dp , q  | + | d^2(phi)/(dp)^2 * p   *  Dp , q  |
          |                      |   |                    (i)           |
           \                    /       \                                /

           /                            \
          |                  n+1        |
       +  | d(phi)/dp * div vs   Dp, q  |
          |                  (i)        |
           \                            /

      */

      LINALG::Matrix<1, nen_> deriv_gridvel;
      deriv_gridvel.MultiplyTN(gridvel_int_, Base::derxy_);

      for (int vi = 0; vi < nen_; ++vi)
      {
        const double u = timefacfacpre * dphi_dp * deriv_gridvel(vi);
        const double v =
            timefacfacpre * Base::funct_(vi) * ((dphi_dJdp * J_ + dphi_dp) * gridvel_div_);
        const double w = Base::fac_ * Base::funct_(vi);
        for (int ui = 0; ui < nen_; ++ui)
        {
          const double funct_ui = Base::funct_(ui);
          ppmat(vi, ui) += u * funct_ui + v * funct_ui + w * funct_ui * dphi_dp +
                           w * dphi_dpp * funct_ui * press_;
        }
      }

      // coupling term on right hand side
      for (int vi = 0; vi < nen_; ++vi)
      {
        preforce(vi) -= rhsfac * porosity_ * deriv_gridvel(vi);
        preforce(vi) -= rhsfac * Base::funct_(vi) * porosity_ * gridvel_div_;
      }
    }
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::GaussPointLoop(Teuchos::ParameterList& params,
    const LINALG::Matrix<nsd_, nen_>& ebofoaf, const LINALG::Matrix<nsd_, nen_>& evelaf,
    const LINALG::Matrix<nsd_, nen_>& evelnp, const LINALG::Matrix<nsd_, nen_>& eveln,
    const LINALG::Matrix<nen_, 1>& epreaf, const LINALG::Matrix<nen_, 1>& eprenp,
    const LINALG::Matrix<nen_, 1>& epren, const LINALG::Matrix<nen_, 1>& epressnp_timederiv,
    const LINALG::Matrix<nen_, 1>& epressam_timederiv,
    const LINALG::Matrix<nen_, 1>& epressn_timederiv, const LINALG::Matrix<nsd_, nen_>& eaccam,
    const LINALG::Matrix<nsd_, nen_>& edispnp, const LINALG::Matrix<nsd_, nen_>& edispn,
    const LINALG::Matrix<nsd_, nen_>& egridv, const LINALG::Matrix<nsd_, nen_>& egridvn,
    const LINALG::Matrix<nen_, 1>& escaaf, const LINALG::Matrix<nsd_, nen_>& emhist,
    const LINALG::Matrix<nen_, 1>& echist, const LINALG::Matrix<nen_, 1>* eporositynp,
    const LINALG::Matrix<nen_, 1>* eporositydot, const LINALG::Matrix<nen_, 1>* eporositydotn,
    LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_u, LINALG::Matrix<nen_ * nsd_, nen_>& estif_p_v,
    LINALG::Matrix<nen_, nen_ * nsd_>& estif_q_u, LINALG::Matrix<nen_, nen_>& ppmat,
    LINALG::Matrix<nen_, 1>& preforce, LINALG::Matrix<nsd_, nen_>& velforce,
    Teuchos::RCP<const MAT::Material> material, const DRT::UTILS::GaussIntegration& intpoints)
{
  // definition of velocity-based momentum residual vectors
  static LINALG::Matrix<nsd_ * nsd_, nen_> lin_resM_Du(true);
  static LINALG::Matrix<nsd_ * nsd_, nen_> lin_resMRea_Du(true);
  static LINALG::Matrix<nsd_, 1> resM_Du(true);
  static LINALG::Matrix<nsd_, nen_> lin_resM_Dp(true);

  // set element area or volume
  const double vol = Base::fac_;

  for (DRT::UTILS::GaussIntegration::const_iterator iquad = intpoints.begin();
       iquad != intpoints.end(); ++iquad)
  {
    // evaluate shape functions and derivatives at integration point
    Base::EvalShapeFuncAndDerivsAtIntPoint(iquad.Point(), iquad.Weight());

    SetupMaterialDerivatives();

    // -------------------------(material) deformation gradient F = d xyze_ / d XYZE = xyze_ *
    // N_XYZ_^T
    static LINALG::Matrix<nsd_, nsd_> defgrd(false);
    ComputeDefGradient(defgrd, N_XYZ_, Base::xyze_);

    // inverse deformation gradient F^-1
    static LINALG::Matrix<nsd_, nsd_> defgrd_inv(false);
    defgrd_inv.Invert(defgrd);

    // volume change (used for porosity law). Same as J in nonlinear theory.
    double volchange = 0.0;

    // compute J and the volume change
    ComputeJacobianDeterminantVolumeChange(J_, volchange, defgrd, N_XYZ_, edispnp);

    EvaluateVariablesAtGaussPoint(params, ebofoaf, evelaf, evelnp, eveln, epreaf, eprenp, epren,
        epressnp_timederiv, epressam_timederiv, epressn_timederiv, eaccam, edispnp, egridv, egridvn,
        escaaf, emhist, echist, eporositynp, eporositydot, eporositydotn);

    //-----------------------------------auxiliary variables for computing the porosity
    double dphi_dp = 0.0;
    double dphi_dJ = 0.0;
    double dphi_dJdp = 0.0;
    double dphi_dpp = 0.0;
    porosity_ = 0.0;

    // compute scalar at n+alpha_F or n+1
    Teuchos::RCP<std::vector<double>> scalars = Teuchos::rcp(new std::vector<double>(0));
    const double scalaraf = Base::funct_.Dot(escaaf);
    scalars->push_back(scalaraf);
    params.set<Teuchos::RCP<std::vector<double>>>("scalar", scalars);

    ComputePorosity(params, press_, volchange, *(iquad), Base::funct_, eporositynp, porosity_,
        &dphi_dp, &dphi_dJ, &dphi_dJdp,
        nullptr,  // dphi_dJJ not needed
        &dphi_dpp, false);

    //--linearization of porosity gradient w.r.t. pressure at gausspoint
    // d(grad(phi))/dp = dphi/(dJdp)* dJ/dx + d^2phi/(dp)^2 * dp/dx + dphi/dp* N,x
    static LINALG::Matrix<nsd_, nen_> dgradphi_dp(false);

    //--------------------------- dJ/dx
    static LINALG::Matrix<nsd_, 1> gradJ(false);

    // dF/dX
    static LINALG::Matrix<nsd_ * nsd_, nsd_> F_X(false);
    {
      //------------------------------------ build F^-1 as vector 9x1
      static LINALG::Matrix<nsd_ * nsd_, 1> defgrd_inv_vec;
      for (int i = 0; i < nsd_; i++)
        for (int j = 0; j < nsd_; j++) defgrd_inv_vec(i * nsd_ + j) = defgrd_inv(i, j);

      //------------------------------------ build F^-T as vector 9x1
      static LINALG::Matrix<nsd_ * nsd_, 1> defgrd_IT_vec;
      for (int i = 0; i < nsd_; i++)
        for (int j = 0; j < nsd_; j++) defgrd_IT_vec(i * nsd_ + j) = defgrd_inv(j, i);

      // dF/dx
      static LINALG::Matrix<nsd_ * nsd_, nsd_> F_x(false);

      ComputeFDerivative(edispnp, defgrd_inv, F_x, F_X);

      // compute gradients if needed
      ComputeGradients(J_, dphi_dp, dphi_dJ, defgrd_IT_vec, F_x, Base::gradp_, eporositynp, gradJ,
          grad_porosity_, refgrad_porosity_);
    }

    ComputeLinearization(dphi_dp, dphi_dpp, dphi_dJdp, gradJ, dgradphi_dp);

    //----------------------------------------------------------------------
    // potential evaluation of material parameters and/or stabilization
    // parameters at integration point
    //----------------------------------------------------------------------
    // get material parameters at integration point
    GetMaterialParamters(material);

    // set viscous term from previous iteration to zero (required for
    // using routine for evaluation of momentum rhs/residual as given)
    Base::visc_old_.Clear();
    Base::viscs2_.Clear();
    // compute viscous term from previous iteration and viscous operator
    if (Base::is_higher_order_ele_ and Base::visceff_) CalcDivEps(evelaf);

    // get reaction tensor and linearisations of material reaction tensor
    ComputeSpatialReactionTerms(material, defgrd_inv);

    // get stabilization parameters at integration point
    ComputeStabilizationParameters(vol);

    // compute old RHS of momentum equation and subgrid scale velocity
    ComputeOldRHSAndSubgridScaleVelocity();

    // compute old RHS of continuity equation
    ComputeOldRHSConti(dphi_dp);

    // compute strong residual of mixture (structural) equation
    if (porofldpara_->StabBiot() and (not porofldpara_->IsStationaryConti()) and
        struct_mat_->PoroLawType() != INPAR::MAT::m_poro_law_constant)
      ComputeMixtureStrongResidual(params, defgrd, edispnp, edispn, F_X, *iquad, false);

    //----------------------------------------------------------------------
    // set time-integration factors for left- and right-hand side
    //----------------------------------------------------------------------
    const double timefacfac = Base::fldparatimint_->TimeFac() * Base::fac_;
    const double timefacfacpre = Base::fldparatimint_->TimeFacPre() * Base::fac_;
    const double rhsfac = Base::fldparatimint_->TimeFacRhs() * Base::fac_;
    // const double rhsfacpre     = Base::fldparatimint_->TimeFacRhsPre() * Base::fac_;

    // set velocity-based momentum residual vectors to zero
    lin_resM_Du.Clear();
    lin_resMRea_Du.Clear();
    resM_Du.Clear();
    lin_resM_Dp.Clear();

    // compute first version of velocity-based momentum residual containing
    // inertia and reaction term
    ComputeLinResMDu(timefacfac, lin_resM_Du, lin_resMRea_Du);

    //----------------------------------------------------------------------
    // computation of standard Galerkin and stabilization contributions to
    // element matrix and right-hand-side vector
    //----------------------------------------------------------------------
    // 1) standard Galerkin inertia and reaction terms
    /* inertia (contribution to mass matrix) if not is_stationary */
    /*
            /              \
           |                |
           |    rho*Du , v  |
           |                |
            \              /
    */
    /*  reaction */
    /*
            /                \
           |                  |
           |    sigma*Du , v  |
           |                  |
            \                /
    */
    /* convection, convective ALE part  */
    /*
              /                             \
             |  /        n+1       \          |
             | |   rho*us   o nabla | Du , v  |
             |  \       (i)        /          |
              \                             /
    */

    for (int ui = 0; ui < nen_; ++ui)
    {
      const int fui = nsd_ * ui;

      for (int vi = 0; vi < nen_; ++vi)
      {
        const int fvi = nsd_ * vi;

        for (int idim = 0; idim < nsd_; ++idim)
        {
          for (int jdim = 0; jdim < nsd_; ++jdim)
            estif_u(fvi + idim, fui + jdim) +=
                Base::funct_(vi) * lin_resM_Du(idim * nsd_ + jdim, ui);
        }
      }
    }

    // inertia terms on the right hand side for instationary fluids
    if (not porofldpara_->IsStationaryMomentum())
    {
      for (int idim = 0; idim < nsd_; ++idim)
      {
        if (Base::fldparatimint_->IsGenalpha())
          resM_Du(idim) += rhsfac * Base::densam_ * Base::accint_(idim);
        else
          resM_Du(idim) += Base::fac_ * Base::densaf_ * Base::velint_(idim);
      }
    }

    if (not Base::fldparatimint_->IsStationary())
    {
      // coupling part RHS
      // reacoeff * phi * v_s
      for (int vi = 0; vi < nen_; ++vi)
      {
        for (int idim = 0; idim < nsd_; ++idim)
          velforce(idim, vi) -= -rhsfac * Base::funct_(vi) * reac_tensor_gridvel_(idim);
      }
    }

    // convective ALE-part
    for (int idim = 0; idim < nsd_; ++idim)
    {
      resM_Du(idim) += rhsfac * Base::densaf_ * Base::conv_old_(idim);
    }

    // reactive part
    for (int idim = 0; idim < nsd_; ++idim)
    {
      resM_Du(idim) += rhsfac * reac_tensor_vel_(idim);
    }

    for (int vi = 0; vi < nen_; ++vi)
    {
      for (int idim = 0; idim < nsd_; ++idim)
      {
        velforce(idim, vi) -= resM_Du(idim) * Base::funct_(vi);
      }
    }

    /* ********************************************************************** */
    /* Brinkman term: viscosity term */
    /*
                     /                        \
                    |       /  \         / \   |
              2 mu  |  eps | Du | , eps | v |  |
                    |       \  /         \ /   |
                     \                        /
    */

    if (Base::visceff_)
    {
      const double visceff_timefacfac = Base::visceff_ * timefacfac;
      const double porosity_inv = 1.0 / porosity_;

      for (int idim = 0; idim < nsd_; ++idim)
      {
        for (int jdim = 0; jdim < nsd_; ++jdim)
        {
          const double grad_porosity_jdim = grad_porosity_(jdim);

          for (int ui = 0; ui < nen_; ++ui)
          {
            const int fui = nsd_ * ui;
            const int fui_p_jdim = fui + jdim;

            const double derxy_idim_ui = Base::derxy_(idim, ui);
            const double derxy_jdim_ui = Base::derxy_(jdim, ui);

            for (int vi = 0; vi < nen_; ++vi)
            {
              const double viscfac1 = visceff_timefacfac * Base::funct_(vi) * porosity_inv;
              const double viscfac2 = visceff_timefacfac * Base::derxy_(jdim, vi);
              const int fvi_p_idim = nsd_ * vi + idim;

              estif_u(fvi_p_idim, fui_p_jdim) +=
                  viscfac2 * derxy_idim_ui - viscfac1 * derxy_idim_ui * grad_porosity_jdim;
              estif_u(fvi_p_idim, fui + idim) +=
                  viscfac2 * derxy_jdim_ui - viscfac1 * derxy_jdim_ui * grad_porosity_jdim;
            }
          }
        }
      }

      static LINALG::Matrix<nsd_, nsd_> viscstress(false);
      for (int jdim = 0; jdim < nsd_; ++jdim)
      {
        for (int idim = 0; idim < nsd_; ++idim)
        {
          viscstress(idim, jdim) =
              Base::visceff_ * (Base::vderxy_(jdim, idim) + Base::vderxy_(idim, jdim));
        }
      }

      static LINALG::Matrix<nsd_, 1> viscstress_gradphi(false);
      viscstress_gradphi.Multiply(viscstress, grad_porosity_);

      // computation of right-hand-side viscosity term
      for (int idim = 0; idim < nsd_; ++idim)
      {
        const double viscstress_gradphi_idim = viscstress_gradphi(idim);
        for (int vi = 0; vi < nen_; ++vi)
        {
          /* viscosity (brinkman) term on right-hand side */
          velforce(idim, vi) -= -rhsfac * porosity_inv * viscstress_gradphi_idim * Base::funct_(vi);
          for (int jdim = 0; jdim < nsd_; ++jdim)
          {
            /* viscosity term on right-hand side */
            velforce(idim, vi) -= rhsfac * viscstress(idim, jdim) * Base::derxy_(jdim, vi);
          }
        }
      }

      static LINALG::Matrix<nsd_, nen_> viscstress_dgradphidp(false);
      viscstress_dgradphidp.Multiply(viscstress, dgradphi_dp);
      for (int idim = 0; idim < nsd_; ++idim)
      {
        const double viscstress_gradphi_idim = viscstress_gradphi(idim);
        for (int ui = 0; ui < nen_; ++ui)
        {
          lin_resM_Dp(idim, ui) +=
              timefacfacpre * porosity_inv *
              (porosity_inv * viscstress_gradphi_idim * dphi_dp * Base::funct_(ui) -
                  viscstress_dgradphidp(idim, ui));
        }
      }
    }

    /* ********************************************************************** */
    // 3) standard Galerkin pressure term + poroelasticity terms
    /* pressure term */
    /*
         /                \
        |                  |
        |  Dp , nabla o v  |
        |                  |
         \                /
    */
    /* poroelasticity pressure term */
    /*
         /                           \      /                            \
        |         n+1                 |     |         n+1                 |
        |  sigma*u  * dphi/dp*Dp , v  |  -  |  sigma*vs  * dphi/dp*Dp , v |
        |         (i)                 |     |         (i)                 |
         \                           /       \                           /
    */

    for (int ui = 0; ui < nen_; ++ui)
    {
      const double v = timefacfacpre * Base::funct_(ui);
      for (int vi = 0; vi < nen_; ++vi)
      {
        const int fvi = nsd_ * vi;
        for (int idim = 0; idim < nsd_; ++idim)
        {
          estif_p_v(fvi + idim, ui) += v * (-1.0 * Base::derxy_(idim, vi));
        }
      }
    }

    ComputeLinResMDp(timefacfacpre, dphi_dp, lin_resM_Dp);

    for (int vi = 0; vi < nen_; ++vi)
    {
      const int fvi = nsd_ * vi;
      const double funct_vi = Base::funct_(vi);
      for (int ui = 0; ui < nen_; ++ui)
      {
        for (int idim = 0; idim < nsd_; ++idim)
        {
          estif_p_v(fvi + idim, ui) += funct_vi * lin_resM_Dp(idim, ui);
        }
      }
    }

    const double pressfac = press_ * rhsfac;

    for (int vi = 0; vi < nen_; ++vi)
    {
      /* pressure term on right-hand side */
      for (int idim = 0; idim < nsd_; ++idim)
      {
        velforce(idim, vi) += pressfac * (Base::derxy_(idim, vi));
      }
    }

    /* ********************************************************************** */
    // 4) standard Galerkin continuity term + poroelasticity terms

    EvaluatePressureEquation(params, timefacfacpre, rhsfac, dphi_dp, dphi_dJ, dphi_dJdp, dphi_dpp,
        eporositydot, eporositydotn, echist, dgradphi_dp, estif_q_u, ppmat, preforce);

    /* *********************************************************************************************************
     */

    // 5) standard Galerkin bodyforce term on right-hand side
    Base::BodyForceRhsTerm(velforce, rhsfac);

    if (Base::fldpara_->PSPG() or Base::fldpara_->RStab() != INPAR::FLUID::reactive_stab_none)
      ComputeLinResMDuStabilization(timefacfac, lin_resM_Du);

    // 6) PSPG term
    if (Base::fldpara_->PSPG())
    {
      PSPG(estif_q_u, ppmat, preforce, lin_resM_Du, lin_resMRea_Du, lin_resM_Dp, dphi_dp, 0.0,
          timefacfac, timefacfacpre, rhsfac);
    }

    // 7) reactive stabilization term
    if (Base::fldpara_->RStab() != INPAR::FLUID::reactive_stab_none)
    {
      ReacStab(estif_u, estif_p_v, velforce, lin_resM_Du, lin_resM_Dp, dphi_dp, timefacfac,
          timefacfacpre, rhsfac, 0.0);
    }

    // 8) Biot stabilization term
    if (porofldpara_->StabBiot() and (not porofldpara_->IsStationaryConti()) and
        struct_mat_->PoroLawType() != INPAR::MAT::m_poro_law_constant)
    {
      StabBiot(estif_q_u, ppmat, preforce, lin_resM_Du, lin_resMRea_Du, lin_resM_Dp, dphi_dp, 0.0,
          timefacfac, timefacfacpre, rhsfac);
    }

    /* ********************************************************************** */
    // 9) stabilization of continuity equation
    if (Base::fldpara_->CStab())
    {
      dserror("continuity stabilization not implemented for poroelasticity");

      // In the case no continuity stabilization and no LOMA:
      // the factors 'conti_stab_and_vol_visc_fac' and 'conti_stab_and_vol_visc_rhs' are zero
      // therefore there is no contribution to the element stiffness matrix and
      // the viscous stress tensor is NOT altered!!
      //
      // ONLY
      // the rhs contribution of the viscous term is added!!

      double conti_stab_and_vol_visc_fac = 0.0;
      double conti_stab_and_vol_visc_rhs = 0.0;

      if (Base::fldpara_->CStab())
      {
        conti_stab_and_vol_visc_fac += timefacfacpre * Base::tau_(2);
        conti_stab_and_vol_visc_rhs -= rhsfac * Base::tau_(2) * Base::conres_old_;
      }

      /* continuity stabilisation on left hand side */
      /*
                  /                        \
                 |                          |
            tauC | nabla o Du  , nabla o v  |
                 |                          |
                  \                        /
      */
      /* viscosity term - subtraction for low-Mach-number flow */
      /*
                 /                             \             /                        \
                |  1                      / \   |     2 mu  |                          |
         - 2 mu |  - (nabla o u) I , eps | v |  | = - ----- | nabla o Du  , nabla o v  |
                |  3                      \ /   |       3   |                          |
                 \                             /             \                        /
      */
      for (int ui = 0; ui < nen_; ++ui)
      {
        const int fui = nsd_ * ui;

        for (int idim = 0; idim < nsd_; ++idim)
        {
          const int fui_p_idim = fui + idim;
          const double v0 = conti_stab_and_vol_visc_fac * Base::derxy_(idim, ui);
          for (int vi = 0; vi < nen_; ++vi)
          {
            const int fvi = nsd_ * vi;

            for (int jdim = 0; jdim < nsd_; ++jdim)
            {
              estif_u(fvi + jdim, fui_p_idim) += v0 * Base::derxy_(jdim, vi);
            }
          }
        }
      }

      // computation of right-hand-side viscosity term
      for (int vi = 0; vi < nen_; ++vi)
      {
        for (int idim = 0; idim < nsd_; ++idim)
        {
          /* viscosity term on right-hand side */
          velforce(idim, vi) += conti_stab_and_vol_visc_rhs * Base::derxy_(idim, vi);
        }
      }
    }
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::GaussPointLoopOD(Teuchos::ParameterList& params,
    const LINALG::Matrix<nsd_, nen_>& ebofoaf, const LINALG::Matrix<nsd_, nen_>& evelaf,
    const LINALG::Matrix<nsd_, nen_>& evelnp, const LINALG::Matrix<nsd_, nen_>& eveln,
    const LINALG::Matrix<nen_, 1>& epreaf, const LINALG::Matrix<nen_, 1>& eprenp,
    const LINALG::Matrix<nen_, 1>& epren, const LINALG::Matrix<nen_, 1>& epressnp_timederiv,
    const LINALG::Matrix<nen_, 1>& epressam_timederiv,
    const LINALG::Matrix<nen_, 1>& epressn_timederiv, const LINALG::Matrix<nsd_, nen_>& eaccam,
    const LINALG::Matrix<nsd_, nen_>& edispnp, const LINALG::Matrix<nsd_, nen_>& edispn,
    const LINALG::Matrix<nsd_, nen_>& egridv, const LINALG::Matrix<nsd_, nen_>& egridvn,
    const LINALG::Matrix<nen_, 1>& escaaf, const LINALG::Matrix<nsd_, nen_>& emhist,
    const LINALG::Matrix<nen_, 1>& echist, const LINALG::Matrix<nen_, 1>* eporositynp,
    LINALG::Matrix<(nsd_ + 1) * nen_, 1>& eforce,
    LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& ecoupl_u, LINALG::Matrix<nen_, nen_ * nsd_>& ecoupl_p,
    Teuchos::RCP<const MAT::Material> material, const DRT::UTILS::GaussIntegration& intpoints)
{
  // definition of velocity-based momentum residual vectors
  static LINALG::Matrix<nsd_, nen_ * nsd_> lin_resM_Dus(true);
  static LINALG::Matrix<nsd_, nen_ * nsd_> lin_resM_Dus_gridvel(true);

  // set element area or volume
  const double vol = Base::fac_;

  for (DRT::UTILS::GaussIntegration::const_iterator iquad = intpoints.begin();
       iquad != intpoints.end(); ++iquad)
  {
    // reset matrix
    lin_resM_Dus.Clear();
    lin_resM_Dus_gridvel.Clear();

    // evaluate shape functions and derivatives at integration point
    Base::EvalShapeFuncAndDerivsAtIntPoint(iquad.Point(), iquad.Weight());

    // evaluate shape function derivatives w.r.t. to material coordinates at integration point
    SetupMaterialDerivatives();

    // -------------------------(material) deformation gradient F = d xyze_ / d XYZE = xyze_ *
    // N_XYZ_^T
    static LINALG::Matrix<nsd_, nsd_> defgrd(false);
    ComputeDefGradient(defgrd, N_XYZ_, Base::xyze_);

    // inverse deformation gradient F^-1
    static LINALG::Matrix<nsd_, nsd_> defgrd_inv(false);
    defgrd_inv.Invert(defgrd);

    // volume change (used for porosity law). Same as J in nonlinear theory.
    double volchange = 0.0;

    // compute J and the volume change
    ComputeJacobianDeterminantVolumeChange(J_, volchange, defgrd, N_XYZ_, edispnp);

    EvaluateVariablesAtGaussPointOD(params, ebofoaf, evelaf, evelnp, eveln, epreaf, eprenp, epren,
        epressnp_timederiv, epressam_timederiv, epressn_timederiv, eaccam, edispnp, edispn, egridv,
        egridvn, escaaf, emhist, echist, eporositynp);

    //************************************************auxilary variables for computing the porosity_

    // compute scalar at n+alpha_F or n+1
    Teuchos::RCP<std::vector<double>> scalars = Teuchos::rcp(new std::vector<double>(0));
    const double scalaraf = Base::funct_.Dot(escaaf);
    scalars->push_back(scalaraf);
    params.set<Teuchos::RCP<std::vector<double>>>("scalar", scalars);

    double dphi_dp = 0.0;
    double dphi_dJ = 0.0;
    double dphi_dJdp = 0.0;
    double dphi_dJJ = 0.0;
    porosity_ = 0.0;

    ComputePorosity(params, press_, volchange, *(iquad), Base::funct_, eporositynp, porosity_,
        &dphi_dp, &dphi_dJ, &dphi_dJdp, &dphi_dJJ,
        nullptr,  // dphi_dpp not needed
        false);

    double refporositydot = struct_mat_->RefPorosityTimeDeriv();

    //---------------------------  dJ/dx = dJ/dF : dF/dx = JF^-T : dF/dx at gausspoint
    static LINALG::Matrix<nsd_, 1> gradJ(false);
    // spatial porosity gradient
    // LINALG::Matrix<nsd_,1>             grad_porosity(true);
    //--------------------- linearization of porosity w.r.t. structure displacements
    static LINALG::Matrix<1, nsd_ * nen_> dphi_dus(false);

    //------------------------------------------------dJ/dus = dJ/dF : dF/dus = J * F^-T . N_X = J *
    // N_x
    static LINALG::Matrix<1, nsd_ * nen_> dJ_dus(false);
    //------------------ d( grad(\phi) ) / du_s = d\phi/(dJ du_s) * dJ/dx+ d\phi/dJ * dJ/(dx*du_s) +
    // d\phi/(dp*du_s) * dp/dx
    static LINALG::Matrix<nsd_, nen_ * nsd_> dgradphi_dus(false);

    //------------------------------------ build F^-T as vector 9x1
    static LINALG::Matrix<nsd_ * nsd_, 1> defgrd_IT_vec;
    for (int i = 0; i < nsd_; i++)
      for (int j = 0; j < nsd_; j++) defgrd_IT_vec(i * nsd_ + j) = defgrd_inv(j, i);

    // dF/dx
    static LINALG::Matrix<nsd_ * nsd_, nsd_> F_x(false);

    // dF/dX
    static LINALG::Matrix<nsd_ * nsd_, nsd_> F_X(false);

    ComputeFDerivative(edispnp, defgrd_inv, F_x, F_X);

    // compute gradients if needed
    ComputeGradients(J_, dphi_dp, dphi_dJ, defgrd_IT_vec, F_x, Base::gradp_, eporositynp, gradJ,
        grad_porosity_, refgrad_porosity_);

    ComputeLinearizationOD(dphi_dJ, dphi_dJJ, dphi_dJdp, defgrd_inv, defgrd_IT_vec, F_x, F_X, gradJ,
        dJ_dus, dphi_dus, dgradphi_dus);

    //----------------------------------------------------------------------
    // potential evaluation of material parameters and/or stabilization
    // parameters at integration point
    //----------------------------------------------------------------------
    // get material parameters at integration point
    GetMaterialParamters(material);

    ComputeSpatialReactionTerms(material, defgrd_inv);

    // compute linearization of spatial reaction tensor w.r.t. structural displacements
    ComputeLinSpatialReactionTerms(material, defgrd_inv, &dJ_dus, &dphi_dus);

    // set viscous term from previous iteration to zero (required for
    // using routine for evaluation of momentum rhs/residual as given)
    Base::visc_old_.Clear();
    Base::viscs2_.Clear();
    // compute viscous term from previous iteration and viscous operator
    if (Base::is_higher_order_ele_ and Base::visceff_) CalcDivEps(evelaf);

    // get stabilization parameters at integration point
    ComputeStabilizationParameters(vol);

    // compute old RHS of momentum equation and subgrid scale velocity
    ComputeOldRHSAndSubgridScaleVelocity();

    // compute old RHS of continuity equation
    ComputeOldRHSConti(dphi_dp);

    // compute strong residual of mixture (structural) equation
    if (porofldpara_->StabBiot() and (not porofldpara_->IsStationaryConti()) and
        struct_mat_->PoroLawType() != INPAR::MAT::m_poro_law_constant)
      ComputeMixtureStrongResidual(params, defgrd, edispnp, edispn, F_X, *iquad, true);

    //----------------------------------------------------------------------
    // set time-integration factors for left- and right-hand side
    //----------------------------------------------------------------------

    const double timefacfac = Base::fldparatimint_->TimeFac() * Base::fac_;
    const double timefacfacpre = Base::fldparatimint_->TimeFacPre() * Base::fac_;

    //**************************************************************************
    // 1) coupling terms in momentum balance

    FillMatrixMomentumOD(timefacfac, evelaf, egridv, epreaf, dgradphi_dus, dphi_dp, dphi_dJ,
        dphi_dus, refporositydot, lin_resM_Dus, lin_resM_Dus_gridvel, ecoupl_u);

    //**************************************************************************
    // 2) coupling terms in continuity equation

    FillMatrixContiOD(timefacfacpre, dphi_dp, dphi_dJ, dphi_dJJ, dphi_dJdp, refporositydot,
        dgradphi_dus, dphi_dus, dJ_dus, egridv, lin_resM_Dus, lin_resM_Dus_gridvel, ecoupl_p);
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::FillMatrixMomentumOD(const double& timefacfac,
    const LINALG::Matrix<nsd_, nen_>& evelaf, const LINALG::Matrix<nsd_, nen_>& egridv,
    const LINALG::Matrix<nen_, 1>& epreaf, const LINALG::Matrix<nsd_, nen_ * nsd_>& dgradphi_dus,
    const double& dphi_dp, const double& dphi_dJ, const LINALG::Matrix<1, nsd_ * nen_>& dphi_dus,
    const double& refporositydot, LINALG::Matrix<nsd_, nen_ * nsd_>& lin_resM_Dus,
    LINALG::Matrix<nsd_, nen_ * nsd_>& lin_resM_Dus_gridvel,
    LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& ecoupl_u)
{
  // stationary
  /*  reaction */
  /*
   /                                      \
   |                    n+1                |
   |    sigma * dphi/dus * u    * Dus , v  |
   |                        (i)            |
   \                                     /
   */

  for (int idim = 0; idim < nsd_; ++idim)
  {
    for (int jdim = 0; jdim < nsd_; ++jdim)
    {
      for (int ui = 0; ui < nen_; ++ui)
      {
        lin_resM_Dus(idim, nsd_ * ui + jdim) +=
            +timefacfac * reac_tensor_linOD_vel_(idim, nsd_ * ui + jdim);
      }
    }
  }

  // transient terms
  /*  reaction */
  /*
    /                           \        /                                           \
   |                             |      |                            n+1              |
-  |    sigma * phi * D(v_s) , v |  -   |    sigma * d(phi)/d(us) * vs *  D(u_s) , v  |
   |                             |      |                             (i)             |
    \                           /        \                                           /
   */
  /*  reactive ALE term */
  /*
   /                                  \
   |                  n+1             |
   |    - rho * grad u     * Dus , v  |
   |                  (i)             |
   \                                 /
   */

  if (not Base::fldparatimint_->IsStationary())
  {
    const double fac_densaf = Base::fac_ * Base::densaf_;
    for (int idim = 0; idim < nsd_; ++idim)
    {
      for (int jdim = 0; jdim < nsd_; ++jdim)
      {
        const double vderxy_idim_jdim = Base::vderxy_(idim, jdim);
        const double reatensor_idim_jdim = reac_tensor_(idim, jdim);
        for (int ui = 0; ui < nen_; ++ui)
        {
          lin_resM_Dus(idim, nsd_ * ui + jdim) +=
              Base::fac_ * (-1.0) * reatensor_idim_jdim * Base::funct_(ui) -
              timefacfac * reac_tensor_linOD_grid_vel_(idim, nsd_ * ui + jdim);
          lin_resM_Dus_gridvel(idim, nsd_ * ui + jdim) +=
              -fac_densaf * vderxy_idim_jdim * Base::funct_(ui);
        }
      }
    }
  }

  // viscous terms (brinkman terms)
  if (Base::visceff_)
  {
    static LINALG::Matrix<nsd_, nsd_> viscstress(false);
    for (int jdim = 0; jdim < nsd_; ++jdim)
    {
      for (int idim = 0; idim < nsd_; ++idim)
      {
        viscstress(idim, jdim) =
            Base::visceff_ * (Base::vderxy_(jdim, idim) + Base::vderxy_(idim, jdim));
      }
    }

    static LINALG::Matrix<nsd_, 1> viscstress_gradphi(false);
    viscstress_gradphi.Multiply(viscstress, grad_porosity_);

    static LINALG::Matrix<nsd_, nen_ * nsd_> viscstress_dgradphidus(false);
    viscstress_dgradphidus.Multiply(viscstress, dgradphi_dus);

    const double porosity_inv = 1.0 / porosity_;

    for (int idim = 0; idim < nsd_; ++idim)
    {
      const double viscstress_gradphi_idim = viscstress_gradphi(idim);
      for (int ui = 0; ui < nen_; ++ui)
      {
        const int fui = nsd_ * ui;
        for (int jdim = 0; jdim < nsd_; ++jdim)
        {
          lin_resM_Dus(idim, fui + jdim) +=
              timefacfac * porosity_inv *
              (porosity_inv * viscstress_gradphi_idim * dphi_dus(fui + jdim) -
                  viscstress_dgradphidus(idim, fui + jdim));
        }
      }
    }
  }

  for (int ui = 0; ui < nen_; ++ui)
  {
    const int fui = nsd_ * ui;

    for (int vi = 0; vi < nen_; ++vi)
    {
      const int fvi = nsd_ * vi;
      const double funct_vi = Base::funct_(vi);

      for (int idim = 0; idim < nsd_; ++idim)
      {
        for (int jdim = 0; jdim < nsd_; ++jdim)
        {
          ecoupl_u(fvi + idim, fui + jdim) +=
              funct_vi * (lin_resM_Dus(idim, fui + jdim) + lin_resM_Dus_gridvel(idim, fui + jdim));
        }
      }
    }
  }

  //**************************************************************************
  if (Base::fldpara_->RStab() != INPAR::FLUID::reactive_stab_none)
  {
    if (Base::fldpara_->Tds() != INPAR::FLUID::subscales_quasistatic)
      dserror("Is this factor correct? Check for bugs!");
    const double reac_tau = Base::fldpara_->ViscReaStabFac() * Base::reacoeff_ * Base::tau_(1);

    for (int vi = 0; vi < nen_; ++vi)
    {
      const double v = reac_tau * Base::funct_(vi);

      for (int idim = 0; idim < nsd_; ++idim)
      {
        const int fvi_p_idim = nsd_ * vi + idim;

        for (int jdim = 0; jdim < nsd_; ++jdim)
        {
          for (int ui = 0; ui < nen_; ++ui)
          {
            const int fui_p_jdim = nsd_ * ui + jdim;

            ecoupl_u(fvi_p_idim, fui_p_jdim) +=
                v * (lin_resM_Dus(idim, fui_p_jdim) + lin_resM_Dus_gridvel(idim, fui_p_jdim));
          }
        }
      }
    }

    {  // linearization of stabilization parameter w.r.t. structure displacement
      const double v =
          timefacfac * Base::fldpara_->ViscReaStabFac() *
          (Base::reacoeff_ * dtau_dphi_(1) / Base::tau_(1) + Base::reacoeff_ / porosity_);
      for (int vi = 0; vi < nen_; ++vi)
      {
        const double w = -1.0 * v * Base::funct_(vi);

        for (int idim = 0; idim < nsd_; ++idim)
        {
          const double w_sgvelint = w * Base::sgvelint_(idim);
          const int fvi_p_idim = nsd_ * vi + idim;

          for (int ui = 0; ui < nen_; ++ui)
          {
            for (int jdim = 0; jdim < nsd_; ++jdim)
            {
              const int fui_p_jdim = nsd_ * ui + jdim;

              ecoupl_u(fvi_p_idim, fui_p_jdim) += w_sgvelint * dphi_dus(fui_p_jdim);
            }
          }
        }
      }
    }
  }

  //**************************************************************************
  // shape derivatives
  if (nsd_ == 3)
  {
    Lin3DMeshMotionOD(
        ecoupl_u, dphi_dp, dphi_dJ, refporositydot, Base::fldparatimint_->TimeFac(), timefacfac);
  }
  else if (nsd_ == 2)
  {
    Lin2DMeshMotionOD(
        ecoupl_u, dphi_dp, dphi_dJ, refporositydot, Base::fldparatimint_->TimeFac(), timefacfac);
  }
  else
  {
    dserror("Linearization of the mesh motion is only available in 2D and 3D");
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::FillMatrixContiOD(const double& timefacfacpre,
    const double& dphi_dp, const double& dphi_dJ, const double& dphi_dJJ, const double& dphi_dJdp,
    const double& refporositydot, const LINALG::Matrix<nsd_, nen_ * nsd_>& dgradphi_dus,
    const LINALG::Matrix<1, nsd_ * nen_>& dphi_dus, const LINALG::Matrix<1, nsd_ * nen_>& dJ_dus,
    const LINALG::Matrix<nsd_, nen_>& egridv, const LINALG::Matrix<nsd_, nen_ * nsd_>& lin_resM_Dus,
    const LINALG::Matrix<nsd_, nen_ * nsd_>& lin_resM_Dus_gridvel,
    LINALG::Matrix<nen_, nen_ * nsd_>& ecoupl_p)
{
  if (!porofldpara_->IsStationaryConti())
  {
    for (int vi = 0; vi < nen_; ++vi)
    {
      const double w = timefacfacpre * Base::funct_(vi);
      for (int ui = 0; ui < nen_; ++ui)
      {
        const int fui = nsd_ * ui;
        for (int idim = 0; idim < nsd_; ++idim)
        {
          ecoupl_p(vi, fui + idim) += w * dphi_dJdp * (-Base::rhscon_) * dJ_dus(fui + idim);
        }
      }
    }
  }

  if (!porofldpara_->PoroContiPartInt())
  {
    // auxiliary variables
    static LINALG::Matrix<1, nen_ * nsd_> grad_porosity_us_velint;
    grad_porosity_us_velint.MultiplyTN(Base::velint_, dgradphi_dus);

    // structure coupling terms on left-hand side
    /*  stationary */
    /*
      /                                 \      /                                    \
     |                 n+1              |     |                        n+1           |
     |   dphi/dus * div u    * Dus , v  |  +  |   d(grad(phi))/dus * u    * Dus , v  |
     |                  (i)             |     |                       (i)            |
      \                                /       \                                    /
     */
    for (int vi = 0; vi < nen_; ++vi)
    {
      const double v = timefacfacpre * Base::funct_(vi);

      for (int ui = 0; ui < nen_; ++ui)
      {
        const int fui = nsd_ * ui;
        for (int idim = 0; idim < nsd_; ++idim)
        {
          ecoupl_p(vi, fui + idim) +=
              v * dphi_dus(fui + idim) * Base::vdiv_ + v * grad_porosity_us_velint(fui + idim);
        }
      }
    }

    // transient coupling terms
    if (!Base::fldparatimint_->IsStationary())
    {
      static LINALG::Matrix<1, nen_ * nsd_> grad_porosity_us_gridvelint;
      grad_porosity_us_gridvelint.MultiplyTN(gridvel_int_, dgradphi_dus);

      /*
        /                            \       / \ |                              |     | n+1 | |
    dphi/dJ * J * div Dus , v  |   + |   d^2(phi)/(dJ)^2 * dJ/dus  * J * div vs    * Dus , v  | | |
    |                                        (i)             | \                            /      \
    /

       /                                            \        / \ |                           n+1 |
    |                           n+1            |
    +  |   dphi/dJ * dJ/dus * div vs    * Dus, v     |    - |   d(grad(phi))/d(us) *  vs    * Dus ,
    v  | |                           (i)               |      |                           (i) | \ /
    \                                        /

          /                       \
         |                         |
       - |    grad phi * Dus , v   |
         |                         |
          \                       /

          /                                          \
         |                             n+1            |
       + |    dphi/(dpdJ) * dJ/dus  * p    * Dus , v  |
         |                             (i)            |
         \                                           /
       */

      for (int vi = 0; vi < nen_; ++vi)
      {
        const double v = Base::fac_ * Base::funct_(vi);
        const double w = timefacfacpre * Base::funct_(vi);
        for (int ui = 0; ui < nen_; ++ui)
        {
          const int fui = nsd_ * ui;
          const double funct_ui = Base::funct_(ui);
          for (int idim = 0; idim < nsd_; ++idim)
          {
            ecoupl_p(vi, fui + idim) +=
                v * (dphi_dJ * J_ * Base::derxy_(idim, ui)) +
                w * (+gridvel_div_ * (dphi_dJJ * J_ + dphi_dJ) * dJ_dus(fui + idim) -
                        grad_porosity_us_gridvelint(fui + idim)) -
                v * grad_porosity_(idim) * funct_ui + v * dphi_dJdp * press_ * dJ_dus(fui + idim);
          }
        }
      }
    }
  }
  else
  {
    static LINALG::Matrix<1, nen_> deriv_vel;
    deriv_vel.MultiplyTN(Base::velint_, Base::derxy_);

    // structure coupling terms on left-hand side
    /*  stationary */
    /*
      /                                    \
     |                 n+1                  |
     |   -dphi/dus *  u    * Dus , grad(v)  |
     |                  (i)                 |
      \                                    /
     */
    for (int vi = 0; vi < nen_; ++vi)
    {
      const double val_deriv_vel_vi = timefacfacpre * (-1.0) * deriv_vel(vi);
      for (int ui = 0; ui < nen_; ++ui)
      {
        const int fui = nsd_ * ui;
        for (int idim = 0; idim < nsd_; ++idim)
        {
          ecoupl_p(vi, fui + idim) += dphi_dus(fui + idim) * val_deriv_vel_vi;
        }
      }
    }

    // transient coupling terms
    if (!Base::fldparatimint_->IsStationary())
    {
      /*
        /                                    \       / \ |                                      | |
    n+1                 | |   (dphi/dJ * J + phi )* div Dus , v  |   + |   d^2(phi)/(dJ)^2 * dJ/dus
    * J * div vs    * Dus , v  | |                                      |     | (i)             | \
    /      \                                                       /

       /                                            \
       |                           n+1               |
    +  |   dphi/dJ * dJ/dus * div vs    * Dus, v     |
       |                           (i)               |
       \                                            /

       /                                   \    /                        \
       |                    n+1            |    |                        |
    +  |   dphi/dus * div vs  * Dus, v     |  + |   phi *  Dus, grad(v)  |
       |                    (i)            |    |                        |
       \                                  /     \                       /

          /                                          \
         |                             n+1            |
       + |    dphi/(dpdJ) * dJ/dus  * p    * Dus , v  |
         |                             (i)            |
         \                                           /
       */

      static LINALG::Matrix<1, nen_> deriv_gridvel;
      deriv_gridvel.MultiplyTN(gridvel_int_, Base::derxy_);

      for (int vi = 0; vi < nen_; ++vi)
      {
        const double v = Base::fac_ * Base::funct_(vi);
        const double w = timefacfacpre * Base::funct_(vi);
        const double deriv_gridvel_vi = deriv_gridvel(vi);
        for (int ui = 0; ui < nen_; ++ui)
        {
          const int fui = nsd_ * ui;
          const double funct_ui = Base::funct_(ui);
          for (int idim = 0; idim < nsd_; ++idim)
          {
            ecoupl_p(vi, fui + idim) +=
                v * ((dphi_dJ * J_ + porosity_) * Base::derxy_(idim, ui)) +
                Base::fac_ * Base::derxy_(idim, vi) * (porosity_ * funct_ui) +
                w * (+gridvel_div_ *
                        ((dphi_dJJ * J_ + dphi_dJ) * dJ_dus(fui + idim) + dphi_dus(fui + idim))) +
                timefacfacpre * deriv_gridvel_vi * dphi_dus(fui + idim) +
                v * dphi_dJdp * press_ * dJ_dus(fui + idim);
          }
        }
      }
    }
  }

  //**************************************************************************
  // PSPG
  if (Base::fldpara_->PSPG())
  {
    double scal_grad_q = 0.0;

    if (Base::fldpara_->Tds() == INPAR::FLUID::subscales_quasistatic)
    {
      scal_grad_q = Base::tau_(1);
    }
    else
    {
      // this is not tested anyway ...
      scal_grad_q = 0.0;
    }

    {
      // linearization of stabilization parameter w.r.t. structure displacement
      const double v1 = -1.0 * timefacfacpre * dtau_dphi_(1) / scal_grad_q;
      static LINALG::Matrix<1, nen_> temp(false);
      temp.MultiplyTN(Base::sgvelint_, Base::derxy_);

      for (int jdim = 0; jdim < nsd_; ++jdim)
      {
        for (int ui = 0; ui < nen_; ++ui)
        {
          const int fui_p_jdim = nsd_ * ui + jdim;
          const double dphi_dus_fui_p_jdim = dphi_dus(fui_p_jdim);
          for (int vi = 0; vi < nen_; ++vi)
          {
            ecoupl_p(vi, fui_p_jdim) += v1 * temp(vi) * dphi_dus_fui_p_jdim;
          }
        }
      }
    }

    for (int vi = 0; vi < nen_; ++vi)
    {
      for (int idim = 0; idim < nsd_; ++idim)
      {
        const double scal_derxy_idim_vi = scal_grad_q * Base::derxy_(idim, vi);
        for (int ui = 0; ui < nen_; ++ui)
        {
          for (int jdim = 0; jdim < nsd_; ++jdim)
          {
            ecoupl_p(vi, ui * nsd_ + jdim) +=
                scal_derxy_idim_vi * (lin_resM_Dus(idim, ui * nsd_ + jdim) +
                                         lin_resM_Dus_gridvel(idim, ui * nsd_ + jdim));
          }
        }
      }
    }
  }

  //**************************************************************************
  // biot stabilization
  if (porofldpara_->StabBiot() and (not porofldpara_->IsStationaryConti()) and
      struct_mat_->PoroLawType() != INPAR::MAT::m_poro_law_constant)
  {
    const double val = tau_struct_ * porosity_;
    double fac_dens = 0.0;
    double fac_dens2 = 0.0;

    const double timefactor = Base::fldparatimint_->Dt() * Base::fldparatimint_->Theta();
    fac_dens = tau_struct_ * Base::fac_ / J_ * struct_mat_->Density() / timefactor;
    fac_dens2 =
        -tau_struct_ * Base::fac_ / (J_ * J_) * struct_mat_->Density() / Base::fldparatimint_->Dt();

    for (int ui = 0; ui < nen_; ++ui)
    {
      const double funct_ui = Base::funct_(ui);
      for (int vi = 0; vi < nen_; ++vi)
      {
        for (int jdim = 0; jdim < nsd_; ++jdim)
        {
          double stiff_val = 0.0;
          double stiff_val2 = 0.0;
          for (int idim = 0; idim < nsd_; ++idim)
          {
            stiff_val += N_XYZ_(idim, vi) * lin_resM_Dus(idim, ui * nsd_ + jdim);
            stiff_val2 += N_XYZ_(idim, vi) * gridvel_int_(idim);
          }
          ecoupl_p(vi, ui * nsd_ + jdim) += fac_dens * N_XYZ_(jdim, vi) * funct_ui +
                                            fac_dens2 * stiff_val2 * dJ_dus(ui * nsd_ + jdim) -
                                            val * stiff_val;
        }
      }
    }

    {
      const double value = -1.0 * timefacfacpre * tau_struct_;
      for (int vi = 0; vi < nen_; ++vi)
      {
        double stiff_val = 0.0;
        for (int idim = 0; idim < nsd_; ++idim)
          stiff_val += value * N_XYZ_(idim, vi) * (Base::gradp_(idim) + reac_tensor_convvel_(idim));
        for (int ui = 0; ui < nen_; ++ui)
          for (int jdim = 0; jdim < nsd_; ++jdim)
            ecoupl_p(vi, ui * nsd_ + jdim) += stiff_val * dphi_dus(ui * nsd_ + jdim);
      }
    }

    if (Base::is_higher_order_ele_)
    {
      const double value = timefacfacpre * tau_struct_ / J_;
      for (int vi = 0; vi < nen_; ++vi)
      {
        for (int ui = 0; ui < nen_; ++ui)
        {
          for (int jdim = 0; jdim < nsd_; ++jdim)
          {
            for (int idim = 0; idim < nsd_; ++idim)
              ecoupl_p(vi, ui * nsd_ + jdim) +=
                  value * N_XYZ_(idim, vi) * mixresLinOD_(idim, ui * nsd_ + jdim);
          }
        }
      }
    }
  }

  //**************************************************************************
  // shape derivatives
  if (nsd_ == 3)
    LinMeshMotion_3D_Pres_OD(ecoupl_p, dphi_dp, dphi_dJ, refporositydot, timefacfacpre);
  else if (nsd_ == 2)
    LinMeshMotion_2D_Pres_OD(ecoupl_p, dphi_dp, dphi_dJ, refporositydot, timefacfacpre);
  else
    dserror("Linearization of the mesh motion is only available in 2D and 3D");
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::Lin3DMeshMotionOD(
    LINALG::Matrix<nsd_ * nen_, nsd_ * nen_>& ecoupl_u, const double& dphi_dp,
    const double& dphi_dJ, const double& refporositydot, const double& timefac,
    const double& timefacfac)
{
  double addstab = 0.0;
  if (Base::fldpara_->RStab() != INPAR::FLUID::reactive_stab_none)
  {
    if (Base::fldpara_->Tds() == INPAR::FLUID::subscales_quasistatic)
      addstab = Base::fldpara_->ViscReaStabFac() * Base::reacoeff_ * Base::tau_(1);
    else
    {
      dserror("Is this factor correct? Check for bugs!");
    }
  }
  //*************************** linearisation of mesh motion in momentum
  // balance**********************************
  // mass
  if (!porofldpara_->IsStationaryMomentum())
  {
    const double fac0 = Base::fac_ * Base::densam_ * (1.0 + addstab) * Base::velint_(0);
    const double fac1 = Base::fac_ * Base::densam_ * (1.0 + addstab) * Base::velint_(1);
    const double fac2 = Base::fac_ * Base::densam_ * (1.0 + addstab) * Base::velint_(2);

    for (int vi = 0; vi < nen_; ++vi)
    {
      const double funct_vi = Base::funct_(vi, 0);
      const double v0 = fac0 * funct_vi;
      for (int ui = 0; ui < nen_; ++ui)
      {
        ecoupl_u(vi * 3, ui * 3) += v0 * Base::derxy_(0, ui);
        ecoupl_u(vi * 3, ui * 3 + 1) += v0 * Base::derxy_(1, ui);
        ecoupl_u(vi * 3, ui * 3 + 2) += v0 * Base::derxy_(2, ui);
      }

      const double v1 = fac1 * funct_vi;
      for (int ui = 0; ui < nen_; ++ui)
      {
        ecoupl_u(vi * 3 + 1, ui * 3) += v1 * Base::derxy_(0, ui);
        ecoupl_u(vi * 3 + 1, ui * 3 + 1) += v1 * Base::derxy_(1, ui);
        ecoupl_u(vi * 3 + 1, ui * 3 + 2) += v1 * Base::derxy_(2, ui);
      }

      const double v2 = fac2 * funct_vi;
      for (int ui = 0; ui < nen_; ++ui)
      {
        ecoupl_u(vi * 3 + 2, ui * 3) += v2 * Base::derxy_(0, ui);
        ecoupl_u(vi * 3 + 2, ui * 3 + 1) += v2 * Base::derxy_(1, ui);
        ecoupl_u(vi * 3 + 2, ui * 3 + 2) += v2 * Base::derxy_(2, ui);
      }
    }
  }

  // rhs
  {
    const double fac_timint =
        Base::fac_ * Base::fldparatimint_->Dt() * Base::fldparatimint_->Theta();
    const double fac0 = fac_timint * (-1.0) * Base::rhsmom_(0);
    const double fac1 = fac_timint * (-1.0) * Base::rhsmom_(1);
    const double fac2 = fac_timint * (-1.0) * Base::rhsmom_(2);

    for (int vi = 0; vi < nen_; ++vi)
    {
      const double funct_vi = Base::funct_(vi, 0);
      const double v0 = fac0 * funct_vi;
      for (int ui = 0; ui < nen_; ++ui)
      {
        ecoupl_u(vi * 3, ui * 3) += v0 * Base::derxy_(0, ui);
        ecoupl_u(vi * 3, ui * 3 + 1) += v0 * Base::derxy_(1, ui);
        ecoupl_u(vi * 3, ui * 3 + 2) += v0 * Base::derxy_(2, ui);
      }

      const double v1 = fac1 * funct_vi;
      for (int ui = 0; ui < nen_; ++ui)
      {
        ecoupl_u(vi * 3 + 1, ui * 3) += v1 * Base::derxy_(0, ui);
        ecoupl_u(vi * 3 + 1, ui * 3 + 1) += v1 * Base::derxy_(1, ui);
        ecoupl_u(vi * 3 + 1, ui * 3 + 2) += v1 * Base::derxy_(2, ui);
      }

      const double v2 = fac2 * funct_vi;
      for (int ui = 0; ui < nen_; ++ui)
      {
        ecoupl_u(vi * 3 + 2, ui * 3) += v2 * Base::derxy_(0, ui);
        ecoupl_u(vi * 3 + 2, ui * 3 + 1) += v2 * Base::derxy_(1, ui);
        ecoupl_u(vi * 3 + 2, ui * 3 + 2) += v2 * Base::derxy_(2, ui);
      }
    }
  }

  //---------reaction term (darcy term)
  {
    const double fac_reac_conv_vel_0 = reac_tensor_convvel_(0) * timefacfac * (1.0 + addstab);
    const double fac_reac_conv_vel_1 = reac_tensor_convvel_(1) * timefacfac * (1.0 + addstab);
    const double fac_reac_conv_vel_2 = reac_tensor_convvel_(2) * timefacfac * (1.0 + addstab);

    for (int vi = 0; vi < nen_; ++vi)
    {
      const double funct_vi = Base::funct_(vi, 0);
      const double v0 = fac_reac_conv_vel_0 * funct_vi;
      for (int ui = 0; ui < nen_; ++ui)
      {
        ecoupl_u(vi * 3, ui * 3) += v0 * Base::derxy_(0, ui);
        ecoupl_u(vi * 3, ui * 3 + 1) += v0 * Base::derxy_(1, ui);
        ecoupl_u(vi * 3, ui * 3 + 2) += v0 * Base::derxy_(2, ui);
      }

      const double v1 = fac_reac_conv_vel_1 * funct_vi;
      for (int ui = 0; ui < nen_; ++ui)
      {
        ecoupl_u(vi * 3 + 1, ui * 3) += v1 * Base::derxy_(0, ui);
        ecoupl_u(vi * 3 + 1, ui * 3 + 1) += v1 * Base::derxy_(1, ui);
        ecoupl_u(vi * 3 + 1, ui * 3 + 2) += v1 * Base::derxy_(2, ui);
      }

      const double v2 = fac_reac_conv_vel_2 * funct_vi;
      for (int ui = 0; ui < nen_; ++ui)
      {
        ecoupl_u(vi * 3 + 2, ui * 3) += v2 * Base::derxy_(0, ui);
        ecoupl_u(vi * 3 + 2, ui * 3 + 1) += v2 * Base::derxy_(1, ui);
        ecoupl_u(vi * 3 + 2, ui * 3 + 2) += v2 * Base::derxy_(2, ui);
      }
    }
  }

  //---------------convective term


#define derxjm_(r, c, d, i) derxjm_##r##c##d(i)

#define derxjm_001(ui) (Base::deriv_(2, ui) * xjm_1_2 - Base::deriv_(1, ui) * xjm_2_2)
#define derxjm_002(ui) (Base::deriv_(1, ui) * xjm_2_1 - Base::deriv_(2, ui) * xjm_1_1)

#define derxjm_100(ui) (Base::deriv_(1, ui) * xjm_2_2 - Base::deriv_(2, ui) * xjm_1_2)
#define derxjm_102(ui) (Base::deriv_(2, ui) * xjm_1_0 - Base::deriv_(1, ui) * xjm_2_0)

#define derxjm_200(ui) (Base::deriv_(2, ui) * xjm_1_1 - Base::deriv_(1, ui) * xjm_2_1)
#define derxjm_201(ui) (Base::deriv_(1, ui) * xjm_2_0 - Base::deriv_(2, ui) * xjm_1_0)

#define derxjm_011(ui) (Base::deriv_(0, ui) * xjm_2_2 - Base::deriv_(2, ui) * xjm_0_2)
#define derxjm_012(ui) (Base::deriv_(2, ui) * xjm_0_1 - Base::deriv_(0, ui) * xjm_2_1)

#define derxjm_110(ui) (Base::deriv_(2, ui) * xjm_0_2 - Base::deriv_(0, ui) * xjm_2_2)
#define derxjm_112(ui) (Base::deriv_(0, ui) * xjm_2_0 - Base::deriv_(2, ui) * xjm_0_0)

#define derxjm_210(ui) (Base::deriv_(0, ui) * xjm_2_1 - Base::deriv_(2, ui) * xjm_0_1)
#define derxjm_211(ui) (Base::deriv_(2, ui) * xjm_0_0 - Base::deriv_(0, ui) * xjm_2_0)

#define derxjm_021(ui) (Base::deriv_(1, ui) * xjm_0_2 - Base::deriv_(0, ui) * xjm_1_2)
#define derxjm_022(ui) (Base::deriv_(0, ui) * xjm_1_1 - Base::deriv_(1, ui) * xjm_0_1)

#define derxjm_120(ui) (Base::deriv_(0, ui) * xjm_1_2 - Base::deriv_(1, ui) * xjm_0_2)
#define derxjm_122(ui) (Base::deriv_(1, ui) * xjm_0_0 - Base::deriv_(0, ui) * xjm_1_0)

#define derxjm_220(ui) (Base::deriv_(1, ui) * xjm_0_1 - Base::deriv_(0, ui) * xjm_1_1)
#define derxjm_221(ui) (Base::deriv_(0, ui) * xjm_1_0 - Base::deriv_(1, ui) * xjm_0_0)

  const double vderiv_0_0 = Base::vderiv_(0, 0);
  const double vderiv_0_1 = Base::vderiv_(0, 1);
  const double vderiv_0_2 = Base::vderiv_(0, 2);
  const double vderiv_1_0 = Base::vderiv_(1, 0);
  const double vderiv_1_1 = Base::vderiv_(1, 1);
  const double vderiv_1_2 = Base::vderiv_(1, 2);
  const double vderiv_2_0 = Base::vderiv_(2, 0);
  const double vderiv_2_1 = Base::vderiv_(2, 1);
  const double vderiv_2_2 = Base::vderiv_(2, 2);

  const double xjm_0_0 = Base::xjm_(0, 0);
  const double xjm_0_1 = Base::xjm_(0, 1);
  const double xjm_0_2 = Base::xjm_(0, 2);
  const double xjm_1_0 = Base::xjm_(1, 0);
  const double xjm_1_1 = Base::xjm_(1, 1);
  const double xjm_1_2 = Base::xjm_(1, 2);
  const double xjm_2_0 = Base::xjm_(2, 0);
  const double xjm_2_1 = Base::xjm_(2, 1);
  const double xjm_2_2 = Base::xjm_(2, 2);

  const double timefacfac_det = timefacfac / Base::det_;

  {
    const double convvelint_0 = Base::convvelint_(0);
    const double convvelint_1 = Base::convvelint_(1);
    const double convvelint_2 = Base::convvelint_(2);

    for (int ui = 0; ui < nen_; ++ui)
    {
      double v00 =
          +convvelint_1 * (vderiv_0_0 * derxjm_(0, 0, 1, ui) + vderiv_0_1 * derxjm_(0, 1, 1, ui) +
                              vderiv_0_2 * derxjm_(0, 2, 1, ui)) +
          convvelint_2 * (vderiv_0_0 * derxjm_(0, 0, 2, ui) + vderiv_0_1 * derxjm_(0, 1, 2, ui) +
                             vderiv_0_2 * derxjm_(0, 2, 2, ui));
      double v01 =
          +convvelint_0 * (vderiv_0_0 * derxjm_(1, 0, 0, ui) + vderiv_0_1 * derxjm_(1, 1, 0, ui) +
                              vderiv_0_2 * derxjm_(1, 2, 0, ui)) +
          convvelint_2 * (vderiv_0_0 * derxjm_(1, 0, 2, ui) + vderiv_0_1 * derxjm_(1, 1, 2, ui) +
                             vderiv_0_2 * derxjm_(1, 2, 2, ui));
      double v02 =
          +convvelint_0 * (vderiv_0_0 * derxjm_(2, 0, 0, ui) + vderiv_0_1 * derxjm_(2, 1, 0, ui) +
                              vderiv_0_2 * derxjm_(2, 2, 0, ui)) +
          convvelint_1 * (vderiv_0_0 * derxjm_(2, 0, 1, ui) + vderiv_0_1 * derxjm_(2, 1, 1, ui) +
                             vderiv_0_2 * derxjm_(2, 2, 1, ui));
      double v10 =
          +convvelint_1 * (vderiv_1_0 * derxjm_(0, 0, 1, ui) + vderiv_1_1 * derxjm_(0, 1, 1, ui) +
                              vderiv_1_2 * derxjm_(0, 2, 1, ui)) +
          convvelint_2 * (vderiv_1_0 * derxjm_(0, 0, 2, ui) + vderiv_1_1 * derxjm_(0, 1, 2, ui) +
                             vderiv_1_2 * derxjm_(0, 2, 2, ui));
      double v11 =
          +convvelint_0 * (vderiv_1_0 * derxjm_(1, 0, 0, ui) + vderiv_1_1 * derxjm_(1, 1, 0, ui) +
                              vderiv_1_2 * derxjm_(1, 2, 0, ui)) +
          convvelint_2 * (vderiv_1_0 * derxjm_(1, 0, 2, ui) + vderiv_1_1 * derxjm_(1, 1, 2, ui) +
                             vderiv_1_2 * derxjm_(1, 2, 2, ui));
      double v12 =
          +convvelint_0 * (vderiv_1_0 * derxjm_(2, 0, 0, ui) + vderiv_1_1 * derxjm_(2, 1, 0, ui) +
                              vderiv_1_2 * derxjm_(2, 2, 0, ui)) +
          convvelint_1 * (vderiv_1_0 * derxjm_(2, 0, 1, ui) + vderiv_1_1 * derxjm_(2, 1, 1, ui) +
                             vderiv_1_2 * derxjm_(2, 2, 1, ui));
      double v20 =
          +convvelint_1 *
              (Base::vderiv_(2, 0) * derxjm_(0, 0, 1, ui) + vderiv_2_1 * derxjm_(0, 1, 1, ui) +
                  vderiv_2_2 * derxjm_(0, 2, 1, ui)) +
          convvelint_2 * (Base::vderiv_(2, 0) * derxjm_(0, 0, 2, ui) +
                             vderiv_2_1 * derxjm_(0, 1, 2, ui) + vderiv_2_2 * derxjm_(0, 2, 2, ui));
      double v21 =
          +convvelint_0 *
              (Base::vderiv_(2, 0) * derxjm_(1, 0, 0, ui) + vderiv_2_1 * derxjm_(1, 1, 0, ui) +
                  vderiv_2_2 * derxjm_(1, 2, 0, ui)) +
          convvelint_2 * (Base::vderiv_(2, 0) * derxjm_(1, 0, 2, ui) +
                             vderiv_2_1 * derxjm_(1, 1, 2, ui) + vderiv_2_2 * derxjm_(1, 2, 2, ui));
      double v22 =
          +convvelint_0 *
              (Base::vderiv_(2, 0) * derxjm_(2, 0, 0, ui) + vderiv_2_1 * derxjm_(2, 1, 0, ui) +
                  vderiv_2_2 * derxjm_(2, 2, 0, ui)) +
          convvelint_1 * (Base::vderiv_(2, 0) * derxjm_(2, 0, 1, ui) +
                             vderiv_2_1 * derxjm_(2, 1, 1, ui) + vderiv_2_2 * derxjm_(2, 2, 1, ui));

      for (int vi = 0; vi < nen_; ++vi)
      {
        double v = Base::densaf_ * timefacfac_det * Base::funct_(vi) * (1.0 + addstab);

        ecoupl_u(vi * 3 + 0, ui * 3 + 0) += v * v00;
        ecoupl_u(vi * 3 + 0, ui * 3 + 1) += v * v01;
        ecoupl_u(vi * 3 + 0, ui * 3 + 2) += v * v02;

        ecoupl_u(vi * 3 + 1, ui * 3 + 0) += v * v10;
        ecoupl_u(vi * 3 + 1, ui * 3 + 1) += v * v11;
        ecoupl_u(vi * 3 + 1, ui * 3 + 2) += v * v12;

        ecoupl_u(vi * 3 + 2, ui * 3 + 0) += v * v20;
        ecoupl_u(vi * 3 + 2, ui * 3 + 1) += v * v21;
        ecoupl_u(vi * 3 + 2, ui * 3 + 2) += v * v22;
      }
    }

    // pressure
    const double v = press_ * timefacfac_det;
    for (int vi = 0; vi < nen_; ++vi)
    {
      const double deriv_vi_0 = Base::deriv_(0, vi);
      const double deriv_vi_1 = Base::deriv_(1, vi);
      const double deriv_vi_2 = Base::deriv_(2, vi);

      for (int ui = 0; ui < nen_; ++ui)
      {
        ecoupl_u(vi * 3, ui * 3 + 1) +=
            v * (deriv_vi_0 * derxjm_(0, 0, 1, ui) + deriv_vi_1 * derxjm_(0, 1, 1, ui) +
                    deriv_vi_2 * derxjm_(0, 2, 1, ui));
        ecoupl_u(vi * 3, ui * 3 + 2) +=
            v * (deriv_vi_0 * derxjm_(0, 0, 2, ui) + deriv_vi_1 * derxjm_(0, 1, 2, ui) +
                    deriv_vi_2 * derxjm_(0, 2, 2, ui));

        ecoupl_u(vi * 3 + 1, ui * 3 + 0) +=
            v * (deriv_vi_0 * derxjm_(1, 0, 0, ui) + deriv_vi_1 * derxjm_(1, 1, 0, ui) +
                    deriv_vi_2 * derxjm_(1, 2, 0, ui));
        ecoupl_u(vi * 3 + 1, ui * 3 + 2) +=
            v * (deriv_vi_0 * derxjm_(1, 0, 2, ui) + deriv_vi_1 * derxjm_(1, 1, 2, ui) +
                    deriv_vi_2 * derxjm_(1, 2, 2, ui));

        ecoupl_u(vi * 3 + 2, ui * 3 + 0) +=
            v * (deriv_vi_0 * derxjm_(2, 0, 0, ui) + deriv_vi_1 * derxjm_(2, 1, 0, ui) +
                    deriv_vi_2 * derxjm_(2, 2, 0, ui));
        ecoupl_u(vi * 3 + 2, ui * 3 + 1) +=
            v * (deriv_vi_0 * derxjm_(2, 0, 1, ui) + deriv_vi_1 * derxjm_(2, 1, 1, ui) +
                    deriv_vi_2 * derxjm_(2, 2, 1, ui));
      }
    }
  }

  // //---------viscous term (brinkman term)
  const double xji_00 = Base::xji_(0, 0);
  const double xji_01 = Base::xji_(0, 1);
  const double xji_02 = Base::xji_(0, 2);
  const double xji_10 = Base::xji_(1, 0);
  const double xji_11 = Base::xji_(1, 1);
  const double xji_12 = Base::xji_(1, 2);
  const double xji_20 = Base::xji_(2, 0);
  const double xji_21 = Base::xji_(2, 1);
  const double xji_22 = Base::xji_(2, 2);

  const double porosity_inv = 1.0 / porosity_;

  if (Base::visceff_)
  {
    const double vderxy_0_0 = 2.0 * Base::vderxy_(0, 0);
    const double vderxy_1_1 = 2.0 * Base::vderxy_(1, 1);
    const double vderxy_2_2 = 2.0 * Base::vderxy_(2, 2);
    const double vderxy_0_1 = Base::vderxy_(0, 1) + Base::vderxy_(1, 0);
    const double vderxy_0_2 = Base::vderxy_(0, 2) + Base::vderxy_(2, 0);
    const double vderxy_1_2 = Base::vderxy_(1, 2) + Base::vderxy_(2, 1);

    const double refgrad_porosity_0 = refgrad_porosity_(0);
    const double refgrad_porosity_1 = refgrad_porosity_(1);
    const double refgrad_porosity_2 = refgrad_porosity_(2);

    // part 1: derivative of 1/det

    {
      const double v = Base::visceff_ * timefac * Base::fac_ * (1.0 + addstab);
      for (int ui = 0; ui < nen_; ++ui)
      {
        const double derinvJ0 = -v * (Base::deriv_(0, ui) * xji_00 + Base::deriv_(1, ui) * xji_01 +
                                         Base::deriv_(2, ui) * xji_02);
        const double derinvJ1 = -v * (Base::deriv_(0, ui) * xji_10 + Base::deriv_(1, ui) * xji_11 +
                                         Base::deriv_(2, ui) * xji_12);
        const double derinvJ2 = -v * (Base::deriv_(0, ui) * xji_20 + Base::deriv_(1, ui) * xji_21 +
                                         Base::deriv_(2, ui) * xji_22);

        for (int vi = 0; vi < nen_; ++vi)
        {
          const double visres0 = Base::derxy_(0, vi) * vderxy_0_0 +
                                 Base::derxy_(1, vi) * vderxy_0_1 +
                                 Base::derxy_(2, vi) * vderxy_0_2;
          const double visres1 = Base::derxy_(0, vi) * vderxy_0_1 +
                                 Base::derxy_(1, vi) * vderxy_1_1 +
                                 Base::derxy_(2, vi) * vderxy_1_2;
          const double visres2 = Base::derxy_(0, vi) * vderxy_0_2 +
                                 Base::derxy_(1, vi) * vderxy_1_2 +
                                 Base::derxy_(2, vi) * vderxy_2_2;
          ecoupl_u(vi * 3 + 0, ui * 3 + 0) += derinvJ0 * visres0;
          ecoupl_u(vi * 3 + 1, ui * 3 + 0) += derinvJ0 * visres1;
          ecoupl_u(vi * 3 + 2, ui * 3 + 0) += derinvJ0 * visres2;

          ecoupl_u(vi * 3 + 0, ui * 3 + 1) += derinvJ1 * visres0;
          ecoupl_u(vi * 3 + 1, ui * 3 + 1) += derinvJ1 * visres1;
          ecoupl_u(vi * 3 + 2, ui * 3 + 1) += derinvJ1 * visres2;

          ecoupl_u(vi * 3 + 0, ui * 3 + 2) += derinvJ2 * visres0;
          ecoupl_u(vi * 3 + 1, ui * 3 + 2) += derinvJ2 * visres1;
          ecoupl_u(vi * 3 + 2, ui * 3 + 2) += derinvJ2 * visres2;

          const double funct_vi = Base::funct_(vi);
          const double visres0_poro = refgrad_porosity_0 * funct_vi * vderxy_0_0 +
                                      refgrad_porosity_1 * funct_vi * vderxy_0_1 +
                                      refgrad_porosity_2 * funct_vi * vderxy_0_2;
          const double visres1_poro = refgrad_porosity_0 * funct_vi * vderxy_0_1 +
                                      refgrad_porosity_1 * funct_vi * vderxy_1_1 +
                                      refgrad_porosity_2 * funct_vi * vderxy_1_2;
          const double visres2_poro = refgrad_porosity_0 * funct_vi * vderxy_0_2 +
                                      refgrad_porosity_1 * funct_vi * vderxy_1_2 +
                                      refgrad_porosity_2 * funct_vi * vderxy_2_2;

          ecoupl_u(vi * 3 + 0, ui * 3 + 0) += -1.0 * derinvJ0 * porosity_inv * visres0_poro;
          ecoupl_u(vi * 3 + 1, ui * 3 + 0) += -1.0 * derinvJ0 * porosity_inv * visres1_poro;
          ecoupl_u(vi * 3 + 2, ui * 3 + 0) += -1.0 * derinvJ0 * porosity_inv * visres2_poro;

          ecoupl_u(vi * 3 + 0, ui * 3 + 1) += -1.0 * derinvJ1 * porosity_inv * visres0_poro;
          ecoupl_u(vi * 3 + 1, ui * 3 + 1) += -1.0 * derinvJ1 * porosity_inv * visres1_poro;
          ecoupl_u(vi * 3 + 2, ui * 3 + 1) += -1.0 * derinvJ1 * porosity_inv * visres2_poro;

          ecoupl_u(vi * 3 + 0, ui * 3 + 2) += -1.0 * derinvJ2 * porosity_inv * visres0_poro;
          ecoupl_u(vi * 3 + 1, ui * 3 + 2) += -1.0 * derinvJ2 * porosity_inv * visres1_poro;
          ecoupl_u(vi * 3 + 2, ui * 3 + 2) += -1.0 * derinvJ2 * porosity_inv * visres2_poro;
        }
      }
    }

    // part 2: derivative of viscosity residual

    {
      const double v = timefacfac_det * Base::visceff_ * (1.0 + addstab);

      for (int ui = 0; ui < nen_; ++ui)
      {
        const double derxjm_100_ui = derxjm_100(ui);
        const double derxjm_110_ui = derxjm_110(ui);
        const double derxjm_120_ui = derxjm_120(ui);
        const double derxjm_200_ui = derxjm_200(ui);
        const double derxjm_210_ui = derxjm_210(ui);
        const double derxjm_220_ui = derxjm_220(ui);
        const double derxjm_201_ui = derxjm_201(ui);
        const double derxjm_001_ui = derxjm_001(ui);
        const double derxjm_002_ui = derxjm_002(ui);
        const double derxjm_011_ui = derxjm_011(ui);
        const double derxjm_021_ui = derxjm_021(ui);
        const double derxjm_012_ui = derxjm_012(ui);
        const double derxjm_022_ui = derxjm_022(ui);
        const double derxjm_221_ui = derxjm_221(ui);
        const double derxjm_102_ui = derxjm_102(ui);
        const double derxjm_112_ui = derxjm_112(ui);
        const double derxjm_122_ui = derxjm_122(ui);
        const double derxjm_211_ui = derxjm_211(ui);

        {
          const double v0 =
              -vderiv_0_0 * (xji_10 * derxjm_100_ui + xji_10 * derxjm_100_ui +
                                xji_20 * derxjm_200_ui + xji_20 * derxjm_200_ui) -
              vderiv_0_1 * (xji_11 * derxjm_100_ui + xji_10 * derxjm_110_ui +
                               xji_21 * derxjm_200_ui + xji_20 * derxjm_210_ui) -
              vderiv_0_2 * (xji_12 * derxjm_100_ui + xji_10 * derxjm_120_ui +
                               xji_22 * derxjm_200_ui + xji_20 * derxjm_220_ui) -
              vderiv_1_0 * (derxjm_100_ui * xji_00) - vderiv_1_1 * (derxjm_100_ui * xji_01) -
              vderiv_1_2 * (derxjm_100_ui * xji_02) - vderiv_2_0 * (derxjm_200_ui * xji_00) -
              vderiv_2_1 * (derxjm_200_ui * xji_01) - vderiv_2_2 * (derxjm_200_ui * xji_02);
          const double v1 =
              -vderiv_0_0 * (xji_10 * derxjm_110_ui + xji_11 * derxjm_100_ui +
                                xji_20 * derxjm_210_ui + xji_21 * derxjm_200_ui) -
              vderiv_0_1 * (xji_11 * derxjm_110_ui + xji_11 * derxjm_110_ui +
                               xji_21 * derxjm_210_ui + xji_21 * derxjm_210_ui) -
              vderiv_0_2 * (xji_12 * derxjm_110_ui + xji_11 * derxjm_120_ui +
                               xji_22 * derxjm_210_ui + xji_21 * derxjm_220_ui) -
              vderiv_1_0 * (derxjm_110_ui * xji_00) - vderiv_1_1 * (derxjm_110_ui * xji_01) -
              vderiv_1_2 * (derxjm_110_ui * xji_02) - vderiv_2_0 * (derxjm_210_ui * xji_00) -
              vderiv_2_1 * (derxjm_210_ui * xji_01) - vderiv_2_2 * (derxjm_210_ui * xji_02);
          const double v2 =
              -vderiv_0_0 * (xji_10 * derxjm_120_ui + xji_12 * derxjm_100_ui +
                                xji_20 * derxjm_220_ui + xji_22 * derxjm_200_ui) -
              vderiv_0_1 * (xji_11 * derxjm_120_ui + xji_12 * derxjm_110_ui +
                               xji_21 * derxjm_220_ui + xji_22 * derxjm_210_ui) -
              vderiv_0_2 * (xji_12 * derxjm_120_ui + xji_12 * derxjm_120_ui +
                               xji_22 * derxjm_220_ui + xji_22 * derxjm_220_ui) -
              vderiv_1_0 * (derxjm_120_ui * xji_00) - vderiv_1_1 * (derxjm_120_ui * xji_01) -
              vderiv_1_2 * (derxjm_120_ui * xji_02) - vderiv_2_0 * (derxjm_220_ui * xji_00) -
              vderiv_2_1 * (derxjm_220_ui * xji_01) - vderiv_2_2 * (derxjm_220_ui * xji_02);

          for (int vi = 0; vi < nen_; ++vi)
          {
            ecoupl_u(vi * 3 + 0, ui * 3 + 0) +=
                v * (Base::deriv_(0, vi) * v0 + Base::deriv_(1, vi) * v1 +
                        Base::deriv_(2, vi) * v2) -
                v * Base::funct_(vi) * porosity_inv *
                    (refgrad_porosity_0 * v0 + refgrad_porosity_1 * v1 + refgrad_porosity_2 * v2);
          }
        }

        // //////////////////////////////////////////////////////////////
        {
          const double v0 = -vderiv_0_0 * (2 * derxjm_001_ui * xji_00 + 2 * derxjm_001_ui * xji_00 +
                                              xji_20 * derxjm_201_ui + xji_20 * derxjm_201_ui) -
                            vderiv_0_1 * (2 * derxjm_011_ui * xji_00 + 2 * derxjm_001_ui * xji_01 +
                                             xji_21 * derxjm_201_ui + xji_20 * derxjm_211_ui) -
                            vderiv_0_2 * (2 * derxjm_021_ui * xji_00 + 2 * derxjm_001_ui * xji_02 +
                                             xji_22 * derxjm_201_ui + xji_20 * derxjm_221_ui) -
                            vderiv_1_0 * (derxjm_001_ui * xji_10) -
                            vderiv_1_1 * (derxjm_011_ui * xji_10) -
                            vderiv_1_2 * (derxjm_021_ui * xji_10) -
                            vderiv_2_0 * (derxjm_201_ui * xji_00 + derxjm_001_ui * xji_20) -
                            vderiv_2_1 * (derxjm_201_ui * xji_01 + derxjm_011_ui * xji_20) -
                            vderiv_2_2 * (derxjm_201_ui * xji_02 + derxjm_021_ui * xji_20);
          const double v1 = -vderiv_0_0 * (2 * derxjm_011_ui * xji_00 + 2 * derxjm_001_ui * xji_01 +
                                              xji_21 * derxjm_201_ui + xji_20 * derxjm_211_ui) -
                            vderiv_0_1 * (2 * derxjm_011_ui * xji_01 + 2 * derxjm_011_ui * xji_01 +
                                             xji_21 * derxjm_211_ui + xji_21 * derxjm_211_ui) -
                            vderiv_0_2 * (2 * derxjm_011_ui * xji_02 + 2 * derxjm_021_ui * xji_01 +
                                             xji_21 * derxjm_221_ui + xji_22 * derxjm_211_ui) -
                            vderiv_1_0 * (derxjm_001_ui * xji_11) -
                            vderiv_1_1 * (derxjm_011_ui * xji_11) -
                            vderiv_1_2 * (derxjm_021_ui * xji_11) -
                            vderiv_2_0 * (derxjm_211_ui * xji_00 + derxjm_001_ui * xji_21) -
                            vderiv_2_1 * (derxjm_211_ui * xji_01 + derxjm_011_ui * xji_21) -
                            vderiv_2_2 * (derxjm_211_ui * xji_02 + derxjm_021_ui * xji_21);
          const double v2 = -vderiv_0_0 * (2 * derxjm_021_ui * xji_00 + 2 * derxjm_001_ui * xji_02 +
                                              xji_22 * derxjm_201_ui + xji_20 * derxjm_221_ui) -
                            vderiv_0_1 * (2 * derxjm_011_ui * xji_02 + 2 * derxjm_021_ui * xji_01 +
                                             xji_21 * derxjm_221_ui + xji_22 * derxjm_211_ui) -
                            vderiv_0_2 * (2 * derxjm_021_ui * xji_02 + 2 * derxjm_021_ui * xji_02 +
                                             xji_22 * derxjm_221_ui + xji_22 * derxjm_221_ui) -
                            vderiv_1_0 * (derxjm_001_ui * xji_12) -
                            vderiv_1_1 * (derxjm_011_ui * xji_12) -
                            vderiv_1_2 * (derxjm_021_ui * xji_12) -
                            vderiv_2_0 * (derxjm_221_ui * xji_00 + derxjm_001_ui * xji_22) -
                            vderiv_2_1 * (derxjm_221_ui * xji_01 + derxjm_011_ui * xji_22) -
                            vderiv_2_2 * (derxjm_221_ui * xji_02 + derxjm_021_ui * xji_22);

          for (int vi = 0; vi < nen_; ++vi)
          {
            ecoupl_u(vi * 3 + 0, ui * 3 + 1) +=
                v * (Base::deriv_(0, vi) * v0 + Base::deriv_(1, vi) * v1 +
                        Base::deriv_(2, vi) * v2) -
                v * Base::funct_(vi) * porosity_inv *
                    (refgrad_porosity_0 * v0 + refgrad_porosity_1 * v1 + refgrad_porosity_2 * v2);
          }
        }

        // //////////////////////////////////////////////////////////////
        {
          const double v0 = -vderiv_0_0 * (2 * derxjm_002_ui * xji_00 + 2 * derxjm_002_ui * xji_00 +
                                              xji_10 * derxjm_102_ui + xji_10 * derxjm_102_ui) -
                            vderiv_0_1 * (2 * derxjm_012_ui * xji_00 + 2 * derxjm_002_ui * xji_01 +
                                             xji_11 * derxjm_102_ui + xji_10 * derxjm_112_ui) -
                            vderiv_0_2 * (2 * derxjm_022_ui * xji_00 + 2 * derxjm_002_ui * xji_02 +
                                             xji_12 * derxjm_102_ui + xji_10 * derxjm_122_ui) -
                            vderiv_1_0 * (derxjm_002_ui * xji_10 + derxjm_102_ui * xji_00) -
                            vderiv_1_1 * (derxjm_012_ui * xji_10 + derxjm_102_ui * xji_01) -
                            vderiv_1_2 * (derxjm_022_ui * xji_10 + derxjm_102_ui * xji_02) -
                            vderiv_2_0 * (derxjm_002_ui * xji_20) -
                            vderiv_2_1 * (derxjm_012_ui * xji_20) -
                            vderiv_2_2 * (derxjm_022_ui * xji_20);
          const double v1 = -vderiv_0_0 * (2 * derxjm_012_ui * xji_00 + 2 * derxjm_002_ui * xji_01 +
                                              xji_11 * derxjm_102_ui + xji_10 * derxjm_112_ui) -
                            vderiv_0_1 * (2 * derxjm_012_ui * xji_01 + 2 * derxjm_012_ui * xji_01 +
                                             xji_11 * derxjm_112_ui + xji_11 * derxjm_112_ui) -
                            vderiv_0_2 * (2 * derxjm_012_ui * xji_02 + 2 * derxjm_022_ui * xji_01 +
                                             xji_11 * derxjm_122_ui + xji_12 * derxjm_112_ui) -
                            vderiv_1_0 * (derxjm_002_ui * xji_11 + derxjm_112_ui * xji_00) -
                            vderiv_1_1 * (derxjm_012_ui * xji_11 + derxjm_112_ui * xji_01) -
                            vderiv_1_2 * (derxjm_022_ui * xji_11 + derxjm_112_ui * xji_02) -
                            vderiv_2_0 * (derxjm_002_ui * xji_21) -
                            vderiv_2_1 * (derxjm_012_ui * xji_21) -
                            vderiv_2_2 * (derxjm_022_ui * xji_21);
          const double v2 = -vderiv_0_0 * (2 * derxjm_022_ui * xji_00 + 2 * derxjm_002_ui * xji_02 +
                                              xji_12 * derxjm_102_ui + xji_10 * derxjm_122_ui) -
                            vderiv_0_1 * (2 * derxjm_012_ui * xji_02 + 2 * derxjm_022_ui * xji_01 +
                                             xji_11 * derxjm_122_ui + xji_12 * derxjm_112_ui) -
                            vderiv_0_2 * (2 * derxjm_022_ui * xji_02 + 2 * derxjm_022_ui * xji_02 +
                                             xji_12 * derxjm_122_ui + xji_12 * derxjm_122_ui) -
                            vderiv_1_0 * (derxjm_002_ui * xji_12 + derxjm_122_ui * xji_00) -
                            vderiv_1_1 * (derxjm_012_ui * xji_12 + derxjm_122_ui * xji_01) -
                            vderiv_1_2 * (derxjm_022_ui * xji_12 + derxjm_122_ui * xji_02) -
                            vderiv_2_0 * (derxjm_002_ui * xji_22) -
                            vderiv_2_1 * (derxjm_012_ui * xji_22) -
                            vderiv_2_2 * (derxjm_022_ui * xji_22);

          for (int vi = 0; vi < nen_; ++vi)
          {
            ecoupl_u(vi * 3 + 0, ui * 3 + 2) +=
                v * (Base::deriv_(0, vi) * v0 + Base::deriv_(1, vi) * v1 +
                        Base::deriv_(2, vi) * v2) -
                v * Base::funct_(vi) * porosity_inv *
                    (refgrad_porosity_0 * v0 + refgrad_porosity_1 * v1 + refgrad_porosity_2 * v2);
          }
        }

        // //////////////////////////////////////////////////////////////
        {
          const double v0 = -vderiv_0_0 * (derxjm_100_ui * xji_00) -
                            vderiv_0_1 * (derxjm_110_ui * xji_00) -
                            vderiv_0_2 * (derxjm_120_ui * xji_00) -
                            vderiv_1_0 * (2 * xji_10 * derxjm_100_ui + 2 * xji_10 * derxjm_100_ui +
                                             xji_20 * derxjm_200_ui + xji_20 * derxjm_200_ui) -
                            vderiv_1_1 * (2 * xji_11 * derxjm_100_ui + 2 * xji_10 * derxjm_110_ui +
                                             xji_21 * derxjm_200_ui + xji_20 * derxjm_210_ui) -
                            vderiv_1_2 * (2 * xji_12 * derxjm_100_ui + 2 * xji_10 * derxjm_120_ui +
                                             xji_22 * derxjm_200_ui + xji_20 * derxjm_220_ui) -
                            vderiv_2_0 * (derxjm_200_ui * xji_10 + derxjm_100_ui * xji_20) -
                            vderiv_2_1 * (derxjm_200_ui * xji_11 + derxjm_110_ui * xji_20) -
                            vderiv_2_2 * (derxjm_200_ui * xji_12 + derxjm_120_ui * xji_20);
          const double v1 = -vderiv_0_0 * (derxjm_100_ui * xji_01) -
                            vderiv_0_1 * (derxjm_110_ui * xji_01) -
                            vderiv_0_2 * (derxjm_120_ui * xji_01) -
                            vderiv_1_0 * (2 * xji_10 * derxjm_110_ui + 2 * xji_11 * derxjm_100_ui +
                                             xji_20 * derxjm_210_ui + xji_21 * derxjm_200_ui) -
                            vderiv_1_1 * (2 * xji_11 * derxjm_110_ui + 2 * xji_11 * derxjm_110_ui +
                                             xji_21 * derxjm_210_ui + xji_21 * derxjm_210_ui) -
                            vderiv_1_2 * (2 * xji_12 * derxjm_110_ui + 2 * xji_11 * derxjm_120_ui +
                                             xji_22 * derxjm_210_ui + xji_21 * derxjm_220_ui) -
                            vderiv_2_0 * (derxjm_210_ui * xji_10 + derxjm_100_ui * xji_21) -
                            vderiv_2_1 * (derxjm_210_ui * xji_11 + derxjm_110_ui * xji_21) -
                            vderiv_2_2 * (derxjm_210_ui * xji_12 + derxjm_120_ui * xji_21);
          const double v2 = -vderiv_0_0 * (derxjm_100_ui * xji_02) -
                            vderiv_0_1 * (derxjm_110_ui * xji_02) -
                            vderiv_0_2 * (derxjm_120_ui * xji_02) -
                            vderiv_1_0 * (2 * xji_10 * derxjm_120_ui + 2 * xji_12 * derxjm_100_ui +
                                             xji_20 * derxjm_220_ui + xji_22 * derxjm_200_ui) -
                            vderiv_1_1 * (2 * xji_11 * derxjm_120_ui + 2 * xji_12 * derxjm_110_ui +
                                             xji_21 * derxjm_220_ui + xji_22 * derxjm_210_ui) -
                            vderiv_1_2 * (2 * xji_12 * derxjm_120_ui + 2 * xji_12 * derxjm_120_ui +
                                             xji_22 * derxjm_220_ui + xji_22 * derxjm_220_ui) -
                            vderiv_2_0 * (derxjm_220_ui * xji_10 + derxjm_100_ui * xji_22) -
                            vderiv_2_1 * (derxjm_220_ui * xji_11 + derxjm_110_ui * xji_22) -
                            vderiv_2_2 * (derxjm_220_ui * xji_12 + derxjm_120_ui * xji_22);

          for (int vi = 0; vi < nen_; ++vi)
          {
            ecoupl_u(vi * 3 + 1, ui * 3 + 0) +=
                v * (Base::deriv_(0, vi) * v0 + Base::deriv_(1, vi) * v1 +
                        Base::deriv_(2, vi) * v2) -
                v * Base::funct_(vi) * porosity_inv *
                    (refgrad_porosity_0 * v0 + refgrad_porosity_1 * v1 + refgrad_porosity_2 * v2);
          }
        }

        // //////////////////////////////////////////////////////////////
        {
          const double v0 =
              -vderiv_0_0 * (derxjm_001_ui * xji_10) - vderiv_0_1 * (derxjm_001_ui * xji_11) -
              vderiv_0_2 * (derxjm_001_ui * xji_12) -
              vderiv_1_0 * (xji_00 * derxjm_001_ui + xji_00 * derxjm_001_ui +
                               xji_20 * derxjm_201_ui + xji_20 * derxjm_201_ui) -
              vderiv_1_1 * (xji_01 * derxjm_001_ui + xji_00 * derxjm_011_ui +
                               xji_21 * derxjm_201_ui + xji_20 * derxjm_211_ui) -
              vderiv_1_2 * (xji_02 * derxjm_001_ui + xji_00 * derxjm_021_ui +
                               xji_22 * derxjm_201_ui + xji_20 * derxjm_221_ui) -
              vderiv_2_0 * (derxjm_201_ui * xji_10) - vderiv_2_1 * (derxjm_201_ui * xji_11) -
              vderiv_2_2 * (derxjm_201_ui * xji_12);
          const double v1 =
              -vderiv_0_0 * (derxjm_011_ui * xji_10) - vderiv_0_1 * (derxjm_011_ui * xji_11) -
              vderiv_0_2 * (derxjm_011_ui * xji_12) -
              vderiv_1_0 * (xji_00 * derxjm_011_ui + xji_01 * derxjm_001_ui +
                               xji_20 * derxjm_211_ui + xji_21 * derxjm_201_ui) -
              vderiv_1_1 * (xji_01 * derxjm_011_ui + xji_01 * derxjm_011_ui +
                               xji_21 * derxjm_211_ui + xji_21 * derxjm_211_ui) -
              vderiv_1_2 * (xji_02 * derxjm_011_ui + xji_01 * derxjm_021_ui +
                               xji_22 * derxjm_211_ui + xji_21 * derxjm_221_ui) -
              vderiv_2_0 * (derxjm_211_ui * xji_10) - vderiv_2_1 * (derxjm_211_ui * xji_11) -
              vderiv_2_2 * (derxjm_211_ui * xji_12);
          const double v2 =
              -vderiv_0_0 * (derxjm_021_ui * xji_10) - vderiv_0_1 * (derxjm_021_ui * xji_11) -
              vderiv_0_2 * (derxjm_021_ui * xji_12) -
              vderiv_1_0 * (xji_00 * derxjm_021_ui + xji_02 * derxjm_001_ui +
                               xji_20 * derxjm_221_ui + xji_22 * derxjm_201_ui) -
              vderiv_1_1 * (xji_01 * derxjm_021_ui + xji_02 * derxjm_011_ui +
                               xji_21 * derxjm_221_ui + xji_22 * derxjm_211_ui) -
              vderiv_1_2 * (xji_02 * derxjm_021_ui + xji_02 * derxjm_021_ui +
                               xji_22 * derxjm_221_ui + xji_22 * derxjm_221_ui) -
              vderiv_2_0 * (derxjm_221_ui * xji_10) - vderiv_2_1 * (derxjm_221_ui * xji_11) -
              vderiv_2_2 * (derxjm_221_ui * xji_12);

          for (int vi = 0; vi < nen_; ++vi)
          {
            ecoupl_u(vi * 3 + 1, ui * 3 + 1) +=
                v * (Base::deriv_(0, vi) * v0 + Base::deriv_(1, vi) * v1 +
                        Base::deriv_(2, vi) * v2) -
                v * Base::funct_(vi) * porosity_inv *
                    (refgrad_porosity_0 * v0 + refgrad_porosity_1 * v1 + refgrad_porosity_2 * v2);
          }
        }

        // //////////////////////////////////////////////////////////////
        {
          const double v0 =
              -vderiv_0_0 * (derxjm_002_ui * xji_10 + derxjm_102_ui * xji_00) -
              vderiv_0_1 * (derxjm_002_ui * xji_11 + derxjm_112_ui * xji_00) -
              vderiv_0_2 * (derxjm_002_ui * xji_12 + derxjm_122_ui * xji_00) -
              vderiv_1_0 * (xji_00 * derxjm_002_ui + xji_00 * derxjm_002_ui +
                               2 * xji_10 * derxjm_102_ui + 2 * xji_10 * derxjm_102_ui) -
              vderiv_1_1 * (xji_01 * derxjm_002_ui + xji_00 * derxjm_012_ui +
                               2 * xji_11 * derxjm_102_ui + 2 * xji_10 * derxjm_112_ui) -
              vderiv_1_2 * (xji_02 * derxjm_002_ui + xji_00 * derxjm_022_ui +
                               2 * xji_12 * derxjm_102_ui + 2 * xji_10 * derxjm_122_ui) -
              vderiv_2_0 * (derxjm_102_ui * xji_20) - vderiv_2_1 * (derxjm_112_ui * xji_20) -
              vderiv_2_2 * (derxjm_122_ui * xji_20);
          const double v1 =
              -vderiv_0_0 * (derxjm_012_ui * xji_10 + derxjm_102_ui * xji_01) -
              vderiv_0_1 * (derxjm_012_ui * xji_11 + derxjm_112_ui * xji_01) -
              vderiv_0_2 * (derxjm_012_ui * xji_12 + derxjm_122_ui * xji_01) -
              vderiv_1_0 * (xji_00 * derxjm_012_ui + xji_01 * derxjm_002_ui +
                               2 * xji_10 * derxjm_112_ui + 2 * xji_11 * derxjm_102_ui) -
              vderiv_1_1 * (xji_01 * derxjm_012_ui + xji_01 * derxjm_012_ui +
                               2 * xji_11 * derxjm_112_ui + 2 * xji_11 * derxjm_112_ui) -
              vderiv_1_2 * (xji_02 * derxjm_012_ui + xji_01 * derxjm_022_ui +
                               2 * xji_12 * derxjm_112_ui + 2 * xji_11 * derxjm_122_ui) -
              vderiv_2_0 * (derxjm_102_ui * xji_21) - vderiv_2_1 * (derxjm_112_ui * xji_21) -
              vderiv_2_2 * (derxjm_122_ui * xji_21);
          const double v2 =
              -vderiv_0_0 * (derxjm_022_ui * xji_10 + derxjm_102_ui * xji_02) -
              vderiv_0_1 * (derxjm_022_ui * xji_11 + derxjm_112_ui * xji_02) -
              vderiv_0_2 * (derxjm_022_ui * xji_12 + derxjm_122_ui * xji_02) -
              vderiv_1_0 * (xji_00 * derxjm_022_ui + xji_02 * derxjm_002_ui +
                               2 * xji_10 * derxjm_122_ui + 2 * xji_12 * derxjm_102_ui) -
              vderiv_1_1 * (xji_01 * derxjm_022_ui + xji_02 * derxjm_012_ui +
                               2 * xji_11 * derxjm_122_ui + 2 * xji_12 * derxjm_112_ui) -
              vderiv_1_2 * (xji_02 * derxjm_022_ui + xji_02 * derxjm_022_ui +
                               2 * xji_12 * derxjm_122_ui + 2 * xji_12 * derxjm_122_ui) -
              vderiv_2_0 * (derxjm_102_ui * xji_22) - vderiv_2_1 * (derxjm_112_ui * xji_22) -
              vderiv_2_2 * (derxjm_122_ui * xji_22);

          for (int vi = 0; vi < nen_; ++vi)
          {
            ecoupl_u(vi * 3 + 1, ui * 3 + 2) +=
                v * (Base::deriv_(0, vi) * v0 + Base::deriv_(1, vi) * v1 +
                        Base::deriv_(2, vi) * v2) -
                v * Base::funct_(vi) * porosity_inv *
                    (refgrad_porosity_0 * v0 + refgrad_porosity_1 * v1 + refgrad_porosity_2 * v2);
          }
        }

        // //////////////////////////////////////////////////////////////
        {
          const double v0 =
              -vderiv_0_0 * (derxjm_200_ui * xji_00) - vderiv_0_1 * (derxjm_210_ui * xji_00) -
              vderiv_0_2 * (derxjm_220_ui * xji_00) -
              vderiv_1_0 * (derxjm_200_ui * xji_10 + derxjm_100_ui * xji_20) -
              vderiv_1_1 * (derxjm_210_ui * xji_10 + derxjm_100_ui * xji_21) -
              vderiv_1_2 * (derxjm_220_ui * xji_10 + derxjm_100_ui * xji_22) -
              vderiv_2_0 * (xji_10 * derxjm_100_ui + xji_10 * derxjm_100_ui +
                               2 * xji_20 * derxjm_200_ui + 2 * xji_20 * derxjm_200_ui) -
              vderiv_2_1 * (xji_11 * derxjm_100_ui + xji_10 * derxjm_110_ui +
                               2 * xji_21 * derxjm_200_ui + 2 * xji_20 * derxjm_210_ui) -
              vderiv_2_2 * (xji_12 * derxjm_100_ui + xji_10 * derxjm_120_ui +
                               2 * xji_22 * derxjm_200_ui + 2 * xji_20 * derxjm_220_ui);
          const double v1 =
              -vderiv_0_0 * (derxjm_200_ui * xji_01) - vderiv_0_1 * (derxjm_210_ui * xji_01) -
              vderiv_0_2 * (derxjm_220_ui * xji_01) -
              vderiv_1_0 * (derxjm_200_ui * xji_11 + derxjm_110_ui * xji_20) -
              vderiv_1_1 * (derxjm_210_ui * xji_11 + derxjm_110_ui * xji_21) -
              vderiv_1_2 * (derxjm_220_ui * xji_11 + derxjm_110_ui * xji_22) -
              vderiv_2_0 * (xji_10 * derxjm_110_ui + xji_11 * derxjm_100_ui +
                               2 * xji_20 * derxjm_210_ui + 2 * xji_21 * derxjm_200_ui) -
              vderiv_2_1 * (xji_11 * derxjm_110_ui + xji_11 * derxjm_110_ui +
                               2 * xji_21 * derxjm_210_ui + 2 * xji_21 * derxjm_210_ui) -
              vderiv_2_2 * (xji_12 * derxjm_110_ui + xji_11 * derxjm_120_ui +
                               2 * xji_22 * derxjm_210_ui + 2 * xji_21 * derxjm_220_ui);
          const double v2 =
              -vderiv_0_0 * (derxjm_200_ui * xji_02) - vderiv_0_1 * (derxjm_210_ui * xji_02) -
              vderiv_0_2 * (derxjm_220_ui * xji_02) -
              vderiv_1_0 * (derxjm_200_ui * xji_12 + derxjm_120_ui * xji_20) -
              vderiv_1_1 * (derxjm_210_ui * xji_12 + derxjm_120_ui * xji_21) -
              vderiv_1_2 * (derxjm_220_ui * xji_12 + derxjm_120_ui * xji_22) -
              vderiv_2_0 * (xji_10 * derxjm_120_ui + xji_12 * derxjm_100_ui +
                               2 * xji_20 * derxjm_220_ui + 2 * xji_22 * derxjm_200_ui) -
              vderiv_2_1 * (xji_11 * derxjm_120_ui + xji_12 * derxjm_110_ui +
                               2 * xji_21 * derxjm_220_ui + 2 * xji_22 * derxjm_210_ui) -
              vderiv_2_2 * (xji_12 * derxjm_120_ui + xji_12 * derxjm_120_ui +
                               2 * xji_22 * derxjm_220_ui + 2 * xji_22 * derxjm_220_ui);

          for (int vi = 0; vi < nen_; ++vi)
          {
            ecoupl_u(vi * 3 + 2, ui * 3 + 0) +=
                v * (Base::deriv_(0, vi) * v0 + Base::deriv_(1, vi) * v1 +
                        Base::deriv_(2, vi) * v2) -
                v * Base::funct_(vi) * porosity_inv *
                    (refgrad_porosity_0 * v0 + refgrad_porosity_1 * v1 + refgrad_porosity_2 * v2);
          }
        }

        // //////////////////////////////////////////////////////////////
        {
          const double v0 =
              -vderiv_0_0 * (derxjm_201_ui * xji_00 + derxjm_001_ui * xji_20) -
              vderiv_0_1 * (derxjm_211_ui * xji_00 + derxjm_001_ui * xji_21) -
              vderiv_0_2 * (derxjm_221_ui * xji_00 + derxjm_001_ui * xji_22) -
              vderiv_1_0 * (derxjm_201_ui * xji_10) - vderiv_1_1 * (derxjm_211_ui * xji_10) -
              vderiv_1_2 * (derxjm_221_ui * xji_10) -
              vderiv_2_0 * (xji_00 * derxjm_001_ui + xji_00 * derxjm_001_ui +
                               2 * xji_20 * derxjm_201_ui + 2 * xji_20 * derxjm_201_ui) -
              vderiv_2_1 * (xji_01 * derxjm_001_ui + xji_00 * derxjm_011_ui +
                               2 * xji_21 * derxjm_201_ui + 2 * xji_20 * derxjm_211_ui) -
              vderiv_2_2 * (xji_02 * derxjm_001_ui + xji_00 * derxjm_021_ui +
                               2 * xji_22 * derxjm_201_ui + 2 * xji_20 * derxjm_221_ui);
          const double v1 =
              -vderiv_0_0 * (derxjm_201_ui * xji_01 + derxjm_011_ui * xji_20) -
              vderiv_0_1 * (derxjm_211_ui * xji_01 + derxjm_011_ui * xji_21) -
              vderiv_0_2 * (derxjm_221_ui * xji_01 + derxjm_011_ui * xji_22) -
              vderiv_1_0 * (derxjm_201_ui * xji_11) - vderiv_1_1 * (derxjm_211_ui * xji_11) -
              vderiv_1_2 * (derxjm_221_ui * xji_11) -
              vderiv_2_0 * (xji_00 * derxjm_011_ui + xji_01 * derxjm_001_ui +
                               2 * xji_20 * derxjm_211_ui + 2 * xji_21 * derxjm_201_ui) -
              vderiv_2_1 * (xji_01 * derxjm_011_ui + xji_01 * derxjm_011_ui +
                               2 * xji_21 * derxjm_211_ui + 2 * xji_21 * derxjm_211_ui) -
              vderiv_2_2 * (xji_02 * derxjm_011_ui + xji_01 * derxjm_021_ui +
                               2 * xji_22 * derxjm_211_ui + 2 * xji_21 * derxjm_221_ui);
          const double v2 =
              -vderiv_0_0 * (derxjm_201_ui * xji_02 + derxjm_021_ui * xji_20) -
              vderiv_0_1 * (derxjm_211_ui * xji_02 + derxjm_021_ui * xji_21) -
              vderiv_0_2 * (derxjm_221_ui * xji_02 + derxjm_021_ui * xji_22) -
              vderiv_1_0 * (derxjm_201_ui * xji_12) - vderiv_1_1 * (derxjm_211_ui * xji_12) -
              vderiv_1_2 * (derxjm_221_ui * xji_12) -
              vderiv_2_0 * (xji_00 * derxjm_021_ui + xji_02 * derxjm_001_ui +
                               2 * xji_20 * derxjm_221_ui + 2 * xji_22 * derxjm_201_ui) -
              vderiv_2_1 * (xji_01 * derxjm_021_ui + xji_02 * derxjm_011_ui +
                               2 * xji_21 * derxjm_221_ui + 2 * xji_22 * derxjm_211_ui) -
              vderiv_2_2 * (xji_02 * derxjm_021_ui + xji_02 * derxjm_021_ui +
                               2 * xji_22 * derxjm_221_ui + 2 * xji_22 * derxjm_221_ui);

          for (int vi = 0; vi < nen_; ++vi)
          {
            ecoupl_u(vi * 3 + 2, ui * 3 + 1) +=
                v * (Base::deriv_(0, vi) * v0 + Base::deriv_(1, vi) * v1 +
                        Base::deriv_(2, vi) * v2) -
                v * Base::funct_(vi) * porosity_inv *
                    (refgrad_porosity_0 * v0 + refgrad_porosity_1 * v1 + refgrad_porosity_2 * v2);
          }
        }

        // //////////////////////////////////////////////////////////////
        {
          const double v0 =
              -vderiv_0_0 * (derxjm_002_ui * xji_20) - vderiv_0_1 * (derxjm_002_ui * xji_21) -
              vderiv_0_2 * (derxjm_002_ui * xji_22) - vderiv_1_0 * (derxjm_102_ui * xji_20) -
              vderiv_1_1 * (derxjm_102_ui * xji_21) - vderiv_1_2 * (derxjm_102_ui * xji_22) -
              vderiv_2_0 * (xji_00 * derxjm_002_ui + xji_00 * derxjm_002_ui +
                               xji_10 * derxjm_102_ui + xji_10 * derxjm_102_ui) -
              vderiv_2_1 * (xji_01 * derxjm_002_ui + xji_00 * derxjm_012_ui +
                               xji_11 * derxjm_102_ui + xji_10 * derxjm_112_ui) -
              vderiv_2_2 * (xji_02 * derxjm_002_ui + xji_00 * derxjm_022_ui +
                               xji_12 * derxjm_102_ui + xji_10 * derxjm_122_ui);
          const double v1 =
              -vderiv_0_0 * (derxjm_012_ui * xji_20) - vderiv_0_1 * (derxjm_012_ui * xji_21) -
              vderiv_0_2 * (derxjm_012_ui * xji_22) - vderiv_1_0 * (derxjm_112_ui * xji_20) -
              vderiv_1_1 * (derxjm_112_ui * xji_21) - vderiv_1_2 * (derxjm_112_ui * xji_22) -
              vderiv_2_0 * (xji_00 * derxjm_012_ui + xji_01 * derxjm_002_ui +
                               xji_10 * derxjm_112_ui + xji_11 * derxjm_102_ui) -
              vderiv_2_1 * (xji_01 * derxjm_012_ui + xji_01 * derxjm_012_ui +
                               xji_11 * derxjm_112_ui + xji_11 * derxjm_112_ui) -
              vderiv_2_2 * (xji_02 * derxjm_012_ui + xji_01 * derxjm_022_ui +
                               xji_12 * derxjm_112_ui + xji_11 * derxjm_122_ui);
          const double v2 =
              -vderiv_0_0 * (derxjm_022_ui * xji_20) - vderiv_0_1 * (derxjm_022_ui * xji_21) -
              vderiv_0_2 * (derxjm_022_ui * xji_22) - vderiv_1_0 * (derxjm_122_ui * xji_20) -
              vderiv_1_1 * (derxjm_122_ui * xji_21) - vderiv_1_2 * (derxjm_122_ui * xji_22) -
              vderiv_2_0 * (xji_00 * derxjm_022_ui + xji_02 * derxjm_002_ui +
                               xji_10 * derxjm_122_ui + xji_12 * derxjm_102_ui) -
              vderiv_2_1 * (xji_01 * derxjm_022_ui + xji_02 * derxjm_012_ui +
                               xji_11 * derxjm_122_ui + xji_12 * derxjm_112_ui) -
              vderiv_2_2 * (xji_02 * derxjm_022_ui + xji_02 * derxjm_022_ui +
                               xji_12 * derxjm_122_ui + xji_12 * derxjm_122_ui);

          for (int vi = 0; vi < nen_; ++vi)
          {
            ecoupl_u(vi * 3 + 2, ui * 3 + 2) +=
                v * (Base::deriv_(0, vi) * v0 + Base::deriv_(1, vi) * v1 +
                        Base::deriv_(2, vi) * v2) -
                v * Base::funct_(vi) * porosity_inv *
                    (refgrad_porosity_0 * v0 + refgrad_porosity_1 * v1 + refgrad_porosity_2 * v2);
          }
        }
      }
    }
  }

  //*************************** ReacStab**********************************
  if (Base::fldpara_->RStab() != INPAR::FLUID::reactive_stab_none)
  {
    const double refgradp_0 = refgrad_press_(0);
    const double refgradp_1 = refgrad_press_(1);
    const double refgradp_2 = refgrad_press_(2);

    // pressure;
    for (int ui = 0; ui < nen_; ++ui)
    {
      const double v10 = refgradp_0 * derxjm_(0, 0, 1, ui) + refgradp_1 * derxjm_(0, 1, 1, ui) +
                         refgradp_2 * derxjm_(0, 2, 1, ui);
      const double v20 = refgradp_0 * derxjm_(0, 0, 2, ui) + refgradp_1 * derxjm_(0, 1, 2, ui) +
                         refgradp_2 * derxjm_(0, 2, 2, ui);

      const double v01 = refgradp_0 * derxjm_(1, 0, 0, ui) + refgradp_1 * derxjm_(1, 1, 0, ui) +
                         refgradp_2 * derxjm_(1, 2, 0, ui);
      const double v21 = refgradp_0 * derxjm_(1, 0, 2, ui) + refgradp_1 * derxjm_(1, 1, 2, ui) +
                         refgradp_2 * derxjm_(1, 2, 2, ui);

      const double v02 = refgradp_0 * derxjm_(2, 0, 0, ui) + refgradp_1 * derxjm_(2, 1, 0, ui) +
                         refgradp_2 * derxjm_(2, 2, 0, ui);
      const double v12 = refgradp_0 * derxjm_(2, 0, 1, ui) + refgradp_1 * derxjm_(2, 1, 1, ui) +
                         refgradp_2 * derxjm_(2, 2, 1, ui);

      for (int vi = 0; vi < nen_; ++vi)
      {
        const double v = Base::funct_(vi) * timefacfac_det * addstab;

        ecoupl_u(vi * 3 + 1, ui * 3) += v * v10;
        ecoupl_u(vi * 3 + 2, ui * 3) += v * v20;

        ecoupl_u(vi * 3 + 0, ui * 3 + 1) += v * v01;
        ecoupl_u(vi * 3 + 2, ui * 3 + 1) += v * v21;

        ecoupl_u(vi * 3 + 0, ui * 3 + 2) += v * v02;
        ecoupl_u(vi * 3 + 1, ui * 3 + 2) += v * v12;
      }
    }
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::LinMeshMotion_3D_Pres_OD(
    LINALG::Matrix<nen_, nsd_ * nen_>& ecoupl_p, const double& dphi_dp, const double& dphi_dJ,
    const double& refporositydot, const double& timefacfac)
{
  const double xjm_0_0 = Base::xjm_(0, 0);
  const double xjm_0_1 = Base::xjm_(0, 1);
  const double xjm_0_2 = Base::xjm_(0, 2);
  const double xjm_1_0 = Base::xjm_(1, 0);
  const double xjm_1_1 = Base::xjm_(1, 1);
  const double xjm_1_2 = Base::xjm_(1, 2);
  const double xjm_2_0 = Base::xjm_(2, 0);
  const double xjm_2_1 = Base::xjm_(2, 1);
  const double xjm_2_2 = Base::xjm_(2, 2);

  const double vderiv_0_0 = Base::vderiv_(0, 0);
  const double vderiv_0_1 = Base::vderiv_(0, 1);
  const double vderiv_0_2 = Base::vderiv_(0, 2);
  const double vderiv_1_0 = Base::vderiv_(1, 0);
  const double vderiv_1_1 = Base::vderiv_(1, 1);
  const double vderiv_1_2 = Base::vderiv_(1, 2);
  const double vderiv_2_0 = Base::vderiv_(2, 0);
  const double vderiv_2_1 = Base::vderiv_(2, 1);
  const double vderiv_2_2 = Base::vderiv_(2, 2);

  const double gridvelderiv_0_0 = gridvel_deriv_(0, 0);
  const double gridvelderiv_0_1 = gridvel_deriv_(0, 1);
  const double gridvelderiv_0_2 = gridvel_deriv_(0, 2);
  const double gridvelderiv_1_0 = gridvel_deriv_(1, 0);
  const double gridvelderiv_1_1 = gridvel_deriv_(1, 1);
  const double gridvelderiv_1_2 = gridvel_deriv_(1, 2);
  const double gridvelderiv_2_0 = gridvel_deriv_(2, 0);
  const double gridvelderiv_2_1 = gridvel_deriv_(2, 1);
  const double gridvelderiv_2_2 = gridvel_deriv_(2, 2);

  const double velint_0 = Base::velint_(0);
  const double velint_1 = Base::velint_(1);
  const double velint_2 = Base::velint_(2);

  const double gridvelint_0 = gridvel_int_(0);
  const double gridvelint_1 = gridvel_int_(1);
  const double gridvelint_2 = gridvel_int_(2);

  const double convvelint_0 = Base::convvelint_(0);
  const double convvelint_1 = Base::convvelint_(1);
  const double convvelint_2 = Base::convvelint_(2);

  //*************************** linearisation of mesh motion in continuity
  // equation**********************************

  const double timefacfac_det = timefacfac / Base::det_;

  if (!porofldpara_->PoroContiPartInt())
  {
    // (porosity_)*div u
    for (int ui = 0; ui < nen_; ++ui)
    {
      const double v0 = vderiv_1_0 * derxjm_(0, 0, 1, ui) + vderiv_1_1 * derxjm_(0, 1, 1, ui) +
                        vderiv_1_2 * derxjm_(0, 2, 1, ui) + vderiv_2_0 * derxjm_(0, 0, 2, ui) +
                        vderiv_2_1 * derxjm_(0, 1, 2, ui) + vderiv_2_2 * derxjm_(0, 2, 2, ui);

      const double v1 = vderiv_0_0 * derxjm_(1, 0, 0, ui) + vderiv_0_1 * derxjm_(1, 1, 0, ui) +
                        vderiv_0_2 * derxjm_(1, 2, 0, ui) + vderiv_2_0 * derxjm_(1, 0, 2, ui) +
                        vderiv_2_1 * derxjm_(1, 1, 2, ui) + vderiv_2_2 * derxjm_(1, 2, 2, ui);

      const double v2 = vderiv_0_0 * derxjm_(2, 0, 0, ui) + vderiv_0_1 * derxjm_(2, 1, 0, ui) +
                        vderiv_0_2 * derxjm_(2, 2, 0, ui) + vderiv_1_0 * derxjm_(2, 0, 1, ui) +
                        vderiv_1_1 * derxjm_(2, 1, 1, ui) + vderiv_1_2 * derxjm_(2, 2, 1, ui);

      for (int vi = 0; vi < nen_; ++vi)
      {
        const double v = timefacfac_det * Base::funct_(vi, 0) * porosity_;

        ecoupl_p(vi, ui * 3) += v * v0;

        ecoupl_p(vi, ui * 3 + 1) += v * v1;

        ecoupl_p(vi, ui * 3 + 2) += v * v2;
      }
    }

    if (!porofldpara_->IsStationaryConti())
    {
      // (dphi_dJ*J)*div vs
      for (int ui = 0; ui < nen_; ++ui)
      {
        const double v0 =
            gridvelderiv_1_0 * derxjm_(0, 0, 1, ui) + gridvelderiv_1_1 * derxjm_(0, 1, 1, ui) +
            gridvelderiv_1_2 * derxjm_(0, 2, 1, ui) + gridvelderiv_2_0 * derxjm_(0, 0, 2, ui) +
            gridvelderiv_2_1 * derxjm_(0, 1, 2, ui) + gridvelderiv_2_2 * derxjm_(0, 2, 2, ui);

        const double v1 =
            gridvelderiv_0_0 * derxjm_(1, 0, 0, ui) + gridvelderiv_0_1 * derxjm_(1, 1, 0, ui) +
            gridvelderiv_0_2 * derxjm_(1, 2, 0, ui) + gridvelderiv_2_0 * derxjm_(1, 0, 2, ui) +
            gridvelderiv_2_1 * derxjm_(1, 1, 2, ui) + gridvelderiv_2_2 * derxjm_(1, 2, 2, ui);

        const double v2 =
            gridvelderiv_0_0 * derxjm_(2, 0, 0, ui) + gridvelderiv_0_1 * derxjm_(2, 1, 0, ui) +
            gridvelderiv_0_2 * derxjm_(2, 2, 0, ui) + gridvelderiv_1_0 * derxjm_(2, 0, 1, ui) +
            gridvelderiv_1_1 * derxjm_(2, 1, 1, ui) + gridvelderiv_1_2 * derxjm_(2, 2, 1, ui);

        for (int vi = 0; vi < nen_; ++vi)
        {
          const double v = timefacfac_det * Base::funct_(vi, 0) * dphi_dJ * J_;

          ecoupl_p(vi, ui * 3 + 0) += v * v0;

          ecoupl_p(vi, ui * 3 + 1) += v * v1;

          ecoupl_p(vi, ui * 3 + 2) += v * v2;
        }
      }
    }

    //-----------(u-vs)grad(phi)

    const double refgrad_porosity_0 = refgrad_porosity_(0);
    const double refgrad_porosity_1 = refgrad_porosity_(1);
    const double refgrad_porosity_2 = refgrad_porosity_(2);

    for (int ui = 0; ui < nen_; ++ui)
    {
      const double v00 =
          +(velint_1 - gridvelint_1) * (refgrad_porosity_0 * derxjm_(0, 0, 1, ui) +
                                           refgrad_porosity_1 * derxjm_(0, 1, 1, ui) +
                                           refgrad_porosity_2 * derxjm_(0, 2, 1, ui)) +
          (velint_2 - gridvelint_2) * (refgrad_porosity_0 * derxjm_(0, 0, 2, ui) +
                                          refgrad_porosity_1 * derxjm_(0, 1, 2, ui) +
                                          refgrad_porosity_2 * derxjm_(0, 2, 2, ui));
      const double v01 =
          +(velint_0 - gridvelint_0) * (refgrad_porosity_0 * derxjm_(1, 0, 0, ui) +
                                           refgrad_porosity_1 * derxjm_(1, 1, 0, ui) +
                                           refgrad_porosity_2 * derxjm_(1, 2, 0, ui)) +
          (velint_2 - gridvelint_2) * (refgrad_porosity_0 * derxjm_(1, 0, 2, ui) +
                                          refgrad_porosity_1 * derxjm_(1, 1, 2, ui) +
                                          refgrad_porosity_2 * derxjm_(1, 2, 2, ui));
      const double v02 =
          +(velint_0 - gridvelint_0) * (refgrad_porosity_0 * derxjm_(2, 0, 0, ui) +
                                           refgrad_porosity_1 * derxjm_(2, 1, 0, ui) +
                                           refgrad_porosity_2 * derxjm_(2, 2, 0, ui)) +
          (velint_1 - gridvelint_1) * (refgrad_porosity_0 * derxjm_(2, 0, 1, ui) +
                                          refgrad_porosity_1 * derxjm_(2, 1, 1, ui) +
                                          refgrad_porosity_2 * derxjm_(2, 2, 1, ui));

      for (int vi = 0; vi < nen_; ++vi)
      {
        const double v = timefacfac_det * Base::funct_(vi);

        ecoupl_p(vi, ui * 3 + 0) += v * v00;
        ecoupl_p(vi, ui * 3 + 1) += v * v01;
        ecoupl_p(vi, ui * 3 + 2) += v * v02;
      }
    }
  }
  else
  {
    if (!porofldpara_->IsStationaryConti())
    {
      // (dphi_dJ*J+phi)*div vs
      for (int ui = 0; ui < nen_; ++ui)
      {
        const double v0 =
            gridvelderiv_1_0 * derxjm_(0, 0, 1, ui) + gridvelderiv_1_1 * derxjm_(0, 1, 1, ui) +
            gridvelderiv_1_2 * derxjm_(0, 2, 1, ui) + gridvelderiv_2_0 * derxjm_(0, 0, 2, ui) +
            gridvelderiv_2_1 * derxjm_(0, 1, 2, ui) + gridvelderiv_2_2 * derxjm_(0, 2, 2, ui);

        const double v1 =
            gridvelderiv_0_0 * derxjm_(1, 0, 0, ui) + gridvelderiv_0_1 * derxjm_(1, 1, 0, ui) +
            gridvelderiv_0_2 * derxjm_(1, 2, 0, ui) + gridvelderiv_2_0 * derxjm_(1, 0, 2, ui) +
            gridvelderiv_2_1 * derxjm_(1, 1, 2, ui) + gridvelderiv_2_2 * derxjm_(1, 2, 2, ui);

        const double v2 =
            gridvelderiv_0_0 * derxjm_(2, 0, 0, ui) + gridvelderiv_0_1 * derxjm_(2, 1, 0, ui) +
            gridvelderiv_0_2 * derxjm_(2, 2, 0, ui) + gridvelderiv_1_0 * derxjm_(2, 0, 1, ui) +
            gridvelderiv_1_1 * derxjm_(2, 1, 1, ui) + gridvelderiv_1_2 * derxjm_(2, 2, 1, ui);

        for (int vi = 0; vi < nen_; ++vi)
        {
          const double v = timefacfac_det * Base::funct_(vi, 0) * (dphi_dJ * J_ + porosity_);

          ecoupl_p(vi, ui * 3 + 0) += v * v0;

          ecoupl_p(vi, ui * 3 + 1) += v * v1;

          ecoupl_p(vi, ui * 3 + 2) += v * v2;
        }
      }
    }

    //----------- phi * (u-vs)grad(vi)
    const double v = -1.0 * timefacfac_det * porosity_;

    for (int vi = 0; vi < nen_; ++vi)
    {
      const double deriv_vi_0 = Base::deriv_(0, vi);
      const double deriv_vi_1 = Base::deriv_(1, vi);
      const double deriv_vi_2 = Base::deriv_(2, vi);

      for (int ui = 0; ui < nen_; ++ui)
      {
        const double v00 = +(velint_1 - gridvelint_1) * (deriv_vi_0 * derxjm_(0, 0, 1, ui) +
                                                            deriv_vi_1 * derxjm_(0, 1, 1, ui) +
                                                            deriv_vi_2 * derxjm_(0, 2, 1, ui)) +
                           (velint_2 - gridvelint_2) * (deriv_vi_0 * derxjm_(0, 0, 2, ui) +
                                                           deriv_vi_1 * derxjm_(0, 1, 2, ui) +
                                                           deriv_vi_2 * derxjm_(0, 2, 2, ui));
        const double v01 = +(velint_0 - gridvelint_0) * (deriv_vi_0 * derxjm_(1, 0, 0, ui) +
                                                            deriv_vi_1 * derxjm_(1, 1, 0, ui) +
                                                            deriv_vi_2 * derxjm_(1, 2, 0, ui)) +
                           (velint_2 - gridvelint_2) * (deriv_vi_0 * derxjm_(1, 0, 2, ui) +
                                                           deriv_vi_1 * derxjm_(1, 1, 2, ui) +
                                                           deriv_vi_2 * derxjm_(1, 2, 2, ui));
        const double v02 = +(velint_0 - gridvelint_0) * (deriv_vi_0 * derxjm_(2, 0, 0, ui) +
                                                            deriv_vi_1 * derxjm_(2, 1, 0, ui) +
                                                            deriv_vi_2 * derxjm_(2, 2, 0, ui)) +
                           (velint_1 - gridvelint_1) * (deriv_vi_0 * derxjm_(2, 0, 1, ui) +
                                                           deriv_vi_1 * derxjm_(2, 1, 1, ui) +
                                                           deriv_vi_2 * derxjm_(2, 2, 1, ui));

        ecoupl_p(vi, ui * 3 + 0) += v * v00;
        ecoupl_p(vi, ui * 3 + 1) += v * v01;
        ecoupl_p(vi, ui * 3 + 2) += v * v02;
      }
    }

  }  // partial integration

  if (!porofldpara_->IsStationaryConti())
  {
    // dphi_dp*dp/dt + rhs
    for (int vi = 0; vi < nen_; ++vi)
    {
      const double v = Base::fac_ * Base::funct_(vi, 0) * dphi_dp * press_ +
                       timefacfac * Base::funct_(vi, 0) * refporositydot;
      for (int ui = 0; ui < nen_; ++ui)
      {
        ecoupl_p(vi, ui * 3) += v * Base::derxy_(0, ui);
        ecoupl_p(vi, ui * 3 + 1) += v * Base::derxy_(1, ui);
        ecoupl_p(vi, ui * 3 + 2) += v * Base::derxy_(2, ui);
      }
    }
    //  rhs
    for (int vi = 0; vi < nen_; ++vi)
    {
      const double v = -1.0 * timefacfac * Base::funct_(vi, 0) * (dphi_dp * Base::rhscon_);
      for (int ui = 0; ui < nen_; ++ui)
      {
        ecoupl_p(vi, ui * 3) += v * Base::derxy_(0, ui);
        ecoupl_p(vi, ui * 3 + 1) += v * Base::derxy_(1, ui);
        ecoupl_p(vi, ui * 3 + 2) += v * Base::derxy_(2, ui);
      }
    }
  }

  //-------------------
  if (Base::fldpara_->PSPG())
  {
    // PSPG rhs
    {
      const double v = -1.0 * timefacfac_det;

      const double sgvelint_0 = Base::sgvelint_(0);
      const double sgvelint_1 = Base::sgvelint_(1);
      const double sgvelint_2 = Base::sgvelint_(2);

      for (int vi = 0; vi < nen_; ++vi)
      {
        const double deriv_vi_0 = Base::deriv_(0, vi);
        const double deriv_vi_1 = Base::deriv_(1, vi);
        const double deriv_vi_2 = Base::deriv_(2, vi);

        for (int ui = 0; ui < nen_; ++ui)
        {
          const double v00 =
              +sgvelint_1 * (deriv_vi_0 * derxjm_(0, 0, 1, ui) + deriv_vi_1 * derxjm_(0, 1, 1, ui) +
                                deriv_vi_2 * derxjm_(0, 2, 1, ui)) +
              sgvelint_2 * (deriv_vi_0 * derxjm_(0, 0, 2, ui) + deriv_vi_1 * derxjm_(0, 1, 2, ui) +
                               deriv_vi_2 * derxjm_(0, 2, 2, ui));
          const double v01 =
              +sgvelint_0 * (deriv_vi_0 * derxjm_(1, 0, 0, ui) + deriv_vi_1 * derxjm_(1, 1, 0, ui) +
                                deriv_vi_2 * derxjm_(1, 2, 0, ui)) +
              sgvelint_2 * (deriv_vi_0 * derxjm_(1, 0, 2, ui) + deriv_vi_1 * derxjm_(1, 1, 2, ui) +
                               deriv_vi_2 * derxjm_(1, 2, 2, ui));
          const double v02 =
              +sgvelint_0 * (deriv_vi_0 * derxjm_(2, 0, 0, ui) + deriv_vi_1 * derxjm_(2, 1, 0, ui) +
                                deriv_vi_2 * derxjm_(2, 2, 0, ui)) +
              sgvelint_1 * (deriv_vi_0 * derxjm_(2, 0, 1, ui) + deriv_vi_1 * derxjm_(2, 1, 1, ui) +
                               deriv_vi_2 * derxjm_(2, 2, 1, ui));

          ecoupl_p(vi, ui * 3 + 0) += v * v00;
          ecoupl_p(vi, ui * 3 + 1) += v * v01;
          ecoupl_p(vi, ui * 3 + 2) += v * v02;
        }
      }
    }

    double scal_grad_q = 0.0;

    if (Base::fldpara_->Tds() == INPAR::FLUID::subscales_quasistatic)
    {
      scal_grad_q = Base::tau_(1);
    }
    else
    {
      scal_grad_q = 0.0;  // Base::fldpara_->AlphaF()*fac3;
    }

    // pressure
    {
      const double v = timefacfac_det * scal_grad_q;

      const double refgradp_0 = refgrad_press_(0);
      const double refgradp_1 = refgrad_press_(1);
      const double refgradp_2 = refgrad_press_(2);

      // shape derivative of pressure gradient
      for (int vi = 0; vi < nen_; ++vi)
      {
        const double derxy_vi_0 = Base::derxy_(0, vi);
        const double derxy_vi_1 = Base::derxy_(1, vi);
        const double derxy_vi_2 = Base::derxy_(2, vi);

        for (int ui = 0; ui < nen_; ++ui)
        {
          const double v00 =
              +derxy_vi_1 * (refgradp_0 * derxjm_(0, 0, 1, ui) + refgradp_1 * derxjm_(0, 1, 1, ui) +
                                refgradp_2 * derxjm_(0, 2, 1, ui)) +
              derxy_vi_2 * (refgradp_0 * derxjm_(0, 0, 2, ui) + refgradp_1 * derxjm_(0, 1, 2, ui) +
                               refgradp_2 * derxjm_(0, 2, 2, ui));
          const double v01 =
              +derxy_vi_0 * (refgradp_0 * derxjm_(1, 0, 0, ui) + refgradp_1 * derxjm_(1, 1, 0, ui) +
                                refgradp_2 * derxjm_(1, 2, 0, ui)) +
              derxy_vi_2 * (refgradp_0 * derxjm_(1, 0, 2, ui) + refgradp_1 * derxjm_(1, 1, 2, ui) +
                               refgradp_2 * derxjm_(1, 2, 2, ui));
          const double v02 =
              +derxy_vi_0 * (refgradp_0 * derxjm_(2, 0, 0, ui) + refgradp_1 * derxjm_(2, 1, 0, ui) +
                                refgradp_2 * derxjm_(2, 2, 0, ui)) +
              derxy_vi_1 * (refgradp_0 * derxjm_(2, 0, 1, ui) + refgradp_1 * derxjm_(2, 1, 1, ui) +
                               refgradp_2 * derxjm_(2, 2, 1, ui));

          ecoupl_p(vi, ui * 3 + 0) += v * v00;
          ecoupl_p(vi, ui * 3 + 1) += v * v01;
          ecoupl_p(vi, ui * 3 + 2) += v * v02;
        }
      }

      // shape derivative of Jacobian
      static LINALG::Matrix<nen_, 1> temp;
      temp.MultiplyTN(Base::derxy_, Base::gradp_);
      for (int vi = 0; vi < nen_; ++vi)
      {
        const double v3 = -1.0 * timefacfac * scal_grad_q * temp(vi);
        for (int ui = 0; ui < nen_; ++ui)
        {
          ecoupl_p(vi, ui * 3) += v3 * Base::derxy_(0, ui);
          ecoupl_p(vi, ui * 3 + 1) += v3 * Base::derxy_(1, ui);
          ecoupl_p(vi, ui * 3 + 2) += v3 * Base::derxy_(2, ui);
        }
      }
    }

    // convective term
    {
      const double v = Base::densaf_ * timefacfac_det * scal_grad_q;
      for (int vi = 0; vi < nen_; ++vi)
      {
        const double derxy_vi_0 = Base::derxy_(0, vi);
        const double derxy_vi_1 = Base::derxy_(1, vi);
        const double derxy_vi_2 = Base::derxy_(2, vi);

        for (int ui = 0; ui < nen_; ++ui)
        {
          const double v00 =
              +derxy_vi_1 * convvelint_1 *
                  (vderiv_0_0 * derxjm_(0, 0, 1, ui) + vderiv_0_1 * derxjm_(0, 1, 1, ui) +
                      vderiv_0_2 * derxjm_(0, 2, 1, ui)) +
              derxy_vi_2 * convvelint_2 *
                  (vderiv_0_0 * derxjm_(0, 0, 2, ui) + vderiv_0_1 * derxjm_(0, 1, 2, ui) +
                      vderiv_0_2 * derxjm_(0, 2, 2, ui));
          const double v10 =
              +derxy_vi_1 * convvelint_1 *
                  (vderiv_1_0 * derxjm_(0, 0, 1, ui) + vderiv_1_1 * derxjm_(0, 1, 1, ui) +
                      vderiv_1_2 * derxjm_(0, 2, 1, ui)) +
              derxy_vi_2 * convvelint_2 *
                  (vderiv_1_0 * derxjm_(0, 0, 2, ui) + vderiv_1_1 * derxjm_(0, 1, 2, ui) +
                      vderiv_1_2 * derxjm_(0, 2, 2, ui));
          const double v20 =
              +derxy_vi_1 * convvelint_1 *
                  (vderiv_2_0 * derxjm_(0, 0, 1, ui) + vderiv_2_1 * derxjm_(0, 1, 1, ui) +
                      vderiv_2_2 * derxjm_(0, 2, 1, ui)) +
              derxy_vi_2 * convvelint_2 *
                  (vderiv_2_0 * derxjm_(0, 0, 2, ui) + vderiv_2_1 * derxjm_(0, 1, 2, ui) +
                      vderiv_2_2 * derxjm_(0, 2, 2, ui));
          const double v01 =
              +derxy_vi_0 * convvelint_0 *
                  (vderiv_0_0 * derxjm_(1, 0, 0, ui) + vderiv_0_1 * derxjm_(1, 1, 0, ui) +
                      vderiv_0_2 * derxjm_(1, 2, 0, ui)) +
              derxy_vi_2 * convvelint_2 *
                  (vderiv_0_0 * derxjm_(1, 0, 2, ui) + vderiv_0_1 * derxjm_(1, 1, 2, ui) +
                      vderiv_0_2 * derxjm_(1, 2, 2, ui));
          const double v11 =
              +derxy_vi_0 * convvelint_0 *
                  (vderiv_1_0 * derxjm_(1, 0, 0, ui) + vderiv_1_1 * derxjm_(1, 1, 0, ui) +
                      vderiv_1_2 * derxjm_(1, 2, 0, ui)) +
              derxy_vi_2 * convvelint_2 *
                  (vderiv_1_0 * derxjm_(1, 0, 2, ui) + vderiv_1_1 * derxjm_(1, 1, 2, ui) +
                      vderiv_1_2 * derxjm_(1, 2, 2, ui));
          const double v21 =
              +derxy_vi_0 * convvelint_0 *
                  (vderiv_2_0 * derxjm_(1, 0, 0, ui) + vderiv_2_1 * derxjm_(1, 1, 0, ui) +
                      vderiv_2_2 * derxjm_(1, 2, 0, ui)) +
              derxy_vi_2 * convvelint_2 *
                  (vderiv_2_0 * derxjm_(1, 0, 2, ui) + vderiv_2_1 * derxjm_(1, 1, 2, ui) +
                      vderiv_2_2 * derxjm_(1, 2, 2, ui));
          const double v02 =
              +derxy_vi_0 * convvelint_0 *
                  (vderiv_0_0 * derxjm_(2, 0, 0, ui) + vderiv_0_1 * derxjm_(2, 1, 0, ui) +
                      vderiv_0_2 * derxjm_(2, 2, 0, ui)) +
              derxy_vi_1 * convvelint_1 *
                  (vderiv_0_0 * derxjm_(2, 0, 1, ui) + vderiv_0_1 * derxjm_(2, 1, 1, ui) +
                      vderiv_0_2 * derxjm_(2, 2, 1, ui));
          const double v12 =
              +derxy_vi_0 * convvelint_0 *
                  (vderiv_1_0 * derxjm_(2, 0, 0, ui) + vderiv_1_1 * derxjm_(2, 1, 0, ui) +
                      vderiv_1_2 * derxjm_(2, 2, 0, ui)) +
              derxy_vi_1 * convvelint_1 *
                  (vderiv_1_0 * derxjm_(2, 0, 1, ui) + vderiv_1_1 * derxjm_(2, 1, 1, ui) +
                      vderiv_1_2 * derxjm_(2, 2, 1, ui));
          const double v22 =
              +derxy_vi_0 * convvelint_0 *
                  (vderiv_2_0 * derxjm_(2, 0, 0, ui) + vderiv_2_1 * derxjm_(2, 1, 0, ui) +
                      vderiv_2_2 * derxjm_(2, 2, 0, ui)) +
              derxy_vi_1 * convvelint_1 *
                  (vderiv_2_0 * derxjm_(2, 0, 1, ui) + vderiv_2_1 * derxjm_(2, 1, 1, ui) +
                      vderiv_2_2 * derxjm_(2, 2, 1, ui));

          ecoupl_p(vi, ui * 3 + 0) += v * (v00 + v10 + v20);
          ecoupl_p(vi, ui * 3 + 1) += v * (v01 + v11 + v21);
          ecoupl_p(vi, ui * 3 + 2) += v * (v02 + v12 + v22);
        }
      }

      // shape derivative of Jacobian
      static LINALG::Matrix<nen_, 1> temp;
      temp.MultiplyTN(Base::derxy_, Base::conv_old_);
      for (int vi = 0; vi < nen_; ++vi)
      {
        const double v3 = -1.0 * Base::densaf_ * timefacfac * scal_grad_q * temp(vi);
        for (int ui = 0; ui < nen_; ++ui)
        {
          ecoupl_p(vi, ui * 3) += v3 * Base::derxy_(0, ui);
          ecoupl_p(vi, ui * 3 + 1) += v3 * Base::derxy_(1, ui);
          ecoupl_p(vi, ui * 3 + 2) += v3 * Base::derxy_(2, ui);
        }
      }
    }
  }

  if (porofldpara_->StabBiot() and (not porofldpara_->IsStationaryConti()) and
      struct_mat_->PoroLawType() != INPAR::MAT::m_poro_law_constant)
  {
    // shape derivative of Jacobian
    {
      const double mixres_0 = mixres_(0);
      const double mixres_1 = mixres_(1);
      const double mixres_2 = mixres_(2);
      const double v = timefacfac * tau_struct_;
      for (int vi = 0; vi < nen_; ++vi)
      {
        double w =
            v * (N_XYZ_(0, vi) * mixres_0 + N_XYZ_(1, vi) * mixres_1 + N_XYZ_(2, vi) * mixres_2);
        for (int ui = 0; ui < nen_; ++ui)
        {
          ecoupl_p(vi, ui * 3) += w * Base::derxy_(0, ui);
          ecoupl_p(vi, ui * 3 + 1) += w * Base::derxy_(1, ui);
          ecoupl_p(vi, ui * 3 + 2) += w * Base::derxy_(2, ui);
        }
      }
    }

    // shape derivative of gradient of test function
    {
      const double gradp_0 = tau_struct_ * timefacfac_det * Base::gradp_(0);
      const double gradp_1 = tau_struct_ * timefacfac_det * Base::gradp_(1);
      const double gradp_2 = tau_struct_ * timefacfac_det * Base::gradp_(2);

      for (int vi = 0; vi < nen_; ++vi)
      {
        const double deriv_vi_0 = Base::deriv_(0, vi);
        const double deriv_vi_1 = Base::deriv_(1, vi);
        const double deriv_vi_2 = Base::deriv_(2, vi);

        for (int ui = 0; ui < nen_; ++ui)
        {
          const double v00 =
              +gradp_1 * (deriv_vi_0 * derxjm_(0, 0, 1, ui) + deriv_vi_1 * derxjm_(0, 1, 1, ui) +
                             deriv_vi_2 * derxjm_(0, 2, 1, ui)) +
              gradp_2 * (deriv_vi_0 * derxjm_(0, 0, 2, ui) + deriv_vi_1 * derxjm_(0, 1, 2, ui) +
                            deriv_vi_2 * derxjm_(0, 2, 2, ui));
          const double v01 =
              +gradp_0 * (deriv_vi_0 * derxjm_(1, 0, 0, ui) + deriv_vi_1 * derxjm_(1, 1, 0, ui) +
                             deriv_vi_2 * derxjm_(1, 2, 0, ui)) +
              gradp_2 * (deriv_vi_0 * derxjm_(1, 0, 2, ui) + deriv_vi_1 * derxjm_(1, 1, 2, ui) +
                            deriv_vi_2 * derxjm_(1, 2, 2, ui));
          const double v02 =
              +gradp_0 * (deriv_vi_0 * derxjm_(2, 0, 0, ui) + deriv_vi_1 * derxjm_(2, 1, 0, ui) +
                             deriv_vi_2 * derxjm_(2, 2, 0, ui)) +
              gradp_1 * (deriv_vi_0 * derxjm_(2, 0, 1, ui) + deriv_vi_1 * derxjm_(2, 1, 1, ui) +
                            deriv_vi_2 * derxjm_(2, 2, 1, ui));

          ecoupl_p(vi, ui * 3 + 0) += v00;
          ecoupl_p(vi, ui * 3 + 1) += v01;
          ecoupl_p(vi, ui * 3 + 2) += v02;
        }
      }
    }

    const double refgradp_0 = refgrad_press_(0);
    const double refgradp_1 = refgrad_press_(1);
    const double refgradp_2 = refgrad_press_(2);

    // pressure
    {
      const double v = timefacfac_det * tau_struct_;

      // shape derivative of pressure gradient
      for (int vi = 0; vi < nen_; ++vi)
      {
        const double derxy_vi_0 = Base::derxy_(0, vi);
        const double derxy_vi_1 = Base::derxy_(1, vi);
        const double derxy_vi_2 = Base::derxy_(2, vi);

        for (int ui = 0; ui < nen_; ++ui)
        {
          const double v00 =
              +derxy_vi_1 * (refgradp_0 * derxjm_(0, 0, 1, ui) + refgradp_1 * derxjm_(0, 1, 1, ui) +
                                refgradp_2 * derxjm_(0, 2, 1, ui)) +
              derxy_vi_2 * (refgradp_0 * derxjm_(0, 0, 2, ui) + refgradp_1 * derxjm_(0, 1, 2, ui) +
                               refgradp_2 * derxjm_(0, 2, 2, ui));
          const double v01 =
              +derxy_vi_0 * (refgradp_0 * derxjm_(1, 0, 0, ui) + refgradp_1 * derxjm_(1, 1, 0, ui) +
                                refgradp_2 * derxjm_(1, 2, 0, ui)) +
              derxy_vi_2 * (refgradp_0 * derxjm_(1, 0, 2, ui) + refgradp_1 * derxjm_(1, 1, 2, ui) +
                               refgradp_2 * derxjm_(1, 2, 2, ui));
          const double v02 =
              +derxy_vi_0 * (refgradp_0 * derxjm_(2, 0, 0, ui) + refgradp_1 * derxjm_(2, 1, 0, ui) +
                                refgradp_2 * derxjm_(2, 2, 0, ui)) +
              derxy_vi_1 * (refgradp_0 * derxjm_(2, 0, 1, ui) + refgradp_1 * derxjm_(2, 1, 1, ui) +
                               refgradp_2 * derxjm_(2, 2, 1, ui));

          ecoupl_p(vi, ui * 3 + 0) += v * v00;
          ecoupl_p(vi, ui * 3 + 1) += v * v01;
          ecoupl_p(vi, ui * 3 + 2) += v * v02;
        }
      }

      // shape derivative of Jacobian
      static LINALG::Matrix<nen_, 1> temp;
      temp.MultiplyTN(Base::derxy_, Base::gradp_);
      for (int vi = 0; vi < nen_; ++vi)
      {
        const double v3 = -1.0 * timefacfac * tau_struct_ * temp(vi);
        for (int ui = 0; ui < nen_; ++ui)
        {
          ecoupl_p(vi, ui * 3) += v3 * Base::derxy_(0, ui);
          ecoupl_p(vi, ui * 3 + 1) += v3 * Base::derxy_(1, ui);
          ecoupl_p(vi, ui * 3 + 2) += v3 * Base::derxy_(2, ui);
        }
      }
    }

    // shape derivative of pressure gradient
    {
      const double v = timefacfac_det * (-1.0) * porosity_ * tau_struct_;

      for (int vi = 0; vi < nen_; ++vi)
      {
        const double N_XYZ__vi_0 = N_XYZ_(0, vi);
        const double N_XYZ__vi_1 = N_XYZ_(1, vi);
        const double N_XYZ__vi_2 = N_XYZ_(2, vi);

        for (int ui = 0; ui < nen_; ++ui)
        {
          const double v00 =
              +N_XYZ__vi_1 *
                  (refgradp_0 * derxjm_(0, 0, 1, ui) + refgradp_1 * derxjm_(0, 1, 1, ui) +
                      refgradp_2 * derxjm_(0, 2, 1, ui)) +
              N_XYZ__vi_2 * (refgradp_0 * derxjm_(0, 0, 2, ui) + refgradp_1 * derxjm_(0, 1, 2, ui) +
                                refgradp_2 * derxjm_(0, 2, 2, ui));
          const double v01 =
              +N_XYZ__vi_0 *
                  (refgradp_0 * derxjm_(1, 0, 0, ui) + refgradp_1 * derxjm_(1, 1, 0, ui) +
                      refgradp_2 * derxjm_(1, 2, 0, ui)) +
              N_XYZ__vi_2 * (refgradp_0 * derxjm_(1, 0, 2, ui) + refgradp_1 * derxjm_(1, 1, 2, ui) +
                                refgradp_2 * derxjm_(1, 2, 2, ui));
          const double v02 =
              +N_XYZ__vi_0 *
                  (refgradp_0 * derxjm_(2, 0, 0, ui) + refgradp_1 * derxjm_(2, 1, 0, ui) +
                      refgradp_2 * derxjm_(2, 2, 0, ui)) +
              N_XYZ__vi_1 * (refgradp_0 * derxjm_(2, 0, 1, ui) + refgradp_1 * derxjm_(2, 1, 1, ui) +
                                refgradp_2 * derxjm_(2, 2, 1, ui));

          ecoupl_p(vi, ui * 3 + 0) += v * v00;
          ecoupl_p(vi, ui * 3 + 1) += v * v01;
          ecoupl_p(vi, ui * 3 + 2) += v * v02;
        }
      }

      // shape derivative of Jacobian
      static LINALG::Matrix<nen_, 1> temp;
      temp.MultiplyTN(N_XYZ_, Base::gradp_);
      for (int vi = 0; vi < nen_; ++vi)
      {
        const double v3 = timefacfac * porosity_ * tau_struct_ * temp(vi);
        for (int ui = 0; ui < nen_; ++ui)
        {
          ecoupl_p(vi, ui * 3) += v3 * Base::derxy_(0, ui);
          ecoupl_p(vi, ui * 3 + 1) += v3 * Base::derxy_(1, ui);
          ecoupl_p(vi, ui * 3 + 2) += v3 * Base::derxy_(2, ui);
        }
      }
    }
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::Lin2DMeshMotionOD(
    LINALG::Matrix<nsd_ * nen_, nsd_ * nen_>& ecoupl_u, const double& dphi_dp,
    const double& dphi_dJ, const double& refporositydot, const double& timefac,
    const double& timefacfac)
{
  double addstab = 0.0;
  if (Base::fldpara_->RStab() != INPAR::FLUID::reactive_stab_none)
  {
    if (Base::fldpara_->Tds() == INPAR::FLUID::subscales_quasistatic)
      addstab = Base::fldpara_->ViscReaStabFac() * Base::reacoeff_ * Base::tau_(1);
    else
    {
      dserror("Is this factor correct? Check for bugs!");
      // addstab = Base::fldpara_->ViscReaStabFac()*Base::reacoeff_*Base::fldpara_->AlphaF()*fac3;
    }
  }

  //*************************** linearisation of mesh motion in momentum
  // balance**********************************
  // mass
  if (!porofldpara_->IsStationaryMomentum())
  {
    const double fac0 = Base::fac_ * Base::densam_ * (1.0 + addstab) * Base::velint_(0);
    const double fac1 = Base::fac_ * Base::densam_ * (1.0 + addstab) * Base::velint_(1);

    for (int vi = 0; vi < nen_; ++vi)
    {
      const double funct_vi = Base::funct_(vi, 0);
      const double v0 = fac0 * funct_vi;
      const double v1 = fac1 * funct_vi;
      for (int ui = 0; ui < nen_; ++ui)
      {
        ecoupl_u(vi * 2, ui * 2) += v0 * Base::derxy_(0, ui);
        ecoupl_u(vi * 2, ui * 2 + 1) += v0 * Base::derxy_(1, ui);

        ecoupl_u(vi * 2 + 1, ui * 2) += v1 * Base::derxy_(0, ui);
        ecoupl_u(vi * 2 + 1, ui * 2 + 1) += v1 * Base::derxy_(1, ui);
      }
    }
  }

  // rhs
  {
    const double fac_timint =
        Base::fac_ * Base::fldparatimint_->Dt() * Base::fldparatimint_->Theta();
    const double fac0 = fac_timint * (-1.0) * Base::rhsmom_(0);
    const double fac1 = fac_timint * (-1.0) * Base::rhsmom_(1);
    for (int vi = 0; vi < nen_; ++vi)
    {
      const double funct_vi = Base::funct_(vi, 0);
      const double v0 = fac0 * funct_vi;
      const double v1 = fac1 * funct_vi;
      for (int ui = 0; ui < nen_; ++ui)
      {
        ecoupl_u(vi * 2, ui * 2) += v0 * Base::derxy_(0, ui);
        ecoupl_u(vi * 2, ui * 2 + 1) += v0 * Base::derxy_(1, ui);

        ecoupl_u(vi * 2 + 1, ui * 2) += v1 * Base::derxy_(0, ui);
        ecoupl_u(vi * 2 + 1, ui * 2 + 1) += v1 * Base::derxy_(1, ui);
      }
    }
  }

  //---------reaction term (darcy term)
  {
    const double fac_reac_conv_vel_0 = reac_tensor_convvel_(0) * timefacfac * (1.0 + addstab);
    const double fac_reac_conv_vel_1 = reac_tensor_convvel_(1) * timefacfac * (1.0 + addstab);
    for (int vi = 0; vi < nen_; ++vi)
    {
      const double funct_vi = Base::funct_(vi, 0);
      const double v0 = fac_reac_conv_vel_0 * funct_vi;
      const double v1 = fac_reac_conv_vel_1 * funct_vi;
      for (int ui = 0; ui < nen_; ++ui)
      {
        ecoupl_u(vi * 2, ui * 2) += v0 * Base::derxy_(0, ui);
        ecoupl_u(vi * 2, ui * 2 + 1) += v0 * Base::derxy_(1, ui);

        ecoupl_u(vi * 2 + 1, ui * 2) += v1 * Base::derxy_(0, ui);
        ecoupl_u(vi * 2 + 1, ui * 2 + 1) += v1 * Base::derxy_(1, ui);
      }
    }
  }

  const double vderiv_0_0 = Base::vderiv_(0, 0);
  const double vderiv_0_1 = Base::vderiv_(0, 1);
  const double vderiv_1_0 = Base::vderiv_(1, 0);
  const double vderiv_1_1 = Base::vderiv_(1, 1);
  //---------------convective term
  {
    const double convvelint_0 = Base::convvelint_(0);
    const double convvelint_1 = Base::convvelint_(1);
    for (int vi = 0; vi < nen_; ++vi)
    {
      const int tvi = 2 * vi;
      const int tvip = tvi + 1;
      const double v = Base::densaf_ * timefacfac / Base::det_ * Base::funct_(vi) * (1.0 + addstab);
      for (int ui = 0; ui < nen_; ++ui)
      {
        const int tui = 2 * ui;
        const int tuip = tui + 1;

        ecoupl_u(tvi, tui) += v * (+convvelint_1 * (-vderiv_0_0 * Base::deriv_(1, ui) +
                                                       vderiv_0_1 * Base::deriv_(0, ui)));

        ecoupl_u(tvi, tuip) +=
            v *
            (+convvelint_0 * (vderiv_0_0 * Base::deriv_(1, ui) - vderiv_0_1 * Base::deriv_(0, ui)));

        ecoupl_u(tvip, tui) += v * (+convvelint_1 * (-vderiv_1_0 * Base::deriv_(1, ui) +
                                                        vderiv_1_1 * Base::deriv_(0, ui)));

        ecoupl_u(tvip, tuip) +=
            v *
            (+convvelint_0 * (vderiv_1_0 * Base::deriv_(1, ui) - vderiv_1_1 * Base::deriv_(0, ui)));
      }
    }
  }

  // pressure
  for (int vi = 0; vi < nen_; ++vi)
  {
    const int tvi = 2 * vi;
    const int tvip = tvi + 1;
    const double v = press_ * timefacfac / Base::det_;
    const double deriv_0_vi = Base::deriv_(0, vi);
    const double deriv_1_vi = Base::deriv_(1, vi);
    for (int ui = 0; ui < nen_; ++ui)
    {
      const int tui = 2 * ui;
      ecoupl_u(tvi, tui + 1) -=
          v * (deriv_0_vi * Base::deriv_(1, ui) - Base::deriv_(0, ui) * deriv_1_vi);
      ecoupl_u(tvip, tui) -=
          v * (-deriv_0_vi * Base::deriv_(1, ui) + Base::deriv_(0, ui) * deriv_1_vi);
    }
  }

  // //---------viscous term (brinkman term)
  const double xji_00 = Base::xji_(0, 0);
  const double xji_01 = Base::xji_(0, 1);
  const double xji_10 = Base::xji_(1, 0);
  const double xji_11 = Base::xji_(1, 1);

  if (Base::visceff_)
  {
    const double vderxy_0_0 = 2.0 * Base::vderxy_(0, 0);
    const double vderxy_1_1 = 2.0 * Base::vderxy_(1, 1);
    const double vderxy_0_1 = Base::vderxy_(0, 1) + Base::vderxy_(1, 0);

    const double refgrad_porosity_0 = refgrad_porosity_(0);
    const double refgrad_porosity_1 = refgrad_porosity_(1);

    // part 1: derivative of det
    {
      const double v = Base::visceff_ * timefac * Base::fac_ * (1.0 + addstab);
      for (int ui = 0; ui < nen_; ++ui)
      {
        const double derinvJ0 = -v * (Base::deriv_(0, ui) * xji_00 + Base::deriv_(1, ui) * xji_01);
        const double derinvJ1 = -v * (Base::deriv_(0, ui) * xji_10 + Base::deriv_(1, ui) * xji_11);
        for (int vi = 0; vi < nen_; ++vi)
        {
          const double visres0 =
              Base::derxy_(0, vi) * vderxy_0_0 + Base::derxy_(1, vi) * vderxy_0_1;
          const double visres1 =
              Base::derxy_(0, vi) * vderxy_0_1 + Base::derxy_(1, vi) * vderxy_1_1;

          ecoupl_u(vi * 2, ui * 2) += derinvJ0 * visres0;
          ecoupl_u(vi * 2 + 1, ui * 2) += derinvJ0 * visres1;

          ecoupl_u(vi * 2, ui * 2 + 1) += derinvJ1 * visres0;
          ecoupl_u(vi * 2 + 1, ui * 2 + 1) += derinvJ1 * visres1;

          const double visres0_poro = refgrad_porosity_0 * Base::funct_(vi) * vderxy_0_0 +
                                      refgrad_porosity_1 * Base::funct_(vi) * vderxy_0_1;
          const double visres1_poro = refgrad_porosity_0 * Base::funct_(vi) * vderxy_0_1 +
                                      refgrad_porosity_1 * Base::funct_(vi) * vderxy_1_1;

          ecoupl_u(vi * 2 + 0, ui * 2 + 0) += -1.0 * derinvJ0 / porosity_ * visres0_poro;
          ecoupl_u(vi * 2 + 1, ui * 2 + 0) += -1.0 * derinvJ0 / porosity_ * visres1_poro;

          ecoupl_u(vi * 2 + 0, ui * 2 + 1) += -1.0 * derinvJ1 / porosity_ * visres0_poro;
          ecoupl_u(vi * 2 + 1, ui * 2 + 1) += -1.0 * derinvJ1 / porosity_ * visres1_poro;
        }
      }
    }

    // part 2: derivative of viscosity residual
    {
      const double v = timefacfac * Base::visceff_ / Base::det_ * (1.0 + addstab);
      for (int ui = 0; ui < nen_; ++ui)
      {
        const double deriv_0_ui = Base::deriv_(0, ui);
        const double deriv_1_ui = Base::deriv_(1, ui);
        {
          const double v0 = -vderiv_0_0 * (xji_10 * deriv_1_ui + xji_10 * deriv_1_ui) -
                            vderiv_0_1 * (xji_11 * deriv_1_ui + xji_10 * deriv_0_ui) -
                            vderiv_1_0 * (deriv_1_ui * xji_00) - vderiv_1_1 * (deriv_1_ui * xji_01);
          const double v1 = -vderiv_0_0 * (xji_10 * deriv_0_ui + xji_11 * deriv_1_ui) -
                            vderiv_0_1 * (xji_11 * deriv_0_ui + xji_11 * deriv_0_ui) -
                            vderiv_1_0 * (deriv_0_ui * xji_00) - vderiv_1_1 * (deriv_0_ui * xji_01);

          for (int vi = 0; vi < nen_; ++vi)
          {
            ecoupl_u(vi * 2 + 0, ui * 2 + 0) +=
                v * (Base::deriv_(0, vi) * v0 + Base::deriv_(1, vi) * v1) -
                v * Base::funct_(vi) / porosity_ *
                    (refgrad_porosity_0 * v0 + refgrad_porosity_1 * v1);
          }
        }

        // //////////////////////////////////////////////////////////////
        {
          const double v0 = -vderiv_0_0 * (2 * deriv_1_ui * xji_00 + 2 * deriv_1_ui * xji_00) -
                            vderiv_0_1 * (2 * deriv_0_ui * xji_00 + 2 * deriv_1_ui * xji_01) -
                            vderiv_1_0 * (deriv_1_ui * xji_10) - vderiv_1_1 * (deriv_0_ui * xji_10);
          const double v1 = -vderiv_0_0 * (2 * deriv_0_ui * xji_00 + 2 * deriv_1_ui * xji_01) -
                            vderiv_0_1 * (2 * deriv_0_ui * xji_01 + 2 * deriv_0_ui * xji_01) -
                            vderiv_1_0 * (deriv_1_ui * xji_11) - vderiv_1_1 * (deriv_0_ui * xji_11);

          for (int vi = 0; vi < nen_; ++vi)
          {
            ecoupl_u(vi * 2 + 0, ui * 2 + 1) +=
                v * (Base::deriv_(0, vi) * v0 + Base::deriv_(1, vi) * v1) -
                v * Base::funct_(vi) / porosity_ *
                    (refgrad_porosity_0 * v0 + refgrad_porosity_1 * v1);
          }
        }

        // //////////////////////////////////////////////////////////////
        {
          const double v0 = -vderiv_0_0 * (deriv_1_ui * xji_00) -
                            vderiv_0_1 * (deriv_0_ui * xji_00) -
                            vderiv_1_0 * (2 * xji_10 * deriv_1_ui + 2 * xji_10 * deriv_1_ui) -
                            vderiv_1_1 * (2 * xji_11 * deriv_1_ui + 2 * xji_10 * deriv_0_ui);
          const double v1 = -vderiv_0_0 * (deriv_1_ui * xji_01) -
                            vderiv_0_1 * (deriv_0_ui * xji_01) -
                            vderiv_1_0 * (2 * xji_10 * deriv_0_ui + 2 * xji_11 * deriv_1_ui) -
                            vderiv_1_1 * (2 * xji_11 * deriv_0_ui + 2 * xji_11 * deriv_0_ui);

          for (int vi = 0; vi < nen_; ++vi)
          {
            ecoupl_u(vi * 2 + 1, ui * 2 + 0) +=
                v * (Base::deriv_(0, vi) * v0 + Base::deriv_(1, vi) * v1) -
                v * Base::funct_(vi) / porosity_ *
                    (refgrad_porosity_0 * v0 + refgrad_porosity_1 * v1);
          }
        }

        // //////////////////////////////////////////////////////////////
        {
          const double v0 = -vderiv_0_0 * (deriv_1_ui * xji_10) -
                            vderiv_0_1 * (deriv_1_ui * xji_11) -
                            vderiv_1_0 * (xji_00 * deriv_1_ui + xji_00 * deriv_1_ui) -
                            vderiv_1_1 * (xji_01 * deriv_1_ui + xji_00 * deriv_0_ui);
          const double v1 = -vderiv_0_0 * (deriv_0_ui * xji_10) -
                            vderiv_0_1 * (deriv_0_ui * xji_11) -
                            vderiv_1_0 * (xji_00 * deriv_0_ui + xji_01 * deriv_1_ui) -
                            vderiv_1_1 * (xji_01 * deriv_0_ui + xji_01 * deriv_0_ui);

          for (int vi = 0; vi < nen_; ++vi)
          {
            ecoupl_u(vi * 2 + 1, ui * 2 + 1) +=
                v * (Base::deriv_(0, vi) * v0 + Base::deriv_(1, vi) * v1) -
                v * Base::funct_(vi) / porosity_ *
                    (refgrad_porosity_0 * v0 + refgrad_porosity_1 * v1);
          }
        }
      }
    }
  }

  //*************************** ReacStab**********************************
  if (Base::fldpara_->RStab() != INPAR::FLUID::reactive_stab_none)
  {
    const double refgradp_0 = refgrad_press_(0);
    const double refgradp_1 = refgrad_press_(1);
    // pressure;
    for (int vi = 0; vi < nen_; ++vi)
    {
      double v = Base::funct_(vi) * timefacfac / Base::det_ * addstab;
      for (int ui = 0; ui < nen_; ++ui)
      {
        ecoupl_u(vi * 2 + 1, ui * 2) +=
            v * (-refgradp_0 * Base::deriv_(1, ui) + refgradp_1 * Base::deriv_(0, ui));

        ecoupl_u(vi * 2 + 0, ui * 2 + 1) +=
            v * (refgradp_0 * Base::deriv_(1, ui) - refgradp_1 * Base::deriv_(0, ui));
      }
    }
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::LinMeshMotion_2D_Pres_OD(
    LINALG::Matrix<nen_, nsd_ * nen_>& ecoupl_p, const double& dphi_dp, const double& dphi_dJ,
    const double& refporositydot, const double& timefacfac)
{
  const double vderiv_0_0 = Base::vderiv_(0, 0);
  const double vderiv_0_1 = Base::vderiv_(0, 1);
  const double vderiv_1_0 = Base::vderiv_(1, 0);
  const double vderiv_1_1 = Base::vderiv_(1, 1);

  const double gridvelderiv_0_0 = gridvel_deriv_(0, 0);
  const double gridvelderiv_0_1 = gridvel_deriv_(0, 1);
  const double gridvelderiv_1_0 = gridvel_deriv_(1, 0);
  const double gridvelderiv_1_1 = gridvel_deriv_(1, 1);

  const double velint_0 = Base::velint_(0);
  const double velint_1 = Base::velint_(1);

  const double gridvelint_0 = gridvel_int_(0);
  const double gridvelint_1 = gridvel_int_(1);

  const double convvelint_0 = Base::convvelint_(0);
  const double convvelint_1 = Base::convvelint_(1);

  if (!porofldpara_->IsStationaryConti())
  {
    // dphi_dp*dp/dt
    for (int vi = 0; vi < nen_; ++vi)
    {
      const double v = Base::fac_ * Base::funct_(vi, 0) * (dphi_dp * press_) +
                       timefacfac * Base::funct_(vi, 0) * refporositydot;
      for (int ui = 0; ui < nen_; ++ui)
      {
        ecoupl_p(vi, ui * 2) += v * Base::derxy_(0, ui);
        ecoupl_p(vi, ui * 2 + 1) += v * Base::derxy_(1, ui);
      }
    }
    // rhs
    for (int vi = 0; vi < nen_; ++vi)
    {
      const double v = -1.0 * timefacfac * Base::funct_(vi, 0) * dphi_dp * Base::rhscon_;
      for (int ui = 0; ui < nen_; ++ui)
      {
        ecoupl_p(vi, ui * 2) += v * Base::derxy_(0, ui);
        ecoupl_p(vi, ui * 2 + 1) += v * Base::derxy_(1, ui);
      }
    }
  }

  const double timefacfac_det = timefacfac / Base::det_;
  if (!porofldpara_->PoroContiPartInt())
  {
    // (porosity)*div u
    for (int vi = 0; vi < nen_; ++vi)
    {
      const double v = timefacfac_det * Base::funct_(vi, 0) * porosity_;
      for (int ui = 0; ui < nen_; ++ui)
      {
        ecoupl_p(vi, ui * 2) +=
            v * (-vderiv_1_0 * Base::deriv_(1, ui) + vderiv_1_1 * Base::deriv_(0, ui));

        ecoupl_p(vi, ui * 2 + 1) +=
            v * (+vderiv_0_0 * Base::deriv_(1, ui) - vderiv_0_1 * Base::deriv_(0, ui));
      }
    }


    if (!porofldpara_->IsStationaryConti())
    {
      // (dphi_dJ*J_)*div vs
      for (int vi = 0; vi < nen_; ++vi)
      {
        const double v = timefacfac_det * Base::funct_(vi, 0) * dphi_dJ * J_;
        for (int ui = 0; ui < nen_; ++ui)
        {
          ecoupl_p(vi, ui * 2) += v * (-gridvelderiv_1_0 * Base::deriv_(1, ui) +
                                          gridvelderiv_1_1 * Base::deriv_(0, ui));

          ecoupl_p(vi, ui * 2 + 1) += v * (+gridvelderiv_0_0 * Base::deriv_(1, ui) -
                                              gridvelderiv_0_1 * Base::deriv_(0, ui));
        }
      }
    }

    //-----------(u-vs)grad(phi)

    for (int ui = 0; ui < nen_; ++ui)
    {
      const double v00 =
          +(velint_1 - gridvelint_1) * (-refgrad_porosity_(0) * Base::deriv_(1, ui) +
                                           refgrad_porosity_(1) * Base::deriv_(0, ui));
      const double v01 =
          +(velint_0 - gridvelint_0) * (+refgrad_porosity_(0) * Base::deriv_(1, ui) -
                                           refgrad_porosity_(1) * Base::deriv_(0, ui));

      for (int vi = 0; vi < nen_; ++vi)
      {
        const double v = timefacfac_det * Base::funct_(vi);

        ecoupl_p(vi, ui * 2) += v * v00;
        ecoupl_p(vi, ui * 2 + 1) += v * v01;
      }
    }
  }
  else
  {
    if (!porofldpara_->IsStationaryConti())
    {
      // (dphi_dJ*J+phi)*div vs
      for (int vi = 0; vi < nen_; ++vi)
      {
        const double v = timefacfac_det * Base::funct_(vi, 0) * (dphi_dJ * J_ + porosity_);
        for (int ui = 0; ui < nen_; ++ui)
        {
          ecoupl_p(vi, ui * 2) += v * (-gridvelderiv_1_0 * Base::deriv_(1, ui) +
                                          gridvelderiv_1_1 * Base::deriv_(0, ui));

          ecoupl_p(vi, ui * 2 + 1) += v * (+gridvelderiv_0_0 * Base::deriv_(1, ui) -
                                              gridvelderiv_0_1 * Base::deriv_(0, ui));
        }
      }
    }

    //----------- phi * (u-vs)grad(vi)

    const double v00 = -1.0 * timefacfac_det * porosity_ * (velint_1 - gridvelint_1);
    const double v01 = -1.0 * timefacfac_det * porosity_ * (velint_0 - gridvelint_0);

    for (int ui = 0; ui < nen_; ++ui)
    {
      const double deriv_0_ui = Base::deriv_(0, ui);
      const double deriv_1_ui = Base::deriv_(1, ui);

      for (int vi = 0; vi < nen_; ++vi)
      {
        ecoupl_p(vi, ui * 2) +=
            v00 * (-Base::deriv_(0, vi) * deriv_1_ui + Base::deriv_(1, vi) * deriv_0_ui);
        ecoupl_p(vi, ui * 2 + 1) +=
            v01 * (+Base::deriv_(0, vi) * deriv_1_ui - Base::deriv_(1, vi) * deriv_0_ui);
      }
    }

  }  // partial integration
  //-------------------
  if (Base::fldpara_->PSPG())
  {
    // shape derivative of gradient of test function
    {
      const double v00 = -1.0 * timefacfac_det * Base::sgvelint_(1);
      const double v01 = -1.0 * timefacfac_det * Base::sgvelint_(0);

      for (int ui = 0; ui < nen_; ++ui)
      {
        const double deriv_0_ui = Base::deriv_(0, ui);
        const double deriv_1_ui = Base::deriv_(1, ui);
        for (int vi = 0; vi < nen_; ++vi)
        {
          ecoupl_p(vi, ui * 2) +=
              v00 * (-Base::deriv_(0, vi) * deriv_1_ui + Base::deriv_(1, vi) * deriv_0_ui);
          ecoupl_p(vi, ui * 2 + 1) +=
              v01 * (+Base::deriv_(0, vi) * deriv_1_ui - Base::deriv_(1, vi) * deriv_0_ui);
        }
      }
    }

    double scal_grad_q = 0.0;

    if (Base::fldpara_->Tds() == INPAR::FLUID::subscales_quasistatic)
    {
      scal_grad_q = Base::tau_(1);
    }
    else
    {
      scal_grad_q = 0.0;  // Base::fldpara_->AlphaF()*fac3;
    }

    // pressure
    {
      // shape derivative of pressure gradient
      const double refgradp_0 = refgrad_press_(0);
      const double refgradp_1 = refgrad_press_(1);
      const double v = timefacfac_det * scal_grad_q;
      for (int ui = 0; ui < nen_; ++ui)
      {
        const double deriv_0_ui = Base::deriv_(0, ui);
        const double deriv_1_ui = Base::deriv_(1, ui);

        for (int vi = 0; vi < nen_; ++vi)
        {
          double v00 = +Base::derxy_(1, vi) * (-refgradp_0 * deriv_1_ui + refgradp_1 * deriv_0_ui);
          double v01 = +Base::derxy_(0, vi) * (refgradp_0 * deriv_1_ui - refgradp_1 * deriv_0_ui);

          ecoupl_p(vi, ui * 2 + 0) += v * v00;
          ecoupl_p(vi, ui * 2 + 1) += v * v01;
        }
      }

      // shape derivative of Jacobian
      static LINALG::Matrix<nen_, 1> temp;
      temp.MultiplyTN(Base::derxy_, Base::gradp_);
      for (int vi = 0; vi < nen_; ++vi)
      {
        const double v1 = -1.0 * timefacfac * scal_grad_q * temp(vi);
        for (int ui = 0; ui < nen_; ++ui)
        {
          ecoupl_p(vi, ui * 2) += v1 * Base::derxy_(0, ui);
          ecoupl_p(vi, ui * 2 + 1) += v1 * Base::derxy_(1, ui);
        }
      }
    }

    // convective term
    {
      const double v = Base::densaf_ * timefacfac_det * scal_grad_q;

      // shape derivative of gradient of velocity
      for (int ui = 0; ui < nen_; ++ui)
      {
        const double deriv_0_ui = Base::deriv_(0, ui);
        const double deriv_1_ui = Base::deriv_(1, ui);

        for (int vi = 0; vi < nen_; ++vi)
        {
          double v00 = +Base::derxy_(1, vi) * convvelint_1 *
                       (-vderiv_0_0 * deriv_1_ui + vderiv_0_1 * deriv_0_ui);
          double v10 = +Base::derxy_(1, vi) * convvelint_1 *
                       (vderiv_1_0 * deriv_1_ui - vderiv_1_1 * deriv_0_ui);
          double v01 = +Base::derxy_(0, vi) * convvelint_0 *
                       (-vderiv_0_0 * deriv_1_ui + vderiv_0_1 * deriv_0_ui);
          double v11 = +Base::derxy_(0, vi) * convvelint_0 *
                       (vderiv_1_0 * deriv_1_ui - vderiv_1_1 * deriv_0_ui);

          ecoupl_p(vi, ui * 2 + 0) += v * (v00 + v10);
          ecoupl_p(vi, ui * 2 + 1) += v * (v01 + v11);
        }
      }

      // shape derivative of Jacobian
      static LINALG::Matrix<nen_, 1> temp;
      temp.MultiplyTN(Base::derxy_, Base::conv_old_);
      for (int vi = 0; vi < nen_; ++vi)
      {
        const double v1 = -1.0 * Base::densaf_ * timefacfac * scal_grad_q * temp(vi);
        for (int ui = 0; ui < nen_; ++ui)
        {
          ecoupl_p(vi, ui * 2) += v1 * Base::derxy_(0, ui);
          ecoupl_p(vi, ui * 2 + 1) += v1 * Base::derxy_(1, ui);
        }
      }
    }
  }

  if (porofldpara_->StabBiot() and (not porofldpara_->IsStationaryConti()) and
      struct_mat_->PoroLawType() != INPAR::MAT::m_poro_law_constant)
  {
    // shape derivative of Jacobian
    {
      const double mixres_0 = mixres_(0);
      const double mixres_1 = mixres_(1);
      const double v = timefacfac * tau_struct_;
      for (int vi = 0; vi < nen_; ++vi)
      {
        double w = v * (N_XYZ_(0, vi) * mixres_0 + N_XYZ_(1, vi) * mixres_1);
        for (int ui = 0; ui < nen_; ++ui)
        {
          ecoupl_p(vi, ui * 2) += w * Base::derxy_(0, ui);
          ecoupl_p(vi, ui * 2 + 1) += w * Base::derxy_(1, ui);
        }
      }
    }

    // shape derivative of gradient of test function
    {
      const double v00 = tau_struct_ * timefacfac_det * Base::gradp_(1);
      const double v01 = tau_struct_ * timefacfac_det * Base::gradp_(0);

      for (int ui = 0; ui < nen_; ++ui)
      {
        const double deriv_0_ui = Base::deriv_(0, ui);
        const double deriv_1_ui = Base::deriv_(1, ui);
        for (int vi = 0; vi < nen_; ++vi)
        {
          ecoupl_p(vi, ui * 2) +=
              v00 * (-Base::deriv_(0, vi) * deriv_1_ui + Base::deriv_(1, vi) * deriv_0_ui);
          ecoupl_p(vi, ui * 2 + 1) +=
              v01 * (+Base::deriv_(0, vi) * deriv_1_ui - Base::deriv_(1, vi) * deriv_0_ui);
        }
      }
    }

    // pressure
    {
      // shape derivative of pressure gradient
      const double refgradp_0 = refgrad_press_(0);
      const double refgradp_1 = refgrad_press_(1);
      const double v = timefacfac_det * tau_struct_;
      for (int ui = 0; ui < nen_; ++ui)
      {
        const double deriv_0_ui = Base::deriv_(0, ui);
        const double deriv_1_ui = Base::deriv_(1, ui);

        for (int vi = 0; vi < nen_; ++vi)
        {
          double v00 = +Base::derxy_(1, vi) * (-refgradp_0 * deriv_1_ui + refgradp_1 * deriv_0_ui);
          double v01 = +Base::derxy_(0, vi) * (refgradp_0 * deriv_1_ui - refgradp_1 * deriv_0_ui);

          ecoupl_p(vi, ui * 2 + 0) += v * v00;
          ecoupl_p(vi, ui * 2 + 1) += v * v01;
        }
      }

      // shape derivative of Jacobian
      static LINALG::Matrix<nen_, 1> temp;
      temp.MultiplyTN(Base::derxy_, Base::gradp_);
      for (int vi = 0; vi < nen_; ++vi)
      {
        const double v1 = -1.0 * timefacfac * tau_struct_ * temp(vi);
        for (int ui = 0; ui < nen_; ++ui)
        {
          ecoupl_p(vi, ui * 2) += v1 * Base::derxy_(0, ui);
          ecoupl_p(vi, ui * 2 + 1) += v1 * Base::derxy_(1, ui);
        }
      }
    }
    {
      // shape derivative of pressure gradient
      const double refgradp_0 = refgrad_press_(0);
      const double refgradp_1 = refgrad_press_(1);
      const double v = timefacfac_det * (-1.0) * porosity_ * tau_struct_;
      for (int ui = 0; ui < nen_; ++ui)
      {
        const double deriv_0_ui = Base::deriv_(0, ui);
        const double deriv_1_ui = Base::deriv_(1, ui);

        for (int vi = 0; vi < nen_; ++vi)
        {
          double v00 = +N_XYZ_(1, vi) * (-refgradp_0 * deriv_1_ui + refgradp_1 * deriv_0_ui);
          double v01 = +N_XYZ_(0, vi) * (refgradp_0 * deriv_1_ui - refgradp_1 * deriv_0_ui);

          ecoupl_p(vi, ui * 2 + 0) += v * v00;
          ecoupl_p(vi, ui * 2 + 1) += v * v01;
        }
      }
    }

    // shape derivative of Jacobian
    {
      const double gradp_0 = Base::gradp_(0);
      const double gradp_1 = Base::gradp_(1);
      const double v = timefacfac * porosity_ * tau_struct_;
      for (int vi = 0; vi < nen_; ++vi)
      {
        double w = v * (N_XYZ_(0, vi) * gradp_0 + N_XYZ_(1, vi) * gradp_1);
        for (int ui = 0; ui < nen_; ++ui)
        {
          ecoupl_p(vi, ui * 2) += w * Base::derxy_(0, ui);
          ecoupl_p(vi, ui * 2 + 1) += w * Base::derxy_(1, ui);
        }
      }
    }
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::PSPG(LINALG::Matrix<nen_, nen_ * nsd_>& estif_q_u,
    LINALG::Matrix<nen_, nen_>& ppmat, LINALG::Matrix<nen_, 1>& preforce,
    const LINALG::Matrix<nsd_ * nsd_, nen_>& lin_resM_Du,
    const LINALG::Matrix<nsd_ * nsd_, nen_>& lin_resMRea_Du,
    const LINALG::Matrix<nsd_, nen_>& lin_resM_Dp, const double& dphi_dp, const double& fac3,
    const double& timefacfac, const double& timefacfacpre, const double& rhsfac)
{
  // conservative, stabilization terms are neglected (Hughes)

  /* pressure stabilisation:                                            */
  /*
              /                 \
             |  ~n+af            |
           - |  u     , nabla q  |
             |                   |
              \                 /
  */

  double scal_grad_q = 0.0;

  if (Base::fldpara_->Tds() == INPAR::FLUID::subscales_quasistatic)
  {
    scal_grad_q = Base::tau_(1);
  }
  else
  {
    scal_grad_q = Base::fldparatimint_->AlphaF() * fac3;
  }

  /* pressure stabilisation: inertia if not stationary*/
  /*
            /                  \
           |                    |
           |  rho*Du , nabla q  |
           |                    |
            \                  /
  */
  /* pressure stabilisation: convection, convective part */
  /*
            /                                   \
           |  /       n+1       \                |
           | |   rho*u   o nabla | Du , nabla q  |
           |  \      (i)        /                |
            \                                   /
  */
  /* pressure stabilisation: convection, reactive part if Newton */
  /*
            /                                   \
           |  /                \   n+1           |
           | |   rho*Du o nabla | u     , grad q |
           |  \                /   (i)           |
            \                                   /
  */
  /* pressure stabilisation: reaction if included */
  /*
            /                     \
           |                      |
           |  sigma*Du , nabla q  |
           |                      |
            \                    /
  */
  /* pressure stabilisation: viscosity (-L_visc_u) */
  /*
            /                              \
           |               /  \             |
       mu  |  nabla o eps | Du | , nabla q  |
           |               \  /             |
            \                              /
  */

  for (int jdim = 0; jdim < nsd_; ++jdim)
  {
    for (int ui = 0; ui < nen_; ++ui)
    {
      const int fui_p_jdim = nsd_ * ui + jdim;

      for (int idim = 0; idim < nsd_; ++idim)
      {
        const int nsd_idim = nsd_ * idim;
        const double lin_resM_Du_idim_jdim_ui = lin_resM_Du(nsd_idim + jdim, ui);

        for (int vi = 0; vi < nen_; ++vi)
        {
          estif_q_u(vi, fui_p_jdim) +=
              lin_resM_Du_idim_jdim_ui * Base::derxy_(idim, vi) * scal_grad_q;
        }
      }
    }
  }


  for (int ui = 0; ui < nen_; ++ui)
  {
    for (int vi = 0; vi < nen_; ++vi)
    {
      /* pressure stabilisation: pressure( L_pres_p) */
      /*
           /                    \
          |                      |
          |  nabla Dp , nabla q  |
          |                      |
           \                    /
      */
      double sum = 0.;
      double sum2 = 0.;
      for (int idim = 0; idim < nsd_; ++idim)
      {
        sum += Base::derxy_(idim, ui) * Base::derxy_(idim, vi);
        sum2 += lin_resM_Dp(idim, ui) * Base::derxy_(idim, vi);
      }

      ppmat(vi, ui) += timefacfacpre * scal_grad_q * sum;
      ppmat(vi, ui) += scal_grad_q * sum2;
    }
  }

  {
    const double v1 = -1.0 * timefacfacpre * dtau_dphi_(1) / scal_grad_q * dphi_dp;
    for (int ui = 0; ui < nen_; ++ui)
    {
      for (int idim = 0; idim < nsd_; ++idim)
      {
        const double v = v1 * Base::sgvelint_(idim) * Base::funct_(ui);

        for (int vi = 0; vi < nen_; ++vi)
        {
          ppmat(vi, ui) += v * Base::derxy_(idim, vi);
        }
      }
    }
  }

  for (int idim = 0; idim < nsd_; ++idim)
  {
    const double temp = rhsfac * Base::sgvelint_(idim);

    for (int vi = 0; vi < nen_; ++vi)
    {
      // pressure stabilisation
      preforce(vi) -= -1.0 * temp * Base::derxy_(idim, vi);
    }
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::StabBiot(
    LINALG::Matrix<nen_, nen_ * nsd_>& estif_q_u, LINALG::Matrix<nen_, nen_>& ppmat,
    LINALG::Matrix<nen_, 1>& preforce, const LINALG::Matrix<nsd_ * nsd_, nen_>& lin_resM_Du,
    const LINALG::Matrix<nsd_ * nsd_, nen_>& lin_resMRea_Du,
    const LINALG::Matrix<nsd_, nen_>& lin_resM_Dp, const double& dphi_dp, const double& fac3,
    const double& timefacfac, const double& timefacfacpre, const double& rhsfac)
{
  for (int idim = 0; idim < nsd_; ++idim)
  {
    const double temp = tau_struct_ * rhsfac / J_ * mixres_(idim);

    for (int vi = 0; vi < nen_; ++vi)
    {
      // pressure stabilisation
      preforce(vi) -= temp * N_XYZ_(idim, vi);
    }
  }

  for (int vi = 0; vi < nen_; ++vi)
  {
    double visc_N_XYZ_ = 0.;
    double gradp_N_XYZ = 0.;

    for (int idim = 0; idim < nsd_; ++idim)
    {
      visc_N_XYZ_ += Base::visc_old_(idim) * N_XYZ_(idim, vi);
      gradp_N_XYZ += Base::gradp_(idim) * N_XYZ_(idim, vi);
    }

    for (int ui = 0; ui < nen_; ++ui)
    {
      double sum = 0.;

      for (int idim = 0; idim < nsd_; ++idim)
      {
        sum += Base::derxy_(idim, ui) * N_XYZ_(idim, vi);
      }

      ppmat(vi, ui) += timefacfacpre * tau_struct_ * (sum - dphi_dp * visc_N_XYZ_) -
                       timefacfacpre * tau_struct_ *
                           (porosity_ * sum + dphi_dp * gradp_N_XYZ * Base::funct_(ui));
    }
  }

  {
    const double val = -1.0 * tau_struct_ * porosity_;
    for (int jdim = 0; jdim < nsd_; ++jdim)
    {
      for (int ui = 0; ui < nen_; ++ui)
      {
        const int fui_p_jdim = nsd_ * ui + jdim;

        for (int idim = 0; idim < nsd_; ++idim)
        {
          const int nsd_idim = nsd_ * idim;
          const double lin_resMRea_Du_idim_jdim_ui = lin_resMRea_Du(nsd_idim + jdim, ui);

          for (int vi = 0; vi < nen_; ++vi)
          {
            estif_q_u(vi, fui_p_jdim) += val * lin_resMRea_Du_idim_jdim_ui * N_XYZ_(idim, vi);
          }
        }
      }
    }
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::ComputeFDerivative(
    const LINALG::Matrix<nsd_, nen_>& edispnp, const LINALG::Matrix<nsd_, nsd_>& defgrd_inv,
    LINALG::Matrix<nsd_ * nsd_, nsd_>& F_x, LINALG::Matrix<nsd_ * nsd_, nsd_>& F_X)
{
  F_X.Clear();

  if (Base::is_higher_order_ele_ and kintype_ != INPAR::STR::kinem_linear)
  {
    for (int i = 0; i < nsd_; i++)
    {
      for (int j = 0; j < nsd_; j++)
      {
        for (int k = 0; k < nsd_; k++)
          for (int n = 0; n < nen_; n++)
            F_X(i * nsd_ + j, k) += N_XYZ2full_(j * nsd_ + k, n) * edispnp(i, n);
      }
    }
  }

  F_x.Multiply(F_X, defgrd_inv);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::ComputeGradients(const double& J,
    const double& dphidp, const double& dphidJ, const LINALG::Matrix<nsd_ * nsd_, 1>& defgrd_IT_vec,
    const LINALG::Matrix<nsd_ * nsd_, nsd_>& F_x, const LINALG::Matrix<nsd_, 1>& gradp,
    const LINALG::Matrix<nen_, 1>* eporositynp, LINALG::Matrix<nsd_, 1>& gradJ,
    LINALG::Matrix<nsd_, 1>& grad_porosity, LINALG::Matrix<nsd_, 1>& refgrad_porosity)
{
  //---------------------------  dJ/dx = dJ/dF : dF/dx = JF^-T : dF/dx
  gradJ.MultiplyTN(J, F_x, defgrd_IT_vec);

  // if(grad_porosity or refgrad_porosity)
  ComputePorosityGradient(
      dphidp, dphidJ, gradJ, gradp, eporositynp, grad_porosity, refgrad_porosity);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::ComputePorosityGradient(const double& dphidp,
    const double& dphidJ, const LINALG::Matrix<nsd_, 1>& gradJ,
    const LINALG::Matrix<nsd_, 1>& gradp, const LINALG::Matrix<nen_, 1>* eporositynp,
    LINALG::Matrix<nsd_, 1>& grad_porosity, LINALG::Matrix<nsd_, 1>& refgrad_porosity)
{
  // if( (Base::fldpara_->PoroContiPartInt() == false) or Base::visceff_)
  {
    //--------------------- current porosity gradient
    for (int idim = 0; idim < nsd_; ++idim)
      (grad_porosity)(idim) = dphidp * gradp(idim) + dphidJ * gradJ(idim);

    refgrad_porosity.Multiply(Base::xjm_, grad_porosity);
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::ComputeLinearization(const double& dphi_dp,
    const double& dphi_dpp, const double& dphi_dJdp, const LINALG::Matrix<nsd_, 1>& gradJ,
    LINALG::Matrix<nsd_, nen_>& dgradphi_dp)
{
  if (!porofldpara_->PoroContiPartInt() or Base::visceff_)
  {
    //--linearization of porosity gradient w.r.t. pressure at gausspoint
    // d(grad(phi))/dp = dphi/(dJdp)* dJ/dx + d^2phi/(dp)^2 * dp/dx + dphi/dp* N,x
    dgradphi_dp.MultiplyNT(dphi_dJdp, gradJ, Base::funct_);
    dgradphi_dp.MultiplyNT(dphi_dpp, Base::gradp_, Base::funct_, 1.0);
    dgradphi_dp.Update(dphi_dp, Base::derxy_, 1.0);
  }
  else
    dgradphi_dp.Clear();
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::ComputeLinearizationOD(const double& dphi_dJ,
    const double& dphi_dJJ, const double& dphi_dJp, const LINALG::Matrix<nsd_, nsd_>& defgrd_inv,
    const LINALG::Matrix<nsd_ * nsd_, 1>& defgrd_IT_vec,
    const LINALG::Matrix<nsd_ * nsd_, nsd_>& F_x, const LINALG::Matrix<nsd_ * nsd_, nsd_>& F_X,
    const LINALG::Matrix<nsd_, 1>& gradJ, LINALG::Matrix<1, nsd_ * nen_>& dJ_dus,
    LINALG::Matrix<1, nsd_ * nen_>& dphi_dus, LINALG::Matrix<nsd_, nen_ * nsd_>& dgradphi_dus)
{
  //------------------------------------------------dJ/dus = dJ/dF : dF/dus = J * F^-T . N_X = J *
  // N_x
  for (int i = 0; i < nen_; i++)
    for (int j = 0; j < nsd_; j++) dJ_dus(j + i * nsd_) = J_ * Base::derxy_(j, i);

  //--------------------- linearization of porosity w.r.t. structure displacements
  dphi_dus.Update(dphi_dJ, dJ_dus);

  if (!porofldpara_->PoroContiPartInt() or Base::visceff_)
  {
    //---------------------d(gradJ)/dus =  dJ/dus * F^-T . : dF/dx + J * dF^-T/dus : dF/dx + J *
    // F^-T : N_X_x

    // dF^-T/dus : dF/dx = - (F^-1. dN/dx . u_s)^T  : dF/dx
    static LINALG::Matrix<nsd_, nsd_ * nen_> dFinvdus_dFdx(false);
    dFinvdus_dFdx.Clear();
    for (int i = 0; i < nsd_; i++)
    {
      for (int n = 0; n < nen_; n++)
      {
        for (int j = 0; j < nsd_; j++)
        {
          const int gid = nsd_ * n + j;
          const double defgrd_inv_ij = defgrd_inv(i, j);
          for (int k = 0; k < nsd_; k++)
          {
            const double derxy_kn = Base::derxy_(k, n);
            for (int p = 0; p < nsd_; p++)
              dFinvdus_dFdx(p, gid) += -defgrd_inv_ij * derxy_kn * F_x(k * nsd_ + i, p);
          }
        }
      }
    }

    // F^-T : d(dF/dx)/dus =  F^-T : (N,XX * F^ -1 + dF/dX * F^-1 * N,x)
    static LINALG::Matrix<nsd_, nsd_ * nen_> FinvT_dFx_dus(false);
    FinvT_dFx_dus.Clear();

    if (Base::is_higher_order_ele_)
    {
      for (int n = 0; n < nen_; n++)
      {
        for (int j = 0; j < nsd_; j++)
        {
          const int gid = nsd_ * n + j;
          for (int p = 0; p < nsd_; p++)
          {
            double val = 0.0;
            const double derxy_p_n = Base::derxy_(p, n);
            for (int k = 0; k < nsd_; k++)
            {
              const double defgrd_inv_kj = defgrd_inv(k, j);
              const double defgrd_inv_kp = defgrd_inv(k, p);
              for (int i = 0; i < nsd_; i++)
              {
                val += defgrd_inv(i, j) * N_XYZ2full_(i * nsd_ + k, n) * defgrd_inv_kp;
                for (int l = 0; l < nsd_; l++)
                  val += -defgrd_inv(i, l) * F_X(i * nsd_ + l, k) * defgrd_inv_kj * derxy_p_n;
              }
            }
            FinvT_dFx_dus(p, gid) += val;
          }
        }
      }
    }

    //----d(gradJ)/dus =  dJ/dus * F^-T . : dF/dx + J * dF^-T/dus : dF/dx + J * F^-T : N_X_x
    static LINALG::Matrix<1, nsd_> temp;
    temp.MultiplyTN(defgrd_IT_vec, F_x);

    //----d(gradJ)/dus =  dJ/dus * F^-T . : dF/dx + J * dF^-T/dus : dF/dx + J * F^-T : N_X_x
    static LINALG::Matrix<nsd_, nen_ * nsd_> dgradJ_dus;

    dgradJ_dus.MultiplyTN(temp, dJ_dus);

    dgradJ_dus.Update(J_, dFinvdus_dFdx, 1.0);

    dgradJ_dus.Update(J_, FinvT_dFx_dus, 1.0);

    //------------------ d( grad(\phi) ) / du_s = d\phi/(dJ du_s) * dJ/dx+ d\phi/dJ * dJ/(dx*du_s) +
    // d\phi/(dp*du_s) * dp/dx
    dgradphi_dus.Multiply(dphi_dJJ, gradJ, dJ_dus);
    dgradphi_dus.Update(dphi_dJ, dgradJ_dus, 1.0);
    dgradphi_dus.Multiply(dphi_dJp, Base::gradp_, dJ_dus, 1.0);
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::ComputePorosity(Teuchos::ParameterList& params,
    const double& press, const double& J, const int& gp, const LINALG::Matrix<nen_, 1>& shapfct,
    const LINALG::Matrix<nen_, 1>* myporosity, double& porosity, double* dphi_dp, double* dphi_dJ,
    double* dphi_dJdp, double* dphi_dJJ, double* dphi_dpp, bool save)
{
  struct_mat_->ComputePorosity(
      params, press, J, gp, porosity, dphi_dp, dphi_dJ, dphi_dJdp, dphi_dJJ, dphi_dpp, save);
}

template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::FluidEleCalcPoro<distype>::SetupMaterialDerivatives()
{
  //------------------------get determinant of Jacobian dX / ds
  // transposed jacobian "dX/ds"
  LINALG::Matrix<nsd_, nsd_> xjm0;
  xjm0.MultiplyNT(Base::deriv_, xyze0_);

  // inverse of transposed jacobian "ds/dX"
  LINALG::Matrix<nsd_, nsd_> xji0(true);
  double det0 = xji0.Invert(xjm0);

  // ----------------------compute derivatives N_XYZ_ at gp w.r.t. material coordinates
  N_XYZ_.Multiply(xji0, Base::deriv_);

  if (Base::is_higher_order_ele_)
  {
    // get the second derivatives of standard element at current GP w.r.t. XYZ
    DRT::UTILS::gder2<distype, nen_>(xjm0, N_XYZ_, Base::deriv2_, xyze0_, N_XYZ2_);

    if (nsd_ == 3)
    {
      for (int n = 0; n < nen_; n++)
      {
        N_XYZ2full_(0, n) = N_XYZ2_(0, n);
        N_XYZ2full_(1, n) = N_XYZ2_(3, n);
        N_XYZ2full_(2, n) = N_XYZ2_(4, n);

        N_XYZ2full_(3, n) = N_XYZ2_(3, n);
        N_XYZ2full_(4, n) = N_XYZ2_(1, n);
        N_XYZ2full_(5, n) = N_XYZ2_(5, n);

        N_XYZ2full_(6, n) = N_XYZ2_(4, n);
        N_XYZ2full_(7, n) = N_XYZ2_(5, n);
        N_XYZ2full_(8, n) = N_XYZ2_(2, n);
      }
    }
    else
    {
      for (int n = 0; n < nen_; n++)
      {
        N_XYZ2full_(0, n) = N_XYZ2_(0, n);
        N_XYZ2full_(1, n) = N_XYZ2_(2, n);

        N_XYZ2full_(2, n) = N_XYZ2_(2, n);
        N_XYZ2full_(3, n) = N_XYZ2_(1, n);
      }
    }
  }
  else
  {
    N_XYZ2_.Clear();
    N_XYZ2full_.Clear();
  }

  return det0;
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::GetStructMaterial(DRT::ELEMENTS::Fluid* ele)
{
  // get fluid material
  {
    // access second material in structure element
    if (ele->NumMaterial() > 1)
    {
      struct_mat_ = Teuchos::rcp_dynamic_cast<MAT::StructPoro>(ele->Material(1));
      if (struct_mat_->MaterialType() != INPAR::MAT::m_structporo and
          struct_mat_->MaterialType() != INPAR::MAT::m_structpororeaction and
          struct_mat_->MaterialType() != INPAR::MAT::m_structpororeactionECM)
        dserror("invalid structure material for poroelasticity");
    }
    else
      dserror("no second material defined for element %i", Base::eid_);
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::ReacStab(
    LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_u, LINALG::Matrix<nen_ * nsd_, nen_>& estif_p_v,
    LINALG::Matrix<nsd_, nen_>& velforce, LINALG::Matrix<nsd_ * nsd_, nen_>& lin_resM_Du,
    const LINALG::Matrix<nsd_, nen_>& lin_resM_Dp, const double& dphi_dp, const double& timefacfac,
    const double& timefacfacpre, const double& rhsfac, const double& fac3)
{
  //  Base::ReacStab(estif_u,
  //           estif_p_v,
  //           velforce,
  //           lin_resM_Du,
  //           timefacfac,
  //           timefacfacpre,
  //           rhsfac,
  //           fac3);

  // do almost the same as the standard implementation
  double reac_tau;
  if (Base::fldpara_->Tds() == INPAR::FLUID::subscales_quasistatic)
    reac_tau = Base::fldpara_->ViscReaStabFac() * Base::reacoeff_ * Base::tau_(1);
  else
  {
    dserror("Is this factor correct? Check for bugs!");
    reac_tau =
        Base::fldpara_->ViscReaStabFac() * Base::reacoeff_ * Base::fldparatimint_->AlphaF() * fac3;
  }


  /* reactive stabilisation, inertia part if not stationary */
  /*
               /                    \
              |                      |
          -/+ |    rho*Du , sigma*v  |
              |                      |
               \                    /
  */
  /* reactive stabilisation, convective part, convective type */
  /*
             /                                  \
            |  /       n+1       \               |
        -/+ | |   rho*u   o nabla | Du , sigma*v |
            |  \       (i)       /               |
             \                                  /
  */
  /* reactive stabilisation, reactive part of convection */
  /*
             /                                   \
            |  /                \   n+1           |
        -/+ | |   rho*Du o nabla | u    , sigma*v |
            |  \                /   (i)           |
             \                                   /
  */
  /* reactive stabilisation, reaction part if included */
  /*
               /                      \
              |                        |
          -/+ |    sigma*Du , sigma*v  |
              |                        |
               \                      /
  */
  /* reactive stabilisation, viscous part (-L_visc_u) */
  /*
             /                             \
            |               /  \            |
       +/-  |  nabla o eps | Du | , sigma*v |
            |               \  /            |
             \                             /
  */
  for (int vi = 0; vi < nen_; ++vi)
  {
    const double v = reac_tau * Base::funct_(vi);

    for (int idim = 0; idim < nsd_; ++idim)
    {
      const int nsd_idim = nsd_ * idim;

      const int fvi_p_idim = nsd_ * vi + idim;

      for (int jdim = 0; jdim < nsd_; ++jdim)
      {
        const int nsd_idim_p_jdim = nsd_idim + jdim;

        for (int ui = 0; ui < nen_; ++ui)
        {
          const int fui_p_jdim = nsd_ * ui + jdim;

          estif_u(fvi_p_idim, fui_p_jdim) += v * lin_resM_Du(nsd_idim_p_jdim, ui);
        }
      }
    }
  }


  /* reactive stabilisation, pressure part ( L_pres_p) */
  /*
             /                    \
            |                      |
       -/+  |  nabla Dp , sigma*v  |
            |                      |
             \                    /
  */
  const double reac_tau_timefacfacpre = reac_tau * timefacfacpre;
  for (int vi = 0; vi < nen_; ++vi)
  {
    const double v = reac_tau_timefacfacpre * Base::funct_(vi);

    for (int idim = 0; idim < nsd_; ++idim)
    {
      const int fvi = nsd_ * vi + idim;

      for (int ui = 0; ui < nen_; ++ui)
      {
        estif_p_v(fvi, ui) += v * Base::derxy_(idim, ui);
      }
    }
  }

  const double reac_fac = Base::fldpara_->ViscReaStabFac() * rhsfac * Base::reacoeff_;
  for (int idim = 0; idim < nsd_; ++idim)
  {
    const double v = reac_fac * Base::sgvelint_(idim);

    for (int vi = 0; vi < nen_; ++vi)
    {
      velforce(idim, vi) += v * Base::funct_(vi);
    }
  }

  // add poro specific linearizations
  for (int vi = 0; vi < nen_; ++vi)
  {
    const double v = reac_tau * Base::funct_(vi);

    for (int idim = 0; idim < nsd_; ++idim)
    {
      const int fvi = nsd_ * vi + idim;

      for (int ui = 0; ui < nen_; ++ui)
      {
        estif_p_v(fvi, ui) += v * lin_resM_Dp(idim, ui);
      }
    }
  }

  {  // linearization of stabilization parameter w.r.t. fluid pressure
    const double v =
        timefacfac * Base::fldpara_->ViscReaStabFac() * dphi_dp *
        (Base::reacoeff_ * dtau_dphi_(1) / Base::tau_(1) + Base::reacoeff_ / porosity_);
    for (int vi = 0; vi < nen_; ++vi)
    {
      const double w = -1.0 * v * Base::funct_(vi);

      for (int idim = 0; idim < nsd_; ++idim)
      {
        const int fvi = nsd_ * vi + idim;

        for (int ui = 0; ui < nen_; ++ui)
        {
          estif_p_v(fvi, ui) += w * Base::sgvelint_(idim) * Base::funct_(ui);
        }
      }
    }
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::GetMaterialParamters(
    Teuchos::RCP<const MAT::Material> material)
{
  if (Base::fldpara_->MatGp())
  {
    Teuchos::RCP<const MAT::FluidPoro> actmat =
        Teuchos::rcp_static_cast<const MAT::FluidPoro>(material);
    if (actmat->MaterialType() != INPAR::MAT::m_fluidporo)
      dserror("invalid fluid material for poroelasticity");

    // set density at n+alpha_F/n+1 and n+alpha_M/n+1
    Base::densaf_ = actmat->Density();
    Base::densam_ = Base::densaf_;
    Base::densn_ = Base::densaf_;

    // calculate reaction coefficient
    Base::reacoeff_ = actmat->ComputeReactionCoeff() * porosity_;

    Base::visceff_ = actmat->EffectiveViscosity();
  }
  else
    dserror("Fluid material parameters have to be evaluated at gauss point for porous flow!");
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::ComputeSpatialReactionTerms(
    Teuchos::RCP<const MAT::Material> material, const LINALG::Matrix<nsd_, nsd_>& invdefgrd)
{
  Teuchos::RCP<const MAT::FluidPoro> actmat =
      Teuchos::rcp_static_cast<const MAT::FluidPoro>(material);

  // Acquire anisotropic permeability coefficients at the GP in case of a nodal orthotropic material
  const std::vector<double> anisotropic_permeability_coeffs = std::invoke(
      [&]()
      {
        if (actmat->IsNodalOrthotropic())
        {
          return ComputeAnisotropicPermeabilityCoeffsAtGP();
        }
        else
        {
          return std::vector<double>();
        }
      });

  // material reaction tensor = inverse material permeability
  actmat->ComputeReactionTensor(mat_reac_tensor_, J_, porosity_,
      anisotropic_permeability_directions_, anisotropic_permeability_coeffs);

  // spatial reaction tensor = J * F^-T * material reaction tensor * F^-1
  static LINALG::Matrix<nsd_, nsd_> temp(false);
  temp.Multiply(J_ * porosity_, mat_reac_tensor_, invdefgrd);
  reac_tensor_.MultiplyTN(invdefgrd, temp);

  reac_tensor_vel_.Multiply(reac_tensor_, Base::velint_);
  reac_tensor_gridvel_.Multiply(reac_tensor_, gridvel_int_);
  reac_tensor_convvel_.Multiply(reac_tensor_, convvel_);

  // linearisations of material reaction tensor
  actmat->ComputeLinMatReactionTensor(
      mat_reac_tensor_linporosity_, mat_reac_tensor_linJ_, J_, porosity_);

  static LINALG::Matrix<nsd_, nsd_> lin_p_tmp_1(false);
  static LINALG::Matrix<nsd_, nsd_> lin_p_tmp_2(false);

  lin_p_tmp_1.MultiplyTN(J_, invdefgrd, mat_reac_tensor_linporosity_);
  lin_p_tmp_2.Multiply(lin_p_tmp_1, invdefgrd);

  lin_p_vel_.Multiply(lin_p_tmp_2, Base::velint_);
  lin_p_vel_grid_.Multiply(lin_p_tmp_2, gridvel_int_);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::ComputeLinSpatialReactionTerms(
    Teuchos::RCP<const MAT::Material> material, const LINALG::Matrix<nsd_, nsd_>& defgrd_inv,
    const LINALG::Matrix<1, nsd_ * nen_>* dJ_dus, const LINALG::Matrix<1, nsd_ * nen_>* dphi_dus)
{
  Teuchos::RCP<const MAT::FluidPoro> actmat =
      Teuchos::rcp_static_cast<const MAT::FluidPoro>(material);
  if (actmat->VaryingPermeability()) dserror("varying material permeability not yet supported!");

  const double porosity_inv = 1.0 / porosity_;
  const double J_inv = 1.0 / J_;

  reac_tensor_linOD_vel_.Clear();
  reac_tensor_linOD_grid_vel_.Clear();

  // check for constant or not given derivatives
  const bool const_phi = (dphi_dus == nullptr);
  const bool const_J = (dJ_dus == nullptr);


  for (int idim = 0; idim < nsd_; ++idim)
  {
    const double reac_vel_idim = reac_tensor_vel_(idim);
    const double reac_grid_vel_idim = reac_tensor_gridvel_(idim);

    for (int jdim = 0; jdim < nsd_; ++jdim)
    {
      for (int inode = 0; inode < nen_; ++inode)
      {
        double val_reatensorlinODvel = 0.0;
        double val_reatensorlinODgridvel = 0.0;

        const int gid = nsd_ * inode + jdim;

        if (!const_J)
        {
          val_reatensorlinODvel += (*dJ_dus)(gid)*J_inv * reac_vel_idim;
          val_reatensorlinODgridvel += (*dJ_dus)(gid)*J_inv * reac_grid_vel_idim;
        }
        if (!const_phi)
        {
          val_reatensorlinODvel += (*dphi_dus)(gid)*porosity_inv * reac_vel_idim;
          val_reatensorlinODgridvel += (*dphi_dus)(gid)*porosity_inv * reac_grid_vel_idim;
        }

        if (kintype_ != INPAR::STR::kinem_linear)
        {
          const double derxy_idim_inode = Base::derxy_(idim, inode);
          for (int ldim = 0; ldim < nsd_; ++ldim)
          {
            const double defgrd_inv_ldim_jdim = defgrd_inv(ldim, jdim);
            const double defgrd_inv_ldim_idim = defgrd_inv(ldim, idim);
            for (int mdim = 0; mdim < nsd_; ++mdim)
            {
              const double matreatensor_ldim_mdim = mat_reac_tensor_(ldim, mdim);
              const double defgrd_inv_mdim_jdim = defgrd_inv(mdim, jdim);
              for (int kdim = 0; kdim < nsd_; ++kdim)
              {
                val_reatensorlinODvel += J_ * porosity_ * Base::velint_(kdim) *
                                         (-defgrd_inv_ldim_jdim * derxy_idim_inode *
                                                 matreatensor_ldim_mdim * defgrd_inv(mdim, kdim) -
                                             defgrd_inv_ldim_idim * matreatensor_ldim_mdim *
                                                 defgrd_inv_mdim_jdim * Base::derxy_(kdim, inode));
                val_reatensorlinODgridvel +=
                    J_ * porosity_ * gridvel_int_(kdim) *
                    (-defgrd_inv_ldim_jdim * derxy_idim_inode * matreatensor_ldim_mdim *
                            defgrd_inv(mdim, kdim) -
                        defgrd_inv_ldim_idim * matreatensor_ldim_mdim * defgrd_inv_mdim_jdim *
                            Base::derxy_(kdim, inode));
              }
            }
          }
        }
        if (!const_permeability_)  // check if derivatives of reaction tensor are zero -->
                                   // significant speed up
        {
          if (!const_phi)
          {
            const double dphi_dus_gid = (*dphi_dus)(gid);
            for (int j = 0; j < nsd_; ++j)
            {
              const double velint_j = Base::velint_(j);
              const double gridvelint_j = Base::gridvelint_(j);
              for (int k = 0; k < nsd_; ++k)
              {
                const double defgrd_inv_k_idim = defgrd_inv(k, idim);
                for (int l = 0; l < nsd_; ++l)
                {
                  val_reatensorlinODvel += J_ * porosity_ * defgrd_inv_k_idim *
                                           mat_reac_tensor_linporosity_(k, l) * defgrd_inv(l, j) *
                                           velint_j * dphi_dus_gid;
                  val_reatensorlinODgridvel += J_ * porosity_ * defgrd_inv_k_idim *
                                               mat_reac_tensor_linporosity_(k, l) *
                                               defgrd_inv(l, j) * gridvelint_j * dphi_dus_gid;
                }
              }
            }
          }

          if (!const_J)
          {
            const double dJ_dus_gid = (*dJ_dus)(gid);
            for (int j = 0; j < nsd_; ++j)
            {
              const double velint_j = Base::velint_(j);
              const double gridvelint_j = Base::gridvelint_(j);
              for (int k = 0; k < nsd_; ++k)
              {
                const double defgrd_inv_k_idim = defgrd_inv(k, idim);
                for (int l = 0; l < nsd_; ++l)
                {
                  val_reatensorlinODvel += J_ * porosity_ * defgrd_inv_k_idim *
                                           mat_reac_tensor_linJ_(k, l) * defgrd_inv(l, j) *
                                           velint_j * dJ_dus_gid;
                  val_reatensorlinODgridvel += J_ * porosity_ * defgrd_inv_k_idim *
                                               mat_reac_tensor_linJ_(k, l) * defgrd_inv(l, j) *
                                               gridvelint_j * dJ_dus_gid;
                }
              }
            }
          }
        }

        reac_tensor_linOD_vel_(idim, gid) += val_reatensorlinODvel;
        reac_tensor_linOD_grid_vel_(idim, gid) += val_reatensorlinODgridvel;
      }
    }
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::ComputeOldRHSAndSubgridScaleVelocity()
{
  //----------------------------------------------------------------------
  // computation of various residuals and residual-based values such as
  // the subgrid-scale velocity
  //----------------------------------------------------------------------
  if (Base::fldparatimint_->IsGenalpha() and (not porofldpara_->IsStationaryMomentum()))
  {
    // rhs of momentum equation: density*bodyforce at n+alpha_F
    Base::rhsmom_.Update(Base::densaf_, Base::bodyforce_, 0.0);

    // evaluate momentum residual once for all stabilization right hand sides
    for (int rr = 0; rr < nsd_; ++rr)
    {
      Base::momres_old_(rr) = Base::densam_ * Base::accint_(rr) +
                              Base::densaf_ * Base::conv_old_(rr) + Base::gradp_(rr) -
                              2.0 * Base::visceff_ * Base::visc_old_(rr) +
                              reac_tensor_convvel_(rr) - Base::densaf_ * Base::bodyforce_(rr);
    }
  }
  else
  {
    if (not porofldpara_->IsStationaryMomentum())
    {
      // rhs of instationary momentum equation:
      // density*theta*bodyforce at n+1 + density*(histmom/dt)
      //                                      f = rho * g
      // Base::rhsmom_.Update((Base::densn_/Base::fldparatimint_->Dt()),Base::histmom_,Base::densaf_*Base::fldparatimint_->Theta(),Base::bodyforce_);
      Base::rhsmom_.Update(
          (Base::densn_ / Base::fldparatimint_->Dt() / Base::fldparatimint_->Theta()),
          Base::histmom_, Base::densaf_, Base::bodyforce_);

      // compute instationary momentum residual:
      // momres_old = u_(n+1)/dt + theta ( ... ) - Base::histmom_/dt - theta*Base::bodyforce_
      for (int rr = 0; rr < nsd_; ++rr)
      {
        /*Base::momres_old_(rr) = Base::densaf_*Base::velint_(rr)/Base::fldparatimint_->Dt()
                           +Base::fldparatimint_->Theta()*(Base::densaf_*conv_old_(rr)+Base::gradp_(rr)
                           -2*Base::visceff_*visc_old_(rr)+Base::reacoeff_*Base::velint_(rr))-Base::rhsmom_(rr);*/
        Base::momres_old_(rr) =
            ((Base::densaf_ * Base::velint_(rr) / Base::fldparatimint_->Dt() +
                 Base::fldparatimint_->Theta() *
                     (Base::densaf_ * Base::conv_old_(rr) + Base::gradp_(rr) -
                         2.0 * Base::visceff_ * Base::visc_old_(rr) + reac_tensor_convvel_(rr))) /
                Base::fldparatimint_->Theta()) -
            Base::rhsmom_(rr);
      }
    }
    else
    {
      // rhs of stationary momentum equation: density*bodyforce
      //                                       f = rho * g
      Base::rhsmom_.Update(Base::densaf_, Base::bodyforce_, 0.0);

      // compute stationary momentum residual:
      for (int rr = 0; rr < nsd_; ++rr)
      {
        Base::momres_old_(rr) = Base::densaf_ * Base::conv_old_(rr) + Base::gradp_(rr) -
                                2.0 * Base::visceff_ * Base::visc_old_(rr) +
                                reac_tensor_convvel_(rr) - Base::rhsmom_(rr);
      }
    }
  }
  //-------------------------------------------------------
  // compute subgrid-scale velocity
  Base::sgvelint_.Update(-Base::tau_(1), Base::momres_old_, 0.0);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::ComputeStabilizationParameters(const double& vol)
{
  // calculate stabilization parameters at integration point
  if (Base::fldpara_->TauGp())
  {
    // check stabilization parameter definition for porous flow
    if (not(Base::fldpara_->WhichTau() ==
                INPAR::FLUID::tau_franca_madureira_valentin_badia_codina or
            Base::fldpara_->WhichTau() ==
                INPAR::FLUID::tau_franca_madureira_valentin_badia_codina_wo_dt or
            Base::fldpara_->WhichTau() == INPAR::FLUID::tau_not_defined))
      dserror("incorrect definition of stabilization parameter for porous flow");

    if (porosity_ < 1e-15) dserror("zero porosity!");

    /*
    This stabilization parameter is only intended to be used for
    (viscous-)reactive problems such as Darcy(-Stokes/Brinkman) problems.

    literature:
    1) L.P. Franca, A.L. Madureira, F. Valentin, Towards multiscale
       functions: enriching finite element spaces with local but not
       bubble-like functions, Comput. Methods Appl. Mech. Engrg. 194
       (2005) 3006-3021.
    2) S. Badia, R. Codina, Stabilized continuous and discontinuous
       Galerkin techniques for Darcy flow, Comput. Methods Appl.
       Mech. Engrg. 199 (2010) 1654-1667.

    */

    // get element-type constant for tau
    const double mk = DRT::ELEMENTS::MK<distype>();

    // total reaction coefficient sigma_tot: sum of "artificial" reaction
    // due to time factor and reaction coefficient
    double sigma_tot = Base::reacoeff_;

    if (not porofldpara_->IsStationaryMomentum())
    {
      sigma_tot += 1.0 / Base::fldparatimint_->TimeFac();
    }

    // calculate characteristic element length
    double h_u = 0.0;
    double h_p = 0.0;
    Base::CalcCharEleLength(vol, 0.0, h_u, h_p);

    // various parameter computations for case with dt:
    // relating viscous to reactive part
    const double re11 = 2.0 * Base::visceff_ / (mk * Base::densaf_ * sigma_tot * ((h_p) * (h_p)));

    // respective "switching" parameter
    const double xi11 = std::max(re11, 1.0);

    // constants c_u and c_p as suggested in Badia and Codina (2010), method A
    const double c_u = 4.0;
    const double c_p = 4.0;

    // tau_Mu not required for porous flow
    Base::tau_(0) = 0.0;
    Base::tau_(1) = ((h_p) * (h_p)) / (c_u * ((h_p) * (h_p)) * Base::densaf_ * sigma_tot * xi11 +
                                          (2.0 * Base::visceff_ / mk));
    Base::tau_(2) = c_p * ((h_p) * (h_p)) * Base::reacoeff_ / porosity_;

    dtau_dphi_(0) = 0.0;
    if (xi11 == 1.0)
    {
      dtau_dphi_(1) =
          -1.0 * Base::tau_(1) * Base::tau_(1) * c_u * Base::densaf_ * Base::reacoeff_ / porosity_;
    }
    else
      dtau_dphi_(1) = 0.0;
    dtau_dphi_(2) = 0.0;

    if (porofldpara_->StabBiot() and (not porofldpara_->IsStationaryConti()) and
        struct_mat_->PoroLawType() != INPAR::MAT::m_poro_law_constant)
    {
      /*
      Stabilization parameter for Biot problems

      literature:
      1) Badia S., Quaini A. Quateroni A. , Coupling Biot and Navier-Stokes equations for
         modelling fluid-poroelastic media interaction, Comput. Physics. 228
         (2009) 7686-8014.
      2) Wan J., Stabilized finite element methods for coupled geomechanics and multiphase flow
         , PHD Thesis, Stanford University, 2002

      */

      const double scaling = porofldpara_->StabBiotScaling();
      const double effective_stiffness = ComputeEffectiveStiffness();
      // note: I do not know if the stabilization parameter should be devided by TimeFac(). It seems
      // to work, though...
      tau_struct_ =
          scaling * ((h_p) * (h_p)) / 12.0 / effective_stiffness / Base::fldparatimint_->TimeFac();
    }
  }
  else
    dserror("Fluid stabilization parameters have to be evaluated at gauss point for porous flow!");
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::ComputeOldRHSConti(double dphi_dp)
{
  double vel_grad_porosity = 0.0;
  for (int idim = 0; idim < nsd_; ++idim)
    vel_grad_porosity += grad_porosity_(idim) * Base::velint_(idim);

  double grad_porosity_gridvelint = 0.0;
  for (int j = 0; j < nsd_; j++) grad_porosity_gridvelint += grad_porosity_(j) * gridvel_int_(j);

  if (not porofldpara_->IsStationaryConti())
  {
    if (Base::fldparatimint_->IsGenalpha())
    {
      // The continuity equation is formulated, as such that the time derivative
      // of the porosity is replaced by the time derivative of the pressure and the Jacobian
      // before discretizing, i.e.
      // $\frac{d\phi}{dt}=\frac{d\phi}{d p}\frac{d p}{d t}+\frac{d\phi}{dJ}\frac{d J}{d t}$

      Base::rhscon_ = 0.0;
    }
    else if (Base::fldparatimint_->IsOneStepTheta())
    {
      // In this case the continuity equation is formulated, as such that the time derivative
      // of the porosity is replaced by the time derivative of the pressure and the Jacobian
      // before discretizing, i.e.
      // $\frac{d\phi}{dt}=\frac{d\phi}{d p}\frac{d p}{d t}+\frac{d\phi}{dJ}\frac{d J}{d t}$

      // rhs of continuity equation
      Base::rhscon_ = 1.0 / Base::fldparatimint_->Dt() / Base::fldparatimint_->Theta() * hist_con_;

      // this is only needed for conti_stab (deactivated for now). If used, it needs to be checked
      // again!!!
      Base::conres_old_ =
          Base::fldparatimint_->Theta() *
              (Base::vdiv_ * porosity_ + vel_grad_porosity - grad_porosity_gridvelint) +
          dphi_dp * press_ / Base::fldparatimint_->Dt() / Base::fldparatimint_->Theta() -
          Base::rhscon_;
    }
  }
  else
  {
    // no time derivatives -> no history
    Base::rhscon_ = 0.0;

    Base::conres_old_ = Base::vdiv_ * porosity_ + vel_grad_porosity;
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::ComputeLinResMDu(const double& timefacfac,
    LINALG::Matrix<nsd_ * nsd_, nen_>& lin_resM_Du,
    LINALG::Matrix<nsd_ * nsd_, nen_>& lin_resMRea_Du)
{
  std::array<int, nsd_> idim_nsd_p_idim;
  for (int idim = 0; idim < nsd_; ++idim)
  {
    idim_nsd_p_idim[idim] = idim * nsd_ + idim;
  }

  // mass
  if (!porofldpara_->IsStationaryMomentum())
  {
    const double fac_densam = Base::fac_ * Base::densam_;

    for (int ui = 0; ui < nen_; ++ui)
    {
      const double v = fac_densam * Base::funct_(ui);

      for (int idim = 0; idim < nsd_; ++idim)
      {
        lin_resM_Du(idim_nsd_p_idim[idim], ui) += v;
      }
    }
  }

  // reactive part
  for (int ui = 0; ui < nen_; ++ui)
  {
    const double v = timefacfac * Base::funct_(ui);

    for (int idim = 0; idim < nsd_; ++idim)
    {
      for (int jdim = 0; jdim < nsd_; ++jdim)
      {
        lin_resMRea_Du(idim * nsd_ + jdim, ui) += v * reac_tensor_(idim, jdim);
        lin_resM_Du(idim * nsd_ + jdim, ui) += v * reac_tensor_(idim, jdim);
      }
    }
  }

  // convective ALE-part
  const double timefacfac_densaf = timefacfac * Base::densaf_;

  for (int ui = 0; ui < nen_; ++ui)
  {
    const double v = timefacfac_densaf * Base::conv_c_(ui);

    for (int idim = 0; idim < nsd_; ++idim)
    {
      lin_resM_Du(idim_nsd_p_idim[idim], ui) += v;
    }
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::ComputeLinResMDuStabilization(
    const double& timefacfac, LINALG::Matrix<nsd_ * nsd_, nen_>& lin_resM_Du)
{
  if (Base::is_higher_order_ele_ and Base::visceff_)
  {
    const double v = -2.0 * Base::visceff_ * timefacfac;
    for (int idim = 0; idim < nsd_; ++idim)
    {
      const int nsd_idim = nsd_ * idim;

      for (int jdim = 0; jdim < nsd_; ++jdim)
      {
        const int nsd_idim_p_jdim = nsd_idim + jdim;

        for (int ui = 0; ui < nen_; ++ui)
        {
          lin_resM_Du(nsd_idim_p_jdim, ui) += v * Base::viscs2_(nsd_idim_p_jdim, ui);
        }
      }
    }
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::CalcDivEps(const LINALG::Matrix<nsd_, nen_>& evelaf)
{
  /*--- viscous term: div(epsilon(u)) --------------------------------*/
  /*   /                                                \
       |  2 N_x,xx + N_x,yy + N_y,xy + N_x,zz + N_z,xz  |
     1 |                                                |
     - |  N_y,xx + N_x,yx + 2 N_y,yy + N_z,yz + N_y,zz  |
     2 |                                                |
       |  N_z,xx + N_x,zx + N_y,zy + N_z,yy + 2 N_z,zz  |
       \                                                /

       with N_x .. x-line of N
       N_y .. y-line of N                                             */

  // set visc_old to zero
  Base::visc_old_.Clear();

  double porosity_inv = 1.0 / porosity_;
  if (nsd_ == 3)
  {
    const double grad_porosity_0 = grad_porosity_(0);
    const double grad_porosity_1 = grad_porosity_(1);
    const double grad_porosity_2 = grad_porosity_(2);

    for (int inode = 0; inode < nen_; ++inode)
    {
      const double derxy_0_inode = Base::derxy_(0, inode);
      const double derxy_1_inode = Base::derxy_(1, inode);
      const double derxy_2_inode = Base::derxy_(2, inode);

      const double derxy2_0_inode = Base::derxy2_(0, inode);
      const double derxy2_1_inode = Base::derxy2_(1, inode);
      const double derxy2_2_inode = Base::derxy2_(2, inode);
      const double derxy2_3_inode = Base::derxy2_(3, inode);
      const double derxy2_4_inode = Base::derxy2_(4, inode);
      const double derxy2_5_inode = Base::derxy2_(5, inode);

      const double sum = (derxy2_0_inode + derxy2_1_inode + derxy2_2_inode);
      Base::viscs2_(0, inode) =
          0.5 * (sum + derxy2_0_inode) +
          0.5 * porosity_inv *
              (2 * derxy_0_inode * grad_porosity_0 + derxy_1_inode * grad_porosity_1 +
                  derxy_2_inode * grad_porosity_2);
      Base::viscs2_(1, inode) =
          0.5 * derxy2_3_inode + 0.5 * porosity_inv * derxy_0_inode * grad_porosity_1;
      Base::viscs2_(2, inode) =
          0.5 * derxy2_4_inode + 0.5 * porosity_inv * derxy_0_inode * grad_porosity_2;
      /* ****************************************************************** */
      Base::viscs2_(3, inode) =
          0.5 * derxy2_3_inode + 0.5 * porosity_inv * derxy_1_inode * grad_porosity_0;
      Base::viscs2_(4, inode) =
          0.5 * (sum + derxy2_1_inode) +
          0.5 * porosity_inv *
              (derxy_0_inode * grad_porosity_0 + 2 * derxy_1_inode * grad_porosity_1 +
                  derxy_2_inode * grad_porosity_2);
      Base::viscs2_(5, inode) =
          0.5 * derxy2_5_inode + 0.5 * porosity_inv * derxy_1_inode * grad_porosity_2;
      /* ****************************************************************** */
      Base::viscs2_(6, inode) =
          0.5 * derxy2_4_inode + 0.5 * porosity_inv * derxy_2_inode * grad_porosity_0;
      Base::viscs2_(7, inode) =
          0.5 * derxy2_5_inode + 0.5 * porosity_inv * derxy_2_inode * grad_porosity_1;
      Base::viscs2_(8, inode) =
          0.5 * (sum + derxy2_2_inode) +
          0.5 * porosity_inv *
              (derxy_0_inode * grad_porosity_0 + derxy_1_inode * grad_porosity_1 +
                  2 * derxy_2_inode * grad_porosity_2);

      for (int idim = 0; idim < nsd_; ++idim)
      {
        const int nsd_idim = idim * nsd_;
        for (int jdim = 0; jdim < nsd_; ++jdim)
        {
          Base::visc_old_(idim) += Base::viscs2_(nsd_idim + jdim, inode) * evelaf(jdim, inode);
        }
      }
    }
  }
  else if (nsd_ == 2)
  {
    const double grad_porosity_0 = grad_porosity_(0);
    const double grad_porosity_1 = grad_porosity_(1);

    for (int inode = 0; inode < nen_; ++inode)
    {
      const double derxy_0_inode = Base::derxy_(0, inode);
      const double derxy_1_inode = Base::derxy_(1, inode);

      const double derxy2_0_inode = Base::derxy2_(0, inode);
      const double derxy2_1_inode = Base::derxy2_(1, inode);
      const double derxy2_2_inode = Base::derxy2_(2, inode);

      const double sum = (derxy2_0_inode + derxy2_1_inode);
      Base::viscs2_(0, inode) =
          0.5 * (sum + derxy2_0_inode) +
          0.5 * porosity_inv *
              (2 * derxy_0_inode * grad_porosity_0 + derxy_1_inode * grad_porosity_1);
      Base::viscs2_(1, inode) =
          0.5 * derxy2_2_inode + 0.5 * porosity_inv * derxy_0_inode * grad_porosity_1;
      /* ****************************************************************** */
      Base::viscs2_(2, inode) =
          0.5 * derxy2_2_inode + 0.5 * porosity_inv * derxy_1_inode * grad_porosity_0;
      Base::viscs2_(3, inode) =
          0.5 * (sum + derxy2_1_inode) +
          0.5 * porosity_inv *
              (derxy_0_inode * grad_porosity_0 + 2 * derxy_1_inode * grad_porosity_1);

      for (int idim = 0; idim < nsd_; ++idim)
      {
        const int nsd_idim = idim * nsd_;
        for (int jdim = 0; jdim < nsd_; ++jdim)
        {
          Base::visc_old_(idim) += Base::viscs2_(nsd_idim + jdim, inode) * evelaf(jdim, inode);
        }
      }
    }
  }
  else
    dserror("Epsilon(N) is not implemented for the 1D case");
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::ComputeLinResMDp(
    const double& timefacfacpre, const double& dphi_dp, LINALG::Matrix<nsd_, nen_>& lin_resM_Dp)
{
  /* poroelasticity pressure term */
  /*
       /                           \      /                            \
      |         n+1                 |     |         n+1                 |
      |  sigma*u  * dphi/dp*Dp , v  |  -  |  sigma*vs  * dphi/dp*Dp , v |
      |         (i)                 |     |         (i)                 |
       \                           /       \                           /
  */

  for (int ui = 0; ui < nen_; ++ui)
  {
    // const double w = Base::funct_(ui)*timefacfacpre*Base::reacoeff_/porosity_*dphi_dp;
    const double w = Base::funct_(ui) * timefacfacpre * dphi_dp / porosity_;
    for (int idim = 0; idim < nsd_; ++idim)
    {
      lin_resM_Dp(idim, ui) += w * reac_tensor_vel_(idim);
    }
  }
  if (!const_permeability_)  // check if derivatives of reaction tensor are zero --> significant
                             // speed up
  {
    for (int ui = 0; ui < nen_; ++ui)
    {
      const double w1 = Base::funct_(ui) * timefacfacpre * dphi_dp * porosity_;
      for (int idim = 0; idim < nsd_; ++idim)
      {
        lin_resM_Dp(idim, ui) += w1 * lin_p_vel_(idim);
      }
    }
  }

  if (not Base::fldparatimint_->IsStationary())
  {
    const double factor = timefacfacpre / porosity_ * dphi_dp;
    for (int ui = 0; ui < nen_; ++ui)
    {
      const double w = Base::funct_(ui) * factor;
      for (int idim = 0; idim < nsd_; ++idim)
        lin_resM_Dp(idim, ui) += w * (-reac_tensor_gridvel_(idim));
    }
    if (!const_permeability_)  // check if derivatives of reaction tensor are zero --> significant
                               // speed up
    {
      const double factor2 = timefacfacpre * dphi_dp * porosity_;
      for (int ui = 0; ui < nen_; ++ui)
      {
        const double w1 = Base::funct_(ui) * factor2;
        for (int idim = 0; idim < nsd_; ++idim)
          lin_resM_Dp(idim, ui) += -w1 * lin_p_vel_grid_(idim);
      }
    }
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::EvaluateVariablesAtGaussPoint(
    Teuchos::ParameterList& params, const LINALG::Matrix<nsd_, nen_>& ebofoaf,
    const LINALG::Matrix<nsd_, nen_>& evelaf, const LINALG::Matrix<nsd_, nen_>& evelnp,
    const LINALG::Matrix<nsd_, nen_>& eveln, const LINALG::Matrix<nen_, 1>& epreaf,
    const LINALG::Matrix<nen_, 1>& eprenp, const LINALG::Matrix<nen_, 1>& epren,
    const LINALG::Matrix<nen_, 1>& epressnp_timederiv,
    const LINALG::Matrix<nen_, 1>& epressam_timederiv,
    const LINALG::Matrix<nen_, 1>& epressn_timederiv, const LINALG::Matrix<nsd_, nen_>& eaccam,
    const LINALG::Matrix<nsd_, nen_>& edispnp, const LINALG::Matrix<nsd_, nen_>& egridv,
    const LINALG::Matrix<nsd_, nen_>& egridvn, const LINALG::Matrix<nen_, 1>& escaaf,
    const LINALG::Matrix<nsd_, nen_>& emhist, const LINALG::Matrix<nen_, 1>& echist,
    const LINALG::Matrix<nen_, 1>* eporositynp, const LINALG::Matrix<nen_, 1>* eporositydot,
    const LINALG::Matrix<nen_, 1>* eporositydotn)
{
  //----------------------------------------------------------------------
  //  evaluation of various values at integration point:
  //  1) velocity (including derivatives and grid velocity)
  //  2) pressure (including derivatives)
  //  3) body-force vector
  //  4) "history" vector for momentum equation
  //----------------------------------------------------------------------
  // get velocity at integration point
  // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
  Base::velint_.Multiply(evelaf, Base::funct_);

  // get velocity derivatives at integration point
  // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
  Base::vderxy_.MultiplyNT(evelaf, Base::derxy_);

  // get velocity at integration point
  // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
  gridvel_int_.Multiply(egridv, Base::funct_);
  // get velocity at integration point
  // (values at n)
  gridvel_n_int_.Multiply(egridvn, Base::funct_);

  // get convective velocity at integration point
  // (ALE case handled implicitly here using the (potential
  //  mesh-movement-dependent) convective velocity, avoiding
  //  various ALE terms used to be calculated before)
  // convmy::velint_.Update(Base::velint_);
  // Base::convvelint_.Multiply(-1.0, egridv, Base::funct_, 0.0);
  Base::convvelint_.Update(-1.0, gridvel_int_, 0.0);

  convvel_.Update(-1.0, gridvel_int_, 1.0, Base::velint_);

  // get pressure at integration point
  // (value at n+alpha_F for generalized-alpha scheme,
  //  value at n+alpha_F for generalized-alpha-NP schemen, n+1 otherwise)
  if (Base::fldparatimint_->IsGenalphaNP())
    press_ = Base::funct_.Dot(eprenp);
  else
    press_ = Base::funct_.Dot(epreaf);

  // get pressure time derivative at integration point
  // (value at n+alpha_M for generalized-alpha scheme, n+1 otherwise)
  press_dot_ = Base::funct_.Dot(epressam_timederiv);

  //  // get pressure time derivative at integration point
  //  // (value at n )
  //  pressdotn_ = Base::funct_.Dot(epressn_timederiv);
  //
  //  // get pressure time derivative at integration point
  //  // (value at n )
  //  pressdotnp_ = Base::funct_.Dot(epressnp_timederiv);

  // get pressure gradient at integration point
  // (value at n+alpha_F for generalized-alpha scheme,
  //  value at n+alpha_F for generalized-alpha-NP schemen, n+1 otherwise)
  if (Base::fldparatimint_->IsGenalphaNP())
    Base::gradp_.Multiply(Base::derxy_, eprenp);
  else
    Base::gradp_.Multiply(Base::derxy_, epreaf);

  // fluid pressure at gradient w.r.t to reference coordinates at gauss point
  refgrad_press_.Multiply(Base::deriv_, epreaf);

  // get bodyforce at integration point
  // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
  Base::bodyforce_.Multiply(ebofoaf, Base::funct_);

  // get momentum history data at integration point
  // (only required for one-step-theta and BDF2 time-integration schemes)
  Base::histmom_.Multiply(emhist, Base::funct_);

  // "history" of continuity equation, i.e. p^n + \Delta t * (1-theta) * \dot{p}^n
  hist_con_ = Base::funct_.Dot(echist);

  // get acceleration at time n+alpha_M at integration point
  Base::accint_.Multiply(eaccam, Base::funct_);

  // get structure velocity derivatives at integration point
  // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
  LINALG::Matrix<nsd_, nsd_> gridvelderxy;
  gridvelderxy.MultiplyNT(egridv, Base::derxy_);

  // structure velocity derivatives w.r.t. reference coordinates at integration point
  gridvel_deriv_.MultiplyNT(egridv, Base::deriv_);

  //----------------------------------------------------------------------
  //  evaluation of various partial operators at integration point
  //  1) convective term from previous iteration (mandatorily set to zero)
  //  2) viscous term from previous iteration and viscous operator
  //  3) divergence of velocity from previous iteration
  //----------------------------------------------------------------------
  // set convective term from previous iteration to zero (required for
  // using routine for evaluation of momentum rhs/residual as given)
  // conv_old_.Clear();

  // set old convective term to ALE-Term only
  Base::conv_old_.Multiply(Base::vderxy_, Base::convvelint_);
  Base::conv_c_.MultiplyTN(Base::derxy_, Base::convvelint_);

  // compute divergence of velocity from previous iteration
  Base::vdiv_ = 0.0;

  gridvel_div_ = 0.0;

  if (not Base::fldparatimint_->IsGenalphaNP())
  {
    for (int idim = 0; idim < nsd_; ++idim)
    {
      Base::vdiv_ += Base::vderxy_(idim, idim);
      gridvel_div_ += gridvelderxy(idim, idim);
    }
  }
  else
  {
    for (int idim = 0; idim < nsd_; ++idim)
    {
      // get vdiv at time n+1 for np_genalpha,
      LINALG::Matrix<nsd_, nsd_> vderxy;
      vderxy.MultiplyNT(evelnp, Base::derxy_);
      Base::vdiv_ += vderxy(idim, idim);

      gridvel_div_ += gridvelderxy(idim, idim);
    }
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::EvaluateVariablesAtGaussPointOD(
    Teuchos::ParameterList& params, const LINALG::Matrix<nsd_, nen_>& ebofoaf,
    const LINALG::Matrix<nsd_, nen_>& evelaf, const LINALG::Matrix<nsd_, nen_>& evelnp,
    const LINALG::Matrix<nsd_, nen_>& eveln, const LINALG::Matrix<nen_, 1>& epreaf,
    const LINALG::Matrix<nen_, 1>& eprenp, const LINALG::Matrix<nen_, 1>& epren,
    const LINALG::Matrix<nen_, 1>& epressnp_timederiv,
    const LINALG::Matrix<nen_, 1>& epressam_timederiv,
    const LINALG::Matrix<nen_, 1>& epressn_timederiv, const LINALG::Matrix<nsd_, nen_>& eaccam,
    const LINALG::Matrix<nsd_, nen_>& edispnp, const LINALG::Matrix<nsd_, nen_>& edispn,
    const LINALG::Matrix<nsd_, nen_>& egridv, const LINALG::Matrix<nsd_, nen_>& egridvn,
    const LINALG::Matrix<nen_, 1>& escaaf, const LINALG::Matrix<nsd_, nen_>& emhist,
    const LINALG::Matrix<nen_, 1>& echist, const LINALG::Matrix<nen_, 1>* eporositynp)
{
  //----------------------------------------------------------------------
  //  evaluation of various values at integration point:
  //  1) velocity (including Base::derivatives and grid velocity)
  //  2) pressure (including Base::derivatives)
  //  3) body-force vector
  //  4) "history" vector for momentum equation
  //  5) and more
  //----------------------------------------------------------------------
  // get velocity at integration point
  // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
  Base::velint_.Multiply(evelaf, Base::funct_);

  // get velocity Base::derivatives at integration point
  // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
  Base::vderxy_.MultiplyNT(evelaf, Base::derxy_);

  Base::vderiv_.MultiplyNT(evelaf, Base::deriv_);

  // get velocity at integration point
  // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
  gridvel_int_.Multiply(egridv, Base::funct_);

  convvel_.Update(-1.0, gridvel_int_, 1.0, Base::velint_);

  // get convective velocity at integration point
  // (ALE case handled implicitly here using the (potential
  //  mesh-movement-dependent) convective velocity, avoiding
  //  various ALE terms used to be calculated before)
  // convvelint_.Update(Base::velint_);
  Base::convvelint_.Update(-1.0, gridvel_int_, 0.0);

  // get pressure at integration point
  // (value at n+alpha_F for generalized-alpha scheme,
  //  value at n+alpha_F for generalized-alpha-NP schemen, n+1 otherwise)
  // double press(true);
  if (Base::fldparatimint_->IsGenalphaNP())
    press_ = Base::funct_.Dot(eprenp);
  else
    press_ = Base::funct_.Dot(epreaf);

  // fluid pressure at gradient w.r.t to paramater space coordinates at gauss point
  refgrad_press_.Multiply(Base::deriv_, epreaf);

  // get pressure time derivative at integration point
  // (value at n+alpha_M for generalized-alpha scheme, n+1 otherwise)
  press_dot_ = Base::funct_.Dot(epressam_timederiv);

  // get pressure time derivative at integration point
  // (value at n )
  // pressdotn_ = Base::funct_.Dot(epressn_timederiv);

  // get pressure time derivative at integration point
  // (value at n )
  // pressdotnp_ = Base::funct_.Dot(epressnp_timederiv);

  // get pressure gradient at integration point
  // (value at n+alpha_F for generalized-alpha scheme,
  //  value at n+alpha_F for generalized-alpha-NP schemen, n+1 otherwise)
  if (Base::fldparatimint_->IsGenalphaNP())
    Base::gradp_.Multiply(Base::derxy_, eprenp);
  else
    Base::gradp_.Multiply(Base::derxy_, epreaf);

  // get displacement Base::derivatives at integration point
  // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
  LINALG::Matrix<nsd_, nsd_> gridvelderxy;
  gridvelderxy.MultiplyNT(egridv, Base::derxy_);

  gridvel_deriv_.MultiplyNT(egridv, Base::deriv_);

  // get bodyforce at integration point
  // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
  Base::bodyforce_.Multiply(ebofoaf, Base::funct_);

  // get momentum history data at integration point
  // (only required for one-step-theta and BDF2 time-integration schemes)
  Base::histmom_.Multiply(emhist, Base::funct_);

  // "history" of continuity equation, i.e. p^n + \Delta t * (1-theta) * \dot{p}^n
  hist_con_ = Base::funct_.Dot(echist);

  // get acceleration at time n+alpha_M at integration point
  Base::accint_.Multiply(eaccam, Base::funct_);

  //----------------------------------------------------------------------
  //  evaluation of various partial operators at integration point
  //  1) convective term from previous iteration (mandatorily set to zero)
  //  2) viscous term from previous iteration and viscous operator
  //  3) divergence of velocity from previous iteration
  //----------------------------------------------------------------------
  // set convective term from previous iteration to zero (required for
  // using routine for evaluation of momentum rhs/residual as given)
  //  conv_old_.Clear();

  // set old convective term to ALE-Term only
  Base::conv_old_.Multiply(Base::vderxy_, Base::convvelint_);
  Base::conv_c_.MultiplyTN(Base::derxy_, Base::convvelint_);

  // compute divergence of velocity from previous iteration
  Base::vdiv_ = 0.0;

  gridvel_div_ = 0.0;
  if (not Base::fldparatimint_->IsGenalphaNP())
  {
    for (int idim = 0; idim < nsd_; ++idim)
    {
      Base::vdiv_ += Base::vderxy_(idim, idim);

      gridvel_div_ += gridvelderxy(idim, idim);
    }
  }
  else
  {
    for (int idim = 0; idim < nsd_; ++idim)
    {
      // get vdiv at time n+1 for np_genalpha,
      LINALG::Matrix<nsd_, nsd_> vderxy;
      vderxy.MultiplyNT(evelnp, Base::derxy_);
      Base::vdiv_ += vderxy(idim, idim);

      gridvel_div_ += gridvelderxy(idim, idim);
    }
  }
}

template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcPoro<distype>::ComputeVolume(Teuchos::ParameterList& params,
    DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1)
{
  //----------------------------------------------------------------
  // Now do the nurbs specific stuff (for isogeometric elements)
  //----------------------------------------------------------------
  if (Base::isNurbs_)
  {
    // access knots and weights for this element
    bool zero_size =
        DRT::NURBS::GetMyNurbsKnotsAndWeights(discretization, ele, Base::myknots_, Base::weights_);

    // if we have a zero sized element due to a interpolated point -> exit here
    if (zero_size) return (0);
  }  // Nurbs specific stuff

  // get node coordinates
  GEO::fillInitialPositionArray<distype, nsd_, LINALG::Matrix<nsd_, nen_>>(ele, Base::xyze_);
  // set element id
  Base::eid_ = ele->Id();

  LINALG::Matrix<nsd_, nen_> edispnp(true);
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, &edispnp, nullptr, "dispnp");

  LINALG::Matrix<nsd_, nen_> evelnp(true);
  LINALG::Matrix<nen_, 1> epressnp(true);
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, &evelnp, &epressnp, "velnp");

  xyze0_ = Base::xyze_;
  // get new node positions of ALE mesh
  Base::xyze_ += edispnp;

  // integration loop
  for (DRT::UTILS::GaussIntegration::iterator iquad = Base::intpoints_.begin();
       iquad != Base::intpoints_.end(); ++iquad)
  {
    // evaluate shape functions and derivatives at integration point
    Base::EvalShapeFuncAndDerivsAtIntPoint(iquad.Point(), iquad.Weight());

    //------------------------get determinant of Jacobian dX / ds
    // transposed jacobian "dX/ds"
    LINALG::Matrix<nsd_, nsd_> xjm0;
    xjm0.MultiplyNT(Base::deriv_, xyze0_);

    // inverse of transposed jacobian "ds/dX"
    const double det0 = xjm0.Determinant();

    // determinant of deformationgradient det F = det ( d x / d X ) = det (dx/ds) * ( det(dX/ds)
    // )^-1
    J_ = Base::det_ / det0;

    // pressure at integration point
    press_ = Base::funct_.Dot(epressnp);

    //-----------------------------------computing the porosity
    porosity_ = 0.0;

    // compute scalar at n+alpha_F or n+1
    // const double scalaraf = Base::funct_.Dot(escaaf);
    // params.set<double>("scalar",scalaraf);

    ComputePorosity(params, press_, J_, *(iquad), Base::funct_, nullptr, porosity_, nullptr,
        nullptr, nullptr, nullptr, nullptr, false);

    elevec1(0) += porosity_ * Base::fac_;
  }

  return 0;
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::ComputeDefGradient(
    LINALG::Matrix<nsd_, nsd_>& defgrd, const LINALG::Matrix<nsd_, nen_>& N_XYZ,
    const LINALG::Matrix<nsd_, nen_>& xcurr)
{
  if (kintype_ == INPAR::STR::kinem_nonlinearTotLag)  // total lagrange (nonlinear)
  {
    // (material) deformation gradient F = d xcurr / d xrefe = xcurr * N_XYZ^T
    defgrd.MultiplyNT(xcurr, N_XYZ);  //  (6.17)
  }
  else if (kintype_ == INPAR::STR::kinem_linear)  // linear kinematics
  {
    defgrd.Clear();
    for (int i = 0; i < nsd_; i++) defgrd(i, i) = 1.0;
  }
  else
    dserror("invalid kinematic type! %d", kintype_);
}

template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcPoro<distype>::ComputeError(DRT::ELEMENTS::Fluid* ele,
    Teuchos::ParameterList& params, Teuchos::RCP<MAT::Material>& mat,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseVector& elevec1)
{
  // integrations points and weights
  // more GP than usual due to (possible) cos/exp fcts in analytical solutions
  // degree 5
  const DRT::UTILS::GaussIntegration intpoints(distype, 5);
  return ComputeError(ele, params, mat, discretization, lm, elevec1, intpoints);
}

template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcPoro<distype>::ComputeError(DRT::ELEMENTS::Fluid* ele,
    Teuchos::ParameterList& params, Teuchos::RCP<MAT::Material>& mat,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseVector& elevec1,
    const DRT::UTILS::GaussIntegration& intpoints)
{
  // analytical solution
  LINALG::Matrix<nsd_, 1> u(true);
  double p = 0.0;

  // error
  LINALG::Matrix<nsd_, 1> deltavel(true);
  double deltap = 0.0;

  const int calcerr = DRT::INPUT::get<INPAR::FLUID::CalcError>(params, "calculate error");

  //----------------------------------------------------------------------------
  //   Extract velocity/pressure from global vectors
  //----------------------------------------------------------------------------

  // fill the local element vector/matrix with the global values
  // af_genalpha: velocity/pressure at time n+alpha_F
  // np_genalpha: velocity at time n+alpha_F, pressure at time n+1
  // ost:         velocity/pressure at time n+1
  LINALG::Matrix<nsd_, nen_> evelaf(true);
  LINALG::Matrix<nen_, 1> epreaf(true);
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, &evelaf, &epreaf, "velaf");

  // np_genalpha: additional vector for velocity at time n+1
  LINALG::Matrix<nsd_, nen_> evelnp(true);
  LINALG::Matrix<nen_, 1> eprenp(true);
  if (Base::fldparatimint_->IsGenalphaNP())
    Base::ExtractValuesFromGlobalVector(
        discretization, lm, *Base::rotsymmpbc_, &evelnp, &eprenp, "velnp");

  //----------------------------------------------------------------------------
  //                         ELEMENT GEOMETRY
  //----------------------------------------------------------------------------

  // get node coordinates
  GEO::fillInitialPositionArray<distype, nsd_, LINALG::Matrix<nsd_, nen_>>(ele, Base::xyze_);
  // set element id
  Base::eid_ = ele->Id();

  LINALG::Matrix<nsd_, nen_> edispnp(true);
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, &edispnp, nullptr, "dispnp");

  // get new node positions for isale
  Base::xyze_ += edispnp;

  //------------------------------------------------------------------
  //                       INTEGRATION LOOP
  //------------------------------------------------------------------

  for (DRT::UTILS::GaussIntegration::iterator iquad = intpoints.begin(); iquad != intpoints.end();
       ++iquad)
  {
    // evaluate shape functions and derivatives at integration point
    Base::EvalShapeFuncAndDerivsAtIntPoint(iquad.Point(), iquad.Weight());

    // get velocity at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    Base::velint_.Multiply(evelaf, Base::funct_);

    // get pressure at integration point
    // (value at n+alpha_F for generalized-alpha scheme,
    //  value at n+alpha_F for generalized-alpha-NP schemen, n+1 otherwise)
    double preint(true);
    if (Base::fldparatimint_->IsGenalphaNP())
      preint = Base::funct_.Dot(eprenp);
    else
      preint = Base::funct_.Dot(epreaf);

    /* H1 -error norm
    // compute first derivative of the velocity
    LINALG::Matrix<nsd_,nsd_> dervelint;
    dervelint.MultiplyNT(evelaf,derxy_);
    */

    // get coordinates at integration point
    LINALG::Matrix<nsd_, 1> xyzint(true);
    xyzint.Multiply(Base::xyze_, Base::funct_);

    //  the error is evaluated at the specific time of the used time integration scheme
    //  n+alpha_F for generalized-alpha scheme
    //  value at n+alpha_F for generalized-alpha-NP schemen, n+1 otherwise)
    const double t = Base::fldparatimint_->Time();

    // Compute analytical solution
    switch (calcerr)
    {
      case INPAR::FLUID::byfunct:
      {
        const int func_no = 1;


        // function evaluation requires a 3D position vector!!
        std::array<double, 3> position = {0.0, 0.0, 0.0};

        if (nsd_ == 2)
        {
          position[0] = xyzint(0);
          position[1] = xyzint(1);
          position[2] = 0.0;
        }
        else if (nsd_ == 3)
        {
          position[0] = xyzint(0);
          position[1] = xyzint(1);
          position[2] = xyzint(2);
        }
        else
          dserror("invalid nsd %d", nsd_);

        if (nsd_ == 2)
        {
          const double u_exact_x = DRT::Problem::Instance()
                                       ->FunctionById<DRT::UTILS::FunctionOfSpaceTime>(func_no - 1)
                                       .Evaluate(position.data(), t, 0);
          const double u_exact_y = DRT::Problem::Instance()
                                       ->FunctionById<DRT::UTILS::FunctionOfSpaceTime>(func_no - 1)
                                       .Evaluate(position.data(), t, 1);
          const double p_exact = DRT::Problem::Instance()
                                     ->FunctionById<DRT::UTILS::FunctionOfSpaceTime>(func_no - 1)
                                     .Evaluate(position.data(), t, 2);

          u(0) = u_exact_x;
          u(1) = u_exact_y;
          p = p_exact;
        }
        else if (nsd_ == 3)
        {
          const double u_exact_x = DRT::Problem::Instance()
                                       ->FunctionById<DRT::UTILS::FunctionOfSpaceTime>(func_no - 1)
                                       .Evaluate(position.data(), t, 0);
          const double u_exact_y = DRT::Problem::Instance()
                                       ->FunctionById<DRT::UTILS::FunctionOfSpaceTime>(func_no - 1)
                                       .Evaluate(position.data(), t, 1);
          const double u_exact_z = DRT::Problem::Instance()
                                       ->FunctionById<DRT::UTILS::FunctionOfSpaceTime>(func_no - 1)
                                       .Evaluate(position.data(), t, 2);
          const double p_exact = DRT::Problem::Instance()
                                     ->FunctionById<DRT::UTILS::FunctionOfSpaceTime>(func_no - 1)
                                     .Evaluate(position.data(), t, 3);

          u(0) = u_exact_x;
          u(1) = u_exact_y;
          u(2) = u_exact_z;
          p = p_exact;
        }
        else
          dserror("invalid dimension");
      }
      break;
      default:
        dserror("analytical solution is not defined");
        break;
    }

    // compute difference between analytical solution and numerical solution
    deltap = preint - p;
    deltavel.Update(1.0, Base::velint_, -1.0, u);

    /* H1 -error norm
    // compute error for first velocity derivative
    for (int i=0;i<nsd_;++i)
      for (int j=0;j<nsd_;++j)
        deltadervel(i,j)= dervelint(i,j) - dervel(i,j);
    */

    // L2 error
    // 0: vel_mag
    // 1: p
    // 2: vel_mag,analytical
    // 3: p_analytic
    // (4: vel_x)
    // (5: vel_y)
    // (6: vel_z)
    for (int isd = 0; isd < nsd_; isd++)
    {
      elevec1[0] += deltavel(isd) * deltavel(isd) * Base::fac_;
      // integrate analytical velocity (computation of relative error)
      elevec1[2] += u(isd) * u(isd) * Base::fac_;
      // velocity components
      // elevec1[isd+4] += deltavel(isd)*deltavel(isd)*fac_;
    }
    elevec1[1] += deltap * deltap * Base::fac_;
    // integrate analytical pressure (computation of relative error)
    elevec1[3] += p * p * Base::fac_;

    /*
    //H1-error norm: first derivative of the velocity
    elevec1[4] += deltadervel.Dot(deltadervel)*fac_;
    */
  }

  return 0;
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::ComputeMixtureStrongResidual(
    Teuchos::ParameterList& params, const LINALG::Matrix<nsd_, nsd_>& defgrd,
    const LINALG::Matrix<nsd_, nen_>& edispnp, const LINALG::Matrix<nsd_, nen_>& edispn,
    const LINALG::Matrix<nsd_ * nsd_, nsd_>& F_X, const int gp, const bool computeLinOD)
{
  const double dens_struct = struct_mat_->Density();
  mixres_.Clear();

  if (not Base::fldparatimint_->IsStationary())
  {
    for (int rr = 0; rr < nsd_; ++rr)
      mixres_(rr) = dens_struct * (gridvel_int_(rr) - gridvel_n_int_(rr)) /
                    Base::fldparatimint_->Theta() / Base::fldparatimint_->Dt();
  }

  for (int rr = 0; rr < nsd_; ++rr)
  {
    mixres_(rr) += -dens_struct * Base::bodyforce_(rr) + J_ * (1.0 - porosity_) * Base::gradp_(rr) -
                   J_ * porosity_ * Base::visc_old_(rr) - J_ * porosity_ * reac_tensor_convvel_(rr);
  }


  if (Base::is_higher_order_ele_)
  {
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    static LINALG::Matrix<6, 1> glstrain(true);
    glstrain.Clear();
    if (kintype_ == INPAR::STR::kinem_nonlinearTotLag)
    {
      // Right Cauchy-Green tensor = F^T * F
      LINALG::Matrix<nsd_, nsd_> cauchygreen;
      cauchygreen.MultiplyTN(defgrd, defgrd);
      // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
      // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
      if (nsd_ == 3)
      {
        glstrain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
        glstrain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
        glstrain(2) = 0.5 * (cauchygreen(2, 2) - 1.0);
        glstrain(3) = cauchygreen(0, 1);
        glstrain(4) = cauchygreen(1, 2);
        glstrain(5) = cauchygreen(2, 0);
      }
      else if (nsd_ == 2)
      {
        glstrain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
        glstrain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
        glstrain(2) = 0.0;
        glstrain(3) = cauchygreen(0, 1);
        glstrain(4) = 0.0;
        glstrain(5) = 0.0;
      }
    }
    else
    {
      LINALG::Matrix<6, nsd_ * nen_> bop(true);
      for (int i = 0; i < nen_; ++i)
      {
        for (int j = 0; j < nsd_; ++j)
          for (int k = 0; k < nsd_; ++k) bop(j, nsd_ * i + k) = defgrd(k, j) * N_XYZ_(j, i);
        /* ~~~ */
        if (nsd_ == 3)
        {
          bop(3, nsd_ * i + 0) = defgrd(0, 0) * N_XYZ_(1, i) + defgrd(0, 1) * N_XYZ_(0, i);
          bop(3, nsd_ * i + 1) = defgrd(1, 0) * N_XYZ_(1, i) + defgrd(1, 1) * N_XYZ_(0, i);
          bop(3, nsd_ * i + 2) = defgrd(2, 0) * N_XYZ_(1, i) + defgrd(2, 1) * N_XYZ_(0, i);
          bop(4, nsd_ * i + 0) = defgrd(0, 1) * N_XYZ_(2, i) + defgrd(0, 2) * N_XYZ_(1, i);
          bop(4, nsd_ * i + 1) = defgrd(1, 1) * N_XYZ_(2, i) + defgrd(1, 2) * N_XYZ_(1, i);
          bop(4, nsd_ * i + 2) = defgrd(2, 1) * N_XYZ_(2, i) + defgrd(2, 2) * N_XYZ_(1, i);
          bop(5, nsd_ * i + 0) = defgrd(0, 2) * N_XYZ_(0, i) + defgrd(0, 0) * N_XYZ_(2, i);
          bop(5, nsd_ * i + 1) = defgrd(1, 2) * N_XYZ_(0, i) + defgrd(1, 0) * N_XYZ_(2, i);
          bop(5, nsd_ * i + 2) = defgrd(2, 2) * N_XYZ_(0, i) + defgrd(2, 0) * N_XYZ_(2, i);
        }
        else if (nsd_ == 2)
        {
          bop(3, nsd_ * i + 0) = 0.5 * (defgrd(0, 0) * N_XYZ_(1, i) + defgrd(0, 1) * N_XYZ_(0, i));
          bop(3, nsd_ * i + 1) = 0.5 * (defgrd(1, 0) * N_XYZ_(1, i) + defgrd(1, 1) * N_XYZ_(0, i));
        }
      }

      // build the linearised strain epsilon = B_L . d
      for (int i = 0; i < nen_; ++i)
      {
        for (int j = 0; j < nsd_; ++j)
        {
          const int dof = i * nsd_ + j;
          glstrain(0) += bop(0, dof) * edispnp(j, i);
          glstrain(1) += bop(1, dof) * edispnp(j, i);
          glstrain(2) += bop(2, dof) * edispnp(j, i);
          glstrain(3) += bop(3, dof) * edispnp(j, i);
          glstrain(4) += bop(4, dof) * edispnp(j, i);
          glstrain(5) += bop(5, dof) * edispnp(j, i);
        }
      }
    }

    static LINALG::Matrix<6, 1> stress_vec(true);
    static LINALG::Matrix<6, 6> cmat(true);
    stress_vec.Clear();
    cmat.Clear();
    struct_mat_->Evaluate(nullptr, &glstrain, params, &stress_vec, &cmat, gp, Base::eid_);

    static LINALG::Matrix<6, nsd_> E_X(true);
    E_X.Clear();

    if (kintype_ == INPAR::STR::kinem_nonlinearTotLag)
    {
      if (nsd_ == 3)
      {
        for (int i = 0; i < nsd_; ++i)
        {
          for (int j = 0; j < nsd_; ++j)
          {
            E_X(0, i) += F_X(j * nsd_ + 0, i) * defgrd(j, 0);
            E_X(1, i) += F_X(j * nsd_ + 1, i) * defgrd(j, 1);
            E_X(2, i) += F_X(j * nsd_ + 2, i) * defgrd(j, 2);

            E_X(3, i) +=
                0.5 * (F_X(j * nsd_ + 0, i) * defgrd(j, 1) + defgrd(j, 0) * F_X(j * nsd_ + 1, i));
            E_X(4, i) +=
                0.5 * (F_X(j * nsd_ + 0, i) * defgrd(j, 2) + defgrd(j, 0) * F_X(j * nsd_ + 2, i));
            E_X(5, i) +=
                0.5 * (F_X(j * nsd_ + 1, i) * defgrd(j, 2) + defgrd(j, 1) * F_X(j * nsd_ + 2, i));
          }
        }
      }
      else if (nsd_ == 2)
      {
        for (int i = 0; i < nsd_; ++i)
        {
          for (int j = 0; j < nsd_; ++j)
          {
            E_X(0, i) += F_X(j * nsd_ + 0, i) * defgrd(j, 0);
            E_X(1, i) += F_X(j * nsd_ + 1, i) * defgrd(j, 1);

            E_X(3, i) +=
                0.5 * (F_X(j * nsd_ + 0, i) * defgrd(j, 1) + defgrd(j, 0) * F_X(j * nsd_ + 1, i));
          }
        }
      }
    }
    else
    {
      if (nsd_ == 3)
      {
        dserror("not implemented");
      }
      else if (nsd_ == 2)
      {
        for (int i = 0; i < nsd_; ++i)
        {
          for (int k = 0; k < nen_; ++k)
          {
            for (int j = 0; j < nsd_; ++j)
            {
              E_X(0, i) += N_XYZ2full_(i * nsd_ + 0, k) * defgrd(j, 0) * edispnp(j, k);
              E_X(1, i) += N_XYZ2full_(i * nsd_ + 1, k) * defgrd(j, 1) * edispnp(j, k);
            }

            E_X(3, i) += 0.5 * (N_XYZ2full_(i * nsd_ + 0, k) * defgrd(1, 1) * edispnp(1, k) +
                                   N_XYZ2full_(i * nsd_ + 1, k) * defgrd(0, 0) * edispnp(0, k));
          }
        }
      }
    }

    static LINALG::Matrix<6, nsd_> cmat_E_X(true);
    cmat_E_X.Multiply(cmat, E_X);
    static LINALG::Matrix<nsd_, 1> cmat_E_X_vec(true);
    if (nsd_ == 3)
    {
      cmat_E_X_vec(0) = cmat_E_X(0, 0) + cmat_E_X(3, 1) + cmat_E_X(4, 2);
      cmat_E_X_vec(1) = cmat_E_X(3, 0) + cmat_E_X(1, 1) + cmat_E_X(5, 2);
      cmat_E_X_vec(2) = cmat_E_X(4, 0) + cmat_E_X(5, 1) + cmat_E_X(2, 2);
    }
    else if (nsd_ == 2)
    {
      cmat_E_X_vec(0) = cmat_E_X(0, 0) + cmat_E_X(3, 1);
      cmat_E_X_vec(1) = cmat_E_X(3, 0) + cmat_E_X(1, 1);
    }

    static LINALG::Matrix<nsd_, nsd_> stress(false);
    if (nsd_ == 3)
    {
      stress(0, 0) = stress_vec(0);
      stress(0, 1) = stress_vec(3);
      stress(0, 2) = stress_vec(4);
      stress(1, 0) = stress_vec(3);
      stress(1, 1) = stress_vec(1);
      stress(1, 2) = stress_vec(5);
      stress(2, 0) = stress_vec(4);
      stress(2, 1) = stress_vec(5);
      stress(2, 2) = stress_vec(2);
    }
    else if (nsd_ == 2)
    {
      stress(0, 0) = stress_vec(0);
      stress(0, 1) = stress_vec(3);
      stress(1, 0) = stress_vec(3);
      stress(1, 1) = stress_vec(1);
    }

    for (int i = 0; i < nsd_; ++i)
      for (int j = 0; j < nsd_; ++j)
        for (int k = 0; k < nsd_; ++k) mixres_(i) -= F_X(i * nsd_ + j, k) * stress(j, k);

    mixres_.MultiplyNN(-1.0, defgrd, cmat_E_X_vec, 1.0);
    if (computeLinOD)
    {
      mixresLinOD_.Clear();
      for (int i = 0; i < nsd_; ++i)
      {
        for (int j = 0; j < nen_; ++j)
          for (int k = 0; k < nsd_; ++k)
            mixresLinOD_(i, j * nsd_ + k) -= N_XYZ_(i, j) * cmat_E_X_vec(k);
      }

      static LINALG::Matrix<6 * nsd_, nsd_ * nen_> E_X_Lin(true);
      E_X_Lin.Clear();
      if (nsd_ == 3)
      {
        for (int i = 0; i < nsd_; ++i)
        {
          for (int j = 0; j < nsd_; ++j)
          {
            for (int k = 0; k < nen_; ++k)
            {
              for (int l = 0; l < nsd_; ++l)
              {
                const int dof = k * nsd_ + l;
                E_X_Lin(i * 6 + 0, dof) += N_XYZ2full_(0 * nsd_ + i, k) * defgrd(l, 0);
                E_X_Lin(i * 6 + 1, dof) += N_XYZ2full_(1 * nsd_ + i, k) * defgrd(l, 1);
                E_X_Lin(i * 6 + 2, dof) += N_XYZ2full_(2 * nsd_ + i, k) * defgrd(l, 2);

                E_X_Lin(i * 6 + 3, dof) += 0.5 * (N_XYZ2full_(0 * nsd_ + i, k) * defgrd(l, 1) +
                                                     defgrd(l, 0) * N_XYZ2full_(1 * nsd_ + i, k));
                E_X_Lin(i * 6 + 4, dof) += 0.5 * (N_XYZ2full_(0 * nsd_ + i, k) * defgrd(l, 2) +
                                                     defgrd(l, 0) * N_XYZ2full_(2 * nsd_ + i, k));
                E_X_Lin(i * 6 + 5, dof) += 0.5 * (N_XYZ2full_(1 * nsd_ + i, k) * defgrd(l, 2) +
                                                     defgrd(l, 1) * N_XYZ2full_(2 * nsd_ + i, k));
              }

              E_X_Lin(i * 6 + 0, k * nsd_ + 0) += F_X(j * nsd_ + 0, i) * N_XYZ_(0, k);
              E_X_Lin(i * 6 + 1, k * nsd_ + 1) += F_X(j * nsd_ + 1, i) * N_XYZ_(1, k);
              E_X_Lin(i * 6 + 2, k * nsd_ + 2) += F_X(j * nsd_ + 2, i) * N_XYZ_(2, k);

              E_X_Lin(i * 6 + 3, k * nsd_ + 0) += 0.5 * (N_XYZ_(0, k) * F_X(j * nsd_ + 1, i));
              E_X_Lin(i * 6 + 3, k * nsd_ + 1) += 0.5 * (F_X(j * nsd_ + 0, i) * N_XYZ_(1, k));
              E_X_Lin(i * 6 + 4, k * nsd_ + 0) += 0.5 * (N_XYZ_(0, k) * F_X(j * nsd_ + 2, i));
              E_X_Lin(i * 6 + 4, k * nsd_ + 2) += 0.5 * (F_X(j * nsd_ + 0, i) * N_XYZ_(2, k));
              E_X_Lin(i * 6 + 5, k * nsd_ + 1) += 0.5 * (N_XYZ_(1, k) * F_X(j * nsd_ + 2, i));
              E_X_Lin(i * 6 + 5, k * nsd_ + 2) += 0.5 * (F_X(j * nsd_ + 1, i) * N_XYZ_(2, k));
            }
          }
        }
      }
      else if (nsd_ == 2)
      {
        for (int i = 0; i < nsd_; ++i)
        {
          for (int j = 0; j < nsd_; ++j)
          {
            for (int k = 0; k < nen_; ++k)
            {
              const int dof = k * nsd_ + j;
              E_X_Lin(i * 6 + 0, dof) += N_XYZ2full_(0 * nsd_ + i, k) * defgrd(j, 0);
              E_X_Lin(i * 6 + 1, dof) += N_XYZ2full_(1 * nsd_ + i, k) * defgrd(j, 1);

              E_X_Lin(i * 6 + 3, dof) += 0.5 * (N_XYZ2full_(0 * nsd_ + i, k) * defgrd(j, 1) +
                                                   defgrd(j, 0) * N_XYZ2full_(1 * nsd_ + i, k));

              E_X_Lin(i * 6 + 0, k * nsd_ + 0) += F_X(j * nsd_ + 0, i) * N_XYZ_(0, k);
              E_X_Lin(i * 6 + 1, k * nsd_ + 1) += F_X(j * nsd_ + 1, i) * N_XYZ_(1, k);

              E_X_Lin(i * 6 + 3, k * nsd_ + 0) += 0.5 * (N_XYZ_(0, k) * F_X(j * nsd_ + 1, i));
              E_X_Lin(i * 6 + 3, k * nsd_ + 1) += 0.5 * (F_X(j * nsd_ + 0, i) * N_XYZ_(1, k));
            }
          }
        }
      }

      static LINALG::Matrix<6 * nsd_, nsd_ * nen_> cmat_E_X_Lin(true);
      cmat_E_X_Lin.Clear();
      for (int i = 0; i < nsd_; ++i)
      {
        for (int j = 0; j < nsd_; ++j)
        {
          for (int k = 0; k < nen_; ++k)
          {
            const int dof = k * nsd_ + j;
            for (int l = 0; l < 6; ++l)
              for (int m = 0; m < 6; ++m)
                cmat_E_X_Lin(i * 6 + m, dof) += cmat(m, l) * E_X_Lin(i * 6 + l, dof);
          }
        }
      }

      static LINALG::Matrix<nsd_, nsd_ * nen_> cmat_E_X_vec_Lin(true);
      if (nsd_ == 3)
      {
        for (int j = 0; j < nsd_; ++j)
        {
          for (int k = 0; k < nen_; ++k)
          {
            const int dof = k * nsd_ + j;
            cmat_E_X_vec_Lin(0, dof) = cmat_E_X_Lin(0 * 6 + 0, dof) + cmat_E_X_Lin(1 * 6 + 3, dof) +
                                       cmat_E_X_Lin(2 * 6 + 4, dof);
            cmat_E_X_vec_Lin(1, dof) = cmat_E_X_Lin(0 * 6 + 3, dof) + cmat_E_X_Lin(1 * 6 + 1, dof) +
                                       cmat_E_X_Lin(2 * 6 + 5, dof);
            cmat_E_X_vec_Lin(2, dof) = cmat_E_X_Lin(0 * 6 + 4, dof) + cmat_E_X_Lin(1 * 6 + 5, dof) +
                                       cmat_E_X_Lin(2 * 6 + 2, dof);
          }
        }
      }
      else if (nsd_ == 2)
      {
        for (int j = 0; j < nsd_; ++j)
        {
          for (int k = 0; k < nen_; ++k)
          {
            const int dof = k * nsd_ + j;
            cmat_E_X_vec_Lin(0, dof) = cmat_E_X_Lin(0 * 6 + 0, dof) + cmat_E_X_Lin(1 * 6 + 3, dof);
            cmat_E_X_vec_Lin(1, dof) = cmat_E_X_Lin(0 * 6 + 3, dof) + cmat_E_X_Lin(1 * 6 + 1, dof);
          }
        }
      }

      mixresLinOD_.MultiplyNN(-1.0, defgrd, cmat_E_X_vec_Lin, 1.0);
    }
  }
}

template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::FluidEleCalcPoro<distype>::ComputeEffectiveStiffness()
{
  Teuchos::RCP<MAT::Material> curmat = struct_mat_->GetMaterial();
  double effective_stiffness = 0.0;

  switch (struct_mat_->GetMaterial()->MaterialType())
  {
    case INPAR::MAT::m_stvenant:
    {
      Teuchos::RCP<MAT::StVenantKirchhoff> stvmat =
          Teuchos::rcp_dynamic_cast<MAT::StVenantKirchhoff>(curmat);
      effective_stiffness = stvmat->ShearMod();
      break;
    }
    case INPAR::MAT::m_elasthyper:
    {
      Teuchos::RCP<MAT::ElastHyper> ehmat = Teuchos::rcp_dynamic_cast<MAT::ElastHyper>(curmat);
      effective_stiffness = ehmat->ShearMod();
      break;
    }
    default:
      dserror(
          "calculation of effective stiffness for biot stabilization not implemented for Material "
          "Type %i",
          struct_mat_->GetMaterial()->MaterialType());
      break;
  }

  if (effective_stiffness < 1e-14)
  {
    dserror(
        "Effective stiffness is very small (%f). ShearMod() method not implemented in material?",
        effective_stiffness);
  }

  return effective_stiffness;
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::ComputeJacobianDeterminantVolumeChange(double& J,
    double& volchange, const LINALG::Matrix<nsd_, nsd_>& defgrd,
    const LINALG::Matrix<nsd_, nen_>& N_XYZ, const LINALG::Matrix<nsd_, nen_>& nodaldisp)
{
  // compute J
  J = defgrd.Determinant();

  if (kintype_ == INPAR::STR::kinem_nonlinearTotLag)  // total lagrange (nonlinear)
  {
    // for nonlinear kinematics the Jacobian of the deformation gradient is the volume change
    volchange = J;
  }
  else if (kintype_ == INPAR::STR::kinem_linear)  // linear kinematics
  {
    // for linear kinematics the volume change is the trace of the linearized strains

    // gradient of displacements
    static LINALG::Matrix<nsd_, nsd_> dispgrad;
    dispgrad.Clear();
    // gradient of displacements
    dispgrad.MultiplyNT(nodaldisp, N_XYZ);

    volchange = 1.0;
    // volchange = 1 + trace of the linearized strains (= trace of displacement gradient)
    for (int i = 0; i < nsd_; ++i) volchange += dispgrad(i, i);
  }
  else
    dserror("invalid kinematic type!");
}

template <DRT::Element::DiscretizationType distype>
std::vector<double>
DRT::ELEMENTS::FluidEleCalcPoro<distype>::ComputeAnisotropicPermeabilityCoeffsAtGP() const
{
  std::vector<double> anisotropic_permeability_coeffs(nsd_, 0.0);

  for (int node = 0; node < nen_; ++node)
  {
    const double shape_val = Base::funct_(node);
    for (int dim = 0; dim < nsd_; ++dim)
    {
      anisotropic_permeability_coeffs[dim] +=
          shape_val * anisotropic_permeability_nodal_coeffs_[dim][node];
    }
  }

  return anisotropic_permeability_coeffs;
}

template class DRT::ELEMENTS::FluidEleCalcPoro<DRT::Element::hex8>;
template class DRT::ELEMENTS::FluidEleCalcPoro<DRT::Element::hex20>;
template class DRT::ELEMENTS::FluidEleCalcPoro<DRT::Element::hex27>;
template class DRT::ELEMENTS::FluidEleCalcPoro<DRT::Element::tet4>;
template class DRT::ELEMENTS::FluidEleCalcPoro<DRT::Element::tet10>;
template class DRT::ELEMENTS::FluidEleCalcPoro<DRT::Element::wedge6>;
template class DRT::ELEMENTS::FluidEleCalcPoro<DRT::Element::wedge15>;
template class DRT::ELEMENTS::FluidEleCalcPoro<DRT::Element::pyramid5>;
template class DRT::ELEMENTS::FluidEleCalcPoro<DRT::Element::quad4>;
template class DRT::ELEMENTS::FluidEleCalcPoro<DRT::Element::quad8>;
template class DRT::ELEMENTS::FluidEleCalcPoro<DRT::Element::quad9>;
template class DRT::ELEMENTS::FluidEleCalcPoro<DRT::Element::tri3>;
template class DRT::ELEMENTS::FluidEleCalcPoro<DRT::Element::tri6>;
template class DRT::ELEMENTS::FluidEleCalcPoro<DRT::Element::nurbs9>;
template class DRT::ELEMENTS::FluidEleCalcPoro<DRT::Element::nurbs27>;