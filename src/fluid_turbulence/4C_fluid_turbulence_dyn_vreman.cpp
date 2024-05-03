/*----------------------------------------------------------------------*/
/*! \file

\brief computes Vreman constant C_v dynamically


\level 2

*/
/*----------------------------------------------------------------------*/

#include "4C_fluid_turbulence_dyn_vreman.hpp"

#include "4C_fluid_ele_action.hpp"
#include "4C_fluid_turbulence_boxfilter.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_newtonianfluid.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_scatra.hpp"
#include "4C_scatra_ele_action.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  Constructor (public)                                     krank 09/13|
 *----------------------------------------------------------------------*/
FLD::Vreman::Vreman(Teuchos::RCP<DRT::Discretization> actdis, Teuchos::ParameterList& params)
    :  // call constructor for "nontrivial" objects
      discret_(actdis),
      params_(params),
      physicaltype_(CORE::UTILS::GetAsEnum<INPAR::FLUID::PhysicalType>(params_, "Physical Type"))
{
  boxf_ = Teuchos::rcp(new FLD::Boxfilter(discret_, params_));
  // Initialize Boxfilter
  boxf_->InitializeVreman();

  return;
}


/*----------------------------------------------------------------------*
 | add some scatra specific parameters                  rasthofer 08/12 |
 * ---------------------------------------------------------------------*/
void FLD::Vreman::AddScatra(Teuchos::RCP<DRT::Discretization> scatradis)
{
  scatradiscret_ = scatradis;

  boxfsc_ = Teuchos::rcp(new FLD::Boxfilter(scatradiscret_, params_));

  // Initialize Boxfilter
  boxfsc_->InitializeVremanScatra(scatradiscret_);

  return;
}

/*----------------------------------------------------------------------*
 | Perform box filter operation, compare filtered quantities            |
 | to solution to get an estimate for Cv, using volume averaging        |
 |                                                             (public) |
 |                                                           krank 09/13|
 *----------------------------------------------------------------------*/
void FLD::Vreman::ApplyFilterForDynamicComputationOfCv(
    const Teuchos::RCP<const Epetra_Vector> velocity,
    const Teuchos::RCP<const Epetra_Vector> scalar, const double thermpress,
    const Teuchos::RCP<const Epetra_Vector> dirichtoggle)
{
  const Epetra_Map* nodecolmap = discret_->NodeColMap();


  col_filtered_strainrate_ = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap, 9, true));
  col_filtered_alphaij_ = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap, 9, true));
  col_filtered_expression_ = Teuchos::rcp(new Epetra_Vector(*nodecolmap, true));
  col_filtered_alpha2_ = Teuchos::rcp(new Epetra_Vector(*nodecolmap, true));


  // perform filtering
  boxf_->ApplyFilter(velocity, scalar, thermpress, dirichtoggle);

  // get fitered fields
  boxf_->GetFilteredVremanStrainrate(col_filtered_strainrate_);
  boxf_->GetFilteredVremanAlphaij(col_filtered_alphaij_);
  boxf_->GetExpression(col_filtered_expression_);
  boxf_->GetAlpha2(col_filtered_alpha2_);

  // compute Cv
  Cv_ = DynVremanComputeCv();

  return;
}

void FLD::Vreman::ApplyFilterForDynamicComputationOfDt(
    const Teuchos::RCP<const Epetra_Vector> scalar, const double thermpress,
    const Teuchos::RCP<const Epetra_Vector> dirichtoggle, Teuchos::ParameterList& extraparams,
    const int ndsvel)
{
  const Epetra_Map* nodecolmap = scatradiscret_->NodeColMap();

  col_filtered_phi_ = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap, 3, true));
  col_filtered_phi2_ = Teuchos::rcp(new Epetra_Vector(*nodecolmap, true));
  col_filtered_phiexpression_ = Teuchos::rcp(new Epetra_Vector(*nodecolmap, true));
  col_filtered_alphaijsc_ = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap, 9, true));


  // perform filtering
  boxfsc_->ApplyFilterScatra(scalar, thermpress, dirichtoggle, ndsvel);
  boxfsc_->GetFilteredPhi(col_filtered_phi_);
  boxfsc_->GetFilteredPhi2(col_filtered_phi2_);
  boxfsc_->GetFilteredPhiexpression(col_filtered_phiexpression_);
  boxfsc_->GetFilteredVremanAlphaijsc(col_filtered_alphaijsc_);
  DynVremanComputeDt(extraparams);
  return;
}



/*----------------------------------------------------------------------*
 | compute Cv from filtered quantities.                       (private) |
 |                                                      krank     09/13 |
 *----------------------------------------------------------------------*/
double FLD::Vreman::DynVremanComputeCv()
{
  double Cv = 0.0;
  double cv_numerator_volumeav = 0.0;
  double cv_denominator_volumeav = 0.0;
  int id = GLOBAL::Problem::Instance()->Materials()->FirstIdByType(CORE::Materials::m_fluid);
  const CORE::MAT::PAR::Parameter* mat =
      GLOBAL::Problem::Instance()->Materials()->ParameterById(id);
  const MAT::PAR::NewtonianFluid* actmat = static_cast<const MAT::PAR::NewtonianFluid*>(mat);
  double dens = actmat->density_;
  double visc = actmat->viscosity_ / dens;
  // action for elements
  // generate a parameterlist for communication and control
  Teuchos::ParameterList calc_vreman_params;
  calc_vreman_params.set<int>("action", FLD::calc_vreman_const);

  // hand filtered global vectors down to the element
  calc_vreman_params.set("col_filtered_strainrate", col_filtered_strainrate_);
  calc_vreman_params.set("col_filtered_alphaij", col_filtered_alphaij_);
  calc_vreman_params.set("col_filtered_alpha2", col_filtered_alpha2_);
  calc_vreman_params.set("col_filtered_expression", col_filtered_expression_);


  // loop all elements on this proc (excluding ghosted ones)
  Teuchos::RCP<CORE::LINALG::SerialDenseVector> Cv_num_denom =
      Teuchos::rcp(new CORE::LINALG::SerialDenseVector(2));


  // call loop over elements (assemble nothing)
  discret_->EvaluateScalars(calc_vreman_params, Cv_num_denom);
  discret_->ClearState();


  cv_numerator_volumeav = (*Cv_num_denom)[0];
  cv_denominator_volumeav = (*Cv_num_denom)[1];
  // multiply with viscosity
  if (sqrt(cv_denominator_volumeav * cv_denominator_volumeav) <
      1.0e-12)  // the denominator might also become negative
    Cv = 0.0;   // constant vreman
  else
    Cv = (-1.0) * visc / 2.0 * cv_numerator_volumeav / cv_denominator_volumeav;
  params_.set<double>("C_vreman_theoretical", Cv);
  if (Cv < 0.0)
  {
    Cv = 0.0;
    if (discret_->Comm().MyPID() == 0)
      std::cout << "!!   Vreman constant negative --> clipping: Cv=0.0   !!" << std::endl;
  }

  // std::cout <<"Vreman constant:   "<< Cv << std::endl;
  params_.set<double>("C_vreman", Cv);



  return Cv;
}  // end FLD::Vreman::DynVremanComputeCv

void FLD::Vreman::DynVremanComputeDt(Teuchos::ParameterList& extraparams)
{
  double Dt = 0.0;
  double dt_numerator_volumeav = 0.0;
  double dt_denominator_volumeav = 0.0;
  int idscatra = GLOBAL::Problem::Instance()->Materials()->FirstIdByType(CORE::Materials::m_scatra);
  const CORE::MAT::PAR::Parameter* matscatra =
      GLOBAL::Problem::Instance()->Materials()->ParameterById(idscatra);
  const MAT::PAR::ScatraMat* actmatscatra = static_cast<const MAT::PAR::ScatraMat*>(matscatra);
  double diffus = MAT::PAR::ScatraMat(*actmatscatra)
                      .GetParameter(actmatscatra->diff, -1);  // actmatscatra->diffusivity_;

  // generate a parameterlist for communication and control
  Teuchos::ParameterList calc_vreman_params_scatra;
  CORE::UTILS::AddEnumClassToParameterList<SCATRA::Action>(
      "action", SCATRA::Action::calc_vreman_scatra, calc_vreman_params_scatra);
  calc_vreman_params_scatra.set("col_filtered_phi", col_filtered_phi_);
  calc_vreman_params_scatra.set("col_filtered_phi2", col_filtered_phi2_);
  calc_vreman_params_scatra.set("col_filtered_phiexpression", col_filtered_phiexpression_);
  calc_vreman_params_scatra.set("col_filtered_alphaijsc", col_filtered_alphaijsc_);
  // loop all elements on this proc (excluding ghosted ones)
  Teuchos::RCP<CORE::LINALG::SerialDenseVector> Dt_num_denom =
      Teuchos::rcp(new CORE::LINALG::SerialDenseVector(2));
  // call loop over elements (assemble nothing)
  scatradiscret_->EvaluateScalars(calc_vreman_params_scatra, Dt_num_denom);
  scatradiscret_->ClearState();
  dt_numerator_volumeav = (*Dt_num_denom)[0];
  dt_denominator_volumeav = (*Dt_num_denom)[1];
  if (sqrt(dt_denominator_volumeav * dt_denominator_volumeav) <
      1.0e-12)  // the denominator might also become negative
    Dt = 0.0;   // constant vreman
  else
    Dt = 1.0 / diffus * dt_numerator_volumeav / dt_denominator_volumeav;

  // remark:
  // Dt does not contain Cv, since Cv cancells out during the calculation of the subgrid diffusivity

  Teuchos::ParameterList* modelparams_scatra = &(extraparams.sublist("TURBULENCE MODEL"));

  modelparams_scatra->set<double>("Dt_vreman", Dt);
  params_.set<double>("Dt_vreman", Dt);
  return;
}  // end FLD::Vreman::DynVremanComputeDt()

FOUR_C_NAMESPACE_CLOSE