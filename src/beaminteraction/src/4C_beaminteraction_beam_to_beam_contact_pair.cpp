// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_beaminteraction_beam_to_beam_contact_pair.hpp"

#include "4C_beam3_euler_bernoulli.hpp"
#include "4C_beam3_kirchhoff.hpp"
#include "4C_beam3_reissner.hpp"
#include "4C_beam3_spatial_discretization_utils.hpp"
#include "4C_beamcontact_input.hpp"
#include "4C_beaminteraction_beam_to_beam_contact_defines.hpp"
#include "4C_beaminteraction_beam_to_beam_contact_params.hpp"
#include "4C_beaminteraction_beam_to_beam_contact_tangentsmoothing.hpp"
#include "4C_beaminteraction_beam_to_beam_contact_utils.hpp"
#include "4C_beaminteraction_beam_to_beam_contact_variables.hpp"
#include "4C_beaminteraction_contact_pair.hpp"
#include "4C_beaminteraction_contact_params.hpp"
#include "4C_contact_input.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_structure_timint_impl.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_TimeMonitor.hpp>

#include <memory>

FOUR_C_NAMESPACE_OPEN

// TODO: Abfangen, dass Kontaktpunkte am Elementuebergang zweimal ausgewertet werden!!!

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::BeamToBeamContactPair()
    : BeamContactPair(),
      r1_(0.0),
      r2_(0.0),
      maxactivegap_(0.0),
      maxsegdist1_(0.0),
      maxsegdist2_(0.0),
      numseg1_(1),
      numseg2_(1),
      boundarynode1_(std::make_pair(0, 0)),
      boundarynode2_(std::make_pair(0, 0))
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::setup()
{
  check_init();

  // call setup of base class first
  BeamContactPair::setup();


  for (unsigned int i = 0; i < 3 * numnodes * numnodalvalues; i++)
  {
    ele1pos_(i) = 0.0;
    ele2pos_(i) = 0.0;
  }

  r1_ = BeamInteraction::calc_ele_radius(element1());
  r2_ = BeamInteraction::calc_ele_radius(element2());

  maxactivegap_ = get_max_active_dist();

  // In case we want to apply a segment-based integration at the endpoints of the physical beam (in
  // order to avoid strong discontinuities in the integrand) we have to check, if a master beam
  // element node coincides with a beams endpoint!
  bool determine_neighbors = false;
  if (params()->beam_to_beam_contact_params()->end_point_penalty()) determine_neighbors = true;

#ifdef ENDPOINTSEGMENTATION
  determine_neighbors = true;
#endif

  if (determine_neighbors)
  {
    neighbors1_ = BeamInteraction::Beam3TangentSmoothing::determine_neighbors(element1());
    neighbors2_ = BeamInteraction::Beam3TangentSmoothing::determine_neighbors(element2());

    bool leftboundarynode1 = false;
    bool rightboundarynode1 = false;

    if (neighbors1_->get_left_neighbor() == nullptr) leftboundarynode1 = true;

    if (neighbors1_->get_right_neighbor() == nullptr) rightboundarynode1 = true;

    boundarynode1_ = std::make_pair(leftboundarynode1, rightboundarynode1);

    bool leftboundarynode2 = false;
    bool rightboundarynode2 = false;

    if (neighbors2_->get_left_neighbor() == nullptr) leftboundarynode2 = true;

    if (neighbors2_->get_right_neighbor() == nullptr) rightboundarynode2 = true;

    boundarynode2_ = std::make_pair(leftboundarynode2, rightboundarynode2);
  }


  // TODO maybe we can even cast the class variables element1_ and element2_ to Beam3Base here in
  // Constructor?! Calculate initial length of beam elements
  const Discret::Elements::Beam3Base* ele1ptr =
      dynamic_cast<const Discret::Elements::Beam3Base*>(element1());
  double l1 = ele1ptr->ref_length();
  const Discret::Elements::Beam3Base* ele2ptr =
      dynamic_cast<const Discret::Elements::Beam3Base*>(element2());
  double l2 = ele2ptr->ref_length();

  if (element1()->element_type() != element2()->element_type())
    FOUR_C_THROW(
        "The class BeamToBeamContactPair only works for contact pairs of the same beam element "
        "type!");

  if (element1()->id() >= element2()->id())
  {
    FOUR_C_THROW(
        "Element 1 has to have the smaller element-ID. Ele1GID: {}, Ele2GID: {}. Adapt your "
        "contact search!",
        element1()->id(), element2()->id());
  }

  cpvariables_.resize(0);
  gpvariables_.resize(0);
  epvariables_.resize(0);

  if (params()->beam_to_beam_contact_params()->gap_shift() != 0.0 and
      params()->beam_to_beam_contact_params()->penalty_law() != BeamContact::pl_lpqp)
    FOUR_C_THROW("BEAMS_GAPSHIFTPARAM only possible for penalty law LinPosQuadPen!");

  double perpshiftangle1 =
      params()->beam_to_beam_contact_params()->beam_to_beam_perp_shifting_angle1();
  double parshiftangle2 =
      params()->beam_to_beam_contact_params()->beam_to_beam_perp_shifting_angle2();

  double segangle = params()->beam_to_beam_contact_params()->segmentation_angle();

  double safetyfac = 1.5;
  // Determine bound for search of large-angle contact pairs
  deltalargeangle_ = perpshiftangle1 - safetyfac * 2 * segangle;

  // In case of a negative value of deltalargeangle_ all pairs have to be evaluate by the
  // large-angle contact formulation
  if (deltalargeangle_ <= 0) deltalargeangle_ = 0;

  // Determine bound for search of small-angle contact pairs
  deltasmallangle_ = parshiftangle2 + safetyfac * 2 * segangle;

  // Check, if we have enough gauss points in order to find every contact point!!
  // Calculate maximal length distance between two gauss points (the factor 1.5 takes into account
  // the not evenly distributed locations of the Gauss points -> this does hold for a number of
  // Gauss points <= 10!!!)
  Core::FE::IntegrationPoints1D gausspoints = Core::FE::IntegrationPoints1D(BEAMCONTACTGAUSSRULE);
  int intintervals = params()->beam_to_beam_contact_params()->num_integration_intervals();

  //  double deltal1=1.5*l1/(intintervals*gausspoints.nquad);

  if (l2 + 1.0e-8 < l1 / intintervals)
    FOUR_C_THROW(
        "Length of second (master) beam has to be larger than length of one integration interval "
        "on first (slave) beam!");

  if (gausspoints.nquad > 10) FOUR_C_THROW("So far, not more than 10 Gauss points are allowed!");

  // TODO We have not considered the factor of 4 occurring in the formula of maximal Gauss point
  // distance, therefore we have an additional safety factor here...
  // TODO!!!!
  //  if((deltal1>R1_/sin(parshiftangle2)) and !beamsdebug)
  //    FOUR_C_THROW("Not enough Gauss points crossing of beams possible!!!");


  issetup_ = true;


  // *********************** DEBUG ************************************************
  //  std::cout << "\nSuccessful Creation&Init&Setup of";
  //  print(std::cout);
  // ******************* END DEBUG ************************************************
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::pre_evaluate()
{
  // do nothing
  return;
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Evaluate the element (public)                             meier 02/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
bool BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::evaluate(
    Core::LinAlg::SerialDenseVector* forcevec1, Core::LinAlg::SerialDenseVector* forcevec2,
    Core::LinAlg::SerialDenseMatrix* stiffmat11, Core::LinAlg::SerialDenseMatrix* stiffmat12,
    Core::LinAlg::SerialDenseMatrix* stiffmat21, Core::LinAlg::SerialDenseMatrix* stiffmat22)
{
  unsigned int dim1 = 3 * numnodes * numnodalvalues;
  unsigned int dim2 = 3 * numnodes * numnodalvalues;

  // resize and initialize variables to zero
  if (forcevec1 != nullptr) forcevec1->size(dim1);
  if (forcevec2 != nullptr) forcevec2->size(dim2);

  if (stiffmat11 != nullptr) stiffmat11->shape(dim1, dim1);
  if (stiffmat12 != nullptr) stiffmat12->shape(dim1, dim2);
  if (stiffmat21 != nullptr) stiffmat21->shape(dim2, dim1);
  if (stiffmat22 != nullptr) stiffmat22->shape(dim2, dim2);


  //**********************************************************************
  // Evaluation of contact forces and stiffness
  //**********************************************************************
  // (1) Closest Point Projection (CPP)
  //     -> find closest point where contact forces are evaluated
  // (2) Compute some auxiliary quantities
  //     -> normal vector, gap, shape functions, contact flag,
  //     -> linearizations of all geometric quantities
  // (3) Compute contact forces and stiffness
  //     -> stiffness terms are directly assembled to global matrix
  //     -> contact forces are only returned as global vector
  // (4) Perform some finite difference checks
  //     -> only if the flag BEAMCONTACTFDCHECKS is defined
  //***************Get some parameters in the beginning*******************

  // All updates that have to be done in every iteration have do be done here,
  // since most of the elements leave directly after the closest point projection!
  clear_class_variables();

  double pp = params()->beam_to_beam_contact_params()->beam_to_beam_point_penalty_param();

  // Subdevide the two elements in segments with linear approximation
  std::vector<Core::LinAlg::Matrix<3, 1, double>> endpoints1(0);
  std::vector<Core::LinAlg::Matrix<3, 1, double>> endpoints2(0);

  // TODO: remove 0 and 1: So far the number 0 and 1 are used in order to distinguish
  // between element 1 and element 2. However, this is only necessary for debugging purposes
  // and can be removed later!
  maxsegdist1_ = create_segments(element1(), endpoints1, numseg1_, 0);
  maxsegdist2_ = create_segments(element2(), endpoints2, numseg2_, 1);

  // Make pairs of close segments: Most of the pairs are already sorted out
  // at this point and don't have to be considered further in the following CPP
  // Additionally, we store the relative orientation of the pairs

  std::map<std::pair<int, int>, Core::LinAlg::Matrix<3, 1, double>> closelargeanglesegments;
  std::map<std::pair<int, int>, Core::LinAlg::Matrix<3, 1, double>> closesmallanglesegments;
  std::vector<std::pair<int, int>> closeendpointsegments(0);
  closelargeanglesegments.clear();
  closesmallanglesegments.clear();

  bool endpoint_penalty = params()->beam_to_beam_contact_params()->end_point_penalty();

// Sub-division of contact elements in search segments or not?
#ifndef NOSEGMENTATION
  get_close_segments(endpoints1, endpoints2, closesmallanglesegments, closelargeanglesegments,
      closeendpointsegments, maxactivegap_);
#else
  Core::LinAlg::Matrix<3, 1, double> segmentdata(true);
  segmentdata(0) = 0.0;  // segment angle
  segmentdata(1) = 0.0;  // eta1_seg
  segmentdata(2) = 0.0;  // eta2_seg
  closesmallanglesegments[std::make_pair(0, 0)] = segmentdata;
  closelargeanglesegments[std::make_pair(0, 0)] = segmentdata;

  if (endpoint_penalty)
  {
    if (boundarynode1_.first or boundarynode1_.second or boundarynode2_.first or
        boundarynode2_.second)
      closeendpointsegments.push_back(std::make_pair(0, 0));
  }
#endif

  //**********************************************************************
  // (1) Closest Point Projection for all close large angle segments(CPP)
  //**********************************************************************

  // Treat large-angle contact pairs if existing
  if (closelargeanglesegments.size() > 0)
  {
    // Get active large angle pairs (valid closest point projections) and create vector of
    // cpvariables_
    get_active_large_angle_pairs(endpoints1, endpoints2, closelargeanglesegments, pp);

    // Evaluate contact contribution of large-angle-contact (residual and stiffness) for all closest
    // points found before
    evaluate_active_large_angle_pairs(
        forcevec1, forcevec2, stiffmat11, stiffmat12, stiffmat21, stiffmat22);
  }

  // Treat small angle contact pairs if existing
  if (closesmallanglesegments.size() > 0)
  {
#ifndef ENDPOINTSEGMENTATION
    // Get active small angle pairs (valid Gauss points) and create vector of gpvariables_
    get_active_small_angle_pairs(closesmallanglesegments);

    // Evaluate contact contribution of small-angle-contact (residual and stiffness) for all closest
    // points found before
    evaluate_active_small_angle_pairs(
        forcevec1, forcevec2, stiffmat11, stiffmat12, stiffmat21, stiffmat22);
#else
    // In case of endpoint segmentation some additional quantities have to be transferred between
    // the methods get_active_small_angle_pairs() and evaluate_active_small_angle_pairs().
    std::pair<int, int> iminmax = std::make_pair(0, 0);
    std::pair<bool, bool> leftrightsolutionwithinsegment = std::make_pair(false, false);
    std::pair<double, double> eta1_leftrightboundary = std::make_pair(0.0, 0.0);

    // Get active small angle pairs (valid Gauss points) and create vector of gpvariables_
    get_active_small_angle_pairs(closesmallanglesegments, &iminmax, &leftrightsolutionwithinsegment,
        &eta1_leftrightboundary);

    // Evaluate contact contribution of small-angle-contact (residual and stiffness) for all closest
    // points found before
    evaluate_active_small_angle_pairs(
        stiffmatrix, fint, &iminmax, &leftrightsolutionwithinsegment, &eta1_leftrightboundary);
#endif
  }

  if (endpoint_penalty)
  {
    // Treat endpoint contact pairs if existing
    if (closeendpointsegments.size() > 0)
    {
      // Get active endpoint pairs and create vector of epvariables_
      get_active_end_point_pairs(closeendpointsegments, pp);

      // Evaluate contact contribution of endpoint-contact (residual and stiffness) for all closest
      // points found before
      evaluate_active_end_point_pairs(
          forcevec1, forcevec2, stiffmat11, stiffmat12, stiffmat21, stiffmat22);
    }
  }

  // Check whether non-zero contributions are returned and hence need assembly in system variables
  // Todo avoid computation of norm here and rather use the maximum of all 'DoNotAssemble' flags
  //      which is set for each cp, gp and ep of this element pair (in evaluate_fc_contact
  //      evaluate_stiffc_contact)
  if ((forcevec1 != nullptr and forcevec1->normInf() > 0) or
      (forcevec2 != nullptr and forcevec2->normInf() > 0) or
      (stiffmat11 != nullptr and stiffmat11->normInf() > 0) or
      (stiffmat12 != nullptr and stiffmat12->normInf() > 0) or
      (stiffmat21 != nullptr and stiffmat21->normInf() > 0) or
      (stiffmat22 != nullptr and stiffmat22->normInf() > 0))
  {
    return true;
  }
  else
  {
    return false;
  }
}
/*----------------------------------------------------------------------*
 |  end: Evaluate the element
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Get active large angle pairs                             meier 10/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::get_active_large_angle_pairs(
    std::vector<Core::LinAlg::Matrix<3, 1, double>>& endpoints1,
    std::vector<Core::LinAlg::Matrix<3, 1, double>>& endpoints2,
    std::map<std::pair<int, int>, Core::LinAlg::Matrix<3, 1, double>>& closelargeanglesegments,
    const double& pp)
{
  std::map<std::pair<int, int>, Core::LinAlg::Matrix<3, 1, double>>::iterator iter;

  for (iter = closelargeanglesegments.begin(); iter != closelargeanglesegments.end(); ++iter)
  {
    Core::LinAlg::Matrix<3, 1, double> segmentdata = iter->second;
    std::pair<int, int> leftpoint_ids = iter->first;
    int nseg1 = endpoints1.size() - 1;
    int nseg2 = endpoints2.size() - 1;
    int segid1 = leftpoint_ids.first;
    int segid2 = leftpoint_ids.second;
    double l1 = 2.0 / nseg1;
    double l2 = 2.0 / nseg2;
    double eta_left1 = -1.0 + segid1 * l1;
    double eta_left2 = -1.0 + segid2 * l2;

    std::pair<TYPE, TYPE> closestpoint(0.0, 0.0);
    bool validpairfound = false;

    // The method closest_point_projection() only delivers a valid solution (validpairfound=true),
    // if eta1 \in [eta_left1,eta_left1+l1], eta2 \in [eta_left2,eta_left2+l2] and
    // gap<maxactivegap_!
    validpairfound = closest_point_projection(
        eta_left1, eta_left2, l1, l2, segmentdata, closestpoint, segid1, segid2);

    // With the following block we sort out identical contact points that occur more than once
    // within this element pair -> this is possible, when the contact point lies on the boundary
    // between two segments!
    bool already_found = false;

    for (unsigned int i = 0; i < cpvariables_.size(); i++)
    {
      double eta1_eval = Core::FADUtils::cast_to_double(cpvariables_[i]->get_cp().first);
      double eta2_eval = Core::FADUtils::cast_to_double(cpvariables_[i]->get_cp().second);

      if (std::fabs(eta1_eval - Core::FADUtils::cast_to_double(closestpoint.first)) <
              XIETARESOLUTIONFAC * XIETAITERATIVEDISPTOL and
          std::fabs(eta2_eval - Core::FADUtils::cast_to_double(closestpoint.second)) <
              XIETARESOLUTIONFAC * XIETAITERATIVEDISPTOL)
        already_found = true;
    }

    if (validpairfound and !already_found)
    {
      std::pair<int, int> integration_ids = std::make_pair(-2, -2);
      cpvariables_.push_back(std::make_shared<BeamToBeamContactVariables<numnodes, numnodalvalues>>(
          closestpoint, leftpoint_ids, integration_ids, pp, 1.0));
    }
  }
}
/*----------------------------------------------------------------------*
 |  end: Get active large angle pairs
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Evaluate active large angle pairs                        meier 10/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes,
    numnodalvalues>::evaluate_active_large_angle_pairs(Core::LinAlg::SerialDenseVector* forcevec1,
    Core::LinAlg::SerialDenseVector* forcevec2, Core::LinAlg::SerialDenseMatrix* stiffmat11,
    Core::LinAlg::SerialDenseMatrix* stiffmat12, Core::LinAlg::SerialDenseMatrix* stiffmat21,
    Core::LinAlg::SerialDenseMatrix* stiffmat22)
{
  for (int numcp = 0; numcp < (int)cpvariables_.size(); numcp++)
  {
    //**********************************************************************
    // (2) Compute some auxiliary quantities
    //**********************************************************************

    // vectors for shape functions and their derivatives
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N1(true);       // = N1
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N2(true);       // = N2
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N1_xi(true);    // = N1,xi
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N2_xi(true);    // = N2,eta
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N1_xixi(true);  // = N1,xixi
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N2_xixi(true);  // = N2,etaeta

    // coords and derivatives of the two contacting points
    Core::LinAlg::Matrix<3, 1, TYPE> r1(true);       // = r1
    Core::LinAlg::Matrix<3, 1, TYPE> r2(true);       // = r2
    Core::LinAlg::Matrix<3, 1, TYPE> r1_xi(true);    // = r1,xi
    Core::LinAlg::Matrix<3, 1, TYPE> r2_xi(true);    // = r2,eta
    Core::LinAlg::Matrix<3, 1, TYPE> r1_xixi(true);  // = r1,xixi
    Core::LinAlg::Matrix<3, 1, TYPE> r2_xixi(true);  // = r2,etaeta
    Core::LinAlg::Matrix<3, 1, TYPE> delta_r(true);  // = r1-r2

    TYPE eta1 = cpvariables_[numcp]->get_cp().first;
    TYPE eta2 = cpvariables_[numcp]->get_cp().second;

#ifdef AUTOMATICDIFF
    BeamContact::SetFADParCoordDofs<numnodes, numnodalvalues>(eta1, eta2);
#endif

    // update shape functions and their derivatives
    get_shape_functions(N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi, eta1, eta2);
    // update coordinates and derivatives of contact points
    compute_coords_and_derivs(
        r1, r2, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi);

    // call function to compute scaled normal and gap of contact point
    compute_normal(r1, r2, r1_xi, r2_xi, cpvariables_[numcp], 0);

    // call function to compute penalty force
    calc_penalty_law(*cpvariables_[numcp]);

    // get shift angles from input file
    double perpshiftangle1 =
        params()->beam_to_beam_contact_params()->beam_to_beam_perp_shifting_angle1();
    double perpshiftangle2 =
        params()->beam_to_beam_contact_params()->beam_to_beam_perp_shifting_angle2();

    // call function to compute scale factor of penalty parameter
    calc_perp_penalty_scale_fac(
        *cpvariables_[numcp], r1_xi, r2_xi, perpshiftangle1, perpshiftangle2);

    // In case of large-angle-contact, the length specific energy and the 'real' energy are
    // identical
    double lengthspec_energy = Core::FADUtils::cast_to_double(cpvariables_[numcp]->get_energy());
    cpvariables_[numcp]->set_integrated_energy(lengthspec_energy);

    //    std::cout << "cpvariables_[numcp]->GetNormal(): " << cpvariables_[numcp]->GetNormal() <<
    //    std::endl; std::cout << "numcp: " << numcp << std::endl; std::cout << "xi: " <<
    //    cpvariables_[numcp]->GetCP().first.val() << std::endl; std::cout << "eta: " <<
    //    cpvariables_[numcp]->GetCP().second.val() << std::endl; std::cout << "gap: " <<
    //    cpvariables_[numcp]->GetGap().val() << std::endl; std::cout << "angle: " <<
    //    cpvariables_[numcp]->get_angle()/M_PI*180.0 << std::endl; std::cout << "r1_xi: " << r1_xi
    //    << std::endl; std::cout << "r2_xi: " << r2_xi << std::endl; std::cout << "|r1_xi|: " <<
    //    r1_xi.Norm2() << std::endl; std::cout << "|r2_xi|: " << r2_xi.Norm2() << std::endl;
    //    std::cout << "r1_xi*r2_xi: " << Core::FADUtils::ScalarProduct(r1_xi,r2_xi) << std::endl;
    //    std::cout << "cpvariables_[numcp]->Getfp(): " << cpvariables_[numcp]->Getfp() <<
    //    std::endl;

    // call function to compute contact contribution to residual vector
    if (forcevec1 != nullptr and forcevec2 != nullptr)
      evaluate_fc_contact(*forcevec1, *forcevec2, r1, r2, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2,
          N1_xi, N2_xi, *cpvariables_[numcp], 1.0, true, false, false, false);

    // call function to compute contact contribution to stiffness matrix
    if (stiffmat11 != nullptr and stiffmat12 != nullptr and stiffmat21 != nullptr and
        stiffmat22 != nullptr)
      evaluate_stiffc_contact(*stiffmat11, *stiffmat12, *stiffmat21, *stiffmat22, r1, r2, r1_xi,
          r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi, *cpvariables_[numcp],
          1.0, true, false, false, false);
  }
}
/*----------------------------------------------------------------------*
 |  end: Evaluate active large angle pairs
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Get active small angle pairs                             meier 10/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::get_active_small_angle_pairs(
    std::map<std::pair<int, int>, Core::LinAlg::Matrix<3, 1, double>>& closesmallanglesegments,
    std::pair<int, int>* iminmax, std::pair<bool, bool>* leftrightsolutionwithinsegment,
    std::pair<double, double>* eta1_leftrightboundary)
{
  // lengths in parameter space of created segments
  double l1 = 2.0 / numseg1_;
  double l2 = 2.0 / numseg2_;

  int numpairs = closesmallanglesegments.size();
  std::vector<std::pair<double, double>> inversepairs(numpairs, std::make_pair(0.0, 0.0));
  int pairiter = 0;
  int intintervals = params()->beam_to_beam_contact_params()->num_integration_intervals();

  std::map<std::pair<int, int>, Core::LinAlg::Matrix<3, 1, double>>::iterator iter;

  for (iter = closesmallanglesegments.begin(); iter != closesmallanglesegments.end(); ++iter)
  {
    std::pair<int, int> segment_ids = iter->first;
    int segid1 = segment_ids.first;
    int segid2 = segment_ids.second;
    double eta1_segleft = -1.0 + segid1 * l1;
    double eta2_segleft = -1.0 + segid2 * l2;

    inversepairs[numpairs - pairiter - 1] = std::make_pair(eta1_segleft, eta2_segleft);
    pairiter++;
  }
  int imin = 0;
  int imax = intintervals - 1;

#ifdef ENDPOINTSEGMENTATION
  double eta1_leftboundary = -1.0;
  double eta1_rightboundary = 1.0;
  double leftintervallength = 2.0 / intintervals;
  double rightintervallength = 2.0 / intintervals;
  bool leftsolutionwithinsegment = false;
  bool rightsolutionwithinsegment = false;

  if (boundarynode2_.first)
  {
    for (iter = closesmallanglesegments.begin(); iter != closesmallanglesegments.end(); ++iter)
    {
      std::pair<int, int> segment_ids = iter->first;
      int segid2 = segment_ids.second;

      if (segid2 == 0)
      {
        int segid1 = segment_ids.first;
        double eta1_segleft = -1.0 + segid1 * l1;
        double eta2_segleft = -1.0;
        double eta1_boundary_trial = 0.0;
        bool dummy;

        bool solutionwithinsegment = false;

        double gap_dummy = 0.0;
        double alpha_dummy = 0.0;

#ifndef CHANGEENDPOINTPROJECTION
        solutionwithinsegment = point_to_line_projection(eta2_segleft, eta1_segleft, l1,
            eta1_boundary_trial, gap_dummy, alpha_dummy, dummy, true, true);
#else
        solutionwithinsegment = point_to_line_projection(eta2_segleft, eta1_segleft, l1,
            eta1_boundary_trial, gap_dummy, alpha_dummy, dummy, true, true, true);
#endif

        if (solutionwithinsegment)
        {
          // Determine if the projection eta1_boundary_trial is a left boundary of the integration
          // segment or a right boundary of the integration segment This is done in the following
          // way: First, we determine the tangent of the master boundary node in a way, such that
          // the tangent points into the elements interior. Then, we determine the tangent on slave
          // beam at the projection point eta1_boundary_trial. This tangent automatically points
          // into positive eta1-direction=integration direction. Thus, if the scalar product of
          // these two tangents is positive, the master element evolves in positive eta1-direction
          // and consequently, eta1_boundary_trial is the left boundary of the integration segment.
          // If the scalar product is negative, eta1_boundary_trial is the right boundary of the
          // integration segment
          Core::LinAlg::Matrix<3, 1, TYPE> inward_tangent_master = r_xi(eta2_segleft, Element2());
          Core::LinAlg::Matrix<3, 1, TYPE> tangent_slave = r_xi(eta1_boundary_trial, Element1());
          double orientation = Core::FADUtils::cast_to_double(
              Core::FADUtils::ScalarProduct(inward_tangent_master, tangent_slave));
          if (orientation > 0)  // left boundary
          {
            leftsolutionwithinsegment = true;
            eta1_leftboundary = eta1_boundary_trial;
            // determine ID of integration interval in which the point eta1_leftboundary lies
            imin = BeamContact::GetIntervalId(eta1_leftboundary, intintervals, true);
            // get length of segmented integration interval
            leftintervallength = -1.0 + (imin + 1) * 2.0 / intintervals - eta1_leftboundary;
            break;
          }
          else if (orientation < 0)  // right boundary
          {
            rightsolutionwithinsegment = true;
            eta1_rightboundary = eta1_boundary_trial;
            // determine ID of integration interval in which the point eta1_leftboundary lies
            imax = BeamContact::GetIntervalId(eta1_rightboundary, intintervals, false);
            // get length of segmented integration interval
            rightintervallength = eta1_rightboundary - (-1.0 + imax * 2.0 / intintervals);
            break;
          }
          else  // This can only happen, if both beams are exactly perpendicular AND the master beam
                // endpoint projects perpendicular on the slave beam!
            FOUR_C_THROW("The very unlikely case orientation==0 is not implemented so far!");
        }
      }
    }
  }
  if (boundarynode2_.second)
  {
    for (iter = closesmallanglesegments.begin(); iter != closesmallanglesegments.end(); ++iter)
    {
      std::pair<int, int> segment_ids = iter->first;
      int segid2 = segment_ids.second;

      if (segid2 == numseg2_ - 1)
      {
        int segid1 = segment_ids.first;
        double eta1_segleft = -1.0 + segid1 * l1;
        double eta2_segright = 1.0;
        double eta1_boundary_trial = 0.0;
        bool dummy;

        bool solutionwithinsegment = false;

        double gap_dummy = 0.0;
        double alpha_dummy = 0.0;

#ifndef CHANGEENDPOINTPROJECTION
        solutionwithinsegment = point_to_line_projection(eta2_segright, eta1_segleft, l1,
            eta1_boundary_trial, gap_dummy, alpha_dummy, dummy, true, true);
#else
        solutionwithinsegment = point_to_line_projection(eta2_segright, eta1_segleft, l1,
            eta1_boundary_trial, gap_dummy, alpha_dummy, dummy, true, true, true);
#endif

        if (solutionwithinsegment)
        {
          // Determine if the projection eta1_boundary_trial is a left boundary of the integration
          // segment or a right boundary of the integration segment This is done in the following
          // way: First, we determine the tangent of the master boundary node in a way, such that
          // the tangent points into the elements interior. Then, we determine the tangent on slave
          // beam at the projection point eta1_boundary_trial. This tangent automatically points
          // into positive eta1-direction=integration direction. Thus, if the scalar product of
          // these two tangents is positive, the master element evolves in positive eta1-direction
          // and consequently, eta1_boundary_trial is the left boundary of the integration segment.
          // If the scalar product is negative, eta1_boundary_trial is the right boundary of the
          // integration segment
          Core::LinAlg::Matrix<3, 1, TYPE> inward_tangent_master = r_xi(eta2_segright, Element2());
          // Scale tangent of right element node (eta2=1.0) in order to get inward tangent!
          inward_tangent_master.scale(-1.0);
          Core::LinAlg::Matrix<3, 1, TYPE> tangent_slave = r_xi(eta1_boundary_trial, Element1());
          double orientation = Core::FADUtils::cast_to_double(
              Core::FADUtils::ScalarProduct(inward_tangent_master, tangent_slave));
          if (orientation > 0)  // left boundary
          {
            if (leftsolutionwithinsegment)
              FOUR_C_THROW(
                  "Something went wrong here: both boundary nodes of the master beam (discretized "
                  "by one finite element?!?) are projected as left boundary of the integration "
                  "segment!");

            leftsolutionwithinsegment = true;
            eta1_leftboundary = eta1_boundary_trial;
            // determine ID of integration interval in which the point eta1_leftboundary lies
            imin = BeamContact::GetIntervalId(eta1_leftboundary, intintervals, true);
            // get length of segmented integration interval
            leftintervallength = -1.0 + (imin + 1) * 2.0 / intintervals - eta1_leftboundary;
            break;
          }
          else if (orientation < 0)  // right boundary
          {
            if (rightsolutionwithinsegment)
              FOUR_C_THROW(
                  "Something went wrong here: both boundary nodes of the master beam (discretized "
                  "by one finite element?!?) are projected as right boundary of the integration "
                  "segment!");

            rightsolutionwithinsegment = true;
            eta1_rightboundary = eta1_boundary_trial;
            // determine ID of integration interval in which the point eta1_leftboundary lies
            imax = BeamContact::GetIntervalId(eta1_rightboundary, intintervals, false);
            // get length of segmented integration interval
            rightintervallength = eta1_rightboundary - (-1.0 + imax * 2.0 / intintervals);
            break;
          }
          else  // This can only happen, if both beams are exactly perpendicular AND the master beam
                // endpoint projects perpendicular on the slave beam!
            FOUR_C_THROW("The very unlikely case orientation==0 is not implemented so far!");
        }
      }
    }
  }
  if (leftsolutionwithinsegment and rightsolutionwithinsegment and imin == imax)
    FOUR_C_THROW(
        "It is not possible to cut an integration interval from both sides, choose a larger value "
        "intintervals!");
#endif

  // gaussian points
  Core::FE::IntegrationPoints1D gausspoints = Core::FE::IntegrationPoints1D(BEAMCONTACTGAUSSRULE);

  // loop over all integration intervals
  for (int interval = imin; interval <= imax; interval++)
  {
    // Calculate parameter bounds of considered integration interval
    double eta1_min = -1.0 + interval * 2.0 / intintervals;
    double eta1_max = -1.0 + (interval + 1) * 2.0 / intintervals;

    // Get jacobi factor of considered interval
    TYPE jacobi_interval = 1.0;

// standard case of equidistant intervals
#ifndef ENDPOINTSEGMENTATION
    // map from segment coordinate xi to element coordinate eta
    jacobi_interval = 1.0 / intintervals;
// case of smaller integration intervals due to segmentation at the beams endpoints
#else
    if (interval == imin and leftsolutionwithinsegment)
    {
      jacobi_interval = leftintervallength / 2.0;
    }
    else if (interval == imax and rightsolutionwithinsegment)
    {
      jacobi_interval = rightintervallength / 2.0;
      // std::cout << "rightintervallength: " << rightintervallength << std::endl;
    }
    else
    {
      jacobi_interval = 1.0 / intintervals;
    }
#endif

    std::vector<std::pair<double, double>> curintsegpairs(0);
    int size = inversepairs.size();

    // All segment pairs for which the segment on the slave beam 1 intersects with the considered
    // integration interval are filtered out and stored in the vector curintsegpairs. These pairs
    // are relevant for the integration procedure on the current interval.
    for (int k = 0; k < size; k++)
    {
      double eta1_segleft = (inversepairs[size - 1 - k]).first;
      double eta1_segright = eta1_segleft + l1;
      // Since the vector inversepairs is sorted with respect to the location of the slave segment
      // (the slave segment with the lowest bounding parameter coordinates eta1_segleft and
      // eta1_segright lie on the last position of the vector inversepairs), it is sufficient to
      // start with the last element and leave the k-loop as soon as we have found the first segment
      // pair without intersection. This procedure only works, if we delete a segment pair as soon
      // as we realize that it will not be relevant for the next integration interval anymore, see
      // comment at (*).
      if (eta1_segleft < eta1_max + 1.0e-10)
      {
        // store relevant pairs in new vector
        curintsegpairs.push_back(inversepairs[size - 1 - k]);

        //(*) In case eta1_segright (the largest parameter coordinate lying within the slave
        // segment) is smaller than eta1_max(the upper bound of the integration interval), the
        // considered segment will not be relevant for the next integration interval at i+1 and can
        // be deleted.
        if (eta1_segright < eta1_max - 1.0e-10)
        {
          inversepairs.pop_back();
        }
      }
      // In case we have no relevant segment pair, we will leave the loop already after the first
      // iteration!
      else
      {
        break;
      }
    }

    // If segments exist, evaluate the corresponding Gauss points
    if (curintsegpairs.size() > 0)
    {
      // loop over Gauss point of considered integration interval
      for (int numgp = 0; numgp < gausspoints.nquad; ++numgp)
      {
        // integration points in parameter space and weights
        const double xi = gausspoints.qxg[numgp][0];

        // Get Gauss point coordinate at slave element
        double eta1_slave = 0.0;

// standard case of equidistant intervals
#ifndef ENDPOINTSEGMENTATION
        // map from segment coordinate xi to element coordinate eta
        eta1_slave = eta1_min + (1.0 + xi) / intintervals;
// case of smaller integration intervals due to segmentation at the beams endpoints
#else
        if (interval == imin and leftsolutionwithinsegment)
        {
          eta1_slave = eta1_leftboundary + (1.0 + xi) / 2.0 * leftintervallength;
        }
        else if (interval == imax and rightsolutionwithinsegment)
        {
          eta1_slave = eta1_min + (1.0 + xi) / 2.0 * rightintervallength;
        }
        else
        {
          eta1_slave = eta1_min + (1.0 + xi) / intintervals;
        }
#endif

        for (int k = 0; k < (int)curintsegpairs.size(); k++)
        {
          double eta1_segleft = (curintsegpairs[k]).first;
          double eta1_segright = eta1_segleft + l1;

          // TODO: This procedure can also be made more efficient by deleting all segments of
          // curintsegpairs which are not relevant for the following Gauss points anymore (see
          // intersection of integration intervals and segment pairs)
          if (BeamInteraction::within_interval(eta1_slave, eta1_segleft, eta1_segright))
          {
            double eta2_segleft = (curintsegpairs[k]).second;
            double eta2_master = 0.0;
            bool pairactive = false;

            double gap_dummy = 0.0;
            double alpha_dummy = 0.0;

            bool solutionwithinsegment = point_to_line_projection(eta1_slave, eta2_segleft, l2,
                eta2_master, gap_dummy, alpha_dummy, pairactive, true);

            if (solutionwithinsegment)
            {
              if (pairactive)
              {
                TYPE eta1 = eta1_slave;
                TYPE eta2 = eta2_master;
                int leftpoint_id1 = BeamInteraction::get_segment_id(eta1_slave, numseg1_);
                int leftpoint_id2 = BeamInteraction::get_segment_id(eta2_master, numseg2_);
                std::pair<TYPE, TYPE> closestpoint(std::make_pair(eta1, eta2));
                std::pair<int, int> integration_ids = std::make_pair(numgp, interval);
                std::pair<int, int> leftpoint_ids = std::make_pair(leftpoint_id1, leftpoint_id2);
                TYPE jacobi = get_jacobi_at_xi(element1(), eta1_slave) * jacobi_interval;

                const double parallel_pp =
                    params()->beam_to_beam_contact_params()->beam_to_beam_line_penalty_param();

                if (parallel_pp < 0.0) FOUR_C_THROW("BEAMS_BTBLINEPENALTYPARAM not set!");

                // Create data container for each Gauss point (in case of small-angle contact the
                // number of the Gauss point [numgp] and the number of the integration interval
                // [interval] are stored in the pair segids_ of the class beamcontactvaribles!)
                gpvariables_.push_back(
                    std::make_shared<BeamToBeamContactVariables<numnodes, numnodalvalues>>(
                        closestpoint, leftpoint_ids, integration_ids, parallel_pp, jacobi));
              }
              // We can leave the k-loop as soon as we have found a valid projection for the given
              // Gauss point eta1_slave
              break;
            }
          }
        }
      }  // for (int numgp=0; numgp<gausspoints.nquad; ++numgp)
    }  // if(curintsegpairs.size()>0)
  }  // for(int interval=imin;interval<=imax;interval++)

#ifdef ENDPOINTSEGMENTATION
  if (iminmax == nullptr or leftrightsolutionwithinsegment == nullptr or
      eta1_leftrightboundary == nullptr)
    FOUR_C_THROW("In case of ENDPOINTSEGMENTATION no NUll pointer should be handed in!!!");

  *iminmax = std::make_pair(imin, imax);
  *leftrightsolutionwithinsegment =
      std::make_pair(leftsolutionwithinsegment, rightsolutionwithinsegment);
  *eta1_leftrightboundary = std::make_pair(eta1_leftboundary, eta1_rightboundary);
#endif
}
/*----------------------------------------------------------------------*
 |  end: Get active small angle pairs
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Evaluate active small angle pairs                        meier 10/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes,
    numnodalvalues>::evaluate_active_small_angle_pairs(Core::LinAlg::SerialDenseVector* forcevec1,
    Core::LinAlg::SerialDenseVector* forcevec2, Core::LinAlg::SerialDenseMatrix* stiffmat11,
    Core::LinAlg::SerialDenseMatrix* stiffmat12, Core::LinAlg::SerialDenseMatrix* stiffmat21,
    Core::LinAlg::SerialDenseMatrix* stiffmat22, std::pair<int, int>* iminmax,
    std::pair<bool, bool>* leftrightsolutionwithinsegment,
    std::pair<double, double>* eta1_leftrightboundary)
{
// Compute linearizations of integration interval boundaries if necessary
#ifdef ENDPOINTSEGMENTATION

  if (iminmax == nullptr or leftrightsolutionwithinsegment == nullptr or
      eta1_leftrightboundary == nullptr)
    FOUR_C_THROW("In case of ENDPOINTSEGMENTATION no NUll pointer should be handed in!!!");

  int imin = (*iminmax).first;
  int imax = (*iminmax).second;

  bool leftsolutionwithinsegment = (*leftrightsolutionwithinsegment).first;
  bool rightsolutionwithinsegment = (*leftrightsolutionwithinsegment).second;

  double eta1_leftboundary = (*eta1_leftrightboundary).first;
  double eta1_rightboundary = (*eta1_leftrightboundary).second;

  Core::LinAlg::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE> delta_xi_R(true);
  Core::LinAlg::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE> delta_xi_L(true);

  if (leftsolutionwithinsegment)
  {
    TYPE eta1_bound = eta1_leftboundary;
    TYPE eta2 = -1.0;
    compute_lin_xi_bound(delta_xi_L, eta1_bound, eta2);
  }
  if (rightsolutionwithinsegment)
  {
    TYPE eta1_bound = eta1_rightboundary;
    TYPE eta2 = 1.0;
    compute_lin_xi_bound(delta_xi_R, eta1_bound, eta2);
  }
#endif

  // gaussian points
  Core::FE::IntegrationPoints1D gausspoints = Core::FE::IntegrationPoints1D(BEAMCONTACTGAUSSRULE);

  // Evaluate all active Gauss points
  for (int numgptot = 0; numgptot < (int)gpvariables_.size(); numgptot++)
  {
    std::pair<TYPE, TYPE> closestpoint = gpvariables_[numgptot]->get_cp();
    TYPE eta1 = closestpoint.first;
    TYPE eta2 = closestpoint.second;

#ifdef AUTOMATICDIFF
    BeamContact::SetFADParCoordDofs<numnodes, numnodalvalues>(eta1, eta2);
#endif

    // vectors for shape functions and their derivatives
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N1(true);       // = N1
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N2(true);       // = N2
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N1_xi(true);    // = N1,xi
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N2_xi(true);    // = N2,eta
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N1_xixi(true);  // = N1,xixi
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N2_xixi(true);  // = N2,etaeta

    // coords and derivatives of the two contacting points
    Core::LinAlg::Matrix<3, 1, TYPE> r1(true);       // = r1
    Core::LinAlg::Matrix<3, 1, TYPE> r2(true);       // = r2
    Core::LinAlg::Matrix<3, 1, TYPE> r1_xi(true);    // = r1,xi
    Core::LinAlg::Matrix<3, 1, TYPE> r2_xi(true);    // = r2,eta
    Core::LinAlg::Matrix<3, 1, TYPE> r1_xixi(true);  // = r1,xixi
    Core::LinAlg::Matrix<3, 1, TYPE> r2_xixi(true);  // = r2,etaeta
    Core::LinAlg::Matrix<3, 1, TYPE> delta_r(true);  // = r1-r2

    // update shape functions and their derivatives
    get_shape_functions(N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi, eta1, eta2);
    // update coordinates and derivatives of contact points
    compute_coords_and_derivs(
        r1, r2, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi);

    // call function to compute scaled normal and gap of contact point
    compute_normal(r1, r2, r1_xi, r2_xi, gpvariables_[numgptot], 1);

    // call function to compute penalty force
    calc_penalty_law(*gpvariables_[numgptot]);

    // get shift angles from input file
    double parshiftangle1 =
        params()->beam_to_beam_contact_params()->beam_to_beam_parallel_shifting_angle1();
    double parshiftangle2 =
        params()->beam_to_beam_contact_params()->beam_to_beam_parallel_shifting_angle2();

    // call function to compute scale factor of penalty parameter
    calc_par_penalty_scale_fac(
        *gpvariables_[numgptot], r1_xi, r2_xi, parshiftangle1, parshiftangle2);

    //    std::cout << "gpvariables_[numgp]->GetNormal(): " << gpvariables_[numgp]->GetNormal() <<
    //    std::endl; std::cout << "numgptot: " << numgptot << std::endl; std::cout << "xi: " <<
    //    gpvariables_[numgptot]->GetCP().first.val() << std::endl; std::cout << "eta: " <<
    //    gpvariables_[numgptot]->GetCP().second.val() << std::endl; std::cout << "gap: " <<
    //    gpvariables_[numgptot]->GetGap().val() << std::endl; std::cout << "angle: " <<
    //    gpvariables_[numgptot]->get_angle()/M_PI*180.0 << std::endl; std::cout << "r1_xi: " <<
    //    r1_xi << std::endl; std::cout << "r2_xi: " << r2_xi << std::endl; std::cout << "|r1_xi|: "
    //    << r1_xi.Norm2() << std::endl; std::cout << "|r2_xi|: " << r2_xi.Norm2() << std::endl;
    //    std::cout << "r1_xi*r2_xi: " << Core::FADUtils::ScalarProduct(r1_xi,r2_xi) << std::endl;
    //    std::cout << "gpvariables_[numgptot]->Getfp(): " << gpvariables_[numgp]->Getfp() <<
    //    std::endl;

    // Determine the integration-segment-local Gauss point-ID of the considered gpvariable
    int numgploc = gpvariables_[numgptot]->get_int_ids().first;

    double weight = gausspoints.qwgt[numgploc];
    TYPE jacobi = gpvariables_[numgptot]->get_jacobi();

#ifdef ENDPOINTSEGMENTATION
    int numinterval = gpvariables_[numgptot]->GetIntIds().second;
#endif

    // The intfac has NOT to be of TYPE FAD in order to deal with non-constant jacobis (in case of
    // ENDPOINTSEGMENTATION) since we explicitly consider the linearization of the jacobi in
    // evaluate_stiffc_contact_int_seg()!
    double intfac = Core::FADUtils::cast_to_double(jacobi) * weight;

    // Convert the length specific energy into a 'real' energy
    // while the length specific energy is used for later calculation, the real (or over the length
    // integrated) energy is a pure output variable and can therefore be of type double!
    double lengthspec_energy = Core::FADUtils::cast_to_double(gpvariables_[numgptot]->get_energy());
    double integrated_energy = lengthspec_energy * intfac;
    gpvariables_[numgptot]->set_integrated_energy(integrated_energy);

    // call function to compute contact contribution to residual vector
    if (forcevec1 != nullptr and forcevec2 != nullptr)
      evaluate_fc_contact(*forcevec1, *forcevec2, r1, r2, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2,
          N1_xi, N2_xi, *gpvariables_[numgptot], intfac, false, true, false, false);

#ifndef ENDPOINTSEGMENTATION
    // call function to compute contact contribution to stiffness matrix
    if (stiffmat11 != nullptr and stiffmat12 != nullptr and stiffmat21 != nullptr and
        stiffmat22 != nullptr)
      evaluate_stiffc_contact(*stiffmat11, *stiffmat12, *stiffmat21, *stiffmat22, r1, r2, r1_xi,
          r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi, *gpvariables_[numgptot],
          intfac, false, true, false, false);
#else
    TYPE jacobi_interval = jacobi / get_jacobi(Element1());
    // In case of segment-based integration, we apply a special FAD linearization technique
    // Case 1: segmentation on left side of integration interval
    if (leftsolutionwithinsegment and numinterval == imin)
    {
      // We need the linearization of the mapping from the element parameter space to the
      // integration interval parameter space:
      // xi_ele=xi_left*(1.0-xi_local)/2.0+xi_right*(1.0+xi_local)/2.0.
      //-> d(xi_ele)/d(xi_left)=(1.0-xi_local)/2.0 and d(xi_ele)/d(xi_right)=(1.0+xi_local)/2.0

      double d_xi_ele_d_xi_left = (1.0 - gausspoints.qxg[numgploc][0]) / 2.0;
      evaluate_stiffc_contact_int_seg(stiffmatrix, delta_xi_L, r1, r2, r1_xi, r2_xi, r1_xixi,
          r2_xixi, N1, N2, N1_xi, N2_xi, gpvariables_[numgptot], intfac, d_xi_ele_d_xi_left,
          -jacobi_interval);
    }
    // Case 2: segmentation on right side of integration interval
    else if (rightsolutionwithinsegment and numinterval == imax)
    {
      double d_xi_ele_d_xi_right = (1.0 + gausspoints.qxg[numgploc][0]) / 2.0;
      evaluate_stiffc_contact_int_seg(stiffmatrix, delta_xi_R, r1, r2, r1_xi, r2_xi, r1_xixi,
          r2_xixi, N1, N2, N1_xi, N2_xi, gpvariables_[numgptot], intfac, d_xi_ele_d_xi_right,
          jacobi_interval);
    }
    // Case 3: No segmentation necessary
    else
    {
      evaluate_stiffc_contact(stiffmatrix, r1, r2, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi,
          N2_xi, N1_xixi, N2_xixi, gpvariables_[numgptot], intfac, false, true, false, false);
    }
#endif
  }
}
/*----------------------------------------------------------------------*
 |  end: Evaluate active small angle pairs
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Get active endpoint pairs                                meier 12/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::get_active_end_point_pairs(
    std::vector<std::pair<int, int>>& closeendpointsegments, const double pp)
{
  for (int i = 0; i < (int)closeendpointsegments.size(); i++)
  {
    int segid1 = closeendpointsegments[i].first;
    int segid2 = closeendpointsegments[i].second;

    // lengths in parameter space of created segments
    double l1 = 2.0 / numseg1_;
    double l2 = 2.0 / numseg2_;
    // element parameter coordinates of left segment boundary
    double eta1_segleft = -1.0 + segid1 * l1;
    double eta2_segleft = -1.0 + segid2 * l2;

    if (segid1 == 0 and boundarynode1_.first)
    {
      // given parameter coordinate
      double eta1 = -1.0;
      // searched parameter coordinate
      double eta2 = 0.0;
      bool pairactive = false;

      double gap_dummy = 0.0;
      double alpha_dummy = 0.0;

      bool solutionwithinsegment = point_to_line_projection(
          eta1, eta2_segleft, l2, eta2, gap_dummy, alpha_dummy, pairactive, false);

      if (solutionwithinsegment and pairactive)
      {
        std::pair<TYPE, TYPE> closestpoint(std::make_pair((TYPE)eta1, (TYPE)eta2));

        std::pair<int, int> integration_ids = std::make_pair(1, 0);
        std::pair<int, int> leftpoint_ids = std::make_pair(segid1, segid2);

        // Create data container for each end point
        // in case of end-point-contact the variable integration_ids contains two bool values (a,b):
        // a \in {0,1} contains the information, if a node of element 1 has been considered as
        // endpoint, b has the same meaning for element 2
        epvariables_.push_back(
            std::make_shared<BeamToBeamContactVariables<numnodes, numnodalvalues>>(
                closestpoint, leftpoint_ids, integration_ids, pp, 1.0));
      }
    }

#ifndef ONLYLEFTENDPOINTCONTACT
    if (segid1 == numseg1_ - 1 and boundarynode1_.second)
    {
      // given parameter coordinate
      double eta1 = 1.0;
      // searched parameter coordinate
      double eta2 = 0.0;
      bool pairactive = false;

      double gap_dummy = 0.0;
      double alpha_dummy = 0.0;

      bool solutionwithinsegment = point_to_line_projection(
          eta1, eta2_segleft, l2, eta2, gap_dummy, alpha_dummy, pairactive, false);

      if (solutionwithinsegment and pairactive)
      {
        std::pair<TYPE, TYPE> closestpoint(std::make_pair((TYPE)eta1, (TYPE)eta2));

        std::pair<int, int> integration_ids = std::make_pair(1, 0);
        std::pair<int, int> leftpoint_ids = std::make_pair(segid1, segid2);

        // Create data container for each end point
        // in case of end-point-contact the variable integration_ids contains two bool values (a,b):
        // a \in {0,1} contains the information, if a node of element 1 has been considered as
        // endpoint, b has the same meaning for element 2
        epvariables_.push_back(
            std::make_shared<BeamToBeamContactVariables<numnodes, numnodalvalues>>(
                closestpoint, leftpoint_ids, integration_ids, pp, 1.0));
      }
    }
#endif

    if (segid2 == 0 and boundarynode2_.first)
    {
      // given parameter coordinate
      double eta2 = -1.0;
      // searched parameter coordinate
      double eta1 = 0.0;
      bool pairactive = false;

      double gap_dummy = 0.0;
      double alpha_dummy = 0.0;

      bool solutionwithinsegment = point_to_line_projection(
          eta2, eta1_segleft, l1, eta1, gap_dummy, alpha_dummy, pairactive, false, true);

      if (solutionwithinsegment and pairactive)
      {
        std::pair<TYPE, TYPE> closestpoint(std::make_pair((TYPE)eta1, (TYPE)eta2));

        std::pair<int, int> integration_ids = std::make_pair(0, 1);
        std::pair<int, int> leftpoint_ids = std::make_pair(segid1, segid2);

        // Create data container for each end point
        // in case of end-point-contact the variable integration_ids contains two bool values (a,b):
        // a \in {0,1} contains the information, if a node of element 1 has been considered as
        // endpoint, b has the same meaning for element 2
        epvariables_.push_back(
            std::make_shared<BeamToBeamContactVariables<numnodes, numnodalvalues>>(
                closestpoint, leftpoint_ids, integration_ids, pp, 1.0));
      }
    }

#ifndef ONLYLEFTENDPOINTCONTACT
    if (segid2 == numseg2_ - 1 and boundarynode2_.second)
    {
      // given parameter coordinate
      double eta2 = 1.0;
      // searched parameter coordinate
      double eta1 = 0.0;
      bool pairactive = false;

      double gap_dummy = 0.0;
      double alpha_dummy = 0.0;

      bool solutionwithinsegment = point_to_line_projection(
          eta2, eta1_segleft, l1, eta1, gap_dummy, alpha_dummy, pairactive, false, true);

      if (solutionwithinsegment and pairactive)
      {
        std::pair<TYPE, TYPE> closestpoint(std::make_pair((TYPE)eta1, (TYPE)eta2));

        std::pair<int, int> integration_ids = std::make_pair(0, 1);
        std::pair<int, int> leftpoint_ids = std::make_pair(segid1, segid2);

        // Create data container for each end point
        // in case of end-point-contact the variable integration_ids contains two bool values (a,b):
        // a \in {0,1} contains the information, if a node of element 1 has been considered as
        // endpoint, b has the same meaning for element 2
        epvariables_.push_back(
            std::make_shared<BeamToBeamContactVariables<numnodes, numnodalvalues>>(
                closestpoint, leftpoint_ids, integration_ids, pp, 1.0));
      }
    }
#endif

    if ((segid1 == 0 and boundarynode1_.first) and (segid2 == 0 and boundarynode2_.first))
    {
      double eta1 = -1.0;
      double eta2 = -1.0;
      Core::LinAlg::Matrix<3, 1> deltanodalpos(true);
      for (int i = 0; i < 3; i++)
      {
        deltanodalpos(i) = Core::FADUtils::cast_to_double(ele2pos_(i) - ele1pos_(i));
      }

      double gap = deltanodalpos.norm2() - r1_ - r2_;
      if (check_contact_status(gap))
      {
        std::pair<TYPE, TYPE> closestpoint(std::make_pair((TYPE)eta1, (TYPE)eta2));

        std::pair<int, int> integration_ids = std::make_pair(1, 1);
        std::pair<int, int> leftpoint_ids = std::make_pair(segid1, segid2);

        // Create data container for each end point
        // in case of end-point-contact the variable integration_ids contains two bool values (a,b):
        // a \in {0,1} contains the information, if a node of element 1 has been considered as
        // endpoint, b has the same meaning for element 2
        epvariables_.push_back(
            std::make_shared<BeamToBeamContactVariables<numnodes, numnodalvalues>>(
                closestpoint, leftpoint_ids, integration_ids, pp, 1.0));
      }
    }

#ifndef ONLYLEFTENDPOINTCONTACT
    if ((segid1 == 0 and boundarynode1_.first) and
        (segid2 == numseg2_ - 1 and boundarynode2_.second))
    {
      double eta1 = -1.0;
      double eta2 = 1.0;
      Core::LinAlg::Matrix<3, 1> deltanodalpos(true);
      for (int i = 0; i < 3; i++)
      {
        deltanodalpos(i) = Core::FADUtils::cast_to_double(ele2pos_(6 + i) - ele1pos_(i));
      }

      double gap = deltanodalpos.norm2() - r1_ - r2_;
      if (check_contact_status(gap))
      {
        std::pair<TYPE, TYPE> closestpoint(std::make_pair((TYPE)eta1, (TYPE)eta2));

        std::pair<int, int> integration_ids = std::make_pair(1, 1);
        std::pair<int, int> leftpoint_ids = std::make_pair(segid1, segid2);

        // Create data container for each end point
        // in case of end-point-contact the variable integration_ids contains two bool values (a,b):
        // a \in {0,1} contains the information, if a node of element 1 has been considered as
        // endpoint, b has the same meaning for element 2
        epvariables_.push_back(
            std::make_shared<BeamToBeamContactVariables<numnodes, numnodalvalues>>(
                closestpoint, leftpoint_ids, integration_ids, pp, 1.0));
      }
    }

    if ((segid1 == numseg1_ - 1 and boundarynode1_.second) and
        (segid2 == 0 and boundarynode2_.first))
    {
      double eta1 = 1.0;
      double eta2 = -1.0;
      Core::LinAlg::Matrix<3, 1> deltanodalpos(true);
      for (int i = 0; i < 3; i++)
      {
        deltanodalpos(i) = Core::FADUtils::cast_to_double(ele2pos_(i) - ele1pos_(6 + i));
      }

      double gap = deltanodalpos.norm2() - r1_ - r2_;
      if (check_contact_status(gap))
      {
        std::pair<TYPE, TYPE> closestpoint(std::make_pair((TYPE)eta1, (TYPE)eta2));

        std::pair<int, int> integration_ids = std::make_pair(1, 1);
        std::pair<int, int> leftpoint_ids = std::make_pair(segid1, segid2);

        // Create data container for each end point
        // in case of end-point-contact the variable integration_ids contains two bool values (a,b):
        // a \in {0,1} contains the information, if a node of element 1 has been considered as
        // endpoint, b has the same meaning for element 2
        epvariables_.push_back(
            std::make_shared<BeamToBeamContactVariables<numnodes, numnodalvalues>>(
                closestpoint, leftpoint_ids, integration_ids, pp, 1.0));
      }
    }

    if ((segid1 == numseg1_ - 1 and boundarynode1_.second) and
        (segid2 == numseg2_ - 1 and boundarynode2_.second))
    {
      double eta1 = 1.0;
      double eta2 = 1.0;
      Core::LinAlg::Matrix<3, 1> deltanodalpos(true);
      for (int i = 0; i < 3; i++)
      {
        deltanodalpos(i) = Core::FADUtils::cast_to_double(ele2pos_(6 + i) - ele1pos_(6 + i));
      }

      double gap = deltanodalpos.norm2() - r1_ - r2_;
      if (check_contact_status(gap))
      {
        std::pair<TYPE, TYPE> closestpoint(std::make_pair((TYPE)eta1, (TYPE)eta2));

        std::pair<int, int> integration_ids = std::make_pair(1, 1);
        std::pair<int, int> leftpoint_ids = std::make_pair(segid1, segid2);

        // Create data container for each end point
        // in case of end-point-contact the variable integration_ids contains two bool values (a,b):
        // a \in {0,1} contains the information, if a node of element 1 has been considered as
        // endpoint, b has the same meaning for element 2
        epvariables_.push_back(
            std::make_shared<BeamToBeamContactVariables<numnodes, numnodalvalues>>(
                closestpoint, leftpoint_ids, integration_ids, pp, 1.0));
      }
    }
#endif
  }  // for(int i=0;i<closeendpointsegments.size();i++)
}
/*----------------------------------------------------------------------*
 |  end: Get active endpoint pairs
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Evaluate active endpoint pairs                           meier 12/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes,
    numnodalvalues>::evaluate_active_end_point_pairs(Core::LinAlg::SerialDenseVector* forcevec1,
    Core::LinAlg::SerialDenseVector* forcevec2, Core::LinAlg::SerialDenseMatrix* stiffmat11,
    Core::LinAlg::SerialDenseMatrix* stiffmat12, Core::LinAlg::SerialDenseMatrix* stiffmat21,
    Core::LinAlg::SerialDenseMatrix* stiffmat22)
{
  for (int numep = 0; numep < (int)epvariables_.size(); numep++)
  {
    //**********************************************************************
    // (2) Compute some auxiliary quantities
    //**********************************************************************

    // vectors for shape functions and their derivatives
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N1(true);       // = N1
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N2(true);       // = N2
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N1_xi(true);    // = N1,xi
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N2_xi(true);    // = N2,eta
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N1_xixi(true);  // = N1,xixi
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N2_xixi(true);  // = N2,etaeta

    // coords and derivatives of the two contacting points
    Core::LinAlg::Matrix<3, 1, TYPE> r1(true);       // = r1
    Core::LinAlg::Matrix<3, 1, TYPE> r2(true);       // = r2
    Core::LinAlg::Matrix<3, 1, TYPE> r1_xi(true);    // = r1,xi
    Core::LinAlg::Matrix<3, 1, TYPE> r2_xi(true);    // = r2,eta
    Core::LinAlg::Matrix<3, 1, TYPE> r1_xixi(true);  // = r1,xixi
    Core::LinAlg::Matrix<3, 1, TYPE> r2_xixi(true);  // = r2,etaeta
    Core::LinAlg::Matrix<3, 1, TYPE> delta_r(true);  // = r1-r2

    TYPE eta1 = epvariables_[numep]->get_cp().first;
    TYPE eta2 = epvariables_[numep]->get_cp().second;

#ifdef AUTOMATICDIFF
    BeamContact::SetFADParCoordDofs<numnodes, numnodalvalues>(eta1, eta2);
#endif

    // update shape functions and their derivatives
    get_shape_functions(N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi, eta1, eta2);
    // update coordinates and derivatives of contact points
    compute_coords_and_derivs(
        r1, r2, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi);

    // call function to compute scaled normal and gap of contact point
    compute_normal(r1, r2, r1_xi, r2_xi, epvariables_[numep], 2);

    // call function to compute penalty force
    calc_penalty_law(*epvariables_[numep]);

    epvariables_[numep]->set_p_pfac(1.0);
    epvariables_[numep]->set_dp_pfac(0.0);

    // In case of endpoint-contact, the length specific energy and the 'real' energy are identical
    double lengthspec_energy = Core::FADUtils::cast_to_double(epvariables_[numep]->get_energy());
    epvariables_[numep]->set_integrated_energy(lengthspec_energy);

    //    std::cout << "epvariables_[numep]->GetNormal(): " << epvariables_[numep]->GetNormal() <<
    //    std::endl; std::cout << "numep: " << numep << std::endl; std::cout << "xi: " <<
    //    epvariables_[numep]->GetCP().first.val() << std::endl; std::cout << "eta: " <<
    //    epvariables_[numep]->GetCP().second.val() << std::endl; std::cout << "gap: " <<
    //    epvariables_[numep]->GetGap().val() << std::endl; std::cout << "angle: " <<
    //    epvariables_[numep]->get_angle()/M_PI*180.0 << std::endl; std::cout << "r1_xi: " << r1_xi
    //    << std::endl; std::cout << "r2_xi: " << r2_xi << std::endl; std::cout << "|r1_xi|: " <<
    //    r1_xi.Norm2() << std::endl; std::cout << "|r2_xi|: " << r2_xi.Norm2() << std::endl;
    //    std::cout << "r1_xi*r2_xi: " << Core::FADUtils::ScalarProduct(r1_xi,r2_xi) << std::endl;
    //    std::cout << "epvariables_[numep]->Getfp(): " << epvariables_[numep]->Getfp() <<
    //    std::endl;

    bool fixedendpointxi = epvariables_[numep]->get_int_ids().first;
    bool fixedendpointeta = epvariables_[numep]->get_int_ids().second;

    // call function to compute contact contribution to residual vector
    if (forcevec1 != nullptr and forcevec2 != nullptr)
      evaluate_fc_contact(*forcevec1, *forcevec2, r1, r2, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2,
          N1_xi, N2_xi, *epvariables_[numep], 1.0, false, false, fixedendpointxi, fixedendpointeta);

    // call function to compute contact contribution to stiffness matrix
    if (stiffmat11 != nullptr and stiffmat12 != nullptr and stiffmat21 != nullptr and
        stiffmat22 != nullptr)
      evaluate_stiffc_contact(*stiffmat11, *stiffmat12, *stiffmat21, *stiffmat22, r1, r2, r1_xi,
          r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi, *epvariables_[numep],
          1.0, false, false, fixedendpointxi, fixedendpointeta);
  }
}
/*----------------------------------------------------------------------*
 |  end: Evaluate active endpoint pairs
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Calculate scalar contact force                           meier 10/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::calc_penalty_law(
    BeamToBeamContactVariables<numnodes, numnodalvalues>& variables)
{
  // First parameter for contact force regularization
  double g0 = params()->beam_to_beam_contact_params()->beam_to_beam_penalty_law_regularization_g0();
  TYPE fp = 0.0;
  TYPE dfp = 0.0;
  TYPE e = 0.0;
  double pp = variables.get_pp();
  TYPE gap = variables.get_gap();

  if (!check_contact_status(Core::FADUtils::cast_to_double(gap))) return;

  switch (params()->beam_to_beam_contact_params()->penalty_law())
  {
    case BeamContact::pl_lp:  // linear penalty force law
    {
      fp = -pp * gap;
      dfp = -pp;
      e = -1.0 / 2.0 * pp * gap * gap;

      break;
    }
    case BeamContact::pl_qp:  // quadratic penalty force law
    {
      fp = pp * gap * gap;
      dfp = 2 * pp * gap;
      e = 1.0 / 3.0 * pp * gap * gap * gap;

      break;
    }
    case BeamContact::pl_lnqp:  // quadratic regularization for negative gaps
    {
      if (g0 == -1.0)
        FOUR_C_THROW("Invalid value of regularization parameter BEAMS_PENREGPARAM_G0!");

      if (gap > -g0)
      {
        // std::cout << "Regularized Penalty!" << std::endl;
        fp = pp / (2.0 * g0) * gap * gap;
        dfp = pp / (g0)*gap;
      }
      else
      {
        // std::cout << "Linear Penalty!" << std::endl;
        fp = -pp * (gap + g0 / 2.0);
        dfp = -pp;
      }

      break;
    }
    case BeamContact::pl_lpqp:  // quadratic regularization for positive gaps
    {
      if (g0 == -1.0)
        FOUR_C_THROW("Invalid value of regularization parameter BEAMS_PENREGPARAM_G0!");

      // Parameter to shift penalty law
      double gbar = params()->beam_to_beam_contact_params()->gap_shift();
      gap = gap + gbar;

      double f0 = g0 * pp / 2.0;
      double factor_a = pp / (g0)-f0 / (g0 * g0);  //=pp/(2*g0)
      double factor_b = -pp;
      double factor_c = f0;
      if (gap > 0)
      {
        // std::cout << "Regularized Penalty!" << std::endl;
        fp = factor_a * gap * gap + factor_b * gap + factor_c;
        dfp = 2 * factor_a * gap + factor_b;
        e = -pp * g0 * g0 / 6.0 +
            (pp / (6.0 * g0) * gap * gap * gap - pp / 2.0 * gap * gap + pp * g0 / 2.0 * gap);
      }
      else
      {
        // std::cout << "Linear Penalty!" << std::endl;
        fp = f0 - pp * gap;
        dfp = -pp;
        e = -(pp * g0 * g0 / 6.0 + pp / 2.0 * gap * gap - pp * g0 / 2.0 * gap);
      }

      break;
    }
    case BeamContact::pl_lpcp:  // cubic regularization for positive gaps
    {
      if (g0 == -1.0)
        FOUR_C_THROW("Invalid value of regularization parameter BEAMS_PENREGPARAM_G0!");

      // Third parameter for contact force regularization
      double c0 =
          params()->beam_to_beam_contact_params()->beam_to_beam_penalty_law_regularization_c0();
      if (c0 == -1.0)
        FOUR_C_THROW("Invalid value of regularization parameter BEAMS_PENREGPARAM_C0!");

      // k \in ~[1;3] delivers sensible results representing a parable without turning point
      // k \in ~[3;6] delivers a parable with turning point and consequently also small negative
      // contact forces ~0.1*f0 k=2.0 is  identical to the quadratic regularization for positive
      // gaps!
      double k = c0;
      double f0 = pp * g0 / k;
      double factor_a = -pp / (g0 * g0) + 2 * f0 / (g0 * g0 * g0);
      double factor_b = 2 * pp / (g0)-3 * f0 / (g0 * g0);
      double factor_c = -pp;
      double factor_d = f0;
      if (gap > 0.0)
      {
        // std::cout << "Regularized Penalty!" << std::endl;
        fp = factor_a * gap * gap * gap + factor_b * gap * gap + factor_c * gap + factor_d;
        dfp = 3 * factor_a * gap * gap + 2 * factor_b * gap + factor_c;
      }
      else
      {
        // std::cout << "Linear Penalty!" << std::endl;
        fp = f0 - pp * gap;
        dfp = -pp;
      }

      break;
    }
    case BeamContact::pl_lpdqp:  // double quadratic regularization for positive gaps
    {
      if (g0 == -1.0)
        FOUR_C_THROW("Invalid value of regularization parameter BEAMS_PENREGPARAM_G0!");

      // Third parameter for contact force regularization
      double c0 =
          params()->beam_to_beam_contact_params()->beam_to_beam_penalty_law_regularization_c0();
      if (c0 == -1.0)
        FOUR_C_THROW("Invalid value of regularization parameter BEAMS_PENREGPARAM_C0!");

      // Second parameter for contact force regularization
      double f0 =
          params()->beam_to_beam_contact_params()->beam_to_beam_penalty_law_regularization_f0();
      if (f0 == -1.0)
        FOUR_C_THROW("Invalid value of regularization parameter BEAMS_PENREGPARAM_F0!");

      double k =
          c0;  // transition between first and second quadratic regularization part: k \in [0;2.0]
      double g1 = k * f0 / pp;
      double c_tilde = f0;
      double b_tilde = -pp;
      double a_bar = (2 * f0 - pp * g1) / (2 * g0 * (g0 - g1));
      double b_bar = -2 * g0 * a_bar;
      double c_bar = -g0 * g0 * a_bar - g0 * b_bar;
      double a_tilde = (2 * g1 * a_bar + b_bar - b_tilde) / (2 * g1);

      if (gap > g1)
      {
        // std::cout << "Regularized Penalty: g1 < gap < g0!" << std::endl;
        fp = a_bar * gap * gap + b_bar * gap + c_bar;
        dfp = 2 * a_bar * gap + b_bar;
      }
      else if (gap > 0)
      {
        // std::cout << "Regularized Penalty: 0 < gap < g1!" << std::endl;
        fp = a_tilde * gap * gap + b_tilde * gap + c_tilde;
        dfp = 2 * a_tilde * gap + b_tilde;
      }
      else
      {
        // std::cout << "Linear Penalty!" << std::endl;
        fp = f0 - pp * gap;
        dfp = -pp;
      }

      break;
    }
    case BeamContact::pl_lpep:  // exponential regularization for positive gaps. Here g0
                                // represents the cut off radius!
    {
      if (g0 == -1.0)
        FOUR_C_THROW("Invalid value of regularization parameter BEAMS_PENREGPARAM_G0!");

      // Second parameter for contact force regularization
      double f0 =
          params()->beam_to_beam_contact_params()->beam_to_beam_penalty_law_regularization_f0();
      if (f0 == -1.0)
        FOUR_C_THROW("Invalid value of regularization parameter BEAMS_PENREGPARAM_F0!");

      if (gap > 0)
      {
        // std::cout << "Regularized Penalty: 0 < gap < g1!" << std::endl;
        fp = f0 * exp(-pp * gap / f0);
        dfp = -pp * exp(-pp * gap / f0);
        if (f0 * exp(-pp * g0 / f0) > 0.01 * f0)
        {
          std::cout << "      Warning - g0: " << g0
                    << " f0*exp(-pp*g0/f0): " << f0 * exp(-pp * g0 / f0)
                    << "-> Choose higher cut-off radius g0!" << std::endl;
        }
      }
      else
      {
        // std::cout << "Linear Penalty!" << std::endl;
        fp = f0 - pp * gap;
        dfp = -pp;
      }

      break;
    }
  }

  variables.setfp(fp);
  variables.setdfp(dfp);
  variables.set_energy(e);

  return;
}
/*----------------------------------------------------------------------*
 |  end: Calculate scalar contact force
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Calculate angle-dependent perp-penalty scale factor       meier 10/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::calc_perp_penalty_scale_fac(
    BeamToBeamContactVariables<numnodes, numnodalvalues>& cpvariables,
    Core::LinAlg::Matrix<3, 1, TYPE>& r1_xi, Core::LinAlg::Matrix<3, 1, TYPE>& r2_xi,
    const double shiftangle1, const double shiftangle2)
{
  // Penalty scale factor that reduces the penalty parameter for small angles
  TYPE ppfac = 1.0;
  TYPE dppfac = 0.0;

  if (shiftangle1 > M_PI / 2.0 and shiftangle2 > M_PI / 2.0)
  {
    ppfac = 0.0;
    dppfac = 0.0;
  }
  else
  {
    TYPE s = fabs(Core::FADUtils::scalar_product(r1_xi, r2_xi) /
                  (Core::FADUtils::vector_norm<3>(r1_xi) * Core::FADUtils::vector_norm<3>(r2_xi)));
    double s1 = cos(shiftangle1);
    double s2 = cos(shiftangle2);

    if (shiftangle1 < 0.0 or shiftangle1 > M_PI / 2.0 or shiftangle2 < 0.0 or
        shiftangle2 > M_PI / 2.0 or shiftangle1 >= shiftangle2)
      FOUR_C_THROW("Invalid choice of shift angles!");

    if (Core::FADUtils::cast_to_double(s) > s1)
      ppfac = 0.0;
    else if (Core::FADUtils::cast_to_double(s) > s2)
    {
#ifndef CONSISTENTTRANSITION
      ppfac = 0.5 * (cos(M_PI * (s - s2) / (s1 - s2)) + 1.0);
      dppfac = -0.5 * M_PI / (s1 - s2) * sin(M_PI * (s - s2) / (s1 - s2));
#else
      if (CONSISTENTTRANSITION == 1)
      {
        TYPE simple_fac = 0.5 * (cos(M_PI * (s - s2) / (s1 - s2)) + 1.0);
        TYPE d_simple_fac = -0.5 * M_PI / (s1 - s2) * sin(M_PI * (s - s2) / (s1 - s2));
        ppfac = simple_fac * simple_fac;
        dppfac = 2 * simple_fac * d_simple_fac;
      }
      else if (CONSISTENTTRANSITION == 2)
      {
        TYPE simple_fac = 0.5 * (-cos(M_PI * (s - s2) / (s1 - s2)) + 1.0);
        TYPE d_simple_fac = 0.5 * M_PI / (s1 - s2) * sin(M_PI * (s - s2) / (s1 - s2));
        ppfac = 1 - simple_fac * simple_fac;
        dppfac = -2 * simple_fac * d_simple_fac;
      }
      else
        FOUR_C_THROW(
            "Inadmissible value of CONSISTENTTRANSITION, only the values 1 and 2 are allowed!");
#endif
    }
  }

  // set class variable
  cpvariables.set_p_pfac(ppfac);
  cpvariables.set_dp_pfac(dppfac);

  return;
}
/*----------------------------------------------------------------------*
 |  end: Calculate angle-dependent perp-penalty scale factor
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Calculate angle-dependent par-penalty scale factor       Meier 10/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::calc_par_penalty_scale_fac(
    BeamToBeamContactVariables<numnodes, numnodalvalues>& gpvariables,
    Core::LinAlg::Matrix<3, 1, TYPE>& r1_xi, Core::LinAlg::Matrix<3, 1, TYPE>& r2_xi,
    const double shiftangle1, const double shiftangle2)
{
  // Penalty scale factor that reduces the penalty parameter for small angles
  TYPE ppfac = 1.0;
  TYPE dppfac = 0.0;

  if (shiftangle1 > M_PI / 2.0 and shiftangle2 > M_PI / 2.0)
  {
    ppfac = 1.0;
    dppfac = 0.0;
  }
  else
  {
    TYPE s = fabs(Core::FADUtils::scalar_product(r1_xi, r2_xi) /
                  (Core::FADUtils::vector_norm<3>(r1_xi) * Core::FADUtils::vector_norm<3>(r2_xi)));
    double s1 = cos(shiftangle1);
    double s2 = cos(shiftangle2);

    if (shiftangle1 < 0.0 or shiftangle1 > M_PI / 2.0 or shiftangle2 < 0.0 or
        shiftangle2 > M_PI / 2.0 or shiftangle1 >= shiftangle2)
      FOUR_C_THROW("Invalid choice of shift angles!");

    if (Core::FADUtils::cast_to_double(s) > s1)
      ppfac = 1.0;
    else if (Core::FADUtils::cast_to_double(s) > s2)
    {
#ifndef CONSISTENTTRANSITION
      ppfac = 0.5 * (-cos(M_PI * (s - s2) / (s1 - s2)) + 1.0);
      dppfac = 0.5 * M_PI / (s1 - s2) * sin(M_PI * (s - s2) / (s1 - s2));
#else
      TYPE simple_fac = 0.5 * (-cos(M_PI * (s - s2) / (s1 - s2)) + 1.0);
      TYPE d_simple_fac = 0.5 * M_PI / (s1 - s2) * sin(M_PI * (s - s2) / (s1 - s2));
      ppfac = simple_fac * simple_fac;
      dppfac = 2 * simple_fac * d_simple_fac;
#endif
    }
    else
      ppfac = 0.0;
  }

  // set class variable
  gpvariables.set_p_pfac(ppfac);
  gpvariables.set_dp_pfac(dppfac);

  return;
}
/*----------------------------------------------------------------------*
 |  end: Calculate angle-dependent par-penalty scale factor
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Subdivide elements into segments for CPP                 meier 10/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
double BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::create_segments(
    const Core::Elements::Element* ele,
    std::vector<Core::LinAlg::Matrix<3, 1, double>>& endpoints_final, int& numsegment, int i)
{
  // endpoints of the segments
  std::vector<Core::LinAlg::Matrix<3, 1, double>> endpoints(
      (int)MAXNUMSEG + 1, Core::LinAlg::Matrix<3, 1, double>(true));
  double segangle = params()->beam_to_beam_contact_params()->segmentation_angle();

  numsegment = 1;
  double deltaxi = 2.0;

  if (i == 0)
  {
    numsegment = INITSEG1;
    deltaxi = 2.0 / INITSEG1;
  }

  if (i == 1)
  {
    numsegment = INITSEG2;
    deltaxi = 2.0 / INITSEG2;
  }

  double xi1 = 0.0;
  double xi2 = 0.0;
  Core::LinAlg::Matrix<3, 1, double> r1(true);
  Core::LinAlg::Matrix<3, 1, double> t1(true);
  Core::LinAlg::Matrix<3, 1, double> r2(true);
  Core::LinAlg::Matrix<3, 1, double> t2(true);
  Core::LinAlg::Matrix<3, 1, double> rm(true);
  double l = 0.0;
  double segdist = 0.0;
  double maxsegdist = 0.0;
  bool moresegments = true;

  while (moresegments)
  {
    // We have to zero maxsegdist for each new segment distribution, otherwise we would get a larger
    // value of a former rougher distribution!
    maxsegdist = 0.0;
    moresegments = false;
    for (int i = 0; i < numsegment; i++)
    {
      if (numsegment > (int)MAXNUMSEG)
        FOUR_C_THROW(
            "Not more segments than MAXNUMSEG per element possible! Increase MAXNUMSEG or apply "
            "finer discretization!");

      xi1 = 0.0;
      xi2 = 0.0;
      xi1 = -1.0 + i / ((double)numsegment) * 2.0;
      xi2 = -1.0 + (i + 1) / ((double)numsegment) *
                       2.0;  // The cast to double is necessary here to avoid integer round-off
      Core::LinAlg::Matrix<3, 1, TYPE> auxmatrix(true);

      auxmatrix = r(xi1, ele);
      r1 = Core::FADUtils::cast_to_double<TYPE, 3, 1>(auxmatrix);
      auxmatrix = r(xi2, ele);
      r2 = Core::FADUtils::cast_to_double<TYPE, 3, 1>(auxmatrix);
      auxmatrix = r_xi(xi1, ele);
      t1 = Core::FADUtils::cast_to_double<TYPE, 3, 1>(auxmatrix);
      auxmatrix = r_xi(xi2, ele);
      t2 = Core::FADUtils::cast_to_double<TYPE, 3, 1>(auxmatrix);
      auxmatrix = r((xi1 + xi2) / 2.0, ele);
      rm = Core::FADUtils::cast_to_double<TYPE, 3, 1>(auxmatrix);

      endpoints[i] = r1;
      endpoints[i + 1] = r2;
      l = Core::FADUtils::vector_norm<3>(Core::FADUtils::diff_vector(r1, r2));
      // TODO: adapt this tolerance if necessary!!!
      segdist = 1.0 * l / 2.0 * std::tan(segangle);

      if (segdist > maxsegdist) maxsegdist = segdist;
      if (!check_segment(r1, t1, r2, t2, rm, segdist)) moresegments = true;
    }

    deltaxi = deltaxi / 2;
    numsegment = numsegment * 2;
  }
  numsegment = numsegment / 2;

#ifdef NOSEGMENTATION
  if (numsegment > 1)
    FOUR_C_THROW(
        "Choose higher SEGANGLE since in the case NOSEGMENTATION only one segment per element is "
        "allowed!");
#endif

  // std::cout << "numsegment: " << numsegment << std::endl;
  endpoints_final.resize(numsegment + 1);

  for (int i = 0; i < numsegment + 1; i++)
  {
    endpoints_final[i] = endpoints[i];
  }

  return maxsegdist;
}
/*----------------------------------------------------------------------*
 |  end: Subdivide elements into segments for CPP
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Max. distance at which a contact force becomes active     meier 10/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
double BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::get_max_active_dist()
{
  double maxactivedist = 0.0;
  int penaltylaw = params()->beam_to_beam_contact_params()->penalty_law();

  switch (penaltylaw)
  {
    case BeamContact::pl_lp:
    case BeamContact::pl_qp:
    case BeamContact::pl_lnqp:
    {
      maxactivedist = 0.0;
      break;
    }
    case BeamContact::pl_lpqp:
    {
      double g0 =
          params()->beam_to_beam_contact_params()->beam_to_beam_penalty_law_regularization_g0();
      if (g0 == -1.0)
        FOUR_C_THROW("Invalid value of regularization parameter BEAMS_PENREGPARAM_G0!");

      // Parameter to shift penalty law
      double gbar = params()->beam_to_beam_contact_params()->gap_shift();

      maxactivedist = g0 - gbar;

      break;
    }
    case BeamContact::pl_lpcp:
    case BeamContact::pl_lpdqp:
    case BeamContact::pl_lpep:
    {
      maxactivedist =
          params()->beam_to_beam_contact_params()->beam_to_beam_penalty_law_regularization_g0();
      if (maxactivedist == -1.0)
        FOUR_C_THROW("Invalid value of regularization parameter BEAMS_PENREGPARAM_G0!");
      break;
    }
  }

  return maxactivedist;
}

/*----------------------------------------------------------------------*
 |  Check, if segments are fine enough                       meier 10/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
bool BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::check_segment(
    Core::LinAlg::Matrix<3, 1, double>& r1, Core::LinAlg::Matrix<3, 1, double>& t1,
    Core::LinAlg::Matrix<3, 1, double>& r2, Core::LinAlg::Matrix<3, 1, double>& t2,
    Core::LinAlg::Matrix<3, 1, double>& rm, double& segdist)
{
  Core::LinAlg::Matrix<3, 1, double> t_lin(true);
  Core::LinAlg::Matrix<3, 1, double> rm_lin(true);
  double angle1(0.0);
  double angle2(0.0);
  double dist(0.0);
  double segangle = params()->beam_to_beam_contact_params()->segmentation_angle();

  // Calculate tangent and midpint of linear nodal interpolation
  for (int i = 0; i < 3; i++)
  {
    t_lin(i) = r2(i) - r1(i);
    rm_lin(i) = (r2(i) + r1(i)) / 2.0;
  }

  Core::LinAlg::Matrix<3, 1, double> diffvec(true);
  diffvec = Core::FADUtils::diff_vector(rm_lin, rm);
  dist = (double)Core::FADUtils::vector_norm<3>(diffvec);
  angle1 = (double)BeamInteraction::calc_angle(t1, t_lin);
  angle2 = (double)BeamInteraction::calc_angle(t2, t_lin);

  if (std::fabs(angle1) < segangle and
      std::fabs(angle2) < segangle)  // segment distribution is fine enough
  {
    if (std::fabs(dist) > segdist)
    {
      this->print(std::cout);

      std::cout << "\n\nangle_left= " << angle1 / M_PI * 180
                << "degrees, angle_right= " << angle2 / M_PI * 180
                << "degrees, segangle= " << segangle / M_PI * 180 << "degrees";
      std::cout << "\n\nsegment midpoint position (linear approx)= " << rm_lin;
      std::cout << "\n\nsegment midpoint position                = " << rm;
      std::cout << "\n\ndist= " << dist << ", segdist= " << segdist;

      FOUR_C_THROW(
          "Segmentation in beam contact search (step 2/2, i.e. 'fine' search) failed!\n"
          "Maximal spatial extension of this segment cannot be approximated by a double cone!\n"
          "Check for extremely deformed beam elements and try to avoid this.\n"
          "Increase value of BEAMS_SEGANGLE ?!");
    }

    return true;
  }
  else  // we still need more segments
    return false;
}

/*----------------------------------------------------------------------*
 |  Find segments close to each other                        meier 10/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::get_close_segments(
    const std::vector<Core::LinAlg::Matrix<3, 1, double>>& endpoints1,
    const std::vector<Core::LinAlg::Matrix<3, 1, double>>& endpoints2,
    std::map<std::pair<int, int>, Core::LinAlg::Matrix<3, 1, double>>& closesmallanglesegments,
    std::map<std::pair<int, int>, Core::LinAlg::Matrix<3, 1, double>>& closelargeanglesegments,
    std::vector<std::pair<int, int>>& closeendpointsegments, double maxactivedist)
{
  Core::LinAlg::Matrix<3, 1, double> t1(true);
  Core::LinAlg::Matrix<3, 1, double> t2(true);
  Core::LinAlg::Matrix<3, 1, double> r1_a(true);
  Core::LinAlg::Matrix<3, 1, double> r1_b(true);
  Core::LinAlg::Matrix<3, 1, double> r2_a(true);
  Core::LinAlg::Matrix<3, 1, double> r2_b(true);
  double angle(0.0);

  bool endpoint_penalty = params()->beam_to_beam_contact_params()->end_point_penalty();

  // Safety factor for determination of close segments
  double safetyfac = 1.1;
  // Distance at which intersection happens
  double distancelimit = safetyfac * (maxsegdist1_ + maxsegdist2_ + maxactivedist + r1_ + r2_);

  int numseg1 = (int)endpoints1.size() - 1;
  int numseg2 = (int)endpoints2.size() - 1;

  // TODO: This check is implemented in a brute force way. However, this should be efficient enough
  // as long as the number of segments per element remains small!
  for (int i = 0; i < numseg1; i++)
  {
    r1_a = endpoints1[i];
    r1_b = endpoints1[i + 1];
    t1 = Core::FADUtils::diff_vector(r1_b, r1_a);
    for (int j = 0; j < numseg2; j++)
    {
      r2_a = endpoints2[j];
      r2_b = endpoints2[j + 1];
      t2 = Core::FADUtils::diff_vector(r2_b, r2_a);

      angle = BeamInteraction::calc_angle(t1, t2);

      //*******1) intersection between two parallel
      // cylinders*********************************************************
      if (std::fabs(angle) < ANGLETOL)
      {
        if (BeamInteraction::intersect_parallel_cylinders(r1_a, r1_b, r2_a, r2_b, distancelimit))
        {
          Core::LinAlg::Matrix<3, 1, double> segmentdata(true);
          segmentdata(0) = angle;   // segment angle
          segmentdata(1) = 1000.0;  // eta1_seg
          segmentdata(2) = 1000.0;  // eta2_seg

          // Add new small angle pair
          if (std::fabs(angle) <= deltasmallangle_)
            closesmallanglesegments[std::make_pair(i, j)] = segmentdata;
          if (std::fabs(angle) >= deltalargeangle_)
            closelargeanglesegments[std::make_pair(i, j)] = segmentdata;

          // If the element lies on the boundary of a physical beam, we sort out the corresponding
          // boundary segments
          if (endpoint_penalty)
          {
            if ((i == 0 and boundarynode1_.first) or (i == numseg1 - 1 and boundarynode1_.second) or
                (j == 0 and boundarynode2_.first) or (j == numseg2 - 1 and boundarynode2_.second))
              closeendpointsegments.push_back(std::make_pair(i, j));
          }
        }
      }
      //*******2) intersection between two arbitrary oriented
      // cylinders************************************************
      else
      {
        std::pair<double, double> closestpoints(std::make_pair(0.0, 0.0));
        bool etaset = false;
        if (BeamInteraction::intersect_arbitrary_cylinders(
                r1_a, r1_b, r2_a, r2_b, distancelimit, closestpoints, etaset))
        {
          Core::LinAlg::Matrix<3, 1, double> segmentdata(true);
          segmentdata(0) = angle;  // segment angle

          if (etaset)
          {
            segmentdata(1) = closestpoints.first;   // eta1_seg
            segmentdata(2) = closestpoints.second;  // eta2_seg
          }
          else
          {
            segmentdata(1) = 1000.0;  // eta1_seg
            segmentdata(2) = 1000.0;  // eta2_seg
          }
          // Add new small angle pair
          if (std::fabs(angle) <= deltasmallangle_)
            closesmallanglesegments[std::make_pair(i, j)] = segmentdata;
          if (std::fabs(angle) >= deltalargeangle_)
            closelargeanglesegments[std::make_pair(i, j)] = segmentdata;

          // If the element lies on the boundary of a physical beam, we sort out the corresponding
          // boundary segments
          if (endpoint_penalty)
          {
            if ((i == 0 and boundarynode1_.first) or (i == numseg1 - 1 and boundarynode1_.second) or
                (j == 0 and boundarynode2_.first) or (j == numseg2 - 1 and boundarynode2_.second))
              closeendpointsegments.push_back(std::make_pair(i, j));
          }
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Closest point projection                                  meier 01/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
bool BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::closest_point_projection(
    double& eta_left1, double& eta_left2, double& l1, double& l2,
    Core::LinAlg::Matrix<3, 1, double>& segmentdata, std::pair<TYPE, TYPE>& solutionpoints,
    int segid1, int segid2)
{
  std::vector<std::pair<double, double>> startingpoints(0);
  bool validpairfound = false;
  double gap = 0.0;
  double eta_right1 = eta_left1 + l1;
  double eta_right2 = eta_left2 + l2;

  double etalocal1 = segmentdata(1);
  double etalocal2 = segmentdata(2);

  if (fabs(etalocal1) <= 1.0 and fabs(etalocal2) <= 1.0)
  {
    startingpoints.push_back(std::make_pair(eta_left1 + 0.5 * l1 * (1 + etalocal1),
        eta_left2 +
            0.5 * l2 * (1 + etalocal2)));  // cp of linear segment approximation as starting point
  }

  startingpoints.push_back(std::make_pair(
      eta_left1 + 0.5 * l1, eta_left2 + 0.5 * l2));  // segment midpoint as starting point

  // Other combinations of (etalocal1, etalocal2 \in {-1;0;1}) for each segment -> 8 additional
  // combinations besides (0,0)
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      if (!(i == 0 and j == 0))  // we already have the segment midpoint combination (0,0)
        startingpoints.push_back(
            std::make_pair(eta_left1 + i * 0.5 * l1, eta_left2 + j * 0.5 * l2));
    }
  }

  for (int numstartpoint = 0; numstartpoint < (int)startingpoints.size(); numstartpoint++)
  {
    // vectors for shape functions and their derivatives
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N1(true);       // = N1
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N2(true);       // = N2
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N1_xi(true);    // = N1,xi
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N2_xi(true);    // = N2,eta
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N1_xixi(true);  // = N1,xixi
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N2_xixi(true);  // = N2,etaeta

    // coords and derivatives of the two contacting points
    Core::LinAlg::Matrix<3, 1, TYPE> r1(true);       // = r1
    Core::LinAlg::Matrix<3, 1, TYPE> r2(true);       // = r2
    Core::LinAlg::Matrix<3, 1, TYPE> r1_xi(true);    // = r1,xi
    Core::LinAlg::Matrix<3, 1, TYPE> r2_xi(true);    // = r2,eta
    Core::LinAlg::Matrix<3, 1, TYPE> r1_xixi(true);  // = r1,xixi
    Core::LinAlg::Matrix<3, 1, TYPE> r2_xixi(true);  // = r2,etaeta
    Core::LinAlg::Matrix<3, 1, TYPE> delta_r(true);  // = r1-r2

    // Tangent and derivatives for tangent field smoothing (only for Reissner beams)
    Core::LinAlg::Matrix<3, 1, TYPE> t1(true);
    Core::LinAlg::Matrix<3, 1, TYPE> t1_xi(true);
    Core::LinAlg::Matrix<3, 1, TYPE> t2(true);
    Core::LinAlg::Matrix<3, 1, TYPE> t2_xi(true);

    // initialize function f and Jacobian df for Newton iteration
    Core::LinAlg::Matrix<2, 1, TYPE> f(true);
    Core::LinAlg::Matrix<2, 2, TYPE> df(true);
    Core::LinAlg::Matrix<2, 2, TYPE> dfinv(true);

    // initial scalar residual (L2-norm of f)
    double residual = 0.0;
    double lastresidual = 0.0;
    double residual0 = 0.0;
    int iter = 0;

    TYPE eta1 = startingpoints[numstartpoint].first;
    TYPE eta2 = startingpoints[numstartpoint].second;
    double eta1_old = Core::FADUtils::cast_to_double(eta1);
    double eta2_old = Core::FADUtils::cast_to_double(eta2);
    bool converged = false;
    bool elementscolinear = false;

#ifdef FADCHECKS
    BeamContact::SetFADParCoordDofs<numnodes, numnodalvalues>(eta1, eta2);
#endif

    // std::cout << "numstartpoint: " << numstartpoint << std::endl;

    //**********************************************************************
    // local Newton iteration
    //**********************************************************************
    for (int i = 0; i < BEAMCONTACTMAXITER; ++i)
    {
      // store residual of last iteration
      lastresidual = residual;
      iter++;

      // reset shape function variables to zero
      N1.clear();
      N2.clear();
      N1_xi.clear();
      N2_xi.clear();
      N1_xixi.clear();
      N2_xixi.clear();

      // update shape functions and their derivatives
      get_shape_functions(N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi, eta1, eta2);
      // update coordinates and derivatives of contact points
      compute_coords_and_derivs(
          r1, r2, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi);

      // use delta_r = r1-r2 as auxiliary quantity
      delta_r = Core::FADUtils::diff_vector(r1, r2);

      // compute norm of difference vector to scale the equations
      // (this yields better conditioning)
      // Note: Even if automatic differentiation via FAD is applied, norm_delta_r has to be of type
      // double since this factor is needed for a pure scaling of the nonlinear CCP and has not to
      // be linearized!
      double norm_delta_r = Core::FADUtils::cast_to_double(Core::FADUtils::vector_norm<3>(delta_r));
      gap = norm_delta_r - r1_ - r2_;

      // The closer the beams get, the smaller is norm_delta_r, but
      // norm_delta_r is not allowed to be too small, else numerical problems occur.
      // It can happen quite often that the centerlines of two beam elements of the same physical
      // beam cross in one point and norm_delta_r = 0. Since in this case |eta1|>1 and |eta2|>1 they
      // will be sorted out later anyways.
      // std::cout << "norm_delta_r: " << norm_delta_r << std::endl;
      if (norm_delta_r < NORMTOL)
      {
        // this excludes pairs with IDs i and i+2, i.e. contact with the next but one element
        if (Core::FADUtils::cast_to_double(Core::FADUtils::norm(eta1)) <= 1.0 and
            Core::FADUtils::cast_to_double(Core::FADUtils::norm(eta2)) <= 1.0)
        {
          this->print(std::cout);
          FOUR_C_THROW("Beam axis identical, choose smaller time step!");
        }
        else
        {
          break;
        }
      }


      // evaluate f at current eta1, eta2
      evaluate_orthogonality_condition(f, delta_r, norm_delta_r, r1_xi, r2_xi, t1, t2);

      TYPE jacobi1 = 1.0;
      TYPE jacobi2 = 1.0;
      jacobi1 = get_jacobi(element1());
      jacobi2 = get_jacobi(element2());

      // compute the scalar residuum
      // The residual is scaled with 1/element_length since an absolute
      // residual norm is used as local CPP convergence criteria and r_xi scales with the
      // element_length
      residual = sqrt((double)Core::FADUtils::cast_to_double(
          (TYPE)(f(0) * f(0) / (jacobi1 * jacobi1) + f(1) * f(1) / (jacobi2 * jacobi2))));

      //      std::cout << "iter: " << iter << std::endl;
      //      std::cout << "residual: " << residual << std::endl;
      //      std::cout << "eta1: " << eta1 << std::endl;
      //      std::cout << "eta2: " << eta2 << std::endl;
      //      std::cout << "angle: " << BeamContact::CalcAngle(r1_xi,r2_xi) << std::endl;
      //      std::cout << "iter: " << iter << std::endl;
      //      std::cout << "residual: " << residual << std::endl;

      if (iter == 1) residual0 = residual;

// check if Newton iteration has converged
#ifndef RELBEAMCONTACTTOL
      if (Core::FADUtils::cast_to_double(residual) < BEAMCONTACTTOL and
          fabs(eta1_old - Core::FADUtils::cast_to_double(eta1)) < XIETAITERATIVEDISPTOL and
          fabs(eta2_old - Core::FADUtils::cast_to_double(eta2)) < XIETAITERATIVEDISPTOL)
      {
        converged = true;
        break;
      }
#else
      if (residual0 > 1.0e-6)
      {
        if (Core::FADUtils::cast_to_double(residual / residual0) < RELBEAMCONTACTTOL and
            fabs(eta1_old - Core::FADUtils::cast_to_double(eta1)) < XIETAITERATIVEDISPTOL and
            fabs(eta2_old - Core::FADUtils::cast_to_double(eta2)) < XIETAITERATIVEDISPTOL)
        {
          converged = true;
          break;
        }
      }
      else
      {
        if (Core::FADUtils::cast_to_double(residual) < BEAMCONTACTTOL and
            fabs(eta1_old - Core::FADUtils::cast_to_double(eta1)) < XIETAITERATIVEDISPTOL and
            fabs(eta2_old - Core::FADUtils::cast_to_double(eta2)) < XIETAITERATIVEDISPTOL)
        {
          converged = true;
          break;
        }
      }
#endif

      // evaluate Jacobian of f at current eta1, eta2
      // Note: Parallel elements can not be handled with this beam contact formulation;
      evaluate_lin_orthogonality_condition(df, dfinv, delta_r, norm_delta_r, r1_xi, r2_xi, r1_xixi,
          r2_xixi, t1, t2, t1_xi, t2_xi, elementscolinear);

#ifdef FADCHECKS
      std::cout << "f: " << f << std::endl;
      std::cout << "df: " << df << std::endl;
      BeamContact::SetFADParCoordDofs<numnodes, numnodalvalues>(eta1, eta2);
      fad_check_lin_orthogonality_condition(delta_r, norm_delta_r, r1_xi, r2_xi, t1, t2);
#endif

      if (elementscolinear)
      {
        break;
      }

      eta1_old = Core::FADUtils::cast_to_double(eta1);
      eta2_old = Core::FADUtils::cast_to_double(eta2);

      // update element coordinates of contact point
      eta1 += -dfinv(0, 0) * f(0) - dfinv(0, 1) * f(1);
      eta2 += -dfinv(1, 0) * f(0) - dfinv(1, 1) * f(1);

#ifdef FADCHECKS
      BeamContact::SetFADParCoordDofs<numnodes, numnodalvalues>(eta1, eta2);
#endif

      // std::cout << "eta1: " << eta1 << std::endl;
      // std::cout << "eta2: " << eta2 << std::endl;
    }  // for (int i=0;i<BEAMCONTACTMAXITER;++i)
       //**********************************************************************

    // Newton iteration unconverged after BEAMCONTACTMAXITER
    if (!converged)
    {
      // Initialize g_min with a very large value, at which no active contact should occur!
      double g_min = 1000 * r2_;
      if (check_contact_status(g_min))
        FOUR_C_THROW("Are sure that contact should be active at such large gaps?");

      double alpha_g_min = 0.0;
      // In case no valid point-to-line solution is found in () (pointtolinesolfound=false)
      // it is assumed, that the distance between the segments is large enough such that no contact
      // can occur.
      bool pointtolinesolfound = false;
      double eta1_min = 0.0;
      double eta2_min = 0.0;

      check_unconverged_segment_pair(eta_left1, eta_left2, l1, l2, eta1_min, eta2_min, g_min,
          alpha_g_min, pointtolinesolfound);

      // Check, if we have found a valid point-to-line projection and if the solution is not a
      // boundary minimum
      if (pointtolinesolfound and eta1_min > (eta_left1 - 1.0e-10) and
          eta1_min < (eta_left1 + l1 + 1.0e-10) and eta_left2 - XIETAITERATIVEDISPTOL <= eta2 and
          eta2 <= eta_right2 + XIETAITERATIVEDISPTOL)
      {
        double perpshiftangle1 =
            params()->beam_to_beam_contact_params()->beam_to_beam_perp_shifting_angle1();
        // Here, we apply the conservative estimate that the closest-point gap is by 0.1*R2_ smaller
        // than g_min
        double g_min_estimate = g_min - 0.1 * r2_;

        // TODO
        if (check_contact_status(g_min_estimate) and fabs(alpha_g_min) >= perpshiftangle1)
        {
          std::cout << std::endl
                    << "Serious Warning!!!!! Local CPP not converged: CP-Approximation applied!"
                    << std::endl;
          std::cout << "Element1()->Id(): " << element1()->id() << std::endl;
          std::cout << "Element2()->Id(): " << element2()->id() << std::endl;
          std::cout << "R2_: " << r2_ << std::endl;
          std::cout << "g_min: " << g_min << std::endl;
          std::cout << "alpha_g_min: " << alpha_g_min / M_PI * 180 << "degrees" << std::endl;
          std::cout << "numstartpoint: " << numstartpoint << std::endl;
          std::cout << "iter: " << iter << std::endl;
          std::cout << "residual0: " << residual0 << std::endl;
          std::cout << "lastresidual: " << lastresidual << std::endl;
          std::cout << "residual: " << residual << std::endl;
          std::cout << "eta1_min: " << Core::FADUtils::cast_to_double(eta1_min) << std::endl;
          std::cout << "eta1: " << Core::FADUtils::cast_to_double(eta1) << std::endl;
          std::cout << "eta1_old: " << eta1_old << std::endl;
          std::cout << "eta2_min: " << Core::FADUtils::cast_to_double(eta2_min) << std::endl;
          std::cout << "eta2: " << Core::FADUtils::cast_to_double(eta2) << std::endl;
          std::cout << "eta2_old: " << eta2_old << std::endl;

          // We need here the original elements of the problem discretization in order to read out
          // time-dependent element quantities (such as the current curvature) since element1_ and
          // element2_ are pure copies of these elements generated once in the beginning of the
          // simulation and which do therefore not contain the current values of such element
          // quantities
          //          Core::Elements::Element* element1 =
          //          pdiscret_.lColElement(pdiscret_.ElementColMap()->LID(element1_->Id()));
          //          Core::Elements::Element* element2 =
          //          pdiscret_.lColElement(pdiscret_.ElementColMap()->LID(element2_->Id()));
          //
          //          const Core::Elements::ElementType & eot = element1->ElementType();
          //          if (eot == Discret::Elements::Beam3ebType::instance())
          //          {
          //            const Discret::Elements::Beam3eb* beam3ebelement1 = dynamic_cast<const
          //            Discret::Elements::Beam3eb*>(element1); double kappamax1 =
          //            beam3ebelement1->GetKappaMax(); const Discret::Elements::Beam3eb*
          //            beam3ebelement2 = dynamic_cast<const Discret::Elements::Beam3eb*>(element2);
          //            double kappamax2 = beam3ebelement2->GetKappaMax();
          //
          //            std::cout << "kappamax1: " << kappamax1 << std::endl;
          //            std::cout << "kappamax2: " << kappamax2 << std::endl << std::endl;
          //
          //          }

// Apply Point-To-Line solution as approximation for CPP or...
#ifdef CPP_APPROX
          eta1 = eta1_min;
          eta2 = eta2_min;
          gap = g_min;
          double angle = alpha_g_min;

          // Fixme
          if (std::fabs(eta1 - 1.0) < 1.1 * XIETAITERATIVEDISPTOL or
              std::fabs(eta1 + 1.0) < 1.1 * XIETAITERATIVEDISPTOL or
              std::fabs(eta2 - 1.0) < 1.1 * XIETAITERATIVEDISPTOL or
              std::fabs(eta2 + 1.0) < 1.1 * XIETAITERATIVEDISPTOL)
          {
            //              this->print(std::cout);
            std::cout << "\n|eta1|=1 or |eta2|=1, danger of multiple gauss point evaluation!\n";
            //              FOUR_C_THROW("|eta1|=1 or |eta2|=1, danger of multiple gauss point
            //              evaluation!");
          }

          //            double perpshiftangle1 = Params()->beam_to_beam_perp_shifting_angle1();

          if (check_contact_status(gap) and angle >= perpshiftangle1) validpairfound = true;

          solutionpoints.first = Core::FADUtils::cast_to_double(eta1);
          solutionpoints.second = Core::FADUtils::cast_to_double(eta2);

          break;

//... or abort simulation
#else
          FOUR_C_THROW(
              "CPP is not converged, eventhough the corresponding closest point is active! "
              "Decrease the value of SEGANGLE or increase your shifting angles!");
#endif
        }
      }

      eta1 = 1e+12;
      eta2 = 1e+12;
    }
    else
    {
      // if we have already found a converged solution with valid closest points eta1 \in
      // [eta_left1;eta_right1] and eta2 \in [eta_left2;eta_right2], we can finish here and don't
      // have to apply more starting points
      if (eta_left1 - XIETAITERATIVEDISPTOL <= eta1 and
          eta1 <= eta_right1 + XIETAITERATIVEDISPTOL and
          eta_left2 - XIETAITERATIVEDISPTOL <= eta2 and eta2 <= eta_right2 + XIETAITERATIVEDISPTOL)
      {
        // Fixme
        if (std::fabs(eta1 - 1.0) < 1.1 * XIETAITERATIVEDISPTOL or
            std::fabs(eta1 + 1.0) < 1.1 * XIETAITERATIVEDISPTOL or
            std::fabs(eta2 - 1.0) < 1.1 * XIETAITERATIVEDISPTOL or
            std::fabs(eta2 + 1.0) < 1.1 * XIETAITERATIVEDISPTOL)
        {
          //          this->print(std::cout);
          std::cout << "\n|eta1|=1 or |eta2|=1, danger of multiple gauss point evaluation!\n";
          //              FOUR_C_THROW("|eta1|=1 or |eta2|=1, danger of multiple gauss point
          //              evaluation!");
        }

        if (Core::FADUtils::cast_to_double(Core::FADUtils::vector_norm<3>(r1_xi)) < 1.0e-8 or
            Core::FADUtils::cast_to_double(Core::FADUtils::vector_norm<3>(r2_xi)) < 1.0e-8)
          FOUR_C_THROW("Tangent vector of zero length, choose smaller time step!");

        double angle =
            fabs(BeamInteraction::calc_angle(Core::FADUtils::cast_to_double<TYPE, 3, 1>(r1_xi),
                Core::FADUtils::cast_to_double<TYPE, 3, 1>(r2_xi)));

        double perpshiftangle1 =
            params()->beam_to_beam_contact_params()->beam_to_beam_perp_shifting_angle1();

        if (check_contact_status(gap) and angle >= perpshiftangle1) validpairfound = true;

        // Here, we perform an additional security check: If a unique CCP solution exists, the
        // Newton scheme should find it with the first starting point. Otherwise, the problem may be
        // ill-conditioned!
        if (validpairfound and numstartpoint != 0)
          FOUR_C_THROW(
              "Valid CCP solution has not been found with the first starting point. Choose smaller "
              "value of SEGANGLE!");

        solutionpoints.first = Core::FADUtils::cast_to_double(eta1);
        solutionpoints.second = Core::FADUtils::cast_to_double(eta2);

        break;
      }
    }
  }  // for (int numstartpoint=0;numstartpoint<startingpoints.size();numstartpoint++)

  return validpairfound;
}
/*----------------------------------------------------------------------*
|  end: Closest point projection
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Closest Point-To-Line Projection                         meier 10/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
bool BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::point_to_line_projection(
    double& eta1_slave, double& eta_left2, double& l2, double& eta2_master, double& gap,
    double& alpha, bool& pairactive, bool smallanglepair, bool invertpairs,
    bool orthogonalprojection)
{
  /* Attention: With the parameters invertpairs and orthogonalprojection, 4 different types of
   * projections can be realized. The parameter invertpairs=true simply changes the meaning of eta1
   * and eta2 within this method, i.e. in this case the given parameter coordinate is eta2 of beam2
   * and the searched parameter is eta1 on beam1 while in the standard case eta1 is given and eta2
   * is searched. The parameter orthogonalprojection changes the projection method: In the standard
   * case (orthogonalprojection=false), the orthogonality condition is satisfied at the projection
   * side of the searched parameter (i.e. on beam2 if invertpairs=false or beam1 if inverpairs=true)
   * while for orthogonalprojection=true the orthogonality condition is satisfied at the projection
   * side of the given parameter (i.e. on beam1 if invertpairs=false or beam2 if inverpairs=true).
   * This leads to the following four possible projections: |invertpairs|orthogonalprojection|given
   * parameter|searched parameter|orthogonality condition satisfied on| |   false   |       false |
   * eta1     |      eta2        |              beam2                 | -> default case |   true |
   * false        |      eta2     |      eta1        |              beam1                 | | false
   * |       true         |      eta1     |      eta2        |              beam1                 |
   * |   true    |       true         |      eta2     |      eta1        |              beam2 |
   */
  std::vector<double> startingpoints(3, 0.0);
  double gap_test = 0.0;
  double eta_right2 = eta_left2 + l2;

  // We start with the segment midpoint, since there it is most likely to find a closest point
  // solution within the segment
  startingpoints[0] = eta_left2 + 0.5 * l2;  // segment midpoint as starting point
  startingpoints[1] = eta_left2;             // left segment point as starting point
  startingpoints[2] = eta_left2 + l2;        // right segment point as starting point

  for (int numstartpoint = 0; numstartpoint < (int)startingpoints.size(); numstartpoint++)
  {
    // vectors for shape functions and their derivatives
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N1(true);       // = N1
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N2(true);       // = N2
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N1_xi(true);    // = N1,xi
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N2_xi(true);    // = N2,eta
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N1_xixi(true);  // = N1,xixi
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N2_xixi(true);  // = N2,etaeta

    // coords and derivatives of the two contacting points
    Core::LinAlg::Matrix<3, 1, TYPE> r1(true);       // = r1
    Core::LinAlg::Matrix<3, 1, TYPE> r2(true);       // = r2
    Core::LinAlg::Matrix<3, 1, TYPE> r1_xi(true);    // = r1,xi
    Core::LinAlg::Matrix<3, 1, TYPE> r2_xi(true);    // = r2,eta
    Core::LinAlg::Matrix<3, 1, TYPE> r1_xixi(true);  // = r1,xixi
    Core::LinAlg::Matrix<3, 1, TYPE> r2_xixi(true);  // = r2,etaeta
    Core::LinAlg::Matrix<3, 1, TYPE> delta_r(true);  // = r1-r2

    // initialize function f and Jacobian df for Newton iteration
    TYPE f = 0.0;
    TYPE df = 0.0;
    Core::LinAlg::Matrix<2, 2, TYPE> dfinv(true);

    // initial scalar residual (L2-norm of f)
    double residual = 0.0;
    double lastresidual = 0.0;
    double residual0 = 0.0;
    int iter = 0;

    TYPE eta1 = eta1_slave;
    TYPE eta2 = startingpoints[numstartpoint];
    double eta2_old = Core::FADUtils::cast_to_double(eta2);

    bool converged = false;

#ifdef FADCHECKS
    BeamContact::SetFADParCoordDofs<numnodes, numnodalvalues>(eta1, eta2);
#endif

    //**********************************************************************
    // local Newton iteration
    //**********************************************************************
    for (int i = 0; i < BEAMCONTACTMAXITER; ++i)
    {
      // store residual of last iteration
      lastresidual = residual;
      iter++;

      // reset shape function variables to zero
      N1.clear();
      N2.clear();
      N1_xi.clear();
      N2_xi.clear();
      N1_xixi.clear();
      N2_xixi.clear();

      bool inversion_possible = false;
      bool endpointpenalty = params()->beam_to_beam_contact_params()->end_point_penalty();
      if (endpointpenalty) inversion_possible = true;

#ifdef ENDPOINTSEGMENTATION
      inversion_possible = true;
#endif

      if (inversion_possible)
      {
        // In the case of ENDPOINTSEGMENTATION or ENDPOINTPENALTY it can be necessary to make an
        // invere projection (from the master beam onto the slave beam). In this case, the local
        // variables (e.g. r1, r1_xi...) inside point_to_line_projection() with index 1 represent
        // the master beam which has the global index 2. In order to get the right nodal positions
        // ele2pos_ for the local variables r1, r1_xi, r1_xixi, we have to invert the arguments of
        // the function call compute_coords_and_derivs()!
        if (!invertpairs)
        {
          // update shape functions and their derivatives
          get_shape_functions(N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi, eta1, eta2);
          // update coordinates and derivatives of contact points
          compute_coords_and_derivs(
              r1, r2, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi);
        }
        else
        {
          // update shape functions and their derivatives
          get_shape_functions(N2, N1, N2_xi, N1_xi, N2_xixi, N1_xixi, eta2, eta1);
          // update coordinates and derivatives of contact points
          compute_coords_and_derivs(
              r2, r1, r2_xi, r1_xi, r2_xixi, r1_xixi, N2, N1, N2_xi, N1_xi, N2_xixi, N1_xixi);
        }
      }
      else
      {
        // update shape functions and their derivatives
        get_shape_functions(N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi, eta1, eta2);
        // update coordinates and derivatives of contact points
        compute_coords_and_derivs(
            r1, r2, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi);
      }

      // use delta_r = r1-r2 as auxiliary quantity
      delta_r = Core::FADUtils::diff_vector(r1, r2);

      // compute norm of difference vector to scale the equations
      // (this yields better conditioning)
      // Note: Even if automatic differentiation via FAD is applied, norm_delta_r has to be of type
      // double since this factor is needed for a pure scaling of the nonlinear CCP and has not to
      // be linearized!
      double norm_delta_r = Core::FADUtils::cast_to_double(Core::FADUtils::vector_norm<3>(delta_r));
      gap_test = norm_delta_r - r1_ - r2_;

      // The closer the beams get, the smaller is norm_delta_r, but
      // norm_delta_r is not allowed to be too small, else numerical problems occur.
      // It can happen quite often that the centerlines of two beam elements of the same physical
      // beam cross in one point and norm_delta_r = 0. Since in this case |eta1|>1 and |eta2|>1 they
      // will be sorted out later anyways.
      // std::cout << "norm_delta_r: " << norm_delta_r << std::endl;
      if (norm_delta_r < NORMTOL)
      {
        // this excludes pairs with IDs i and i+2, i.e. contact with the next but one element
        if (Core::FADUtils::cast_to_double(Core::FADUtils::norm(eta1)) +
                Core::FADUtils::cast_to_double(Core::FADUtils::norm(eta2)) <
            NEIGHBORTOL)
        {
          FOUR_C_THROW("Beam axis identical, choose smaller time step!");
        }
      }

      // evaluate f at current eta1, eta2
      evaluate_ptl_orthogonality_condition(
          f, delta_r, norm_delta_r, r1_xi, r2_xi, orthogonalprojection);

      // The residual will be scaled with the length of the element whose tangent appears in the
      // orthogonality conditions. Which of the two elements is relevant in this context depends of
      // the parameters orthogonalprojection and invertpairs:
      TYPE jacobi = 1.0;

      if (invertpairs == false and orthogonalprojection == false)  // default case
      {
        jacobi = get_jacobi(element2());
      }
      else if (invertpairs == false and orthogonalprojection == true)
      {
        jacobi = get_jacobi(element1());
      }
      else if (invertpairs == true and orthogonalprojection == false)
      {
        jacobi = get_jacobi(element1());
      }
      else if (invertpairs == true and orthogonalprojection == true)
      {
        jacobi = get_jacobi(element2());
      }
      jacobi = get_jacobi(element2());

      // compute the scalar residuum
      // The residual is scaled with 1/element_length since r_xi scales with the element_length
      residual = fabs((double)Core::FADUtils::cast_to_double((TYPE)(f / jacobi)));

      //      std::cout << "iter: " << iter << std::endl;
      //      std::cout << "residual: " << residual << std::endl;
      //      std::cout << "eta1: " << eta1.val() << std::endl;
      //      std::cout << "eta2: " << eta2.val() << std::endl;
      //      std::cout << "r1: " << r1 << std::endl;
      //      std::cout << "r2: " << r2 << std::endl;
      //      std::cout << "r1_xi: " << r1_xi << std::endl;
      //      std::cout << "r2_xi: " << r2_xi << std::endl;
      //      std::cout << "r1_xixi: " << r1_xixi << std::endl;
      //      std::cout << "r2_xixi: " << r2_xixi << std::endl;
      //      std::cout << "ele1pos_: " << ele1pos_ << std::endl;
      //      std::cout << "ele2pos_: " << ele2pos_ << std::endl;
      //      std::cout << "angle: " <<
      //      BeamContact::CalcAngle(Core::FADUtils::cast_to_double<TYPE,3,1>(r1_xi),Core::FADUtils::cast_to_double<TYPE,3,1>(r2_xi))/M_PI*180.0
      //      << std::endl;

      if (iter == 1) residual0 = residual;

// check if Newton iteration has converged
#ifndef RELBEAMCONTACTTOL
      if (Core::FADUtils::cast_to_double(residual) < BEAMCONTACTTOL and
          fabs(eta2_old - Core::FADUtils::cast_to_double(eta2)) < XIETAITERATIVEDISPTOL)
      {
        converged = true;
        break;
      }
#else
      if (residual0 > 1.0e-6)
      {
        if (Core::FADUtils::cast_to_double(residual / residual0) < RELBEAMCONTACTTOL and
            fabs(eta2_old - Core::FADUtils::cast_to_double(eta2)) < XIETAITERATIVEDISPTOL)
        {
          converged = true;
          break;
        }
      }
      else
      {
        if (Core::FADUtils::cast_to_double(residual) < BEAMCONTACTTOL and
            fabs(eta2_old - Core::FADUtils::cast_to_double(eta2)) < XIETAITERATIVEDISPTOL)
        {
          converged = true;
          break;
        }
      }
#endif

      // evaluate Jacobian of f at current eta1, eta2
      // Note: It has to be checked, if the linearization is equal to zero;
      bool validlinearization = evaluate_lin_ptl_orthogonality_condition(
          df, delta_r, norm_delta_r, r1_xi, r2_xi, r2_xixi, orthogonalprojection);

      if (!validlinearization)
        FOUR_C_THROW(
            "Linearization of point to line projection is zero, choose tighter search boxes!");

#ifdef FADCHECKS
      std::cout << "f: " << f << std::endl;
      std::cout << "df: " << df << std::endl;
      BeamContact::SetFADParCoordDofs<numnodes, numnodalvalues>(eta1, eta2);
      fad_check_lin_orthogonality_condition(delta_r, norm_delta_r, r1_xi, r2_xi, t1, t2);
#endif

      eta2_old = Core::FADUtils::cast_to_double(eta2);

      // update master element coordinate of contact point
      eta2 += -f / df;

#ifdef FADCHECKS
      BeamContact::SetFADParCoordDofs<numnodes, numnodalvalues>(eta1, eta2);
#endif

    }  // for (int i=0;i<BEAMCONTACTMAXITER;++i)
       //**********************************************************************

    // Newton iteration unconverged after BEAMCONTACTMAXITER
    if (!converged)
    {
      std::cout << "iter: " << iter << std::endl;
      std::cout << "residual0: " << residual0 << std::endl;
      std::cout << "lastresidual: " << lastresidual << std::endl;
      std::cout << "residual: " << residual << std::endl;
      std::cout << "eta2: " << Core::FADUtils::cast_to_double(eta2) << std::endl;
      std::cout << "eta2_old: " << eta2_old << std::endl;


      // TODO:
      this->print(std::cout);
      FOUR_C_THROW(
          "Local Newton loop unconverged. Adapt segangle or the shift angles for small-angle "
          "contact!");

      eta1 = 1e+12;
      eta2 = 1e+12;
    }
    else
    {
      // if we have already found a converged solution with valid closest point eta2 \in
      // [eta_left2;eta_right2], we can finish here and don't have to apply more starting points
      if (eta_left2 - XIETAITERATIVEDISPTOL <= eta2 and eta2 <= eta_right2 + XIETAITERATIVEDISPTOL)
      {
        if (fabs(eta2 - 1.0) < 1.1 * XIETAITERATIVEDISPTOL)
        {
          bool throw_error = true;

          // There is no danger of multiple evaluation, if the considered node is a boundary node of
          // the beam:
          if ((invertpairs == false and boundarynode2_.second) or
              (invertpairs == true and boundarynode1_.second))
            throw_error = false;

          if (throw_error)
          {
            std::cout << "ID1: " << element1()->id() << std::endl;
            std::cout << "ID2: " << element2()->id() << std::endl;
            std::cout << "eta1: " << eta1 << std::endl;
            std::cout << "eta2: " << eta2 << std::endl;
            // TODO: In some cases a warning is sufficient, but in general we need the
            // FOUR_C_THROW("");
            std::cout << "Serious Warning!!!!! eta2=1, danger of multiple gauss point evaluation! "
                      << std::endl;
            // FOUR_C_THROW("eta2=1, danger of multiple gauss point evaluation!");
          }
        }

        if (fabs(eta2 + 1.0) < 1.1 * XIETAITERATIVEDISPTOL)
        {
          bool throw_error = true;

          // There is no danger of multiple evaluation, if the considered node is a boundary node of
          // the beam:
          if ((invertpairs == false and boundarynode2_.first) or
              (invertpairs == true and boundarynode1_.first))
            throw_error = false;

          if (throw_error)
          {
            this->print(std::cout);
            FOUR_C_THROW("eta2=-1, danger of multiple gauss point evaluation!");
          }
        }

        if (Core::FADUtils::cast_to_double(Core::FADUtils::vector_norm<3>(r1_xi)) < 1.0e-8 or
            Core::FADUtils::cast_to_double(Core::FADUtils::vector_norm<3>(r2_xi)) < 1.0e-8)
          FOUR_C_THROW("Tangent vector of zero length, choose smaller time step!");

        bool relevant_angle = true;
        double angle =
            fabs(BeamInteraction::calc_angle(Core::FADUtils::cast_to_double<TYPE, 3, 1>(r1_xi),
                Core::FADUtils::cast_to_double<TYPE, 3, 1>(r2_xi)));
        if (smallanglepair)
        {
          double parshiftangle2 =
              params()->beam_to_beam_contact_params()->beam_to_beam_parallel_shifting_angle2();

          if (angle > parshiftangle2) relevant_angle = false;
        }
        if (check_contact_status(gap_test) and relevant_angle) pairactive = true;

        eta2_master = Core::FADUtils::cast_to_double(eta2);

        // Here, we perform an additional security check: If a unique CCP solution exists, the
        // Newton scheme should find it with the first starting point. Otherwise, the problem may be
        // ill-conditioned!
        if (pairactive and numstartpoint != 0)
          FOUR_C_THROW(
              "Valid Point-To-Line solution has not been found with the first starting point. "
              "Choose smaller value of SEGANGLE!");

        gap = gap_test;
        alpha = angle;

        // eta2 \in [eta_left2;eta_right2] --> true
        return true;
      }
    }
  }  // for (int numstartpoint=0;numstartpoint<startingpoints.size();numstartpoint++)

  // no eta2 \in [eta_left2;eta_right2] --> false
  return false;
}
/*----------------------------------------------------------------------*
|  end: Point-To-Line Projection
 *----------------------------------------------------------------------*/

/*------------------------------------------------------------------------------------------*
|  Determine minimal distance and contact angle for unconverged segment pair     meier 05/15|
 *------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes,
    numnodalvalues>::check_unconverged_segment_pair(double& eta_left1, double& eta_left2,
    double& l1, double& l2, double& eta1_min, double& eta2_min, double& g_min, double& alpha_g_min,
    bool& pointtolinesolfound)
{
  // Calculate initial length of slave element
  Core::LinAlg::Matrix<3, 1, double> lengthvec1(true);
  for (int i = 0; i < 3; i++)
  {
    lengthvec1(i) = (element1()->nodes())[0]->x()[i] - (element1()->nodes())[1]->x()[i];
  }
  // length1 = physical length; l1=length in parameter space
  double length1 = lengthvec1.norm2();

  int n = 1;
  // subdivide the slave segment by n+1 test points until the distance between the
  // test points is smaller than half of the cross-section radius
  while (l1 / 2 * length1 / n > r2_ / 2)
  {
    n = 2 * n;
  }

  double eta1_closestpoint = 0.0;
  double eta2_closestpoint = 0.0;

  for (int i = 0; i < n + 1; i++)
  {
    double eta1_slave = eta_left1 + i * l1 / n;
    double eta2_segleft = eta_left2;
    double eta2_master = 0.0;
    bool pairactive = false;
    double gap = 0.0;
    double alpha = 0.0;

    bool solutionwithinsegment = point_to_line_projection(
        eta1_slave, eta2_segleft, l2, eta2_master, gap, alpha, pairactive, true);

    if (solutionwithinsegment)
    {
      pointtolinesolfound = true;
      if (gap < g_min)
      {
        g_min = gap;
        alpha_g_min = alpha;
        eta1_closestpoint = eta1_slave;
        eta2_closestpoint = eta2_master;
      }  // search also for the second-smallest gap
    }
  }

  bool cp_at_right_neighbor = false;

  // if we have a boundary minimum on the left, we also investigate the left neighbor point
  if (fabs(eta1_closestpoint - eta_left1) < 1.0e-10)
  {
    double eta1_slave = eta_left1 - l1 / n;
    double eta2_segleft = eta_left2;
    double eta2_master = 0.0;
    bool pairactive = false;
    double gap = 0.0;
    double alpha = 0.0;

    bool solutionwithinsegment = point_to_line_projection(
        eta1_slave, eta2_segleft, l2, eta2_master, gap, alpha, pairactive, true);

    if (solutionwithinsegment)
    {
      if (gap < g_min)
      {
        cp_at_right_neighbor = true;
        g_min = gap;
        alpha_g_min = alpha;
        eta1_closestpoint = eta1_slave;
        eta2_closestpoint = eta2_master;
      }  // search also for the second-smallest gap
    }
  }

  // if we have a boundary minimum on the right, we also investigate the right neighbor point
  if (fabs(eta1_closestpoint - (eta_left1 + l1)) < 1.0e-10)
  {
    double eta1_slave = eta_left1 + (n + 1) * l1 / n;
    double eta2_segleft = eta_left2;
    double eta2_master = 0.0;
    bool pairactive = false;
    double gap = 0.0;
    double alpha = 0.0;

    bool solutionwithinsegment = point_to_line_projection(
        eta1_slave, eta2_segleft, l2, eta2_master, gap, alpha, pairactive, true);

    if (solutionwithinsegment)
    {
      if (gap < g_min)
      {
        if (cp_at_right_neighbor)
          FOUR_C_THROW(
              "This should not happen, that we have a local minimum on the right and on the left "
              "neighbor!");

        g_min = gap;
        alpha_g_min = alpha;
        eta1_closestpoint = eta1_slave;
        eta2_closestpoint = eta2_master;
      }  // search also for the second-smallest gap
    }
  }


  eta1_min = eta1_closestpoint;
  eta2_min = eta2_closestpoint;

  return;
}
/*------------------------------------------------------------------------------------------*
|  end: Determine minimal distance and contact angle for unconverged segment pair
 *------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Compute contact forces                                   meier 10/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::evaluate_fc_contact(
    Core::LinAlg::SerialDenseVector& forcevec1, Core::LinAlg::SerialDenseVector& forcevec2,
    const Core::LinAlg::Matrix<3, 1, TYPE>& r1, const Core::LinAlg::Matrix<3, 1, TYPE>& r2,
    const Core::LinAlg::Matrix<3, 1, TYPE>& r1_xi, const Core::LinAlg::Matrix<3, 1, TYPE>& r2_xi,
    const Core::LinAlg::Matrix<3, 1, TYPE>& r1_xixi,
    const Core::LinAlg::Matrix<3, 1, TYPE>& r2_xixi,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xi,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xi,
    BeamToBeamContactVariables<numnodes, numnodalvalues>& variables, const double& intfac, bool cpp,
    bool gp, bool fixedendpointxi, bool fixedendpointeta,
    Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, TYPE>* fc1_FAD,
    Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, TYPE>* fc2_FAD)
{
  // Check for sensible combinations:
  if ((cpp and (gp or fixedendpointxi or fixedendpointeta)) or
      (gp and (fixedendpointxi or fixedendpointeta)))
    FOUR_C_THROW(
        "This is no possible combination of the parameters cpp, gp, fixedendpointxi and "
        "fixedendpointeta!");

  // get dimensions for vectors fc1 and fc2
  const unsigned int dim1 = 3 * numnodes * numnodalvalues;
  const unsigned int dim2 = 3 * numnodes * numnodalvalues;

  // temporary vectors for contact forces, DOF-GIDs and owning procs
  Core::LinAlg::Matrix<dim1, 1, TYPE> fc1(true);
  Core::LinAlg::Matrix<dim2, 1, TYPE> fc2(true);

  // flag indicating assembly
  bool DoNotAssemble = true;

  TYPE gap = variables.get_gap();
  Core::LinAlg::Matrix<3, 1, TYPE> normal = variables.get_normal();
  TYPE fp = variables.getfp();
  // The factor ppfac reduces the penalty parameter for the large-angle and small-angle formulation
  // in dependence of the current contact angle
  TYPE ppfac = variables.get_p_pfac();

  //**********************************************************************
  // evaluate contact forces for active pairs
  //**********************************************************************
  if (check_contact_status(Core::FADUtils::cast_to_double(gap)))
  {
    DoNotAssemble = false;
#ifndef CONSISTENTTRANSITION
    //********************************************************************
    // Compute Fc1 (force acting on first element)
    //********************************************************************
    // The variable intfac represents the integration factor containing the
    // Gauss weight and the jacobian. This factor is only necessary for the
    // small-angle formulation and is set to 1.0 otherwise!
    for (unsigned int i = 0; i < dim1; ++i)
    {
      for (unsigned int j = 0; j < 3; ++j)
      {
        fc1(i) += N1(j, i) * normal(j) * fp * ppfac * intfac;
      }
    }

    //********************************************************************
    // Compute Fc2 (force acting on second element)
    //********************************************************************
    // The variable intfac represents the integration factor containing the
    // Gauss weight and the jacobian. This factor is only necessary for the
    // small-angle formulation and is set to 1.0 otherwise!
    for (unsigned int i = 0; i < dim2; ++i)
    {
      for (unsigned int j = 0; j < 3; ++j)
      {
        fc2(i) += -N2(j, i) * normal(j) * fp * ppfac * intfac;
      }
    }
#else
    // initialize storage for linearizations
    Core::LinAlg::Matrix<dim1 + dim2, 1, TYPE> delta_xi(true);
    Core::LinAlg::Matrix<dim1 + dim2, 1, TYPE> delta_eta(true);
    Core::LinAlg::Matrix<dim1 + dim2, 1, TYPE> delta_coscontactangle(true);

    Core::LinAlg::Matrix<3, 1, TYPE> delta_r = Core::FADUtils::DiffVector(r1, r2);
    TYPE norm_delta_r = Core::FADUtils::vector_norm<3>(delta_r);
    Core::LinAlg::Matrix<3, 1, TYPE> normal = variables->GetNormal();
    TYPE fp = variables->Getfp();
    TYPE dfp = variables->Getdfp();
    TYPE dppfac = variables->GetDPPfac();
    TYPE e = variables->get_energy();

    // linearization of contact point
    if (cpp)  // in case of large-angle-contact (standard closest-point-projection), we need
              // delta_xi and delta_eta.
    {
      compute_lin_xi_and_lin_eta(
          delta_xi, delta_eta, delta_r, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi);
    }
    else if (gp or (fixedendpointxi and
                       !fixedendpointeta))  // in case of small-angle-contact (xi remains fixed), we
                                            // only need delta_eta, delta_xi remains zero (this does
                                            // not hold in case of ENDPOINTSEGMENTATION)
    {  // this also holds in case of ENDPOINTPENALTY when the endpoint xi is fixed and the endpoint
       // eta not!
      compute_lin_eta_fix_xi(delta_eta, delta_r, r2_xi, r2_xixi, N1, N2, N2_xi);
#ifdef ENDPOINTSEGMENTATION
      FOUR_C_THROW(
          "The combination of ENDPOINTSEGMENTATION and CONSISTENTTRANSITION is not possible!");
#endif
    }
    else if (fixedendpointeta and !fixedendpointxi)  // In case of ENDPOINTPENALTY when the endpoint
                                                     // eta is fixed and the endpoint xi not...
    {
      compute_lin_xi_fix_eta(delta_xi, delta_r, r1_xi, r1_xixi, N2, N1, N1_xi);
    }
    else if (fixedendpointeta and fixedendpointxi)  // In case of ENDPOINTPENALTY when the endpoint
                                                    // eta is fixed and the endpoint xi is fixed...
    {
      //..we need to do nothing since delta_xi and delta_eta have to remain zero!
    }

    // linearization of large-angle/small-angle scale factor
    compute_lin_cos_contact_angle(
        delta_coscontactangle, delta_xi, delta_eta, r1_xi, r2_xi, r1_xixi, r2_xixi, N1_xi, N2_xi);
    //********************************************************************
    // Compute Fc1 (force acting on first element)
    //********************************************************************
    // The variable intfac represents the integration factor containing the
    // Gauss weight and the jacobian. This factor is only necessary for the
    // small-angle formulation and is set to 1.0 otherwise!
    for (unsigned int i = 0; i < dim1; ++i)
    {
      for (unsigned int j = 0; j < 3; ++j)
      {
        fc1(i) += (N1(j, i) * normal(j) * fp * ppfac) * intfac;
      }
      fc1(i) += (e * dppfac * delta_coscontactangle(i)) * intfac;
    }

    //********************************************************************
    // Compute Fc2 (force acting on second element)
    //********************************************************************
    // The variable intfac represents the integration factor containing the
    // Gauss weight and the jacobian. This factor is only necessary for the
    // small-angle formulation and is set to 1.0 otherwise!
    for (unsigned int i = 0; i < dim2; ++i)
    {
      for (unsigned int j = 0; j < 3; ++j)
      {
        fc2(i) += (-N2(j, i) * normal(j) * fp * ppfac) * intfac;
      }
      fc2(i) += (e * dppfac * delta_coscontactangle(dim1 + i)) * intfac;
    }
#endif
  }

// Quantities necessary for automatic differentiation
#ifdef AUTOMATICDIFF
  if (fc1_FAD != nullptr and fc2_FAD != nullptr)
  {
    for (unsigned int i = 0; i < dim1; ++i)
    {
      (*fc1_FAD)(i) = fc1(i);
    }
    for (unsigned int i = 0; i < dim2; ++i)
    {
      (*fc2_FAD)(i) = fc2(i);
    }
  }
#endif

  //**********************************************************************
  // assemble contact forces
  //**********************************************************************
  if (not DoNotAssemble)
  {
    for (unsigned int i = 0; i < dim1; ++i) forcevec1(i) += Core::FADUtils::cast_to_double(fc1(i));

    for (unsigned int i = 0; i < dim2; ++i) forcevec2(i) += Core::FADUtils::cast_to_double(fc2(i));
  }
}
/*----------------------------------------------------------------------*
 |  end: Compute contact forces
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Evaluate contact stiffness                               meier 10/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::evaluate_stiffc_contact(
    Core::LinAlg::SerialDenseMatrix& stiffmat11, Core::LinAlg::SerialDenseMatrix& stiffmat12,
    Core::LinAlg::SerialDenseMatrix& stiffmat21, Core::LinAlg::SerialDenseMatrix& stiffmat22,
    const Core::LinAlg::Matrix<3, 1, TYPE>& r1, const Core::LinAlg::Matrix<3, 1, TYPE>& r2,
    const Core::LinAlg::Matrix<3, 1, TYPE>& r1_xi, const Core::LinAlg::Matrix<3, 1, TYPE>& r2_xi,
    const Core::LinAlg::Matrix<3, 1, TYPE>& r1_xixi,
    const Core::LinAlg::Matrix<3, 1, TYPE>& r2_xixi,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xi,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xi,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xixi,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xixi,
    BeamToBeamContactVariables<numnodes, numnodalvalues>& variables, const double& intfac, bool cpp,
    bool gp, bool fixedendpointxi, bool fixedendpointeta)
{
  // get dimensions for vectors fc1 and fc2
  const unsigned int dim1 = 3 * numnodes * numnodalvalues;
  const unsigned int dim2 = 3 * numnodes * numnodalvalues;

  // temporary matrices for stiffness and vectors for DOF-GIDs and owning procs
  Core::LinAlg::Matrix<dim1, dim1 + dim2, TYPE> stiffc1(true);
  Core::LinAlg::Matrix<dim2, dim1 + dim2, TYPE> stiffc2(true);
  Core::LinAlg::Matrix<dim1, dim1 + dim2, TYPE> stiffc1_FAD(true);
  Core::LinAlg::Matrix<dim2, dim1 + dim2, TYPE> stiffc2_FAD(true);

  // flag indicating assembly
  bool DoNotAssemble = true;
  TYPE gap = variables.get_gap();
  // The factor ppfac reduces the penalty parameter for the large-angle and small-angle formulation
  // in dependence of the current contact angle
  TYPE ppfac = variables.get_p_pfac();
  TYPE dppfac = variables.get_dp_pfac();

  // In order to accelerate convergence, we only apply the basic stiffness part in case of very
  // large gaps!
  double basicstiffgap = params()->beam_to_beam_contact_params()->basic_stiff_gap();
  bool completestiff = true;
  if (basicstiffgap != -1.0)
  {
    if (basicstiffgap < 0.0)
      FOUR_C_THROW("The parameter BEAMS_BASICSTIFFGAP has to be positive!");
    else if (gap < -1.0 * basicstiffgap)
    {
      completestiff = false;
    }
  }

  // Apply additional weighting of the basic stiffness term e.g. in the first iterations or when
  // the Newton scheme oscillates (no convergence after a certain number of iterations)
  double basicstiffweightfac = 1.0;
#ifdef BASICSTIFFWEIGHT
  if (iter_ < 5)
  {
    basicstiffweightfac = BASICSTIFFWEIGHT;
  }
#endif

  //**********************************************************************
  // evaluate contact stiffness for active pairs
  //**********************************************************************
  if (check_contact_status(Core::FADUtils::cast_to_double(gap)))
  {
    DoNotAssemble = false;

    // initialize storage for linearizations
    Core::LinAlg::Matrix<dim1 + dim2, 1, TYPE> delta_xi(true);
    Core::LinAlg::Matrix<dim1 + dim2, 1, TYPE> delta_eta(true);
    Core::LinAlg::Matrix<dim1 + dim2, 1, TYPE> delta_gap(true);
    Core::LinAlg::Matrix<dim1 + dim2, 1, TYPE> delta_gap_t(true);
    Core::LinAlg::Matrix<3, dim1 + dim2, TYPE> delta_x1_minus_x2(true);
    Core::LinAlg::Matrix<3, dim1 + dim2, TYPE> delta_n(true);
    Core::LinAlg::Matrix<dim1 + dim2, 1, TYPE> delta_coscontactangle(true);

    Core::LinAlg::Matrix<3, 1, TYPE> delta_r = Core::FADUtils::diff_vector(r1, r2);
    TYPE norm_delta_r = Core::FADUtils::vector_norm<3>(delta_r);
    Core::LinAlg::Matrix<3, 1, TYPE> normal = variables.get_normal();
    TYPE fp = variables.getfp();
    TYPE dfp = variables.getdfp();

    //********************************************************************
    // evaluate linearizations and distance
    //********************************************************************

    // linearization of contact point
    if (cpp)  // in case of large-angle-contact (standard closest-point-projection), we need
              // delta_xi and delta_eta.
    {
      compute_lin_xi_and_lin_eta(
          delta_xi, delta_eta, delta_r, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi);
    }
    else if (gp or (fixedendpointxi and
                       !fixedendpointeta))  // in case of small-angle-contact (xi remains fixed), we
                                            // only need delta_eta, delta_xi remains zero (this does
                                            // not hold in case of ENDPOINTSEGMENTATION)
    {  // this also holds in case of ENDPOINTPENALTY when the endpoint xi is fixed and the endpoint
       // eta not!
      compute_lin_eta_fix_xi(delta_eta, delta_r, r2_xi, r2_xixi, N1, N2, N2_xi);
    }
    else if (fixedendpointeta and !fixedendpointxi)  // In case of ENDPOINTPENALTY when the endpoint
                                                     // eta is fixed and the endpoint xi not...
    {
      compute_lin_xi_fix_eta(delta_xi, delta_r, r1_xi, r1_xixi, N2, N1, N1_xi);
    }
    else if (fixedendpointeta and fixedendpointxi)  // In case of ENDPOINTPENALTY when the endpoint
                                                    // eta is fixed and the endpoint xi is fixed...
    {
      //..we need to do nothing since delta_xi and delta_eta have to remain zero!
    }

    // linearization of gap function which is equal to delta d
    compute_lin_gap(delta_gap, delta_xi, delta_eta, delta_r, norm_delta_r, r1_xi, r2_xi, N1, N2);

    // linearization of normal vector
    compute_lin_normal(delta_n, delta_xi, delta_eta, delta_r, r1_xi, r2_xi, N1, N2);

    // linearization of large-angle/small-angle scale factor
    compute_lin_cos_contact_angle(
        delta_coscontactangle, delta_xi, delta_eta, r1_xi, r2_xi, r1_xixi, r2_xixi, N1_xi, N2_xi);

#ifdef FADCHECKS
    std::cout << "delta_xi: " << std::endl;
    for (unsigned int i = 0; i < dim1 + dim2; i++) std::cout << delta_xi(i).val() << "  ";
    std::cout << std::endl << "delta_eta: " << std::endl;
    for (unsigned int i = 0; i < dim1 + dim2; i++) std::cout << delta_eta(i).val() << "  ";
    std::cout << std::endl;
    fad_check_lin_xi_and_lin_eta(delta_r, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi);
#endif

    //*************Begin of standard linearization of penalty contact forces**************
    // The full contact stiffness is only applied if the contact flag is true
    // and gap_ > -BEAMS_BASICSTIFFGAP. If gap_ < -BEAMS_BASICSTIFFGAP, only
    // the basic stiffness is applied.

    //********************************************************************
    // evaluate contact stiffness
    // (1) stiffc1 of first element
    //********************************************************************

    //********************************************************************
    // part I - basic stiffness
    //********************************************************************

    Core::LinAlg::Matrix<dim1, 1, TYPE> N1T_normal(true);
    for (unsigned int i = 0; i < 3; i++)
    {
      for (unsigned int j = 0; j < dim1; j++)
      {
        N1T_normal(j) += N1(i, j) * normal(i);
      }
    }
    for (unsigned int i = 0; i < dim1; i++)
    {
      for (unsigned int j = 0; j < dim1 + dim2; j++)
      {
        stiffc1(i, j) += basicstiffweightfac * N1T_normal(i) *
                         (ppfac * dfp * delta_gap(j) + dppfac * fp * delta_coscontactangle(j));
      }
    }

    // The geometric part is only applied for gap_ < -BEAMS_BASICSTIFFGAP
    if (completestiff)
    {
      //********************************************************************
      // part II - geometric stiffness 1
      //********************************************************************
      for (unsigned int i = 0; i < 3; i++)
      {
        for (unsigned int j = 0; j < dim1; j++)
        {
          for (unsigned int k = 0; k < dim1 + dim2; k++)
          {
            stiffc1(j, k) += ppfac * fp * N1(i, j) * delta_n(i, k);
          }
        }
      }
      //********************************************************************
      // part III - geometric stiffness 2
      //********************************************************************
      Core::LinAlg::Matrix<dim1, 1, TYPE> N1xiT_normal(true);
      for (unsigned int i = 0; i < 3; i++)
      {
        for (unsigned int j = 0; j < dim1; j++)
        {
          N1xiT_normal(j) += N1_xi(i, j) * normal(i);
        }
      }

      for (unsigned int i = 0; i < dim1; i++)
      {
        for (unsigned int j = 0; j < dim1 + dim2; j++)
        {
          stiffc1(i, j) += ppfac * fp * N1xiT_normal(i) * delta_xi(j);
        }
      }
    }
    //********************************************************************
    // evaluate contact stiffness
    // (2) stiffc2 of second element
    //********************************************************************

    //********************************************************************
    // part I
    //********************************************************************
    Core::LinAlg::Matrix<dim2, 1, TYPE> N2T_normal(true);
    for (unsigned int i = 0; i < 3; i++)
    {
      for (unsigned int j = 0; j < dim2; j++)
      {
        N2T_normal(j) += N2(i, j) * normal(i);
      }
    }
    for (unsigned int i = 0; i < dim2; i++)
    {
      for (unsigned int j = 0; j < dim1 + dim2; j++)
      {
        stiffc2(i, j) += -basicstiffweightfac * N2T_normal(i) *
                         (ppfac * dfp * delta_gap(j) + dppfac * fp * delta_coscontactangle(j));
      }
    }

    if (completestiff)
    {
      //********************************************************************
      // part II
      //********************************************************************
      for (unsigned int i = 0; i < 3; i++)
      {
        for (unsigned int j = 0; j < dim2; j++)
        {
          for (unsigned int k = 0; k < dim1 + dim2; k++)
          {
            stiffc2(j, k) += -ppfac * fp * N2(i, j) * delta_n(i, k);
          }
        }
      }
      //********************************************************************
      // part III
      //********************************************************************
      Core::LinAlg::Matrix<dim1, 1, TYPE> N2xiT_normal(true);
      for (unsigned int i = 0; i < 3; i++)
      {
        for (unsigned int j = 0; j < dim2; j++)
        {
          N2xiT_normal(j) += N2_xi(i, j) * normal(i);
        }
      }

      for (unsigned int i = 0; i < dim2; i++)
      {
        for (unsigned int j = 0; j < dim1 + dim2; j++)
        {
          stiffc2(i, j) += -ppfac * fp * N2xiT_normal(i) * delta_eta(j);
        }
      }
    }
    //*************End of standard linearization of penalty contact forces****************

    stiffc1.scale(intfac);
    stiffc2.scale(intfac);

// automatic differentiation for debugging
#ifdef AUTOMATICDIFF
    FOUR_C_THROW(
        "check implementation of AUTOMATICDIFF for BeamToBeamContactPair before using it!");

    Core::LinAlg::Matrix<dim1, 1, TYPE> fc1_FAD(true);
    Core::LinAlg::Matrix<dim2, 1, TYPE> fc2_FAD(true);
    evaluate_fc_contact(nullptr, r1, r2, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi,
        variables, intfac, cpp, gp, fixedendpointxi, fixedendpointeta, &fc1_FAD, &fc2_FAD);

    if (cpp)  // in case of large-angle-contact (standard closest-point-projection), we need
              // delta_xi and delta_eta.
    {
      for (unsigned int j = 0; j < dim1 + dim2; j++)
      {
        for (unsigned int i = 0; i < dim1; i++)
          stiffc1_FAD(i, j) = (fc1_FAD(i).dx(j) + fc1_FAD(i).dx(dim1 + dim2) * delta_xi(j) +
                               fc1_FAD(i).dx(dim1 + dim2 + 1) * delta_eta(j));
        for (unsigned int i = 0; i < dim2; i++)
          stiffc2_FAD(i, j) = (fc2_FAD(i).dx(j) + fc2_FAD(i).dx(dim1 + dim2) * delta_xi(j) +
                               fc2_FAD(i).dx(dim1 + dim2 + 1) * delta_eta(j));
      }
    }
    else if (gp or (fixedendpointxi and
                       !fixedendpointeta))  // in case of small-angle-contact (xi remains fixed), we
                                            // only need delta_eta, delta_xi remains zero (this does
                                            // not hold in case of ENDPOINTSEGMENTATION)
    {  // this also holds in case of ENDPOINTPENALTY when the endpoint xi is fixed and the endpoint
       // eta not!
      for (unsigned int j = 0; j < dim1 + dim2; j++)
      {
        for (unsigned int i = 0; i < dim1; i++)
          stiffc1_FAD(i, j) = (fc1_FAD(i).dx(j) + fc1_FAD(i).dx(dim1 + dim2 + 1) * delta_eta(j));
        for (unsigned int i = 0; i < dim2; i++)
          stiffc2_FAD(i, j) = (fc2_FAD(i).dx(j) + fc2_FAD(i).dx(dim1 + dim2 + 1) * delta_eta(j));
      }
    }
    else if (fixedendpointeta and !fixedendpointxi)  // In case of ENDPOINTPENALTY when the endpoint
                                                     // eta is fixed and the endpoint xi not...
    {
      for (unsigned int j = 0; j < dim1 + dim2; j++)
      {
        for (unsigned int i = 0; i < dim1; i++)
          stiffc1_FAD(i, j) = (fc1_FAD(i).dx(j) + fc1_FAD(i).dx(dim1 + dim2) * delta_xi(j));
        for (unsigned int i = 0; i < dim2; i++)
          stiffc2_FAD(i, j) = (fc2_FAD(i).dx(j) + fc2_FAD(i).dx(dim1 + dim2) * delta_xi(j));
      }
    }
    else if (fixedendpointeta and fixedendpointxi)  // In case of ENDPOINTPENALTY when the endpoint
                                                    // eta is fixed and the endpoint xi is fixed...
    {
      for (unsigned int j = 0; j < dim1 + dim2; j++)
      {
        for (unsigned int i = 0; i < dim1; i++) stiffc1_FAD(i, j) = (fc1_FAD(i).dx(j));
        for (unsigned int i = 0; i < dim2; i++) stiffc2_FAD(i, j) = (fc2_FAD(i).dx(j));
      }
    }

    //      std::cout << "Pair: " << Element1()->Id() << " / " << Element2()->Id() << std::endl;

    //      std::cout << "stiffc1: " << std::endl;
    for (unsigned int i = 0; i < dim1; i++)
    {
      for (unsigned int j = 0; j < dim1 + dim2; j++)
      {
        //          std::cout << stiffc1(i,j).val() << " ";
        if (fabs(stiffc1(i, j).val()) > 1.0e-7 and
            fabs((stiffc1(i, j).val() - stiffc1_FAD(i, j).val()) / stiffc1(i, j).val()) > 1.0e-7)
        {
          // std::cout << std::endl << std::endl << "stiffc1(i,j).val(): " << stiffc1(i,j).val() <<
          // "   stiffc1_FAD(i,j).val(): " << stiffc1_FAD(i,j).val() << std::endl;
          // FOUR_C_THROW("Error in linearization!");
        }
      }
      //        std::cout << std::endl;
    }
    //      std::cout << std::endl;
    //      std::cout << "stiffc1_FAD: " << std::endl;
    //      for (unsigned int i=0;i<dim1;i++)
    //      {
    //        for (unsigned int j=0;j<dim1+dim2;j++)
    //        {
    //          std::cout << stiffc1_FAD(i,j).val() << " ";
    //        }
    //        std::cout << std::endl;
    //      }
    //      std::cout << std::endl;
    //      std::cout << "stiffc2: " << std::endl;
    for (unsigned int i = 0; i < dim1; i++)
    {
      for (unsigned int j = 0; j < dim1 + dim2; j++)
      {
        //          std::cout << stiffc2(i,j).val() << " ";
        if (fabs(stiffc2(i, j).val()) > 1.0e-7 and
            fabs((stiffc2(i, j).val() - stiffc2_FAD(i, j).val()) / stiffc2(i, j).val()) > 1.0e-7)
        {
          // std::cout << std::endl << std::endl <<"stiffc2(i,j).val(): " << stiffc2(i,j).val() << "
          // stiffc2_FAD(i,j).val(): " << stiffc2_FAD(i,j).val() << std::endl; FOUR_C_THROW("Error
          // in linearization!");
        }
      }
      //        std::cout << std::endl;
    }
    //      std::cout << std::endl;
    //      std::cout << "stiffc2_FAD: " << std::endl;
    //      for (unsigned int i=0;i<dim1;i++)
    //      {
    //        for (unsigned int j=0;j<dim1+dim2;j++)
    //        {
    //          std::cout << stiffc2_FAD(i,j).val() << " ";
    //        }
    //        std::cout << std::endl;
    //      }
    //      std::cout << std::endl;
#endif
  }  // if (check_contact_status(gap))

  //**********************************************************************
  // assemble contact stiffness
  //**********************************************************************
  // change sign of stiffc1 and stiffc2 due to time integration.
  // according to analytical derivation there is no minus sign, but for
  // our time integration methods the negative stiffness must be assembled.

  // Fixme check the minus sign !!!

  // now finally assemble stiffc1 and stiffc2
  if (not DoNotAssemble)
  {
#ifndef AUTOMATICDIFF
    for (unsigned int j = 0; j < dim1; j++)
    {
      for (unsigned int i = 0; i < dim1; i++)
        stiffmat11(i, j) += -Core::FADUtils::cast_to_double(stiffc1(i, j));
      for (unsigned int i = 0; i < dim2; i++)
        stiffmat21(i, j) += -Core::FADUtils::cast_to_double(stiffc2(i, j));
    }
    for (unsigned int j = 0; j < dim2; j++)
    {
      for (unsigned int i = 0; i < dim1; i++)
        stiffmat12(i, j) += -Core::FADUtils::cast_to_double(stiffc1(i, dim1 + j));
      for (unsigned int i = 0; i < dim2; i++)
        stiffmat22(i, j) += -Core::FADUtils::cast_to_double(stiffc2(i, dim1 + j));
    }
#else
    for (unsigned int j = 0; j < dim1; j++)
    {
      for (unsigned int i = 0; i < dim1; i++)
        stiffmat11(i, j) += -Core::FADUtils::cast_to_double(stiffc1_FAD(i, j));
      for (unsigned int i = 0; i < dim2; i++)
        stiffmat21(i, j) += -Core::FADUtils::cast_to_double(stiffc2_FAD(i, j));
    }
    for (unsigned int j = 0; j < dim2; j++)
    {
      for (unsigned int i = 0; i < dim1; i++)
        stiffmat12(i, j) += -Core::FADUtils::cast_to_double(stiffc1_FAD(i, dim1 + j));
      for (unsigned int i = 0; i < dim2; i++)
        stiffmat22(i, j) += -Core::FADUtils::cast_to_double(stiffc2_FAD(i, dim1 + j));
    }
#endif
  }

  return;
}
/*----------------------------------------------------------------------*
 |  end: Evaluate contact stiffness
 *----------------------------------------------------------------------*/

#ifdef ENDPOINTSEGMENTATION
/*------------------------------------------------------------------------------------------*
 |  FAD-based Evaluation of contact stiffness in case of ENDPOINTSEGMENTATION    meier 10/14|
 *------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes,
    numnodalvalues>::evaluate_stiffc_contact_int_seg(Core::LinAlg::SparseMatrix& stiffmatrix,
    const Core::LinAlg::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_xi_bound,
    const Core::LinAlg::Matrix<3, 1, TYPE>& r1, const Core::LinAlg::Matrix<3, 1, TYPE>& r2,
    const Core::LinAlg::Matrix<3, 1, TYPE>& r1_xi, const Core::LinAlg::Matrix<3, 1, TYPE>& r2_xi,
    const Core::LinAlg::Matrix<3, 1, TYPE>& r1_xixi,
    const Core::LinAlg::Matrix<3, 1, TYPE>& r2_xixi,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xi,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xi,
    BeamToBeamContactVariables<numnodes, numnodalvalues>& cpvariables, const double& intfac,
    const double& d_xi_ele_d_xi_bound, TYPE signed_jacobi_interval)
{
#ifndef AUTOMATICDIFF
  FOUR_C_THROW("This method only works with automatic differentiation!");
#endif

  // get dimensions for vectors fc1 and fc2
  const unsigned int dim1 = 3 * numnodes * numnodalvalues;
  const unsigned int dim2 = 3 * numnodes * numnodalvalues;

  // temporary matrices for stiffness and vectors for DOF-GIDs and owning procs
  Core::LinAlg::Matrix<dim1, dim1 + dim2, TYPE> stiffc1_FAD(true);
  Core::LinAlg::Matrix<dim2, dim1 + dim2, TYPE> stiffc2_FAD(true);
  Core::LinAlg::SerialDenseMatrix stiffcontact1(dim1, dim1 + dim2);
  Core::LinAlg::SerialDenseMatrix stiffcontact2(dim2, dim1 + dim2);
  std::vector<int> lmrow1(dim1);
  std::vector<int> lmrow2(dim2);
  std::vector<int> lmrowowner1(dim1);
  std::vector<int> lmrowowner2(dim2);
  std::vector<int> lmcol1(dim1 + dim2);
  std::vector<int> lmcol2(dim1 + dim2);

  // flag indicating assembly
  bool DoNotAssemble = true;
  TYPE gap = cpvariables->GetGap();

  //**********************************************************************
  // evaluate contact stiffness for active pairs
  //**********************************************************************
  if (check_contact_status(Core::FADUtils::cast_to_double(gap)))
  {
    DoNotAssemble = false;

    // node ids of both elements
    //    const int* node_ids1 = Element1()->NodeIds();
    //    const int* node_ids2 = Element2()->NodeIds();
    //
    // TODO: Introduce this quantities as class variables?
    //********************************************************************
    // prepare assembly
    //********************************************************************
    //    // fill lmrow1 and lmrowowner1
    //    for (unsigned int i=0;i<numnodes;++i)
    //    {
    //      // get pointer and dof ids
    //      Core::Nodes::Node* node = ContactDiscret().gNode(node_ids1[i]);
    //      std::vector<int> NodeDofGIDs =  get_global_dofs(node);
    //
    //      for (unsigned int j=0;j<3*numnodalvalues;++j)
    //      {
    //        lmrow1[3*numnodalvalues*i+j]=NodeDofGIDs[j];
    //        lmrowowner1[3*numnodalvalues*i+j]=node->Owner();
    //      }
    //    }
    //
    //    // fill lmrow2 and lmrowowner2
    //    for (unsigned int i=0;i<numnodes;++i)
    //    {
    //      // get pointer and node ids
    //      Core::Nodes::Node* node = ContactDiscret().gNode(node_ids2[i]);
    //      std::vector<int> NodeDofGIDs =  get_global_dofs(node);
    //
    //      for (unsigned int j=0;j<3*numnodalvalues;++j)
    //      {
    //        lmrow2[3*numnodalvalues*i+j]=NodeDofGIDs[j];
    //        lmrowowner2[3*numnodalvalues*i+j]=node->Owner();
    //      }
    //    }
    //
    //    // fill lmcol1 and lmcol2
    //    for (unsigned int i=0;i<numnodes;++i)
    //    {
    //      // get pointer and node ids
    //      Core::Nodes::Node* node = ContactDiscret().gNode(node_ids1[i]);
    //      std::vector<int> NodeDofGIDs =  get_global_dofs(node);
    //
    //      for (unsigned int j=0;j<3*numnodalvalues;++j)
    //      {
    //        lmcol1[3*numnodalvalues*i+j] = NodeDofGIDs[j];
    //        lmcol2[3*numnodalvalues*i+j] = NodeDofGIDs[j];
    //      }
    //    }
    //
    //    // fill lmcol1 and lmcol2
    //    for (unsigned int i=0;i<numnodes;++i)
    //    {
    //      // get pointer and node ids
    //      Core::Nodes::Node* node = ContactDiscret().gNode(node_ids2[i]);
    //      std::vector<int> NodeDofGIDs =  get_global_dofs(node);
    //
    //      for (unsigned int j=0;j<3*numnodalvalues;++j)
    //      {
    //        lmcol1[3*numnodalvalues*numnodes+3*numnodalvalues*i+j] = NodeDofGIDs[j];
    //        lmcol2[3*numnodalvalues*numnodes+3*numnodalvalues*i+j] = NodeDofGIDs[j];
    //      }
    //    }

    // initialize storage for linearizations
    Core::LinAlg::Matrix<dim1 + dim2, 1, TYPE> delta_eta(true);
    Core::LinAlg::Matrix<3, 1, TYPE> delta_r = Core::FADUtils::DiffVector(r1, r2);

    compute_lin_eta_fix_xi(delta_eta, delta_r, r2_xi, r2_xixi, N1, N2, N2_xi);

    Core::LinAlg::Matrix<dim1, 1, TYPE> fc1_FAD(true);
    Core::LinAlg::Matrix<dim2, 1, TYPE> fc2_FAD(true);
    evaluate_fc_contact(nullptr, r1, r2, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi,
        cpvariables, intfac, false, true, false, false, &fc1_FAD, &fc2_FAD);

#ifdef AUTOMATICDIFF
    TYPE fac1(0.0);
    TYPE fac2(0.0);
    fac2 = -Core::FADUtils::ScalarProduct(r2_xi, r2_xi) +
           Core::FADUtils::ScalarProduct(delta_r, r2_xixi);
    fac1 = Core::FADUtils::ScalarProduct(r2_xi, r1_xi);
    for (unsigned int j = 0; j < dim1 + dim2; j++)
    {
      for (unsigned int i = 0; i < dim1; i++)
        stiffc1_FAD(i, j) =
            fc1_FAD(i).dx(j) +
            fc1_FAD(i).dx(dim1 + dim2) * d_xi_ele_d_xi_bound * delta_xi_bound(j) +
            fc1_FAD(i).dx(dim1 + dim2 + 1) *
                (delta_eta(j) - fac1 / fac2 * d_xi_ele_d_xi_bound * delta_xi_bound(j)) +
            fc1_FAD(i).val() / (2.0 * signed_jacobi_interval) * delta_xi_bound(j);
      // d(f)/d(disp)    +d(f)/d(xi,GP)           *d(xi,GP)/d(disp) +d(f)/d(eta,GP)
      // *d(eta,GP)/d(disp) +d(f)/d(xi,Bound)*d(xi,Bound)/d(disp)

      for (unsigned int i = 0; i < dim2; i++)
        stiffc2_FAD(i, j) =
            fc2_FAD(i).dx(j) +
            fc2_FAD(i).dx(dim1 + dim2) * d_xi_ele_d_xi_bound * delta_xi_bound(j) +
            fc2_FAD(i).dx(dim1 + dim2 + 1) *
                (delta_eta(j) - fac1 / fac2 * d_xi_ele_d_xi_bound * delta_xi_bound(j)) +
            fc2_FAD(i).val() / (2.0 * signed_jacobi_interval) * delta_xi_bound(j);
    }
#endif

  }  // if (check_contact_status(gap))

  //**********************************************************************
  // assemble contact stiffness
  //**********************************************************************
  // change sign of stiffc1 and stiffc2 due to time integration.
  // according to analytical derivation there is no minus sign, but for
  // our time integration methods the negative stiffness must be assembled.

  // now finally assemble stiffc1 and stiffc2
  if (!DoNotAssemble)
  {
#ifndef AUTOMATICDIFF
    FOUR_C_THROW("This method only works with AUTOMATICDIFF");
#else
    for (unsigned int j = 0; j < dim1 + dim2; j++)
    {
      for (unsigned int i = 0; i < dim1; i++)
        stiffcontact1(i, j) = -Core::FADUtils::cast_to_double(stiffc1_FAD(i, j));
      for (unsigned int i = 0; i < dim2; i++)
        stiffcontact2(i, j) = -Core::FADUtils::cast_to_double(stiffc2_FAD(i, j));
    }
#endif

    stiffmatrix.Assemble(0, stiffcontact1, lmrow1, lmrowowner1, lmcol1);
    stiffmatrix.Assemble(0, stiffcontact2, lmrow2, lmrowowner2, lmcol2);
  }

  return;
}
/*------------------------------------------------------------------------------------------*
 |  end: FAD-based Evaluation of contact stiffness in case of ENDPOINTSEGMENTATION
 *------------------------------------------------------------------------------------------*/
#endif

/*----------------------------------------------------------------------*
 |  Linearizations of contact point                          meier 10/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::compute_lin_xi_and_lin_eta(
    Core::LinAlg::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_xi,
    Core::LinAlg::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_eta,
    const Core::LinAlg::Matrix<3, 1, TYPE>& delta_r, const Core::LinAlg::Matrix<3, 1, TYPE>& r1_xi,
    const Core::LinAlg::Matrix<3, 1, TYPE>& r2_xi, const Core::LinAlg::Matrix<3, 1, TYPE>& r1_xixi,
    const Core::LinAlg::Matrix<3, 1, TYPE>& r2_xixi,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xi,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xi)
{
  //**********************************************************************
  // we have to solve the following system of equations:
  //  _              _       _      _       _              _      _       _
  // | L(1,1)  L(1,2) |    | Lin_Xi  |    |  B(1,1)  B(1,2) |   | Lin_d1 |
  // |                | *  |         | =  |                 | * |        |
  // |_L(2,1)  L(2,2)_|    |_Lin_Eta_|    |_B(2,1)  B(2,2)_ |   |_Lin_d2_|
  //
  // this can be done easily because it is a linear 2x2-system.
  // we obtain the solution by inverting matrix L:
  //
  // [Lin_Xi; Lin_Eta] = L^-1 * B * [Lin_d1; Lin_d2] = D * [Lin_d1; Lin_d2]
  //
  //**********************************************************************

  const unsigned int dim1 = 3 * numnodes * numnodalvalues;
  const unsigned int dim2 = 3 * numnodes * numnodalvalues;

  // matrices to compute Lin_Xi and Lin_Eta
  Core::LinAlg::Matrix<2, 2, TYPE> L(true);
  Core::LinAlg::Matrix<2, 2, TYPE> L_inv(true);
  Core::LinAlg::Matrix<2, dim1 + dim2, TYPE> B(true);
  Core::LinAlg::Matrix<2, dim1 + dim2, TYPE> D(true);

  // compute L elementwise
  L(0, 0) = Core::FADUtils::scalar_product(r1_xi, r1_xi) +
            Core::FADUtils::scalar_product(delta_r, r1_xixi);
  L(1, 1) = -Core::FADUtils::scalar_product(r2_xi, r2_xi) +
            Core::FADUtils::scalar_product(delta_r, r2_xixi);
  L(0, 1) = -Core::FADUtils::scalar_product(r2_xi, r1_xi);
  L(1, 0) = -L(0, 1);

  // invert L by hand
  TYPE det_L = L(0, 0) * L(1, 1) - L(0, 1) * L(1, 0);
  if (Core::FADUtils::cast_to_double(Core::FADUtils::norm(det_L)) < DETERMINANTTOL)
    FOUR_C_THROW("ERROR: determinant of L = 0");
  L_inv(0, 0) = L(1, 1) / det_L;
  L_inv(0, 1) = -L(0, 1) / det_L;
  L_inv(1, 0) = -L(1, 0) / det_L;
  L_inv(1, 1) = L(0, 0) / det_L;

  for (unsigned int i = 0; i < 3; i++)
  {
    for (unsigned int j = 0; j < dim1; j++)
    {
      B(0, j) += -delta_r(i) * N1_xi(i, j) - r1_xi(i) * N1(i, j);
      B(1, j) += -r2_xi(i) * N1(i, j);
    }
  }

  for (unsigned int i = 0; i < 3; i++)
  {
    for (unsigned int j = 0; j < dim2; j++)
    {
      B(0, j + dim1) += r1_xi(i) * N2(i, j);
      B(1, j + dim1) += -delta_r(i) * N2_xi(i, j) + r2_xi(i) * N2(i, j);
    }
  }

  // compute D = L^-1 * B
  D.multiply(L_inv, B);

  // finally the linearizations / directional derivatives
  for (unsigned int i = 0; i < dim1 + dim2; i++)
  {
    delta_xi(i) = D(0, i);
    delta_eta(i) = D(1, i);
  }

  return;
}
/*----------------------------------------------------------------------*
 |  end: Linearizations of contact point
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Lin. of contact point coordinate eta with fixed xi        meier 10/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::compute_lin_eta_fix_xi(
    Core::LinAlg::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_eta,
    const Core::LinAlg::Matrix<3, 1, TYPE>& delta_r, const Core::LinAlg::Matrix<3, 1, TYPE>& r2_xi,
    const Core::LinAlg::Matrix<3, 1, TYPE>& r2_xixi,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xi)
{
  const unsigned int dim1 = 3 * numnodes * numnodalvalues;
  const unsigned int dim2 = 3 * numnodes * numnodalvalues;

  // matrices to compute Lin_Xi and Lin_Eta
  TYPE L = 0.0;
  Core::LinAlg::Matrix<1, dim1 + dim2, TYPE> B(true);

  // compute L elementwise
  L = -Core::FADUtils::scalar_product(r2_xi, r2_xi) +
      Core::FADUtils::scalar_product(delta_r, r2_xixi);

  //  std::cout << "r2_xi: " << r2_xi << std::endl;
  //  std::cout << "r2_xixi: " << r2_xixi << std::endl;

  if (fabs(Core::FADUtils::cast_to_double(L)) < COLINEARTOL)
    FOUR_C_THROW("Linearization of point to line projection is zero, choose tighter search boxes!");

  for (unsigned int i = 0; i < 3; i++)
  {
    for (unsigned int j = 0; j < dim1; j++)
    {
      B(0, j) += -r2_xi(i) * N1(i, j);
    }
  }

  for (unsigned int i = 0; i < 3; i++)
  {
    for (unsigned int j = 0; j < dim2; j++)
    {
      B(0, j + dim1) += -delta_r(i) * N2_xi(i, j) + r2_xi(i) * N2(i, j);
    }
  }

  // finally the linearizations / directional derivatives
  for (unsigned int i = 0; i < dim1 + dim2; i++)
  {
    delta_eta(i) = 1.0 / L * B(0, i);
  }

  return;
}
/*----------------------------------------------------------------------*
 |  end: Lin. of contact point coordinate eta with fixed xi
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Lin. of contact point coordinate xi with fixed eta        meier 12/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::compute_lin_xi_fix_eta(
    Core::LinAlg::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_xi,
    const Core::LinAlg::Matrix<3, 1, TYPE>& delta_r, const Core::LinAlg::Matrix<3, 1, TYPE>& r1_xi,
    const Core::LinAlg::Matrix<3, 1, TYPE>& r1_xixi,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xi)
{
  const unsigned int dim1 = 3 * numnodes * numnodalvalues;
  const unsigned int dim2 = 3 * numnodes * numnodalvalues;

  // matrices to compute Lin_Xi and Lin_Eta
  TYPE L = 0.0;
  Core::LinAlg::Matrix<1, dim1 + dim2, TYPE> B(true);

  // compute L elementwise
  L = Core::FADUtils::scalar_product(r1_xi, r1_xi) +
      Core::FADUtils::scalar_product(delta_r, r1_xixi);

  if (fabs(Core::FADUtils::cast_to_double(L)) < COLINEARTOL)
    FOUR_C_THROW("Linearization of point to line projection is zero, choose tighter search boxes!");

  for (unsigned int i = 0; i < 3; i++)
  {
    for (unsigned int j = 0; j < dim1; j++)
    {
      B(0, j) += -delta_r(i) * N1_xi(i, j) - r1_xi(i) * N1(i, j);
    }
  }

  for (unsigned int i = 0; i < 3; i++)
  {
    for (unsigned int j = 0; j < dim2; j++)
    {
      B(0, j + dim1) += r1_xi(i) * N2(i, j);
    }
  }

  // finally the linearizations / directional derivatives
  for (unsigned int i = 0; i < dim1 + dim2; i++)
  {
    delta_xi(i) = 1.0 / L * B(0, i);
  }

  return;
}
/*----------------------------------------------------------------------*
 |  end: Lin. of contact point coordinate xi with fixed eta
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Compute linearization of integration interval bounds      meier 10/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::compute_lin_xi_bound(
    Core::LinAlg::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_xi_bound,
    TYPE& eta1_bound, TYPE eta2)
{
  // vectors for shape functions and their derivatives
  Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N1(true);       // = N1
  Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N2(true);       // = N2
  Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N1_xi(true);    // = N1,xi
  Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N2_xi(true);    // = N2,eta
  Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N1_xixi(true);  // = N1,xixi
  Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N2_xixi(true);  // = N2,etaeta

  // coords and derivatives of the two contacting points
  Core::LinAlg::Matrix<3, 1, TYPE> r1(true);       // = r1
  Core::LinAlg::Matrix<3, 1, TYPE> r2(true);       // = r2
  Core::LinAlg::Matrix<3, 1, TYPE> r1_xi(true);    // = r1,xi
  Core::LinAlg::Matrix<3, 1, TYPE> r2_xi(true);    // = r2,eta
  Core::LinAlg::Matrix<3, 1, TYPE> r1_xixi(true);  // = r1,xixi
  Core::LinAlg::Matrix<3, 1, TYPE> r2_xixi(true);  // = r2,etaeta
  Core::LinAlg::Matrix<3, 1, TYPE> delta_r(true);  // = r1-r2

  // update shape functions and their derivatives
  get_shape_functions(N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi, eta1_bound, eta2);
  // update coordinates and derivatives of contact points
  compute_coords_and_derivs(
      r1, r2, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi);

  delta_r = Core::FADUtils::diff_vector(r1, r2);

  const unsigned int dim1 = 3 * numnodes * numnodalvalues;
  const unsigned int dim2 = 3 * numnodes * numnodalvalues;

  // matrices to compute Lin_Xi and Lin_Eta
  TYPE a_11(0.0);
  Core::LinAlg::Matrix<2, dim1 + dim2, TYPE> B(true);

  a_11 = Core::FADUtils::scalar_product(r1_xi, r1_xi) +
         Core::FADUtils::scalar_product(delta_r, r1_xixi);

#ifdef CHANGEENDPOINTPROJECTION
  TYPE a_21(0.0);
  a_21 = Core::FADUtils::ScalarProduct(r1_xi, r2_xi);
#endif


  for (unsigned int i = 0; i < 3; i++)
  {
    for (unsigned int j = 0; j < dim1; j++)
    {
      B(0, j) += -delta_r(i) * N1_xi(i, j) - r1_xi(i) * N1(i, j);
      B(1, j) += -r2_xi(i) * N1(i, j);
    }
  }

  for (unsigned int i = 0; i < 3; i++)
  {
    for (unsigned int j = 0; j < dim2; j++)
    {
      B(0, j + dim1) += r1_xi(i) * N2(i, j);
      B(1, j + dim1) += -delta_r(i) * N2_xi(i, j) + r2_xi(i) * N2(i, j);
    }
  }

#ifndef CHANGEENDPOINTPROJECTION
  // finally the linearizations / directional derivatives in case the orthogonality condition is
  // fulfilled on beam1
  for (unsigned int i = 0; i < dim1 + dim2; i++)
  {
    delta_xi_bound(i) = B(0, i) / a_11;
  }
#else
  // finally the linearizations / directional derivatives in case the orthogonality condition is
  // fulfilled on beam2
  for (unsigned int i = 0; i < dim1 + dim2; i++)
  {
    delta_xi_bound(i) = B(1, i) / a_21;
  }
#endif

  return;
}
/*----------------------------------------------------------------------*
 |  end: Compute linearization of integration interval bounds
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Compute linearization of gap                              meier 10/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::compute_lin_gap(
    Core::LinAlg::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_gap,
    const Core::LinAlg::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_xi,
    const Core::LinAlg::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_eta,
    const Core::LinAlg::Matrix<3, 1, TYPE>& delta_r, const TYPE& norm_delta_r,
    const Core::LinAlg::Matrix<3, 1, TYPE>& r1_xi, const Core::LinAlg::Matrix<3, 1, TYPE>& r2_xi,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2)
{
  const unsigned int dim1 = 3 * numnodes * numnodalvalues;
  const unsigned int dim2 = 3 * numnodes * numnodalvalues;

  // delta g := delta_r/||delta_r||*auxiliary_matri1 delta d, with auxiliary_matri1 =
  // (r1_xi*delta_xi-r2_xi*delta_eta + (N1, -N2))

  Core::LinAlg::Matrix<3, dim1 + dim2, TYPE> auxiliary_matrix1(true);

  for (unsigned int i = 0; i < 3; i++)
  {
    for (unsigned int j = 0; j < dim1 + dim2; j++)
    {
      auxiliary_matrix1(i, j) += r1_xi(i) * delta_xi(j) - r2_xi(i) * delta_eta(j);
    }
  }

  for (unsigned int i = 0; i < 3; i++)
  {
    for (unsigned int j = 0; j < dim1; j++)
    {
      auxiliary_matrix1(i, j) += N1(i, j);
    }
  }

  for (unsigned int i = 0; i < 3; i++)
  {
    for (unsigned int j = 0; j < dim2; j++)
    {
      auxiliary_matrix1(i, j + dim1) += -N2(i, j);
    }
  }

  // compute linearization of gap
  for (unsigned int i = 0; i < 3; i++)
  {
    for (unsigned int j = 0; j < dim1 + dim2; j++)
    {
      delta_gap(j) += delta_r(i) * auxiliary_matrix1(i, j) / norm_delta_r;
    }
  }

  return;
}
/*----------------------------------------------------------------------*
 | end: Compute linearization of gap
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Compute linearization of cosine of contact angle          meier 10/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::
    compute_lin_cos_contact_angle(
        Core::LinAlg::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_coscontactangle,
        Core::LinAlg::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_xi,
        Core::LinAlg::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_eta,
        const Core::LinAlg::Matrix<3, 1, TYPE>& r1_xi,
        const Core::LinAlg::Matrix<3, 1, TYPE>& r2_xi,
        const Core::LinAlg::Matrix<3, 1, TYPE>& r1_xixi,
        const Core::LinAlg::Matrix<3, 1, TYPE>& r2_xixi,
        const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xi,
        const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xi)
{
  const unsigned int dim1 = 3 * numnodes * numnodalvalues;
  const unsigned int dim2 = 3 * numnodes * numnodalvalues;

  TYPE norm_r1xi = Core::FADUtils::vector_norm<3>(r1_xi);
  TYPE norm_r2xi = Core::FADUtils::vector_norm<3>(r2_xi);
  Core::LinAlg::Matrix<3, 1, TYPE> r1_xi_unit(r1_xi);
  Core::LinAlg::Matrix<3, 1, TYPE> r2_xi_unit(r2_xi);
  r1_xi_unit.scale(1.0 / norm_r1xi);
  r2_xi_unit.scale(1.0 / norm_r2xi);
  TYPE r1xi_unit_r2xi_unit = Core::FADUtils::scalar_product(r1_xi_unit, r2_xi_unit);

  // Pre-factor representing the modulus, since s=|r1xi_unit_r2xi_unit|
  double modulus_factor = 1.0;

  if (r1xi_unit_r2xi_unit < 0.0) modulus_factor = -1.0;

  Core::LinAlg::Matrix<3, 1, TYPE> v1(r2_xi_unit);
  Core::LinAlg::Matrix<3, 1, TYPE> v2(r1_xi_unit);
  v1.update(-r1xi_unit_r2xi_unit, r1_xi_unit, 1.0);
  v2.update(-r1xi_unit_r2xi_unit, r2_xi_unit, 1.0);
  v1.scale(1.0 / norm_r1xi);
  v2.scale(1.0 / norm_r2xi);

  Core::LinAlg::Matrix<3, dim1 + dim2, TYPE> delta_r1_xi(true);
  Core::LinAlg::Matrix<3, dim1 + dim2, TYPE> delta_r2_xi(true);

  for (unsigned int i = 0; i < 3; i++)
  {
    for (unsigned int j = 0; j < dim1 + dim2; j++)
    {
      delta_r1_xi(i, j) = r1_xixi(i) * delta_xi(j);

      if (j < dim1) delta_r1_xi(i, j) += N1_xi(i, j);
    }
  }

  for (unsigned int i = 0; i < 3; i++)
  {
    for (unsigned int j = 0; j < dim1 + dim2; j++)
    {
      delta_r2_xi(i, j) = r2_xixi(i) * delta_eta(j);

      if (j >= dim1) delta_r2_xi(i, j) += N2_xi(i, j - dim1);
    }
  }

  Core::LinAlg::Matrix<1, dim1 + dim2, TYPE> v1_delta_r1_xi(true);
  Core::LinAlg::Matrix<1, dim1 + dim2, TYPE> v2_delta_r2_xi(true);
  v1_delta_r1_xi.multiply_tn(v1, delta_r1_xi);
  v2_delta_r2_xi.multiply_tn(v2, delta_r2_xi);

  for (unsigned int j = 0; j < dim1 + dim2; j++)
  {
    delta_coscontactangle(j) = modulus_factor * (v1_delta_r1_xi(j) + v2_delta_r2_xi(j));
  }

  return;
}
/*----------------------------------------------------------------------*
 | end: Compute linearization of cosine of contact angle
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Compute linearization of normal vector                    meier 10/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::compute_lin_normal(
    Core::LinAlg::Matrix<3, 2 * 3 * numnodes * numnodalvalues, TYPE>& delta_normal,
    const Core::LinAlg::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_xi,
    const Core::LinAlg::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_eta,
    const Core::LinAlg::Matrix<3, 1, TYPE>& delta_r, const Core::LinAlg::Matrix<3, 1, TYPE>& r1_xi,
    const Core::LinAlg::Matrix<3, 1, TYPE>& r2_xi,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2)
{
  const unsigned int dim1 = 3 * numnodes * numnodalvalues;
  const unsigned int dim2 = 3 * numnodes * numnodalvalues;

  // delta n := auxiliary_matri2*auxiliary_matrix1* delta d, with auxiliary_matri2 =
  // (I-nxn)/||r1-r2|| and auxiliary_matri1 = (r1_xi*delta_xi-r2_xi*delta_eta + (N1, -N2))

  TYPE norm_delta_r = Core::FADUtils::vector_norm<3>(delta_r);
  Core::LinAlg::Matrix<3, 1, TYPE> normal(delta_r);
  normal.scale(1.0 / norm_delta_r);

  Core::LinAlg::Matrix<3, dim1 + dim2, TYPE> auxiliary_matrix1(true);
  Core::LinAlg::Matrix<3, 3, TYPE> auxiliary_matrix2(true);

  // compute auxiliary_matrix1
  for (unsigned int i = 0; i < 3; i++)
  {
    for (unsigned int j = 0; j < dim1 + dim2; j++)
    {
      auxiliary_matrix1(i, j) += r1_xi(i) * delta_xi(j) - r2_xi(i) * delta_eta(j);
    }
  }

  for (unsigned int i = 0; i < 3; i++)
  {
    for (unsigned int j = 0; j < dim1; j++)
    {
      auxiliary_matrix1(i, j) += N1(i, j);
    }
  }

  for (unsigned int i = 0; i < 3; i++)
  {
    for (unsigned int j = 0; j < dim2; j++)
    {
      auxiliary_matrix1(i, j + dim1) += -N2(i, j);
    }
  }

  // compute auxiliary_matrix2
  for (unsigned int i = 0; i < 3; i++)
  {
    auxiliary_matrix2(i, i) += 1.0 / norm_delta_r;
    for (unsigned int j = 0; j < 3; j++)
    {
      auxiliary_matrix2(i, j) += -normal(i) * normal(j) / norm_delta_r;
    }
  }

  // compute linearization of normal vector
  for (unsigned int i = 0; i < 3; i++)
    for (unsigned int j = 0; j < 3; j++)
      for (unsigned int k = 0; k < dim1 + dim2; k++)
        delta_normal(i, k) += auxiliary_matrix2(i, j) * auxiliary_matrix1(j, k);

  return;
}
/*----------------------------------------------------------------------*
 | end: Compute linearization of normal vector
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  evaluate shape functions and derivatives                 meier 01/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::get_shape_functions(
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xi,
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xi,
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xixi,
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xixi, const TYPE& eta1,
    const TYPE& eta2)
{
  // get both discretization types
  const Core::FE::CellType distype1 = element1()->shape();
  const Core::FE::CellType distype2 = element2()->shape();

  Core::LinAlg::Matrix<1, numnodes * numnodalvalues, TYPE> N1_i(true);
  Core::LinAlg::Matrix<1, numnodes * numnodalvalues, TYPE> N1_i_xi(true);
  Core::LinAlg::Matrix<1, numnodes * numnodalvalues, TYPE> N1_i_xixi(true);
  Core::LinAlg::Matrix<1, numnodes * numnodalvalues, TYPE> N2_i(true);
  Core::LinAlg::Matrix<1, numnodes * numnodalvalues, TYPE> N2_i_xi(true);
  Core::LinAlg::Matrix<1, numnodes * numnodalvalues, TYPE> N2_i_xixi(true);

  if (numnodalvalues == 1)
  {
    // get values and derivatives of shape functions
    Core::FE::shape_function_1d(N1_i, eta1, distype1);
    Core::FE::shape_function_1d(N2_i, eta2, distype2);
    Core::FE::shape_function_1d_deriv1(N1_i_xi, eta1, distype1);
    Core::FE::shape_function_1d_deriv1(N2_i_xi, eta2, distype2);
    Core::FE::shape_function_1d_deriv2(N1_i_xixi, eta1, distype1);
    Core::FE::shape_function_1d_deriv2(N2_i_xixi, eta2, distype2);
  }
  else if (numnodalvalues == 2)
  {
    // TODO maybe cast class variables to Beam3Base upon construction ?!
    double length1 = (dynamic_cast<const Discret::Elements::Beam3Base*>(element1()))->ref_length();
    double length2 = (dynamic_cast<const Discret::Elements::Beam3Base*>(element2()))->ref_length();

    /* TODO hard set distype to line2 in case of numnodalvalues_=2 because
     *  only 3rd order Hermite interpolation is used (always 2 nodes) */
    const Core::FE::CellType distype1herm = Core::FE::CellType::line2;
    const Core::FE::CellType distype2herm = Core::FE::CellType::line2;

    // get values and derivatives of shape functions
    Core::FE::shape_function_hermite_1d(N1_i, eta1, length1, distype1herm);
    Core::FE::shape_function_hermite_1d(N2_i, eta2, length2, distype2herm);
    Core::FE::shape_function_hermite_1d_deriv1(N1_i_xi, eta1, length1, distype1herm);
    Core::FE::shape_function_hermite_1d_deriv1(N2_i_xi, eta2, length2, distype2herm);
    Core::FE::shape_function_hermite_1d_deriv2(N1_i_xixi, eta1, length1, distype1herm);
    Core::FE::shape_function_hermite_1d_deriv2(N2_i_xixi, eta2, length2, distype2herm);
  }
  else
    FOUR_C_THROW(
        "Only beam elements with one (nodal positions) or two (nodal positions + nodal tangents) "
        "values are valid!");

  // Assemble the individual shape functions in matrices, such that: r1=N1*d1, r1_xi=N1_xi*d1,
  // r1_xixi=N1_xixi*d1, r2=N2*d2, r2_xi=N2_xi*d2, r2_xixi=N2_xixi*d2
  Discret::Utils::Beam::assemble_shape_functions_and_derivs_and2nd_derivs<numnodes, numnodalvalues,
      TYPE>(N1_i, N1_i_xi, N1_i_xixi, N1, N1_xi, N1_xixi);
  Discret::Utils::Beam::assemble_shape_functions_and_derivs_and2nd_derivs<numnodes, numnodalvalues,
      TYPE>(N2_i, N2_i_xi, N2_i_xixi, N2, N2_xi, N2_xixi);

  return;
}
/*----------------------------------------------------------------------*
 |  end: evaluate shape functions and derivatives
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  evaluate shape functions and derivatives                 meier 10/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::get_shape_functions(
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N, const TYPE& eta, int deriv,
    const Core::Elements::Element* ele) const
{
  // get both discretization types
  const Core::FE::CellType distype = ele->shape();
  Core::LinAlg::Matrix<1, numnodes * numnodalvalues, TYPE> N_i(true);

  if (numnodalvalues == 1)
  {
    // get values and derivatives of shape functions
    switch (deriv)
    {
      case 0:
      {
        Core::FE::shape_function_1d(N_i, eta, distype);
        break;
      }
      case 1:
      {
        Core::FE::shape_function_1d_deriv1(N_i, eta, distype);
        break;
      }
      case 2:
      {
        Core::FE::shape_function_1d_deriv2(N_i, eta, distype);
        break;
      }
    }
  }
  else if (numnodalvalues == 2)
  {
    double length = (dynamic_cast<const Discret::Elements::Beam3Base*>(ele))->ref_length();

    /* TODO hard set distype to line2 in case of numnodalvalues_=2 because
     *  only 3rd order Hermite interpolation is used (always 2 nodes) */
    const Core::FE::CellType distypeherm = Core::FE::CellType::line2;

    // get values and derivatives of shape functions
    switch (deriv)
    {
      case 0:
      {
        Core::FE::shape_function_hermite_1d(N_i, eta, length, distypeherm);
        break;
      }
      case 1:
      {
        Core::FE::shape_function_hermite_1d_deriv1(N_i, eta, length, distypeherm);
        break;
      }
      case 2:
      {
        Core::FE::shape_function_hermite_1d_deriv2(N_i, eta, length, distypeherm);
        break;
      }
    }
  }
  else
    FOUR_C_THROW(
        "Only beam elements with one (nodal positions) or two (nodal positions + nodal tangents) "
        "values are valid!");

  // Assemble the individual shape functions in matrices, such that: r1=N1*d1, r1_xi=N1_xi*d1,
  // r1_xixi=N1_xixi*d1, r2=N2*d2, r2_xi=N2_xi*d2, r2_xixi=N2_xixi*d2
  Discret::Utils::Beam::assemble_shape_functions<numnodes, numnodalvalues, TYPE>(N_i, N);

  return;
}
/*----------------------------------------------------------------------*
 |  end: evaluate shape functions and derivatives
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | compute position at given curve point                  meier 10/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
Core::LinAlg::Matrix<3, 1, TYPE>
BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::r(
    const TYPE& eta, const Core::Elements::Element* ele) const
{
  Core::LinAlg::Matrix<3, 1, TYPE> r(true);
  Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N(true);
  get_shape_functions(N, eta, 0, ele);

  if (ele->id() == element1()->id())
  {
    // compute output variable
    for (unsigned int i = 0; i < 3; i++)
    {
      for (unsigned int j = 0; j < 3 * numnodes * numnodalvalues; j++)
      {
        r(i) += N(i, j) * ele1pos_(j);
      }
    }
  }
  else if (ele->id() == element2()->id())
  {
    // compute output variable
    for (unsigned int i = 0; i < 3; i++)
    {
      for (unsigned int j = 0; j < 3 * numnodes * numnodalvalues; j++)
      {
        r(i) += N(i, j) * ele2pos_(j);
      }
    }
  }
  else
    FOUR_C_THROW("This method can only applied to element1_ and element2_!");

  return r;
}
/*----------------------------------------------------------------------*
 | end: compute position at given curve point              meier 02/14|
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | compute tangent at given curve point                  meier 10/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
Core::LinAlg::Matrix<3, 1, TYPE>
BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::r_xi(
    const TYPE& eta, const Core::Elements::Element* ele)
{
  Core::LinAlg::Matrix<3, 1, TYPE> r_xi(true);
  Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE> N_xi(true);
  get_shape_functions(N_xi, eta, 1, ele);

  if (ele->id() == element1()->id())
  {
    // compute output variable
    for (unsigned int i = 0; i < 3; i++)
    {
      for (unsigned int j = 0; j < 3 * numnodes * numnodalvalues; j++)
      {
        r_xi(i) += N_xi(i, j) * ele1pos_(j);
      }
    }
  }
  else if (ele->id() == element2()->id())
  {
    // compute output variable
    for (unsigned int i = 0; i < 3; i++)
    {
      for (unsigned int j = 0; j < 3 * numnodes * numnodalvalues; j++)
      {
        r_xi(i) += N_xi(i, j) * ele2pos_(j);
      }
    }
  }
  else
    FOUR_C_THROW("This method can only applied to element1_ and element2_!");

  return r_xi;
}
/*----------------------------------------------------------------------*
 | end: compute tangent at given curve point              meier 02/14|
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | compute contact point coordinates and their derivatives   meier 02/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::compute_coords_and_derivs(
    Core::LinAlg::Matrix<3, 1, TYPE>& r1, Core::LinAlg::Matrix<3, 1, TYPE>& r2,
    Core::LinAlg::Matrix<3, 1, TYPE>& r1_xi, Core::LinAlg::Matrix<3, 1, TYPE>& r2_xi,
    Core::LinAlg::Matrix<3, 1, TYPE>& r1_xixi, Core::LinAlg::Matrix<3, 1, TYPE>& r2_xixi,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xi,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xi,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xixi,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xixi)
{
  r1.clear();
  r2.clear();
  r1_xi.clear();
  r2_xi.clear();
  r1_xixi.clear();
  r2_xixi.clear();

#ifdef AUTOMATICDIFF
  BeamContact::SetFADDispDofs<numnodes, numnodalvalues>(ele1pos_, ele2pos_);
#endif

  // compute output variable
  for (unsigned int i = 0; i < 3; i++)
  {
    for (unsigned int j = 0; j < 3 * numnodes * numnodalvalues; j++)
    {
      r1(i) += N1(i, j) * ele1pos_(j);
      r2(i) += N2(i, j) * ele2pos_(j);
      r1_xi(i) += N1_xi(i, j) * ele1pos_(j);
      r2_xi(i) += N2_xi(i, j) * ele2pos_(j);
      r1_xixi(i) += N1_xixi(i, j) * ele1pos_(j);
      r2_xixi(i) += N2_xixi(i, j) * ele2pos_(j);
    }
  }

  return;
}
/*----------------------------------------------------------------------*
 | end: compute contact point coordinates and their derivatives         |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Evaluate function f in CPP                               meier 02/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes,
    numnodalvalues>::evaluate_orthogonality_condition(Core::LinAlg::Matrix<2, 1, TYPE>& f,
    const Core::LinAlg::Matrix<3, 1, TYPE>& delta_r, const double norm_delta_r,
    const Core::LinAlg::Matrix<3, 1, TYPE>& r1_xi, const Core::LinAlg::Matrix<3, 1, TYPE>& r2_xi,
    const Core::LinAlg::Matrix<3, 1, TYPE>& t1, const Core::LinAlg::Matrix<3, 1, TYPE>& t2)
{
  // reset f
  f.clear();

  for (unsigned int i = 0; i < 3; i++)
  {
    f(0) += delta_r(i) * r1_xi(i) / norm_delta_r;
    f(1) += -delta_r(i) * r2_xi(i) / norm_delta_r;
  }

  return;
}
/*----------------------------------------------------------------------*
 |  end: Evaluate function f in CPP
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Evaluate Jacobian df in CPP                              meier 02/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes,
    numnodalvalues>::evaluate_lin_orthogonality_condition(Core::LinAlg::Matrix<2, 2, TYPE>& df,
    Core::LinAlg::Matrix<2, 2, TYPE>& dfinv, const Core::LinAlg::Matrix<3, 1, TYPE>& delta_r,
    const double norm_delta_r, const Core::LinAlg::Matrix<3, 1, TYPE>& r1_xi,
    const Core::LinAlg::Matrix<3, 1, TYPE>& r2_xi, const Core::LinAlg::Matrix<3, 1, TYPE>& r1_xixi,
    const Core::LinAlg::Matrix<3, 1, TYPE>& r2_xixi, const Core::LinAlg::Matrix<3, 1, TYPE>& t1,
    const Core::LinAlg::Matrix<3, 1, TYPE>& t2, const Core::LinAlg::Matrix<3, 1, TYPE>& t1_xi,
    const Core::LinAlg::Matrix<3, 1, TYPE>& t2_xi, bool& elementscolinear)

{
  // reset df and dfinv
  df.clear();
  dfinv.clear();

  // evaluate df
  // see Wriggers, Computational Contact Mechanics, equation (12.7)
  for (unsigned int i = 0; i < 3; i++)
  {
    df(0, 0) += (r1_xi(i) * r1_xi(i) + delta_r(i) * r1_xixi(i)) / norm_delta_r;
    df(0, 1) += -r1_xi(i) * r2_xi(i) / norm_delta_r;
    df(1, 0) += -r2_xi(i) * r1_xi(i) / norm_delta_r;
    df(1, 1) += (r2_xi(i) * r2_xi(i) - delta_r(i) * r2_xixi(i)) / norm_delta_r;
  }


  // Inverting (2x2) matrix df by hard coded formula, so that it is
  // possible to handle colinear vectors, because they lead to det(df) =0
  TYPE det_df = df(0, 0) * df(1, 1) - df(1, 0) * df(0, 1);

  //********************************************************************
  // ASSUMPTION:
  // If det_df=0 we assume, that the two elements have an identical
  // neutral axis. These contact objects will be rejected. The outcome
  // of this physically rare phenomenon is that handling of line contact
  // is not possible with this approach.
  //********************************************************************

  // singular df
  if (Core::FADUtils::cast_to_double(Core::FADUtils::norm(det_df)) < COLINEARTOL)
  {
    // sort out
    elementscolinear = true;
  }
  // regular df (inversion possible)
  else
  {
    // do not sort out
    elementscolinear = false;

    // invert df
    dfinv(0, 0) = df(1, 1) / det_df;
    dfinv(0, 1) = -df(0, 1) / det_df;
    dfinv(1, 0) = -df(1, 0) / det_df;
    dfinv(1, 1) = df(0, 0) / det_df;
  }

  return;
}
/*----------------------------------------------------------------------*
 |  end: Evaluate Jacobian df in CPP
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Evaluate orthogonality cond. of point to line projection  meier 10/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes,
    numnodalvalues>::evaluate_ptl_orthogonality_condition(TYPE& f,
    const Core::LinAlg::Matrix<3, 1, TYPE>& delta_r, const double norm_delta_r,
    const Core::LinAlg::Matrix<3, 1, TYPE>& r1_xi, const Core::LinAlg::Matrix<3, 1, TYPE>& r2_xi,
    bool orthogonalprojection)
{
  // reset f
  f = 0.0;

  // evaluate f
  if (orthogonalprojection == false)  // standard case
  {
    for (unsigned int i = 0; i < 3; i++)
    {
      f += -delta_r(i) * r2_xi(i) / norm_delta_r;
    }
  }
  else
  {
    for (unsigned int i = 0; i < 3; i++)
    {
      f += -delta_r(i) * r1_xi(i) / norm_delta_r;
    }
  }


  return;
}
/*----------------------------------------------------------------------*
 |  end: Evaluate orthogonality cond. of point to line projection
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Evaluate Jacobian df of PTLOrthogonalityCondition        meier 10/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
bool BeamInteraction::BeamToBeamContactPair<numnodes,
    numnodalvalues>::evaluate_lin_ptl_orthogonality_condition(TYPE& df,
    const Core::LinAlg::Matrix<3, 1, TYPE>& delta_r, const double norm_delta_r,
    const Core::LinAlg::Matrix<3, 1, TYPE>& r1_xi, const Core::LinAlg::Matrix<3, 1, TYPE>& r2_xi,
    const Core::LinAlg::Matrix<3, 1, TYPE>& r2_xixi, bool orthogonalprojection)

{
  // reset df
  df = 0.0;

  // evaluate df
  if (orthogonalprojection == false)  // standard case
  {
    for (unsigned int i = 0; i < 3; i++)
    {
      df += (r2_xi(i) * r2_xi(i) - delta_r(i) * r2_xixi(i)) / norm_delta_r;
    }
  }
  else
  {
    for (unsigned int i = 0; i < 3; i++)
    {
      df += r1_xi(i) * r2_xi(i) / norm_delta_r;
    }
  }

  // check, if df=0: This can happen e.g. when the master beam 2 describes a circle geometry and the
  // projectiong slave point coincides with the cetern of the circle

  if (fabs(Core::FADUtils::cast_to_double(df)) < COLINEARTOL)
    return false;
  else
    return true;
}
/*----------------------------------------------------------------------*
 |  end: Evaluate Jacobian df in CPP
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Compute normal vector in contact point                   meier 02/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::compute_normal(
    Core::LinAlg::Matrix<3, 1, TYPE>& r1, Core::LinAlg::Matrix<3, 1, TYPE>& r2,
    Core::LinAlg::Matrix<3, 1, TYPE>& r1_xi, Core::LinAlg::Matrix<3, 1, TYPE>& r2_xi,
    std::shared_ptr<BeamToBeamContactVariables<numnodes, numnodalvalues>> variables,
    int contacttype)
{
  // compute non-unit normal
  Core::LinAlg::Matrix<3, 1, TYPE> delta_r = Core::FADUtils::diff_vector(r1, r2);

  // compute length of normal
  TYPE norm_delta_r = Core::FADUtils::vector_norm<3>(delta_r);

  if (Core::FADUtils::cast_to_double(norm_delta_r) < NORMTOL)
    FOUR_C_THROW("ERROR: Normal of length zero! --> change time step!");

  // unit normal
  Core::LinAlg::Matrix<3, 1, TYPE> normal(true);
  normal.update(1.0 / norm_delta_r, delta_r, 0.0);

  TYPE gap = norm_delta_r - r1_ - r2_;

  variables->set_gap(gap);
  variables->set_normal(normal);
  variables->set_angle(
      BeamInteraction::calc_angle(Core::FADUtils::cast_to_double<TYPE, 3, 1>(r1_xi),
          Core::FADUtils::cast_to_double<TYPE, 3, 1>(r2_xi)));

  // Fixme
  //  if (Core::FADUtils::cast_to_double(gap)<-MAXPENETRATIONSAFETYFAC*(R1_+R2_) and numstep_>0)
  if (Core::FADUtils::cast_to_double(gap) < -MAXPENETRATIONSAFETYFAC * (r1_ + r2_))
  {
    this->print(std::cout);
    FOUR_C_THROW(
        "Gap too small, danger of penetration. Choose smaller time step or higher penalty!");
  }

  return;
}
/*----------------------------------------------------------------------*
 |  end: Compute normal vector in contact point
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Check if contact is active or inactive                    meier 02/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
bool BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::check_contact_status(
    const double& gap)
{
  // First parameter for contact force regularization
  double g0 = params()->beam_to_beam_contact_params()->beam_to_beam_penalty_law_regularization_g0();
  bool contactflag = false;

  int penaltylaw = params()->beam_to_beam_contact_params()->penalty_law();

  if (penaltylaw == BeamContact::pl_lp)
  {
    // linear penalty force law
    if (gap < 0)
    {
      contactflag = true;
    }
    else
      contactflag = false;
  }
  else if (penaltylaw == BeamContact::pl_qp)
  {
    // quadratic penalty force law
    if (gap < 0)
    {
      contactflag = true;
    }
    else
      contactflag = false;
  }
  else if (penaltylaw == BeamContact::pl_lpqp)
  {
    // penalty laws with regularization for positive gaps
    if (g0 == -1.0) FOUR_C_THROW("Invalid value of regularization parameter BEAMS_PENREGPARAM_G0!");

    // Parameter to shift penalty law
    double gbar = params()->beam_to_beam_contact_params()->gap_shift();
    TYPE g = gap + gbar;

    if (g < g0)
    {
      contactflag = true;
    }
    else
      contactflag = false;
  }
  else if (penaltylaw == BeamContact::pl_lpcp or penaltylaw == BeamContact::pl_lpdqp or
           penaltylaw == BeamContact::pl_lpep)
  {
    // penalty laws with regularization for positive gaps
    if (g0 == -1.0) FOUR_C_THROW("Invalid value of regularization parameter BEAMS_PENREGPARAM_G0!");

    if (gap < g0)
    {
      contactflag = true;
    }
    else
      contactflag = false;
  }
  else if (penaltylaw == BeamContact::pl_lnqp)
  {
    // penalty law with quadratic regularization for negative gaps
    if (gap < 0)
    {
      contactflag = true;
    }
    else
      contactflag = false;
  }

  return contactflag;
}
/*----------------------------------------------------------------------*
 |  end: Check if contact is active or inactive
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::clear_class_variables()
{
  cpvariables_.clear();
  gpvariables_.clear();
  epvariables_.clear();

  //  double kappa_max = timeintparams.get<double>("kappa_max",-1.0);
  //  if(kappa_max<0)
  //    FOUR_C_THROW("Maximal curvature should be a positive value!");

  // TODO
  //  //Check, if maximal curvature bound is exceeded:
  //  double crosssection_to_curvature_ratio = max(R1_,R2_)*kappa_max;
  //  if(crosssection_to_curvature_ratio>MAXCROSSSECTIONTOCURVATURE)
  //    FOUR_C_THROW("Curvature too large. Choose larger value MAXCROSSSECTIONTOCURVATURE and adapt
  //    shifting angles!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::reset_state(
    const std::vector<double>& centerline_dofvec_ele1,
    const std::vector<double>& centerline_dofvec_ele2)
{
  for (unsigned int i = 0; i < 3 * numnodalvalues; i++)
  {
    for (unsigned int j = 0; j < numnodes; j++)
    {
      ele1pos_(3 * numnodalvalues * j + i) = centerline_dofvec_ele1[3 * numnodalvalues * j + i];
      ele2pos_(3 * numnodalvalues * j + i) = centerline_dofvec_ele2[3 * numnodalvalues * j + i];
    }
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
double BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::get_jacobi(
    const Core::Elements::Element* element1)
{
  double jacobi = 1.0;
  const Core::Elements::ElementType& eot1 = element1->element_type();

  // The jacobi factor is only needed in order to scale the CPP condition. Therefore, we only use
  // the jacobi_ factor corresponding to the first gauss point of the beam element
  if (eot1 == Discret::Elements::Beam3ebType::instance())
  {
    jacobi = (static_cast<const Discret::Elements::Beam3eb*>(element1))->get_jacobi();
  }
  else if (eot1 == Discret::Elements::Beam3rType::instance())
  {
    jacobi = (static_cast<const Discret::Elements::Beam3r*>(element1))->get_jacobi();
  }
  else if (eot1 == Discret::Elements::Beam3kType::instance())
  {
    jacobi = (static_cast<const Discret::Elements::Beam3k*>(element1))->get_jacobi();
  }
  else
  {
    std::cout << "      Warning: No valid jacobi weight in CPP supported by applied beam element!!!"
              << std::endl;
  }

  return jacobi;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::print(
    std::ostream& out) const
{
  check_init_setup();

  out << "\nInstance of BeamToBeamContactPair (EleGIDs " << element1()->id() << " & "
      << element2()->id() << "):";
  out << "\nele1 dofvec: " << ele1pos_;
  out << "\nele2 dofvec: " << ele2pos_;

  out << "\nnumber of active closest-point pairs: " << cpvariables_.size();
  out << "\nnumber of active Gauss-point pairs: " << gpvariables_.size();
  out << "\nnumber of active end-point pairs: " << epvariables_.size();

  if (cpvariables_.size())
  {
    out << "\n\nactive closest-point pairs:";
    for (unsigned int i = 0; i < cpvariables_.size(); ++i) cpvariables_[i]->print(out);
  }
  if (gpvariables_.size())
  {
    out << "\n\nactive gauss-point pairs:";
    for (unsigned int i = 0; i < gpvariables_.size(); ++i) gpvariables_[i]->print(out);
  }
  if (epvariables_.size())
  {
    out << "\n\nactive end-point pairs:";
    for (unsigned int i = 0; i < epvariables_.size(); ++i) epvariables_[i]->print(out);
  }

  out << "\n";
  // Todo add more relevant information here: segments, contact states, angles ...
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes,
    numnodalvalues>::print_summary_one_line_per_active_segment_pair(std::ostream& out) const
{
  check_init_setup();

  unsigned int numactivelargeanglesegmentpairs = cpvariables_.size();
  unsigned int numactivesmallanglesegmentpairs = gpvariables_.size();
  unsigned int numactiveendpointsegmentpairs = epvariables_.size();

  for (unsigned int i = 0; i < numactivelargeanglesegmentpairs; ++i)
  {
    out << "    " << std::setw(5) << std::left << element1()->id() << "(" << std::setw(3)
        << std::right << cpvariables_[i]->get_seg_ids().first + 1 << "/" << std::setw(3) << numseg1_
        << ")"
        << " " << std::setw(5) << std::left << element2()->id() << "(" << std::setw(3) << std::right
        << cpvariables_[i]->get_seg_ids().second + 1 << "/" << std::setw(3) << numseg2_ << ")"
        << "  CP    ";
    cpvariables_[i]->print_summary_in_one_line(out);
    out << "\n";
  }
  for (unsigned int i = 0; i < numactivesmallanglesegmentpairs; ++i)
  {
    out << "    " << std::setw(5) << std::left << element1()->id() << "(" << std::setw(3)
        << std::right << gpvariables_[i]->get_seg_ids().first + 1 << "/" << std::setw(3) << numseg1_
        << ")"
        << " " << std::setw(5) << std::left << element2()->id() << "(" << std::setw(3) << std::right
        << gpvariables_[i]->get_seg_ids().second + 1 << "/" << std::setw(3) << numseg2_ << ")"
        << "  GP    ";
    gpvariables_[i]->print_summary_in_one_line(out);
    out << "\n";
  }
  for (unsigned int i = 0; i < numactiveendpointsegmentpairs; ++i)
  {
    out << "    " << std::setw(5) << std::left << element1()->id() << "(" << std::setw(3)
        << std::right << epvariables_[i]->get_seg_ids().first + 1 << "/" << std::setw(3) << numseg1_
        << ")"
        << " " << std::setw(5) << std::left << element2()->id() << "(" << std::setw(3) << std::right
        << epvariables_[i]->get_seg_ids().second + 1 << "/" << std::setw(3) << numseg2_ << ")"
        << "  EP    ";
    epvariables_[i]->print_summary_in_one_line(out);
    out << "\n";
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::
    get_all_active_contact_point_coords_element1(
        std::vector<Core::LinAlg::Matrix<3, 1, double>>& coords) const
{
  int size1 = cpvariables_.size();
  int size2 = gpvariables_.size();
  int size3 = epvariables_.size();

  coords.resize(size1 + size2 + size3, Core::LinAlg::Matrix<3, 1, double>(true));

  for (int i = 0; i < size1; ++i)
  {
    TYPE eta1 = cpvariables_[i]->get_cp().first;
    for (int j = 0; j < 3; j++)
      coords[i](j) = Core::FADUtils::cast_to_double(r(eta1, element1())(j));
    /* Todo: Element1 might no longer be available on this proc after re-distribution!
     * take care that this method is called BEFORE any re-distribution
     * OR reformulate this such that Element1() is no longer required;
     * e.g. compute and store coord vector in instance of contact variables
     * (just like for gap and contact force) */
  }

  for (int i = size1; i < size2 + size1; ++i)
  {
    TYPE eta1 = gpvariables_[i - size1]->get_cp().first;
    for (int j = 0; j < 3; j++)
      coords[i](j) = Core::FADUtils::cast_to_double(r(eta1, element1())(j));
  }

  for (int i = size1 + size2; i < size1 + size2 + size3; ++i)
  {
    TYPE eta1 = epvariables_[i - size1 - size2]->get_cp().first;
    for (int j = 0; j < 3; j++)
      coords[i](j) = Core::FADUtils::cast_to_double(r(eta1, element1())(j));
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::
    get_all_active_contact_point_coords_element2(
        std::vector<Core::LinAlg::Matrix<3, 1, double>>& coords) const
{
  int size1 = cpvariables_.size();
  int size2 = gpvariables_.size();
  int size3 = epvariables_.size();

  coords.resize(size1 + size2 + size3, Core::LinAlg::Matrix<3, 1, double>(true));

  for (int i = 0; i < size1; ++i)
  {
    TYPE eta2 = cpvariables_[i]->get_cp().second;
    for (int j = 0; j < 3; j++)
      coords[i](j) = Core::FADUtils::cast_to_double(r(eta2, element2())(j));
  }

  for (int i = size1; i < size2 + size1; ++i)
  {
    TYPE eta2 = gpvariables_[i - size1]->get_cp().second;
    for (int j = 0; j < 3; j++)
      coords[i](j) = Core::FADUtils::cast_to_double(r(eta2, element2())(j));
  }

  for (int i = size1 + size2; i < size1 + size2 + size3; ++i)
  {
    TYPE eta2 = epvariables_[i - size1 - size2]->get_cp().second;
    for (int j = 0; j < 3; j++)
      coords[i](j) = Core::FADUtils::cast_to_double(r(eta2, element2())(j));
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes,
    numnodalvalues>::get_all_active_contact_forces(std::vector<double>& forces) const
{
  int size1 = cpvariables_.size();
  int size2 = gpvariables_.size();
  int size3 = epvariables_.size();

  forces.resize(size1 + size2 + size3, 0.0);

  for (int i = 0; i < size1; i++)
  {
    forces[i] =
        Core::FADUtils::cast_to_double(cpvariables_[i]->getfp() * cpvariables_[i]->get_p_pfac());
  }

  for (int i = size1; i < size2 + size1; i++)
  {
    forces[i] = Core::FADUtils::cast_to_double(
        gpvariables_[i - size1]->getfp() * gpvariables_[i - size1]->get_p_pfac());
  }

  for (int i = size1 + size2; i < size1 + size2 + size3; i++)
  {
    forces[i] = Core::FADUtils::cast_to_double(
        epvariables_[i - size1 - size2]->getfp() * epvariables_[i - size1 - size2]->get_p_pfac());
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::get_all_active_contact_gaps(
    std::vector<double>& gaps) const
{
  int size1 = cpvariables_.size();
  int size2 = gpvariables_.size();
  int size3 = epvariables_.size();

  gaps.resize(size1 + size2 + size3, 0.0);

  for (int i = 0; i < size1; i++)
  {
    gaps[i] = Core::FADUtils::cast_to_double(cpvariables_[i]->get_gap());
  }

  for (int i = size1; i < size2 + size1; i++)
  {
    gaps[i] = Core::FADUtils::cast_to_double(gpvariables_[i - size1]->get_gap());
  }

  for (int i = size1 + size2; i < size1 + size2 + size3; i++)
  {
    gaps[i] = Core::FADUtils::cast_to_double(epvariables_[i - size1 - size2]->get_gap());
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
double BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::get_energy() const
{
  if (params()->beam_to_beam_contact_params()->penalty_law() != BeamContact::pl_lp and
      params()->beam_to_beam_contact_params()->penalty_law() != BeamContact::pl_qp and
      params()->beam_to_beam_contact_params()->penalty_law() != BeamContact::pl_lpqp)
    FOUR_C_THROW("Contact Energy calculation not implemented for the chosen penalty law!");

  double energy = 0.0;

  for (unsigned int i = 0; i < cpvariables_.size(); ++i)
  {
    double ppfac = Core::FADUtils::cast_to_double(cpvariables_[i]->get_p_pfac());
    double e = -cpvariables_[i]->get_integrated_energy();
    energy += ppfac * e;
  }

  for (unsigned int i = 0; i < gpvariables_.size(); ++i)
  {
    double ppfac = Core::FADUtils::cast_to_double(gpvariables_[i]->get_p_pfac());
    double e = -gpvariables_[i]->get_integrated_energy();
    energy += ppfac * e;
  }

  for (unsigned int i = 0; i < epvariables_.size(); ++i)
  {
    double ppfac = Core::FADUtils::cast_to_double(epvariables_[i]->get_p_pfac());
    double e = -epvariables_[i]->get_integrated_energy();
    energy += ppfac * e;
  }

  return energy;
}


#ifdef FADCHECKS
/*----------------------------------------------------------------------*
 |  FAD-Check for Linearizations of contact point            meier 02/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes, numnodalvalues>::fad_check_lin_xi_and_lin_eta(
    const Core::LinAlg::Matrix<3, 1, TYPE>& delta_r, const Core::LinAlg::Matrix<3, 1, TYPE>& r1_xi,
    const Core::LinAlg::Matrix<3, 1, TYPE>& r2_xi, const Core::LinAlg::Matrix<3, 1, TYPE>& r1_xixi,
    const Core::LinAlg::Matrix<3, 1, TYPE>& r2_xixi,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xi,
    const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xi)
{
  Core::LinAlg::Matrix<2, 1, TYPE> f(true);
  Core::LinAlg::Matrix<3, 1, TYPE> t1_dummy(true);
  Core::LinAlg::Matrix<3, 1, TYPE> t2_dummy(true);

  // compute norm of difference vector to scale the equations
  // (this yields better conditioning)
  // Note: Even if automatic differentiation via FAD is applied, norm_delta_r has to be of type
  // double since this factor is needed for a pure scaling of the nonlinear CCP and has not to be
  // linearized!
  double norm_delta_r = Core::FADUtils::cast_to_double(Core::FADUtils::vector_norm<3>(delta_r));

  evaluate_orthogonality_condition(f, delta_r, norm_delta_r, r1_xi, r2_xi, t1_dummy, t2_dummy);

  //**********************************************************************
  // we have to solve the following system of equations:
  //  _              _       _      _       _              _      _       _
  // | L(1,1)  L(1,2) |    | Lin_Xi  |    |  B(1,1)  B(1,2) |   | Lin_d1 |
  // |                | *  |         | =  |                 | * |        |
  // |_L(2,1)  L(2,2)_|    |_Lin_Eta_|    |_B(2,1)  B(2,2)_ |   |_Lin_d2_|
  //
  // this can be done easily because it is a linear 2x2-system.
  // we obtain the solution by inverting matrix L:
  //
  // [Lin_Xi; Lin_Eta] = L^-1 * B * [Lin_d1; Lin_d2] = D * [Lin_d1; Lin_d2]
  //
  //**********************************************************************

  const unsigned int dim1 = 3 * numnodes * numnodalvalues;
  const unsigned int dim2 = 3 * numnodes * numnodalvalues;

  // matrices to compute Lin_Xi and Lin_Eta
  Core::LinAlg::Matrix<2, 2, TYPE> L(true);
  Core::LinAlg::Matrix<2, 2, TYPE> L_inv(true);
  Core::LinAlg::Matrix<2, dim1 + dim2, TYPE> B(true);
  Core::LinAlg::Matrix<2, dim1 + dim2, TYPE> D(true);

  // compute L elementwise
  L(0, 0) = f(0).dx(2 * 3 * numnodes * numnodalvalues);
  L(0, 1) = f(0).dx(2 * 3 * numnodes * numnodalvalues + 1);
  L(1, 0) = f(1).dx(2 * 3 * numnodes * numnodalvalues);
  L(1, 1) = f(1).dx(2 * 3 * numnodes * numnodalvalues + 1);

  // invert L by hand
  TYPE det_L = L(0, 0) * L(1, 1) - L(0, 1) * L(1, 0);
  if (Core::FADUtils::cast_to_double(Core::FADUtils::Norm(det_L)) < DETERMINANTTOL)
    FOUR_C_THROW("ERROR: Determinant of L = 0");
  L_inv(0, 0) = L(1, 1) / det_L;
  L_inv(0, 1) = -L(0, 1) / det_L;
  L_inv(1, 0) = -L(1, 0) / det_L;
  L_inv(1, 1) = L(0, 0) / det_L;

  for (unsigned int j = 0; j < dim1 + dim2; j++)
  {
    B(0, j) = -f(0).dx(j);
    B(1, j) = -f(1).dx(j);
  }

  // compute D = L^-1 * B
  D.multiply(L_inv, B);

  std::cout << "linxi and lineta: " << std::endl;

  std::cout << D << std::endl;

  return;
}
/*----------------------------------------------------------------------*
 |  End: FAD-Check for Linearizations of contact point
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  FAD-Check for Linearizations of CCP                      meier 02/14|
 *----------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BeamInteraction::BeamToBeamContactPair<numnodes,
    numnodalvalues>::fad_check_lin_orthogonality_condition(const Core::LinAlg::Matrix<3, 1, TYPE>&
                                                               delta_r,
    const double& norm_delta_r, const Core::LinAlg::Matrix<3, 1, TYPE>& r1_xi,
    const Core::LinAlg::Matrix<3, 1, TYPE>& r2_xi, const Core::LinAlg::Matrix<3, 1, TYPE>& t1,
    const Core::LinAlg::Matrix<3, 1, TYPE>& t2)
{
  Core::LinAlg::Matrix<2, 1, TYPE> f(true);

  evaluate_orthogonality_condition(f, delta_r, norm_delta_r, r1_xi, r2_xi, t1, t2);

  Core::LinAlg::Matrix<2, 2, TYPE> df(true);

  for (unsigned int i = 0; i < 2; i++)
  {
    for (unsigned int j = 0; j < 2; j++)
    {
      df(i, j) = f(i).dx(2 * 3 * numnodes * numnodalvalues + j);
    }
  }

  std::cout << "df_FAD: " << std::endl;

  std::cout << df << std::endl;

  return;
}
/*----------------------------------------------------------------------*
 |  End: FAD-Check for Linearizations of CPP
 *----------------------------------------------------------------------*/
#endif

/*----------------------------------------------------------------------*
|  FD-Check of stiffness matrix                              meier 11/14|
*-----------------------------------------------------------------------*/
// template<unsigned int numnodes, unsigned int numnodalvalues>
// void BeamInteraction::BeamToBeamContactPair<numnodes,
// numnodalvalues>::fd_check(Core::LinAlg::SparseMatrix& stiffmatrix,
//                                                              Core::LinAlg::Vector<double>& fint,
//                                                              const double& pp,
//                                                              std::map<std::pair<int,int>,
//                                                              std::shared_ptr<Beam3contactinterface
//                                                              >
//                                                              >& contactpairmap,
//                                                              Teuchos::ParameterList&
//                                                              timeintparams, bool fdcheck)
//{
//  //This FD-Check is very general, since it applies the complete method "Evaluate" recursively.
//  //Therefore, all changes within this class are automatically considered and have not to be
//  //adapted in this finite difference check!
//  if(fint.GlobalLength()>2*3*numnodes*numnodalvalues)
//    FOUR_C_THROW("So far, this fd_check only works for simulations with two elements!!!");
//
//  Core::LinAlg::Vector<double> fint1(fint);
//  fint1.PutScalar(0.0);
//  Core::LinAlg::Vector<double> fint2(fint);
//  fint2.PutScalar(0.0);
//
//  Core::LinAlg::SparseMatrix stiffmatrix_analyt(stiffmatrix);
//  stiffmatrix_analyt.PutScalar(0.0);
//
//  Core::LinAlg::SparseMatrix stiffmatrix_dummy(stiffmatrix);
//  stiffmatrix_dummy.PutScalar(0.0);
//
//  Core::LinAlg::Matrix<2*3*numnodes*numnodalvalues,2*3*numnodes*numnodalvalues>
//  stiffmatrix_fd(true);
//
//  Core::LinAlg::Matrix<2*3*numnodes*numnodalvalues,2*3*numnodes*numnodalvalues> fint2_mat(true);
//
//  double delta = 1.0e-10;
//
//  std::cout << "undisturbed configuration: " << std::endl;
//
//  this->evaluate(stiffmatrix_analyt,fint1,pp,contactpairmap,timeintparams,true);
//
////  std::cout << std::setprecision(25) << "fint1: " << std::endl;
//
////  fint1.print(std::cout);
//
//  std::vector<double> xi1((int)gpvariables_.size(),0.0);
//  std::vector<double> eta1((int)gpvariables_.size(),0.0);
//
//  std::vector<Core::LinAlg::Matrix<2*3*numnodes*numnodalvalues,1> > xi2((int)gpvariables_.size());
//  std::vector<Core::LinAlg::Matrix<2*3*numnodes*numnodalvalues,1> >
//  eta2((int)gpvariables_.size());
//
//
//  for (int i=0;i<(int)gpvariables_.size();i++)
//  {
//    xi1[i]=Core::FADUtils::cast_to_double(gpvariables_[i]->GetCP().first);
//    eta1[i]=Core::FADUtils::cast_to_double(gpvariables_[i]->GetCP().second);
//  }
//
//  for(int dof=0;dof<2*3*numnodes*numnodalvalues;dof++)
//  {
//    std::cout << "disturbed configuration: " << std::endl;
//
//    //Add delta
//    if(dof<3*numnodes*numnodalvalues)
//      ele1pos_(dof)+=delta;
//    else
//      ele2pos_(dof-3*numnodes*numnodalvalues)+=delta;
//
//    fint2.PutScalar(0.0);
//    stiffmatrix_dummy.PutScalar(0.0);
//
//    this->evaluate(stiffmatrix_dummy,fint2,pp,contactpairmap,timeintparams,true);
//
////    std::cout << std::setprecision(25) << "fint2: " << std::endl;
////
////    fint2.print(std::cout);
//
//    for(int i=0;i<(int)gpvariables_.size();i++)
//    {
//      xi2[i](dof)=Core::FADUtils::cast_to_double(gpvariables_[i]->GetCP().first);
//      eta2[i](dof)=Core::FADUtils::cast_to_double(gpvariables_[i]->GetCP().second);
//    }
//
//    for(unsigned int i=0;i<2*3*numnodalvalues*numnodes;i++)
//    {
//      fint2_mat(i,dof)=fint2[i];
//
//      if(fabs(fint2[i])<1.0e-10 and fabs(fint1[i])<1.0e-10)
//        stiffmatrix_fd(i,dof) = 999999999999;
//      else
//        stiffmatrix_fd(i,dof) = -(fint2[i]-fint1[i])/delta;
//    }
//
//    //restore original displacements
//    if(dof<3*numnodes*numnodalvalues)
//      ele1pos_(dof)-=delta;
//    else
//      ele2pos_(dof-3*numnodes*numnodalvalues)-=delta;
//  }
//
//  std::cout << "FD_LIN: " << std::endl;
//
//  for(unsigned int i=0;i<2*3*numnodes*numnodalvalues;i++)
//  {
//    for(unsigned int j=0;j<2*3*numnodes*numnodalvalues;j++)
//    {
//      std::cout << "row: " << i << "   " << "col: " << j << "   " << stiffmatrix_fd(i,j) << "
//      fint2: " << fint2_mat(i,j) << "   fint1: " << fint1[i] << std::endl;
//    }
//  }
//
//  std::cout << "ANALYT_LIN: " << std::endl;
//
//  std::cout << (*(stiffmatrix_analyt.EpetraMatrix())) << std::endl;
//}
/*----------------------------------------------------------------------*
|  end: FD-Check of stiffness matrix
*-----------------------------------------------------------------------*/

// Possible template cases: this is necessary for the compiler
template class BeamInteraction::BeamToBeamContactPair<2, 1>;
template class BeamInteraction::BeamToBeamContactPair<3, 1>;
template class BeamInteraction::BeamToBeamContactPair<4, 1>;
template class BeamInteraction::BeamToBeamContactPair<5, 1>;
template class BeamInteraction::BeamToBeamContactPair<2, 2>;

FOUR_C_NAMESPACE_CLOSE
