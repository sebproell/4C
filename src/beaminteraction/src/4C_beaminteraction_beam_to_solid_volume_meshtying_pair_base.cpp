// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_beaminteraction_beam_to_solid_volume_meshtying_pair_base.hpp"

#include "4C_beaminteraction_beam_to_beam_contact_defines.hpp"
#include "4C_beaminteraction_beam_to_beam_contact_utils.hpp"
#include "4C_beaminteraction_beam_to_solid_visualization_output_writer_base.hpp"
#include "4C_beaminteraction_beam_to_solid_visualization_output_writer_visualization.hpp"
#include "4C_beaminteraction_beam_to_solid_volume_meshtying_params.hpp"
#include "4C_beaminteraction_beam_to_solid_volume_meshtying_visualization_output_params.hpp"
#include "4C_beaminteraction_contact_pair.hpp"
#include "4C_beaminteraction_contact_params.hpp"
#include "4C_beaminteraction_geometry_pair_access_traits.hpp"
#include "4C_geometry_pair_element_evaluation_functions.hpp"
#include "4C_geometry_pair_factory.hpp"
#include "4C_geometry_pair_line_to_volume.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_utils_densematrix_inverse.hpp"

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
template <typename ScalarType, typename Beam, typename Solid>
BeamInteraction::BeamToSolidVolumeMeshtyingPairBase<ScalarType, Beam,
    Solid>::BeamToSolidVolumeMeshtyingPairBase()
    : base_class(), meshtying_is_evaluated_(false)
{
  // Empty constructor.
}

/**
 *
 */
template <typename ScalarType, typename Beam, typename Solid>
void BeamInteraction::BeamToSolidVolumeMeshtyingPairBase<ScalarType, Beam, Solid>::setup()
{
  // Call setup of base class first.
  base_class::setup();

  // Get the solid element data container
  ele2posref_ = GeometryPair::InitializeElementData<Solid, double>::initialize(this->element2());
  ele2pos_ = GeometryPair::InitializeElementData<Solid, ScalarType>::initialize(this->element2());

  // Set reference nodal positions for the solid element
  for (unsigned int n = 0; n < Solid::n_nodes_; ++n)
  {
    const Core::Nodes::Node* node = this->element2()->nodes()[n];
    for (int d = 0; d < 3; ++d) ele2posref_.element_position_(3 * n + d) = node->x()[d];
  }

  // Initialize current nodal positions for the solid element
  for (unsigned int i = 0; i < Solid::n_dof_; i++) ele2pos_.element_position_(i) = 0.0;
}

/**
 *
 */
template <typename ScalarType, typename Beam, typename Solid>
void BeamInteraction::BeamToSolidVolumeMeshtyingPairBase<ScalarType, Beam,
    Solid>::create_geometry_pair(const Core::Elements::Element* element1,
    const Core::Elements::Element* element2,
    const std::shared_ptr<GeometryPair::GeometryEvaluationDataBase>& geometry_evaluation_data_ptr)
{
  this->geometry_pair_ = GeometryPair::geometry_pair_line_to_volume_factory<double, Beam, Solid>(
      element1, element2, geometry_evaluation_data_ptr);
}


/**
 *
 */
template <typename ScalarType, typename Beam, typename Solid>
void BeamInteraction::BeamToSolidVolumeMeshtyingPairBase<ScalarType, Beam, Solid>::pre_evaluate()
{
  // Call pre_evaluate on the geometry Pair.
  if (!meshtying_is_evaluated_)
  {
    GeometryPair::ElementData<Beam, double> beam_coupling_ref;
    GeometryPair::ElementData<Solid, double> solid_coupling_ref;
    this->get_coupling_reference_position(beam_coupling_ref, solid_coupling_ref);
    cast_geometry_pair()->pre_evaluate(
        beam_coupling_ref, solid_coupling_ref, this->line_to_3D_segments_);
  }
}

/**
 *
 */
template <typename ScalarType, typename Beam, typename Solid>
void BeamInteraction::BeamToSolidVolumeMeshtyingPairBase<ScalarType, Beam, Solid>::reset_state(
    const std::vector<double>& beam_centerline_dofvec,
    const std::vector<double>& solid_nodal_dofvec)
{
  // Call the method in the parent class.
  base_class::reset_state(beam_centerline_dofvec, solid_nodal_dofvec);

  // Solid element.
  for (unsigned int i = 0; i < Solid::n_dof_; i++)
  {
    ele2pos_.element_position_(i) = Core::FADUtils::HigherOrderFadValue<ScalarType>::apply(
        Beam::n_dof_ + Solid::n_dof_, Beam::n_dof_ + i, solid_nodal_dofvec[i]);
  }
}

/**
 *
 */
template <typename ScalarType, typename Beam, typename Solid>
void BeamInteraction::BeamToSolidVolumeMeshtyingPairBase<ScalarType, Beam, Solid>::
    set_restart_displacement(const std::vector<std::vector<double>>& centerline_restart_vec_)
{
  // Call the parent method.
  base_class::set_restart_displacement(centerline_restart_vec_);

  // We only set the restart displacement, if the current section has the restart coupling flag.
  if (this->params()->beam_to_solid_volume_meshtying_params()->get_couple_restart_state())
  {
    for (unsigned int i_dof = 0; i_dof < Beam::n_dof_; i_dof++)
      ele1posref_offset_(i_dof) = centerline_restart_vec_[0][i_dof];

    // Add the displacement at the restart step to the solid reference position.
    for (unsigned int i_dof = 0; i_dof < Solid::n_dof_; i_dof++)
      ele2posref_offset_(i_dof) = centerline_restart_vec_[1][i_dof];
  }
}


/**
 *
 */
template <typename ScalarType, typename Beam, typename Solid>
void BeamInteraction::BeamToSolidVolumeMeshtyingPairBase<ScalarType, Beam,
    Solid>::get_pair_visualization(std::shared_ptr<BeamToSolidVisualizationOutputWriterBase>
                                       visualization_writer,
    Teuchos::ParameterList& visualization_params) const
{
  // Get visualization of base class.
  base_class::get_pair_visualization(visualization_writer, visualization_params);

  // Get the writers.
  std::shared_ptr<BeamInteraction::BeamToSolidOutputWriterVisualization>
      visualization_segmentation =
          visualization_writer->get_visualization_writer("btsv-segmentation");
  std::shared_ptr<BeamInteraction::BeamToSolidOutputWriterVisualization>
      visualization_integration_points =
          visualization_writer->get_visualization_writer("btsv-integration-points");
  if (!visualization_segmentation and visualization_integration_points == nullptr) return;

  const std::shared_ptr<const BeamToSolidVolumeMeshtyingVisualizationOutputParams>&
      output_params_ptr =
          visualization_params
              .get<std::shared_ptr<const BeamToSolidVolumeMeshtyingVisualizationOutputParams>>(
                  "btsv-output_params_ptr");
  const bool write_unique_ids = output_params_ptr->get_write_unique_ids_flag();

  if (visualization_segmentation != nullptr)
  {
    // Setup variables.
    Core::LinAlg::Matrix<3, 1, ScalarType> X;
    Core::LinAlg::Matrix<3, 1, ScalarType> u;
    Core::LinAlg::Matrix<3, 1, ScalarType> r;

    // Get the visualization vectors.
    auto& visualization_data = visualization_segmentation->get_visualization_data();
    std::vector<double>& point_coordinates = visualization_data.get_point_coordinates();
    std::vector<double>& displacement = visualization_data.get_point_data<double>("displacement");

    std::vector<int>* pair_beam_id = nullptr;
    std::vector<int>* pair_solid_id = nullptr;
    if (write_unique_ids)
    {
      pair_beam_id = &(visualization_data.get_point_data<int>("uid_0_pair_beam_id"));
      pair_solid_id = &(visualization_data.get_point_data<int>("uid_1_pair_solid_id"));
    }

    // Loop over the segments on the beam.
    for (const auto& segment : this->line_to_3D_segments_)
    {
      // Add the left and right boundary point of the segment.
      for (const auto& segmentation_point : {segment.get_eta_a(), segment.get_eta_b()})
      {
        GeometryPair::evaluate_position<Beam>(segmentation_point, this->ele1posref_, X);
        GeometryPair::evaluate_position<Beam>(segmentation_point, this->ele1pos_, r);
        u = r;
        u -= X;
        for (unsigned int dim = 0; dim < 3; dim++)
        {
          point_coordinates.push_back(Core::FADUtils::cast_to_double(X(dim)));
          displacement.push_back(Core::FADUtils::cast_to_double(u(dim)));
        }

        if (write_unique_ids)
        {
          pair_beam_id->push_back(this->element1()->id());
          pair_solid_id->push_back(this->element2()->id());
        }
      }
    }
  }

  // If a writer exists for integration point data, add the integration point data.
  if (visualization_integration_points != nullptr)
  {
    // Setup variables.
    Core::LinAlg::Matrix<3, 1, double> X;
    Core::LinAlg::Matrix<3, 1, double> u;
    Core::LinAlg::Matrix<3, 1, double> r;
    Core::LinAlg::Matrix<3, 1, double> r_solid;
    Core::LinAlg::Matrix<3, 1, double> force_integration_point;

    // Get the visualization vectors.
    auto& visualization_data = visualization_integration_points->get_visualization_data();
    std::vector<double>& point_coordinates = visualization_data.get_point_coordinates();
    std::vector<double>& displacement = visualization_data.get_point_data<double>("displacement");
    std::vector<double>& force = visualization_data.get_point_data<double>("force");

    std::vector<int>* pair_beam_id = nullptr;
    std::vector<int>* pair_solid_id = nullptr;
    if (write_unique_ids)
    {
      pair_beam_id = &(visualization_data.get_point_data<int>("uid_0_pair_beam_id"));
      pair_solid_id = &(visualization_data.get_point_data<int>("uid_1_pair_solid_id"));
    }

    // Loop over the segments on the beam.
    for (const auto& segment : this->line_to_3D_segments_)
    {
      // Add the integration points.
      for (const auto& projection_point : segment.get_projection_points())
      {
        this->evaluate_beam_position_double(projection_point, X, true);
        this->evaluate_beam_position_double(projection_point, r, false);
        u = r;
        u -= X;
        GeometryPair::evaluate_position<Solid>(projection_point.get_xi(),
            GeometryPair::ElementDataToDouble<Solid>::to_double(this->ele2pos_), r_solid);
        evaluate_penalty_force_double(r, r_solid, force_integration_point);
        for (unsigned int dim = 0; dim < 3; dim++)
        {
          point_coordinates.push_back(X(dim));
          displacement.push_back(u(dim));
          force.push_back(force_integration_point(dim));
        }

        if (write_unique_ids)
        {
          pair_beam_id->push_back(this->element1()->id());
          pair_solid_id->push_back(this->element2()->id());
        }
      }
    }
  }
}

/**
 *
 */
template <typename ScalarType, typename Beam, typename Solid>
void BeamInteraction::BeamToSolidVolumeMeshtyingPairBase<ScalarType, Beam,
    Solid>::evaluate_penalty_force_double(const Core::LinAlg::Matrix<3, 1, double>& r_beam,
    const Core::LinAlg::Matrix<3, 1, double>& r_solid,
    Core::LinAlg::Matrix<3, 1, double>& force) const
{
  // The base implementation of the force is a simple linear penalty law.
  force = r_solid;
  force -= r_beam;
  force.scale(this->params()->beam_to_solid_volume_meshtying_params()->get_penalty_parameter());
}

/**
 *
 */
template <typename ScalarType, typename Beam, typename Solid>
void BeamInteraction::BeamToSolidVolumeMeshtyingPairBase<ScalarType, Beam,
    Solid>::get_coupling_reference_position(GeometryPair::ElementData<Beam, double>&
                                                beam_coupling_ref,
    GeometryPair::ElementData<Solid, double>& solid_coupling_ref) const
{
  // Add the offset to the reference position.
  beam_coupling_ref = this->ele1posref_;
  beam_coupling_ref.element_position_ += this->ele1posref_offset_;
  solid_coupling_ref = ele2posref_;
  solid_coupling_ref.element_position_ += ele2posref_offset_;
}


/**
 * Explicit template initialization of template class.
 */
namespace BeamInteraction
{
  using namespace GeometryPair;

  template class BeamToSolidVolumeMeshtyingPairBase<double, t_hermite, t_hex8>;
  template class BeamToSolidVolumeMeshtyingPairBase<double, t_hermite, t_hex20>;
  template class BeamToSolidVolumeMeshtyingPairBase<double, t_hermite, t_hex27>;
  template class BeamToSolidVolumeMeshtyingPairBase<double, t_hermite, t_tet4>;
  template class BeamToSolidVolumeMeshtyingPairBase<double, t_hermite, t_tet10>;
  template class BeamToSolidVolumeMeshtyingPairBase<double, t_hermite, t_nurbs27>;

  template class BeamToSolidVolumeMeshtyingPairBase<line_to_volume_scalar_type<t_hermite, t_hex8>,
      t_hermite, t_hex8>;
  template class BeamToSolidVolumeMeshtyingPairBase<line_to_volume_scalar_type<t_hermite, t_hex20>,
      t_hermite, t_hex20>;
  template class BeamToSolidVolumeMeshtyingPairBase<line_to_volume_scalar_type<t_hermite, t_hex27>,
      t_hermite, t_hex27>;
  template class BeamToSolidVolumeMeshtyingPairBase<line_to_volume_scalar_type<t_hermite, t_tet4>,
      t_hermite, t_tet4>;
  template class BeamToSolidVolumeMeshtyingPairBase<line_to_volume_scalar_type<t_hermite, t_tet10>,
      t_hermite, t_tet10>;
  template class BeamToSolidVolumeMeshtyingPairBase<
      line_to_volume_scalar_type<t_hermite, t_nurbs27>, t_hermite, t_nurbs27>;
}  // namespace BeamInteraction

FOUR_C_NAMESPACE_CLOSE
