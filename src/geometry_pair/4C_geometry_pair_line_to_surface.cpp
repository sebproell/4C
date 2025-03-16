// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_geometry_pair_line_to_surface.hpp"

#include "4C_beam3_base.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_geometry_pair_constants.hpp"
#include "4C_geometry_pair_element.hpp"
#include "4C_geometry_pair_element_evaluation_functions.hpp"
#include "4C_geometry_pair_scalar_types.hpp"
#include "4C_geometry_pair_utility_classes.hpp"
#include "4C_linalg_utils_densematrix_inverse.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_fad.hpp"

FOUR_C_NAMESPACE_OPEN

/**
 *
 */
template <typename ScalarType, typename Line, typename Surface>
GEOMETRYPAIR::GeometryPairLineToSurface<ScalarType, Line, Surface>::GeometryPairLineToSurface(
    const Core::Elements::Element* element1, const Core::Elements::Element* element2,
    const std::shared_ptr<GEOMETRYPAIR::LineToSurfaceEvaluationData>&
        line_to_surface_evaluation_data)
    : GeometryPair(element1, element2),
      line_to_surface_evaluation_data_(line_to_surface_evaluation_data),
      is_unit_test_(false)
{
  // For the current implementation, the line element has to be on the same processor as the pair
  // object. This is because the tracking vector in LineTo3DEvaluationData is only local and we
  // need this vector for segmentation e.t.c.
  int myrank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (element1->owner() != myrank)
    FOUR_C_THROW(
        "The GeometryPairLineToSurface pair has to be on the same processor as the line element! "
        "Currently the pair is on rank {}, the line element on {}!",
        myrank, element1->owner());
}

/**
 *
 */
template <typename ScalarType, typename Line, typename Surface>
void GEOMETRYPAIR::GeometryPairLineToSurface<ScalarType, Line, Surface>::project_point_to_other(
    const Core::LinAlg::Matrix<3, 1, ScalarType>& point,
    const ElementData<Surface, ScalarType>& element_data_surface,
    Core::LinAlg::Matrix<3, 1, ScalarType>& xi, ProjectionResult& projection_result,
    const bool min_one_iteration) const
{
  project_point_to_surface(point, element_data_surface, xi, projection_result,
      get_surface_normal_influence_direction(element_data_surface), min_one_iteration);
}


/**
 *
 */
template <typename ScalarType, typename Line, typename Surface>
void GEOMETRYPAIR::GeometryPairLineToSurface<ScalarType, Line, Surface>::intersect_line_with_other(
    const ElementData<Line, ScalarType>& element_data_line,
    const ElementData<Surface, ScalarType>& element_data_surface,
    std::vector<ProjectionPoint1DTo3D<ScalarType>>& intersection_points,
    const ScalarType& eta_start, const Core::LinAlg::Matrix<3, 1, ScalarType>& xi_start) const
{
  unsigned int n_faces;
  std::vector<unsigned int> face_fixed_parameters;
  std::vector<double> face_fixed_values;
  get_face_fixed_parameters(n_faces, face_fixed_parameters, face_fixed_values);

  // Clear the input vector.
  intersection_points.clear();
  intersection_points.reserve(n_faces);

  // Create variables.
  ScalarType eta;
  Core::LinAlg::Matrix<3, 1, ScalarType> xi;
  ProjectionResult intersection_found;

  // Try to intersect the beam with each face.
  for (unsigned int i = 0; i < n_faces; i++)
  {
    // Set starting values.
    xi = xi_start;
    eta = eta_start;

    // Intersect the line with the surface.
    intersect_line_with_surface_edge(element_data_line, element_data_surface,
        face_fixed_parameters[i], face_fixed_values[i], eta, xi, intersection_found);

    // If a valid intersection is found, add it to the output vector.
    if (intersection_found == ProjectionResult::projection_found_valid)
    {
      intersection_points.push_back(ProjectionPoint1DTo3D<ScalarType>(eta, xi));
      intersection_points.back().set_intersection_face(i);
    }
  }
}

/**
 *
 */
template <typename ScalarType, typename Line, typename Surface>
void GEOMETRYPAIR::GeometryPairLineToSurface<ScalarType, Line,
    Surface>::intersect_line_with_surface_edge(const ElementData<Line, ScalarType>&
                                                   element_data_line,
    const ElementData<Surface, ScalarType>& element_data_surface,
    const unsigned int& fixed_parameter, const double& fixed_value, ScalarType& eta,
    Core::LinAlg::Matrix<3, 1, ScalarType>& xi, ProjectionResult& projection_result,
    const bool min_one_iteration) const
{
  // Check the input parameters.
  {
    if (Surface::geometry_type_ == DiscretizationTypeGeometry::quad && fixed_parameter > 1)
      FOUR_C_THROW(
          "Fixed_parameter in intersect_line_with_surface_edge has to be smaller than 2 with a "
          "quad element.");
    else if (Surface::geometry_type_ == DiscretizationTypeGeometry::none)
      FOUR_C_THROW("Wrong DiscretizationTypeGeometry type given.");
    else if (fixed_parameter > 2)
      FOUR_C_THROW("fixed_parameter in intersect_line_with_surface_edge can be 2 at maximum.");
  }

  // Approximated influence size of the surface.
  const double normal_influence_direction =
      get_surface_normal_influence_direction(element_data_surface);

  // Initialize data structures.
  // Point on line.
  Core::LinAlg::Matrix<3, 1, ScalarType> r_line;
  Core::LinAlg::Matrix<3, 1, ScalarType> dr_line;

  // Point on surface.
  Core::LinAlg::Matrix<3, 1, ScalarType> r_surface;
  Core::LinAlg::Matrix<3, 3, ScalarType> dr_surface;

  // Residuum.
  Core::LinAlg::Matrix<4, 1, ScalarType> residuum;
  Core::LinAlg::Matrix<4, 1, ScalarType> delta_xi;
  // Initialize the increment with a value that will not pass the first convergence check.
  delta_xi.put_scalar(10 * Constants::projection_xi_eta_tol);

  // Jacobian / inverse.
  Core::LinAlg::Matrix<4, 4, ScalarType> J_J_inv;

  // Reset the projection result flag.
  projection_result = ProjectionResult::projection_not_found;

  {
    // Local Newton iteration.
    unsigned int counter = 0;
    while (counter < Constants::local_newton_iter_max)
    {
      // Evaluate the position and its derivative on the line.
      evaluate_position<Line>(eta, element_data_line, r_line);
      evaluate_position_derivative1<Line>(eta, element_data_line, dr_line);

      // Evaluate the position and its derivative on the surface.
      evaluate_surface_position_and_derivative(element_data_surface, xi, r_surface, dr_surface);

      // Evaluate the residuum $r_{surface} - r_{line} = R_{pos}$ and $xi(i) - value = R_{edge}$
      J_J_inv.put_scalar(0.);
      residuum.put_scalar(0.);
      for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
        residuum(i_dir) = r_surface(i_dir) - r_line(i_dir);

      if (fixed_parameter < 2)
      {
        residuum(3) = xi(fixed_parameter) - fixed_value;
        J_J_inv(3, fixed_parameter) = 1.;
      }
      else
      {
        for (unsigned int i = 0; i < 2; i++)
        {
          residuum(3) += xi(i);
          J_J_inv(3, i) = 1.;
        }
        residuum(3) -= fixed_value;
      }

      if (counter == 0 and min_one_iteration)
      {
        // if the min_one_iteration flag is set we run at least one iteration, so the dependency on
        // FAD variables is calculated correctly.
      }
      else if (Core::FADUtils::vector_norm(residuum) < Constants::local_newton_res_tol &&
               Core::FADUtils::vector_norm(delta_xi) < Constants::projection_xi_eta_tol)
      {
        // System is solved, now check if the parameter coordinates are valid.
        if (valid_parameter_1d(eta) &&
            valid_parameter_surface<ScalarType, Surface>(xi, normal_influence_direction))
          projection_result = ProjectionResult::projection_found_valid;
        else
          projection_result = ProjectionResult::projection_found_not_valid;
        break;
      }

      // Check if residuum is in a sensible range where we still expect to find a solution.
      if (Core::FADUtils::vector_norm(residuum) > Constants::local_newton_res_max) break;

      // Fill up the jacobian.
      for (unsigned int i = 0; i < 3; i++)
      {
        for (unsigned int j = 0; j < 3; j++)
        {
          J_J_inv(i, j) = dr_surface(i, j);
        }
        J_J_inv(i, 3) = -dr_line(i);
      }

      // Solve the linearized system.
      if (Core::LinAlg::solve_linear_system_do_not_throw_error_on_zero_determinant_scaled(
              J_J_inv, residuum, delta_xi, Constants::local_newton_det_tol))
      {
        // Set the new parameter coordinates.
        eta -= delta_xi(3);
        for (unsigned int i = 0; i < 3; i++) xi(i) -= delta_xi(i);

        // Advance Newton iteration counter.
        counter++;
      }
      else
        break;
    }
  }
}

/**
 *
 */
template <typename ScalarType, typename Line, typename Surface>
double GEOMETRYPAIR::GeometryPairLineToSurface<ScalarType, Line,
    Surface>::get_surface_normal_influence_direction(const ElementData<Surface, ScalarType>&
        element_data_surface) const
{
  if (is_unit_test_)
    return -1.0;
  else
  {
    double surface_size = get_surface_size(element_data_surface);
    double line_tube_size_radius = (dynamic_cast<const Discret::Elements::Beam3Base*>(element1()))
                                       ->get_circular_cross_section_radius_for_interactions();
    return std::max(surface_size, 3.0 * line_tube_size_radius);
  }
}

/**
 *
 */
template <typename ScalarType, typename Line, typename Surface>
double GEOMETRYPAIR::GeometryPairLineToSurface<ScalarType, Line, Surface>::get_surface_size(
    const ElementData<Surface, ScalarType>& element_data_surface) const
{
  // Get the position of the first 3 nodes of the surface.
  Core::LinAlg::Matrix<2, 1, double> xi_corner_node;
  Core::LinAlg::Matrix<3, 1, Core::LinAlg::Matrix<3, 1, double>> corner_nodes;
  Core::LinAlg::SerialDenseMatrix nodal_coordinates =
      Core::FE::get_ele_node_numbering_nodes_paramspace(Surface::discretization_);
  const auto element_data_surface_double =
      ElementDataToDouble<Surface>::to_double(element_data_surface);
  for (unsigned int i_node = 0; i_node < 3; i_node++)
  {
    for (unsigned int i_dim = 0; i_dim < 2; i_dim++)
      xi_corner_node(i_dim) = nodal_coordinates(i_dim, i_node);
    evaluate_position<Surface>(xi_corner_node, element_data_surface_double, corner_nodes(i_node));
  }

  // Calculate the maximum distance between the three points.
  double max_distance = 0.0;
  double distance = 0.0;
  Core::LinAlg::Matrix<3, 1, double> diff;
  for (unsigned int i_node = 0; i_node < 3; i_node++)
  {
    for (unsigned int j_node = 0; j_node < 3; j_node++)
    {
      if (i_node == j_node) continue;

      diff = corner_nodes(j_node);
      diff -= corner_nodes(i_node);
      distance = Core::FADUtils::vector_norm(diff);
      if (distance > max_distance) max_distance = distance;
    }
  }

  return max_distance;
}

/**
 *
 */
template <typename ScalarType, typename Line, typename Surface>
void GEOMETRYPAIR::GeometryPairLineToSurface<ScalarType, Line, Surface>::get_face_fixed_parameters(
    unsigned int& n_faces, std::vector<unsigned int>& face_fixed_parameters,
    std::vector<double>& face_fixed_values) const
{
  if (Surface::geometry_type_ == DiscretizationTypeGeometry::quad)
  {
    n_faces = 4;
    face_fixed_parameters = {0, 0, 1, 1};
    face_fixed_values = {-1., 1., -1., 1.};
  }
  else if (Surface::geometry_type_ == DiscretizationTypeGeometry::triangle)
  {
    n_faces = 3;
    face_fixed_parameters = {0, 1, 2};
    face_fixed_values = {0., 0., 1.};
  }
  else
  {
    FOUR_C_THROW("Wrong DiscretizationTypeGeometry given!");
  }
}


/**
 *
 */
template <typename ScalarType, typename Line, typename Surface>
void GEOMETRYPAIR::GeometryPairLineToSurfaceFADWrapper<ScalarType, Line, Surface>::pre_evaluate(
    const ElementData<Line, ScalarType>& element_data_line,
    const ElementData<Surface, ScalarType>& element_data_surface,
    std::vector<LineSegment<ScalarType>>& segments) const
{
  // Call pre_evaluate on the double pair.
  std::vector<LineSegment<double>> segments_double;
  geometry_pair_double_->pre_evaluate(ElementDataToDouble<Line>::to_double(element_data_line),
      ElementDataToDouble<Surface>::to_double(element_data_surface), segments_double);

  // Convert the created double segments to a segment of scalar type.
  segments.clear();
  for (auto& segment_double : segments_double)
  {
    // Create the segment with the scalar FAD type.
    segments.push_back(LineSegment<ScalarType>());
    copy_segment(segment_double, segments.back());
  }
}

/**
 *
 */
template <typename ScalarType, typename Line, typename Surface>
void GEOMETRYPAIR::GeometryPairLineToSurfaceFADWrapper<ScalarType, Line, Surface>::evaluate(
    const ElementData<Line, ScalarType>& element_data_line,
    const ElementData<Surface, ScalarType>& element_data_surface,
    std::vector<LineSegment<ScalarType>>& segments) const
{
  // Convert the input segments to a segment of scalar type double.
  std::vector<LineSegment<double>> segments_double;
  for (auto& segment : segments)
  {
    // Create the segment with the scalar FAD type.
    segments_double.push_back(LineSegment<double>());
    copy_segment(segment, segments_double.back());
  }

  // Call Evaluate on the double pair.
  geometry_pair_double_->evaluate(ElementDataToDouble<Line>::to_double(element_data_line),
      ElementDataToDouble<Surface>::to_double(element_data_surface), segments_double);

  // Get the face parameters.
  unsigned int n_faces;
  std::vector<unsigned int> face_fixed_parameters;
  std::vector<double> face_fixed_values;
  this->get_face_fixed_parameters(n_faces, face_fixed_parameters, face_fixed_values);

  // Initialize variables for the projections.
  ProjectionResult projection_result = ProjectionResult::none;
  Core::LinAlg::Matrix<3, 1, ScalarType> point_in_space;

  // If segments are found, convert them to FAD segments.
  segments.clear();
  for (auto& segment_double : segments_double)
  {
    // Create the segment with the scalar FAD type.
    segments.push_back(LineSegment<ScalarType>());
    LineSegment<ScalarType>& new_segment = segments.back();

    // Add the projection point to an array.
    std::array<std::reference_wrapper<ProjectionPoint1DTo3D<ScalarType>>, 2>
        segment_start_end_points = {new_segment.get_start_point(), new_segment.get_end_point()};
    segment_start_end_points[0].get().set_from_other_point_double(segment_double.get_start_point());
    segment_start_end_points[1].get().set_from_other_point_double(segment_double.get_end_point());

    // If the start or end points are intersection points, the intersections have to be reevaluated.
    for (auto& point : segment_start_end_points)
    {
      const int intersection_face = point.get().get_intersection_face();
      if (intersection_face >= 0)
      {
        this->intersect_line_with_surface_edge(element_data_line, element_data_surface,
            face_fixed_parameters[intersection_face], face_fixed_values[intersection_face],
            point.get().get_eta(), point.get().get_xi(), projection_result, true);
      }
    }

    // Reevaluate the integration points along the segment.
    std::vector<ProjectionPoint1DTo3D<ScalarType>>& projection_points =
        new_segment.get_projection_points();
    projection_points.resize(segment_double.get_number_of_projection_points());
    for (unsigned int i_point = 0; i_point < segment_double.get_number_of_projection_points();
        i_point++)
    {
      // Position of the projection point within the segment.
      auto& projection_point_double = segment_double.get_projection_points()[i_point];
      const double factor = (projection_point_double.get_eta() - segment_double.get_eta_a()) /
                            segment_double.get_segment_length();

      // Calculate spatial point.
      auto& projection_point = projection_points[i_point];
      projection_point.set_from_other_point_double(projection_point_double);
      projection_point.set_eta(new_segment.get_eta_a() + new_segment.get_segment_length() * factor);

      evaluate_position<Line>(projection_point.get_eta(), element_data_line, point_in_space);

      // Calculate the projection.
      this->project_point_to_other(
          point_in_space, element_data_surface, projection_point.get_xi(), projection_result, true);
    }
  }
}

/**
 *
 */
template <typename ScalarType, typename Surface>
void GEOMETRYPAIR::project_point_to_surface(const Core::LinAlg::Matrix<3, 1, ScalarType>& point,
    const ElementData<Surface, ScalarType>& element_data_surface,
    Core::LinAlg::Matrix<3, 1, ScalarType>& xi, ProjectionResult& projection_result,
    const double normal_influence_direction, const bool min_one_iteration)
{
  // Vectors in 3D.
  Core::LinAlg::Matrix<3, 1, ScalarType> r_surface;
  Core::LinAlg::Matrix<3, 1, ScalarType> delta_xi;
  // Initialize the increment with a value that will not pass the first convergence check.
  delta_xi.put_scalar(10 * Constants::projection_xi_eta_tol);
  Core::LinAlg::Matrix<3, 1, ScalarType> residuum;

  // Jacobian / inverse.
  Core::LinAlg::Matrix<3, 3, ScalarType> J_J_inv;

  // Reset the projection result flag.
  projection_result = ProjectionResult::projection_not_found;

  // Local Newton iteration.
  {
    unsigned int counter = 0;
    while (counter < Constants::local_newton_iter_max)
    {
      // Evaluate the position and its derivative on the surface.
      evaluate_surface_position_and_derivative(element_data_surface, xi, r_surface, J_J_inv);

      // Evaluate the residuum $r_{solid} - r_{point} = R_{pos}$.
      residuum = r_surface;
      residuum -= point;

      if (counter == 0 and min_one_iteration)
      {
        // if the min_one_iteration flag is set we run at least one iteration, so the dependency on
        // FAD variables is calculated correctly.
      }
      else if (Core::FADUtils::vector_norm(residuum) < Constants::local_newton_res_tol &&
               Core::FADUtils::vector_norm(delta_xi) < Constants::projection_xi_eta_tol)
      {
        if (valid_parameter_surface<ScalarType, Surface>(xi, normal_influence_direction))
          projection_result = ProjectionResult::projection_found_valid;
        else
          projection_result = ProjectionResult::projection_found_not_valid;
        break;
      }

      // Check if residuum is in a sensible range where we still expect to find a solution.
      if (Core::FADUtils::vector_norm(residuum) > Constants::local_newton_res_max) break;

      // Solve the linearized system.
      if (Core::LinAlg::solve_linear_system_do_not_throw_error_on_zero_determinant_scaled(
              J_J_inv, residuum, delta_xi, Constants::local_newton_det_tol))
      {
        // Set the new parameter coordinates.
        xi -= delta_xi;

        // Advance Newton iteration counter.
        counter++;
      }
      else
        break;
    }
  }
}


/**
 * Explicit template initialization of template class.
 */
namespace GEOMETRYPAIR
{
  template class GeometryPairLineToSurface<double, t_line2, t_tri3>;
  template class GeometryPairLineToSurface<double, t_line2, t_tri6>;
  template class GeometryPairLineToSurface<double, t_line2, t_quad4>;
  template class GeometryPairLineToSurface<double, t_line2, t_quad8>;
  template class GeometryPairLineToSurface<double, t_line2, t_quad9>;
  template class GeometryPairLineToSurface<double, t_line2, t_nurbs9>;

  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type_1st_order, t_line2,
      t_tri3>;
  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type_1st_order, t_line2,
      t_tri6>;
  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type_1st_order, t_line2,
      t_quad4>;
  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type_1st_order, t_line2,
      t_quad8>;
  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type_1st_order, t_line2,
      t_quad9>;
  template class GeometryPairLineToSurface<
      line_to_surface_patch_scalar_type_fixed_size_1st_order<t_line2, t_nurbs9>, t_line2, t_nurbs9>;

  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type, t_line2, t_tri3>;
  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type, t_line2, t_tri6>;
  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type, t_line2, t_quad4>;
  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type, t_line2, t_quad8>;
  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type, t_line2, t_quad9>;
  template class GeometryPairLineToSurface<
      line_to_surface_patch_scalar_type_fixed_size<t_line2, t_nurbs9>, t_line2, t_nurbs9>;

  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type_1st_order,
      t_line2, t_tri3>;
  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type_1st_order,
      t_line2, t_tri6>;
  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type_1st_order,
      t_line2, t_quad4>;
  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type_1st_order,
      t_line2, t_quad8>;
  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type_1st_order,
      t_line2, t_quad9>;
  template class GeometryPairLineToSurfaceFADWrapper<
      line_to_surface_patch_scalar_type_fixed_size_1st_order<t_line2, t_nurbs9>, t_line2, t_nurbs9>;

  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type, t_line2,
      t_tri3>;
  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type, t_line2,
      t_tri6>;
  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type, t_line2,
      t_quad4>;
  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type, t_line2,
      t_quad8>;
  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type, t_line2,
      t_quad9>;
  template class GeometryPairLineToSurfaceFADWrapper<
      line_to_surface_patch_scalar_type_fixed_size<t_line2, t_nurbs9>, t_line2, t_nurbs9>;

  template class GeometryPairLineToSurface<double, t_hermite, t_tri3>;
  template class GeometryPairLineToSurface<double, t_hermite, t_tri6>;
  template class GeometryPairLineToSurface<double, t_hermite, t_quad4>;
  template class GeometryPairLineToSurface<double, t_hermite, t_quad8>;
  template class GeometryPairLineToSurface<double, t_hermite, t_quad9>;
  template class GeometryPairLineToSurface<double, t_hermite, t_nurbs9>;

  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type_1st_order, t_hermite,
      t_tri3>;
  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type_1st_order, t_hermite,
      t_tri6>;
  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type_1st_order, t_hermite,
      t_quad4>;
  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type_1st_order, t_hermite,
      t_quad8>;
  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type_1st_order, t_hermite,
      t_quad9>;
  template class GeometryPairLineToSurface<
      line_to_surface_patch_scalar_type_fixed_size_1st_order<t_hermite, t_nurbs9>, t_hermite,
      t_nurbs9>;

  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type, t_hermite, t_tri3>;
  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type, t_hermite, t_tri6>;
  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type, t_hermite, t_quad4>;
  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type, t_hermite, t_quad8>;
  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type, t_hermite, t_quad9>;
  template class GeometryPairLineToSurface<
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_nurbs9>, t_hermite, t_nurbs9>;

  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type_1st_order,
      t_hermite, t_tri3>;
  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type_1st_order,
      t_hermite, t_tri6>;
  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type_1st_order,
      t_hermite, t_quad4>;
  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type_1st_order,
      t_hermite, t_quad8>;
  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type_1st_order,
      t_hermite, t_quad9>;
  template class GeometryPairLineToSurfaceFADWrapper<
      line_to_surface_patch_scalar_type_fixed_size_1st_order<t_hermite, t_nurbs9>, t_hermite,
      t_nurbs9>;

  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type, t_hermite,
      t_tri3>;
  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type, t_hermite,
      t_tri6>;
  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type, t_hermite,
      t_quad4>;
  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type, t_hermite,
      t_quad8>;
  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type, t_hermite,
      t_quad9>;
  template class GeometryPairLineToSurfaceFADWrapper<
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_nurbs9>, t_hermite, t_nurbs9>;
}  // namespace GEOMETRYPAIR

FOUR_C_NAMESPACE_CLOSE
