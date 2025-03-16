// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLID_3D_ELE_CALC_LIB_HPP
#define FOUR_C_SOLID_3D_ELE_CALC_LIB_HPP

#include "4C_config.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_cell_type.hpp"
#include "4C_fem_general_cell_type_traits.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_fem_general_element_dof_matrix.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_fem_general_fiber_node_holder.hpp"
#include "4C_fem_general_fiber_node_utils.hpp"
#include "4C_fem_general_utils_gauss_point_postprocess.hpp"
#include "4C_fem_general_utils_gausspoints.hpp"
#include "4C_fem_general_utils_nurbs_shapefunctions.hpp"
#include "4C_fem_nurbs_discretization_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_linalg_fixedsizematrix_generators.hpp"
#include "4C_linalg_fixedsizematrix_solver.hpp"
#include "4C_linalg_utils_densematrix_eigen.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_solid_3D_ele_calc_lib_integration.hpp"
#include "4C_utils_function.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Discret::Elements
{
  template <Core::FE::CellType celltype>
  struct ElementNodes;
}  // namespace Discret::Elements

namespace Discret::Elements::Internal
{
  template <Core::FE::CellType celltype>
  inline static constexpr int num_nodes = Core::FE::num_nodes<celltype>;

  template <Core::FE::CellType celltype>
  inline static constexpr int num_dim = Core::FE::dim<celltype>;

  template <Core::FE::CellType celltype>
  inline static constexpr int num_str = num_dim<celltype> * (num_dim<celltype> + 1) / 2;

  template <Core::FE::CellType celltype>
  inline static constexpr int num_dof_per_ele = num_nodes<celltype> * num_dim<celltype>;
}  // namespace Discret::Elements::Internal

namespace Discret::Elements
{

  /*!
   * @brief Calculate the lumped mass matrix
   */
  inline void lump_matrix(Core::LinAlg::SerialDenseMatrix& matrix)
  {
    FOUR_C_ASSERT(
        matrix.numRows() == matrix.numCols(), "The provided mass matrix is not a square matrix!");

    // we assume mass is a square matrix
    for (int c = 0; c < matrix.numCols(); ++c)  // parse columns
    {
      double d = 0.0;
      for (int r = 0; r < matrix.numRows(); ++r)  // parse rows
      {
        d += matrix(r, c);  // accumulate row entries
        matrix(r, c) = 0.0;
      }
      matrix(c, c) = d;  // apply sum of row entries on diagonal
    }
  }

  /*!
   * @brief A type holding information of the nodes of an element,
   * such as their reference and current position. Additional information
   * is stored for NURBS elements
   *
   * @tparam celltype
   */
  template <Core::FE::CellType celltype>
  struct ElementNodes
  {
    /*!
     * @brief Position of nodes in the reference configuration
     */
    Core::LinAlg::Matrix<Internal::num_dim<celltype>, Internal::num_nodes<celltype>>
        reference_coordinates;

    /*!
     * @brief Position of nodes in the current configuration
     */
    Core::LinAlg::Matrix<Internal::num_dim<celltype>, Internal::num_nodes<celltype>>
        current_coordinates;

    /*!
     * @brief Displacements of the element nodes
     */
    Core::LinAlg::Matrix<Internal::num_dim<celltype>, Internal::num_nodes<celltype>> displacements;
  };

  template <Core::FE::CellType celltype>
    requires(Core::FE::is_nurbs<celltype>)
  struct ElementNodes<celltype>
  {
    /*!
     * @brief Position of nodes in the reference configuration
     */
    Core::LinAlg::Matrix<Internal::num_dim<celltype>, Internal::num_nodes<celltype>>
        reference_coordinates;

    /*!
     * @brief Position of nodes in the current configuration
     */
    Core::LinAlg::Matrix<Internal::num_dim<celltype>, Internal::num_nodes<celltype>>
        current_coordinates;

    /*!
     * @brief Displacements of the element nodes
     */
    Core::LinAlg::Matrix<Internal::num_dim<celltype>, Internal::num_nodes<celltype>> displacements;

    /*!
     * @brief Knot span of a NURBS element
     */
    std::vector<Core::LinAlg::SerialDenseVector> knots;

    /*!
     * @brief Weights of control points
     */
    Core::LinAlg::Matrix<Internal::num_nodes<celltype>, 1, double> weights;
  };

  /*!
   * @brief Evaluate element node information given the element and a displacement vector
   *
   * @tparam celltype
   * @param ele (in): Element
   * @param disp (in) : Vector of nodal displacements of the element
   * @return ElementNodes<celltype>
   */
  template <Core::FE::CellType celltype, std::ranges::contiguous_range R>
    requires std::ranges::sized_range<R>
  ElementNodes<celltype> evaluate_element_nodes(const Core::Elements::Element& ele, const R& disp)
  {
    Discret::Elements::ElementNodes<celltype> element_nodes;
    for (int i = 0; i < Internal::num_nodes<celltype>; ++i)
    {
      for (int d = 0; d < Internal::num_dim<celltype>; ++d)
      {
        element_nodes.reference_coordinates(d, i) = ele.nodes()[i]->x()[d];
      }
    }
    element_nodes.displacements =
        Core::FE::get_element_dof_matrix<celltype, Core::FE::dim<celltype>>(disp);

    element_nodes.current_coordinates = element_nodes.reference_coordinates;
    element_nodes.current_coordinates.update(1.0, element_nodes.displacements, 1.0);

    return element_nodes;
  }

  /*!
   * @brief Evaluates the nodal coordinates from this iteration
   *
   * @param ele (in) : Reference to the element
   * @param discretization (in) : discretization
   * @param lm (in) : Location vector of the element, i.e., global dof numbers of elemental dofs
   */
  template <Core::FE::CellType celltype>
  ElementNodes<celltype> evaluate_element_nodes(const Core::Elements::Element& ele,
      const Core::FE::Discretization& discretization, const std::vector<int>& lm)
  {
    const Core::LinAlg::Vector<double>& displacements = *discretization.get_state("displacement");

    constexpr unsigned num_dofs = Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>;
    std::array<double, num_dofs> mydisp =
        Core::FE::extract_values_as_array<num_dofs>(displacements, lm);

    Discret::Elements::ElementNodes<celltype> element_nodes =
        evaluate_element_nodes<celltype>(ele, mydisp);

    if constexpr (Core::FE::is_nurbs<celltype>)
    {
      // Obtain the information required for a NURBS element
      bool zero_size = Core::FE::Nurbs::get_my_nurbs_knots_and_weights(
          discretization, &ele, element_nodes.knots, element_nodes.weights);
      if (zero_size)
        FOUR_C_THROW("get_my_nurbs_knots_and_weights has to return a non zero size NURBS element.");
    }

    return element_nodes;
  }

  /*!
   * @brief Evaluates the parameter coordinate of the Gauss point according the the Gauss rule
   *
   * @tparam celltype : Cell type
   * @param intpoints (in) : Gauss integration points
   * @param gp (in) : id of the Gauss point
   * @return Core::LinAlg::Matrix<num_dim<celltype>, 1> : Coordinates of the Gauss Point in the
   * parameter space
   */
  template <Core::FE::CellType celltype>
  Core::LinAlg::Matrix<Internal::num_dim<celltype>, 1> evaluate_parameter_coordinate(
      const Core::FE::GaussIntegration& intpoints, const int gp)
  {
    return Core::LinAlg::Matrix<Internal::num_dim<celltype>, 1>(intpoints.point(gp), false);
  }

  /*!
   * @brief Evaluates the parameter coordinate of the element centroid for Hexes
   *
   * Returns xi = [0 0 0].
   *
   * @tparam celltype : Cell type
   * @return Core::LinAlg::Matrix<Internal::num_dim<celltype>, 1> : Coordinates of the centroid in
   * the parameter space
   */
  template <Core::FE::CellType celltype>
  Core::LinAlg::Matrix<Internal::num_dim<celltype>, 1> evaluate_parameter_coordinate_centroid()
    requires(Core::FE::is_hex<celltype> || Core::FE::is_nurbs<celltype>)
  {
    Core::LinAlg::Matrix<Internal::num_dim<celltype>, 1> xi(true);
    return xi;
  }

  /*!
   * @brief Evaluates the parameter coordinate of the element centroid for Tets
   *
   * Returns xi = [0.25 0.25 0.25].
   *
   * @tparam celltype : Cell type
   * @return Core::LinAlg::Matrix<Internal::num_dim<celltype>, 1> : Coordinates of the centroid in
   * the parameter space
   */
  template <Core::FE::CellType celltype>
  Core::LinAlg::Matrix<Internal::num_dim<celltype>, 1> evaluate_parameter_coordinate_centroid()
    requires(Core::FE::is_tet<celltype>)
  {
    Core::LinAlg::Matrix<Internal::num_dim<celltype>, 1> xi;
    xi.put_scalar(0.25);

    return xi;
  }

  /*!
   * @brief Evaluates the parameter coordinate of the element centroid for Pyramids
   *
   * Returns xi = [0 0 0.25].
   *
   * @tparam celltype : Cell type
   * @return Core::LinAlg::Matrix<Internal::num_dim<celltype>, 1> : Coordinates of the centroid in
   * the parameter space
   */
  template <Core::FE::CellType celltype>
  Core::LinAlg::Matrix<Internal::num_dim<celltype>, 1> evaluate_parameter_coordinate_centroid()
    requires(Core::FE::is_pyramid<celltype>)
  {
    Core::LinAlg::Matrix<Internal::num_dim<celltype>, 1> xi(true);
    xi(2) = 0.25;

    return xi;
  }

  /*!
   * @brief Evaluates the parameter coordinate of the element centroid for Wedges
   *
   * Returns xi = [1/3 1/3 0].
   *
   * @tparam celltype : Cell type
   * @return Core::LinAlg::Matrix<Internal::num_dim<celltype>, 1> : Coordinates of the centroid in
   * the parameter space
   */
  template <Core::FE::CellType celltype>
  Core::LinAlg::Matrix<Internal::num_dim<celltype>, 1> evaluate_parameter_coordinate_centroid()
    requires(Core::FE::is_wedge<celltype>)
  {
    Core::LinAlg::Matrix<Internal::num_dim<celltype>, 1> xi(true);
    xi(0) = 1.0 / 3.0;
    xi(1) = 1.0 / 3.0;

    return xi;
  }

  /*!
   * @brief Evaluates a point's coordinates in reference configuration
   *
   * @tparam celltype : Cell type
   * @param nodal_coordinates_reference (in) : Reference coordinates of the nodes of the element
   * @param shape_functions_point (in) : Shape functions evaluated at the specific point
   * @return Core::LinAlg::Matrix<Core::FE::dim<celltype>, 1> : point's reference coordinates
   */
  template <Core::FE::CellType celltype>
  Core::LinAlg::Matrix<Core::FE::dim<celltype>, 1> evaluate_reference_coordinate(
      const Core::LinAlg::Matrix<Internal::num_dim<celltype>, Internal::num_nodes<celltype>>&
          nodal_coordinates_reference,
      const Core::LinAlg::Matrix<Internal::num_nodes<celltype>, 1>& shape_functions_point)
  {
    Core::LinAlg::Matrix<Internal::num_dim<celltype>, 1> coordinates_reference(true);
    coordinates_reference.multiply_nn(nodal_coordinates_reference, shape_functions_point);

    return coordinates_reference;
  }

  /*!
   * @brief Evaluates the element centroid's coordinates in reference configuration
   *
   * @tparam celltype : Cell type
   * @param nodal_coordinates (in) : Reference coordinates of the nodes of the element
   * @return Core::LinAlg::Matrix<Core::FE::dim<celltype>, 1> : Element centroid's coordinates in
   * reference configuration
   */
  template <Core::FE::CellType celltype>
  Core::LinAlg::Matrix<Core::FE::dim<celltype>, 1> evaluate_reference_coordinate_centroid(
      const ElementNodes<celltype>& nodal_coordinates)
    requires(Core::FE::use_lagrange_shapefnct<celltype>)
  {
    const Core::LinAlg::Matrix<Internal::num_dim<celltype>, 1> xi_centroid =
        evaluate_parameter_coordinate_centroid<celltype>();

    Core::LinAlg::Matrix<Internal::num_nodes<celltype>, 1> shape_functions_centroid(true);
    Core::FE::shape_function<celltype>(xi_centroid, shape_functions_centroid);

    return evaluate_reference_coordinate<celltype>(
        nodal_coordinates.reference_coordinates, shape_functions_centroid);
  }

  template <Core::FE::CellType celltype>
  Core::LinAlg::Matrix<Core::FE::dim<celltype>, 1> evaluate_reference_coordinate_centroid(
      const ElementNodes<celltype>& nodal_coordinates)
    requires(Core::FE::is_nurbs<celltype>)
  {
    const Core::LinAlg::Matrix<Internal::num_dim<celltype>, 1> xi_centroid =
        evaluate_parameter_coordinate_centroid<celltype>();

    Core::LinAlg::Matrix<Internal::num_nodes<celltype>, 1> shape_functions_centroid(true);
    Core::FE::Nurbs::nurbs_shape_function_dim(shape_functions_centroid, xi_centroid,
        nodal_coordinates.knots, nodal_coordinates.weights, celltype);

    return evaluate_reference_coordinate<celltype>(
        nodal_coordinates.reference_coordinates, shape_functions_centroid);
  }

  /*!
   * @brief Type holding the shape functions and it's first derivatives evaluated at a specific
   * point.
   *
   * @tparam celltype
   */
  template <Core::FE::CellType celltype>
  struct ShapeFunctionsAndDerivatives
  {
    Core::LinAlg::Matrix<Internal::num_nodes<celltype>, 1> shapefunctions_{};
    Core::LinAlg::Matrix<Internal::num_dim<celltype>, Internal::num_nodes<celltype>> derivatives_{};
  };

  /*!
   * @brief Evaluates the shape functions and their derivatives at the specified point in the
   * parameter space
   *
   * @tparam celltype : Discretizationt type
   * @param xi (in) : Coordinate in the parameter space
   * @return ShapeFunctionsAndDerivatives<celltype> : An object holding the shape functions and the
   * first derivatives evaluated at the respective point in the parameter space
   */
  template <Core::FE::CellType celltype>
  ShapeFunctionsAndDerivatives<celltype> evaluate_shape_functions_and_derivs(
      const Core::LinAlg::Matrix<Internal::num_dim<celltype>, 1>& xi,
      const ElementNodes<celltype>& nodal_coordinates)
    requires(Core::FE::use_lagrange_shapefnct<celltype>)
  {
    ShapeFunctionsAndDerivatives<celltype> shapefcns;
    Core::FE::shape_function<celltype>(xi, shapefcns.shapefunctions_);
    Core::FE::shape_function_deriv1<celltype>(xi, shapefcns.derivatives_);

    return shapefcns;
  }

  template <Core::FE::CellType celltype>
  ShapeFunctionsAndDerivatives<celltype> evaluate_shape_functions_and_derivs(
      const Core::LinAlg::Matrix<Internal::num_dim<celltype>, 1>& xi,
      const ElementNodes<celltype>& nodal_coordinates)
    requires(Core::FE::is_nurbs<celltype>)
  {
    ShapeFunctionsAndDerivatives<celltype> shapefcns;
    Core::FE::Nurbs::nurbs_get_funct_deriv(shapefcns.shapefunctions_, shapefcns.derivatives_, xi,
        nodal_coordinates.knots, nodal_coordinates.weights, celltype);

    return shapefcns;
  }

  template <Core::FE::CellType celltype>
  struct JacobianMapping
  {
    /// determinant of the jacobian
    double determinant_;

    /// Jacobian matrix at a specific point
    Core::LinAlg::Matrix<Internal::num_dim<celltype>, Internal::num_dim<celltype>> jacobian_;

    /// Inverse jacobian matrix at a specific point
    Core::LinAlg::Matrix<Internal::num_dim<celltype>, Internal::num_dim<celltype>>
        inverse_jacobian_;

    /// Derivative of the shape functions w.r.t. the reference coordinates
    Core::LinAlg::Matrix<Internal::num_dim<celltype>, Internal::num_nodes<celltype>> N_XYZ_;
  };

  /*!
   * @brief Evaluates the jacobian mapping of the element
   *
   * @tparam celltype : Cell type
   * @param shapefcns (in) : Shape functions and derivatives evaluated at the respective point in
   * the parameter space
   * @param nodal_coordinates (in) : Reference and current coordinates of the nodes of the element
   * @param gp (in) : Id of the Gauss point
   * @return JacobianMapping<celltype> : An object holding quantities of the jacobian mapping
   * (inverse Jacobian, determinant, derivatives of the shape functions w.r.t. XYZ, integration
   * factor)
   */
  template <Core::FE::CellType celltype>
  JacobianMapping<celltype> evaluate_jacobian_mapping(
      const ShapeFunctionsAndDerivatives<celltype>& shapefcns,
      const ElementNodes<celltype>& nodal_coordinates)
  {
    JacobianMapping<celltype> jacobian;

    jacobian.jacobian_.multiply_nt(shapefcns.derivatives_, nodal_coordinates.reference_coordinates);
    jacobian.inverse_jacobian_ = jacobian.jacobian_;
    jacobian.determinant_ = jacobian.inverse_jacobian_.invert();
    jacobian.N_XYZ_.multiply(jacobian.inverse_jacobian_, shapefcns.derivatives_);

    return jacobian;
  }

  /*!
   * @brief Evaluates the jacobian determinant of the element
   *
   * @tparam celltype : Cell type
   * @param shapefcns (in) : Shape functions and derivatives evaluated at the respective point in
   * the parameter space
   * @param nodal_coordinates (in) : Reference and current coordinates of the nodes of the element
   * @return double : Jacobian determinant
   */
  template <Core::FE::CellType celltype>
  double evaluate_jacobian_determinant(const ShapeFunctionsAndDerivatives<celltype>& shapefcns,
      const ElementNodes<celltype>& nodal_coordinates)
  {
    Core::LinAlg::Matrix<Internal::num_dim<celltype>, Internal::num_dim<celltype>> jacobian;
    jacobian.multiply_nt(shapefcns.derivatives_, nodal_coordinates.reference_coordinates);

    return jacobian.determinant();
  }

  /*!
   * @brief Check for negative Jacobian determinants
   *
   * @tparam celltype : Cell type
   * @param element_nodes (in) : Element node information
   */
  template <Core::FE::CellType celltype>
  void ensure_positive_jacobian_determinant_at_element_nodes(
      const ElementNodes<celltype>& element_nodes)
  {
    const Core::LinAlg::SerialDenseMatrix rst =
        Core::FE::get_ele_node_numbering_nodes_paramspace(celltype);

    for (const auto& xi : Core::FE::get_element_nodes_in_parameter_space<celltype>())
    {
      Core::LinAlg::Matrix<3, 1> xi_mat(xi.data(), true);

      ShapeFunctionsAndDerivatives<celltype> shape_functions =
          evaluate_shape_functions_and_derivs(xi_mat, element_nodes);

      JacobianMapping<celltype> jacobian_mapping =
          evaluate_jacobian_mapping(shape_functions, element_nodes);

      FOUR_C_ASSERT_ALWAYS(jacobian_mapping.determinant_ > 0,
          "determinant of jacobian is {} <= 0 at one node of the element.",
          jacobian_mapping.determinant_);
    }
  }

  /*!
   * @brief Evaluate the jacobian mapping at the element centroid
   *
   * @tparam celltype : Cell type
   * @param nodal_coordinates (in) : Reference and current coordinates of the nodes of the element
   * @return JacobianMapping<celltype> : jacobian mapping at the element centroid
   */
  template <Core::FE::CellType celltype>
  JacobianMapping<celltype> evaluate_jacobian_mapping_centroid(
      const ElementNodes<celltype>& nodal_coordinates)
    requires(Internal::num_dim<celltype> == 3)
  {
    // set coordinates in parameter space at centroid as zero -> xi = [0; 0; 0]
    Core::LinAlg::Matrix<Internal::num_dim<celltype>, 1> xi_centroid =
        evaluate_parameter_coordinate_centroid<celltype>();

    // shape functions and derivatives evaluated at element centroid
    const ShapeFunctionsAndDerivatives<celltype> shape_functions_centroid =
        evaluate_shape_functions_and_derivs<celltype>(xi_centroid, nodal_coordinates);

    // jacobian mapping evaluated at centroid
    const JacobianMapping<celltype> jacobian_mapping_centroid =
        evaluate_jacobian_mapping(shape_functions_centroid, nodal_coordinates);

    return jacobian_mapping_centroid;
  }

  template <Core::FE::CellType celltype>
  struct SpatialMaterialMapping
  {
    double determinant_deformation_gradient_{};
    Core::LinAlg::Matrix<Internal::num_dim<celltype>, Internal::num_dim<celltype>>
        deformation_gradient_{};
    Core::LinAlg::Matrix<Internal::num_dim<celltype>, Internal::num_dim<celltype>>
        inverse_deformation_gradient_{};
  };

  /*!
   * @brief Evaluates the mapping between the spatial and material configuration of the element
   *
   * @tparam celltype : Cell type
   * @param jacobian_mapping (in) : An object holding quantities of the jacobian mapping
   * (inverse Jacobian, determinant, derivatives of the shape functions w.r.t. XYZ)
   * @param nodal_coordinates (in) : Reference and current coordinates of the nodes of the element
   * @param scale_defgrd (in) : scaling for deformation gradient
   * @param kinematictype (in) : kinematic type of element
   * @return SpatialMaterialMapping<celltype> : An object holding quantities of the spatial material
   * mapping (deformation_gradient, inverse_deformation_gradient,
   * determinant_deformation_gradient)
   */
  template <Core::FE::CellType celltype>
  SpatialMaterialMapping<celltype> evaluate_spatial_material_mapping(
      const JacobianMapping<celltype>& jacobian_mapping,
      const ElementNodes<celltype>& nodal_coordinates, const double scale_defgrd = 1.0,
      const Inpar::Solid::KinemType& kinematictype = Inpar::Solid::KinemType::nonlinearTotLag)
  {
    SpatialMaterialMapping<celltype> spatial_material_mapping;

    if (kinematictype == Inpar::Solid::KinemType::nonlinearTotLag)
    {
      spatial_material_mapping.deformation_gradient_ =
          evaluate_deformation_gradient(jacobian_mapping, nodal_coordinates, scale_defgrd);
    }

    spatial_material_mapping.inverse_deformation_gradient_.invert(
        spatial_material_mapping.deformation_gradient_);
    spatial_material_mapping.determinant_deformation_gradient_ =
        spatial_material_mapping.deformation_gradient_.determinant();

    return spatial_material_mapping;
  }

  /*!
   * @brief Evaluate the deformation gradient given the displacements and the jacobian mapping
   *
   * @tparam celltype
   * @param jacobian_mapping
   * @param element_nodes
   * @param scale_defgrd
   * @return Core::LinAlg::Matrix<3, 3>
   */
  template <Core::FE::CellType celltype>
  Core::LinAlg::Matrix<3, 3> evaluate_deformation_gradient(
      const JacobianMapping<celltype>& jacobian_mapping,
      const ElementNodes<celltype>& element_nodes, const double scale_defgrd = 1.0)
  {
    Core::LinAlg::Matrix<3, 3> defgrd(false);
    if constexpr (celltype == Core::FE::CellType::hex8)
    {
      // For some reason, some contact tests with hex8 discretization don't like the computation
      // of the deformation gradient based on the nodal displacements. They only converge if the
      // deformation gradient is computed based on the current coordinates. Similarly, there are
      // ssi and ssti tests with an tet4 discretization that only converge when using the
      // displacements. Until we found the problem, we compute the deformation gradient based on
      // the current coordinates (F=(X+u)^T  dN/dX^T) for hex8 and based on the displacement (F=I
      // + u^T dN/dX^T) for the other celltypes.
      defgrd.multiply_nt(scale_defgrd, element_nodes.current_coordinates, jacobian_mapping.N_XYZ_);
    }
    else
    {
      defgrd = Core::LinAlg::identity_matrix<Core::FE::dim<celltype>>();

      defgrd.multiply_nt(
          scale_defgrd, element_nodes.displacements, jacobian_mapping.N_XYZ_, scale_defgrd);
    }

    return defgrd;
  }

  /*!
   * @brief Evaluate the determinant of the deformation gradient at the element centroid
   *
   * @tparam celltype : Cell type
   * @param nodal_coordinates (in) : Reference and current coordinates of the nodes of the element
   * @return double : determinant of the deformation gradient at the centroid
   */
  template <Core::FE::CellType celltype>
  double evaluate_deformation_gradient_determinant_centroid(
      const ElementNodes<celltype>& nodal_coordinates)
    requires(Internal::num_dim<celltype> == 3)
  {
    // jacobian mapping at centroid of element
    const JacobianMapping<celltype> jacobian_mapping_centroid =
        evaluate_jacobian_mapping_centroid(nodal_coordinates);

    // deformation gradient and strains at centroid of element
    const Discret::Elements::SpatialMaterialMapping<celltype> spatial_material_mapping_centroid =
        evaluate_spatial_material_mapping(jacobian_mapping_centroid, nodal_coordinates);

    return spatial_material_mapping_centroid.determinant_deformation_gradient_;
  }

  /*!
   * @brief Evaluates Green-Lagrange strain from right Cauchy-Green tensor
   *
   * GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
   *
   * @param cauchygreen (in) : Right Cauchy-Green deformation tensor
   * @return Core::LinAlg::Matrix<Internal::num_str<celltype>, 1> : Green-Lagrange strain tensor in
   * strain-like Voigt notation
   */
  inline Core::LinAlg::Matrix<6, 1> evaluate_green_lagrange_strain(
      const Core::LinAlg::Matrix<3, 3>& cauchygreen)
  {
    Core::LinAlg::Matrix<6, 1> gl_strain;

    gl_strain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
    gl_strain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
    gl_strain(2) = 0.5 * (cauchygreen(2, 2) - 1.0);
    gl_strain(3) = cauchygreen(0, 1);
    gl_strain(4) = cauchygreen(1, 2);
    gl_strain(5) = cauchygreen(2, 0);

    return gl_strain;
  }

  /*!
   * @brief Evaluates right Cauchy-Green deformation tensor
   *
   * @tparam celltype: Cell type
   * @param spatial_material_mapping (in) : An object holding quantities of the spatial material
   * mapping (deformation_gradient, inverse_deformation_gradient,
   * determinant_deformation_gradient)
   * @return Core::LinAlg::Matrix<Internal::num_dim<celltype>, Internal::num_dim<celltype>> : Right
   * Cauchy-Green deformation tensor
   */
  template <Core::FE::CellType celltype>
  Core::LinAlg::Matrix<Internal::num_dim<celltype>, Internal::num_dim<celltype>>
  evaluate_cauchy_green(const SpatialMaterialMapping<celltype>& spatial_material_mapping)
    requires(Internal::num_dim<celltype> == 3)
  {
    Core::LinAlg::Matrix<Internal::num_dim<celltype>, Internal::num_dim<celltype>> cauchygreen(
        false);

    cauchygreen.multiply_tn(spatial_material_mapping.deformation_gradient_,
        spatial_material_mapping.deformation_gradient_);

    return cauchygreen;
  }


  /*!
   * @brief Evaluate the Green-Lagrange strain tensor from small displacement assumptions
   */
  template <Core::FE::CellType celltype>
  Core::LinAlg::Matrix<Internal::num_str<celltype>, 1> evaluate_linear_gl_strain(
      const ElementNodes<celltype>& nodal_coordinates,
      const Core::LinAlg::Matrix<Internal::num_str<celltype>,
          Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>>& linear_b_operator)
  {
    Core::LinAlg::Matrix<Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>, 1u> nodal_displs =
        Core::FE::get_element_dof_vector_view<celltype>(nodal_coordinates.displacements);

    Core::LinAlg::Matrix<Internal::num_str<celltype>, 1> gl_strain;
    gl_strain.multiply(linear_b_operator, nodal_displs);

    return gl_strain;
  }

  /*!
   * @brief Evaluates the strain gradient (B-Operator) of the specified element
   *
   * @tparam celltype : Cell type
   * @param jacobian_mapping (in) : Quantities of the jacobian mapping
   * @param spatial_material_mapping (in) :An object holding quantities of the spatial material
   * mapping (deformation_gradient, inverse_deformation_gradient,
   * determinant_deformation_gradient)
   * @return Core::LinAlg::Matrix<num_str<celltype>, num_dim<celltype> * num_nodes<celltype>> :
   * B-Operator
   */
  template <Core::FE::CellType celltype>
  Core::LinAlg::Matrix<Internal::num_str<celltype>,
      Internal::num_dim<celltype> * Internal::num_nodes<celltype>>
  evaluate_strain_gradient(const JacobianMapping<celltype>& jacobian_mapping,
      const SpatialMaterialMapping<celltype>& spatial_material_mapping)
    requires(Internal::num_dim<celltype> == 3)
  {
    // B-operator
    Core::LinAlg::Matrix<Internal::num_str<celltype>,
        Internal::num_dim<celltype> * Internal::num_nodes<celltype>>
        Bop;
    for (int i = 0; i < Internal::num_nodes<celltype>; ++i)
    {
      for (int d = 0; d < Internal::num_dim<celltype>; ++d)
      {
        for (int e = 0; e < Internal::num_dim<celltype>; ++e)
        {
          Bop(d, Internal::num_dim<celltype> * i + e) =
              spatial_material_mapping.deformation_gradient_(e, d) * jacobian_mapping.N_XYZ_(d, i);
        }
      }

      Bop(3, Internal::num_dim<celltype> * i + 0) =
          spatial_material_mapping.deformation_gradient_(0, 0) * jacobian_mapping.N_XYZ_(1, i) +
          spatial_material_mapping.deformation_gradient_(0, 1) * jacobian_mapping.N_XYZ_(0, i);
      Bop(3, Internal::num_dim<celltype> * i + 1) =
          spatial_material_mapping.deformation_gradient_(1, 0) * jacobian_mapping.N_XYZ_(1, i) +
          spatial_material_mapping.deformation_gradient_(1, 1) * jacobian_mapping.N_XYZ_(0, i);
      Bop(3, Internal::num_dim<celltype> * i + 2) =
          spatial_material_mapping.deformation_gradient_(2, 0) * jacobian_mapping.N_XYZ_(1, i) +
          spatial_material_mapping.deformation_gradient_(2, 1) * jacobian_mapping.N_XYZ_(0, i);
      Bop(4, Internal::num_dim<celltype> * i + 0) =
          spatial_material_mapping.deformation_gradient_(0, 1) * jacobian_mapping.N_XYZ_(2, i) +
          spatial_material_mapping.deformation_gradient_(0, 2) * jacobian_mapping.N_XYZ_(1, i);
      Bop(4, Internal::num_dim<celltype> * i + 1) =
          spatial_material_mapping.deformation_gradient_(1, 1) * jacobian_mapping.N_XYZ_(2, i) +
          spatial_material_mapping.deformation_gradient_(1, 2) * jacobian_mapping.N_XYZ_(1, i);
      Bop(4, Internal::num_dim<celltype> * i + 2) =
          spatial_material_mapping.deformation_gradient_(2, 1) * jacobian_mapping.N_XYZ_(2, i) +
          spatial_material_mapping.deformation_gradient_(2, 2) * jacobian_mapping.N_XYZ_(1, i);
      Bop(5, Internal::num_dim<celltype> * i + 0) =
          spatial_material_mapping.deformation_gradient_(0, 2) * jacobian_mapping.N_XYZ_(0, i) +
          spatial_material_mapping.deformation_gradient_(0, 0) * jacobian_mapping.N_XYZ_(2, i);
      Bop(5, Internal::num_dim<celltype> * i + 1) =
          spatial_material_mapping.deformation_gradient_(1, 2) * jacobian_mapping.N_XYZ_(0, i) +
          spatial_material_mapping.deformation_gradient_(1, 0) * jacobian_mapping.N_XYZ_(2, i);
      Bop(5, Internal::num_dim<celltype> * i + 2) =
          spatial_material_mapping.deformation_gradient_(2, 2) * jacobian_mapping.N_XYZ_(0, i) +
          spatial_material_mapping.deformation_gradient_(2, 0) * jacobian_mapping.N_XYZ_(2, i);
    }

    return Bop;
  }

  /*!
   * @brief Evaluates the linear strain gradient (B-operator) for small displacements
   *
   * @tparam celltype
   * @param jacobian_mapping (in) : Quantities of the jacobian mapping
   * @return Core::LinAlg::Matrix<Internal::num_str<celltype>,
   * Internal::num_dim<celltype> * Internal::num_nodes<celltype>>
   */
  template <Core::FE::CellType celltype>
  Core::LinAlg::Matrix<Internal::num_str<celltype>,
      Internal::num_dim<celltype> * Internal::num_nodes<celltype>>
  evaluate_linear_strain_gradient(const JacobianMapping<celltype>& jacobian_mapping)
    requires(Internal::num_dim<celltype> == 3)
  {
    // B-operator
    Core::LinAlg::Matrix<Internal::num_str<celltype>,
        Internal::num_dim<celltype> * Internal::num_nodes<celltype>>
        Bop;
    for (int i = 0; i < Internal::num_nodes<celltype>; ++i)
    {
      for (int d = 0; d < Internal::num_dim<celltype>; ++d)
      {
        Bop(d, Internal::num_dim<celltype> * i + d) = jacobian_mapping.N_XYZ_(d, i);
      }

      Bop(3, Internal::num_dim<celltype> * i + 0) = jacobian_mapping.N_XYZ_(1, i);
      Bop(3, Internal::num_dim<celltype> * i + 1) = jacobian_mapping.N_XYZ_(0, i);
      Bop(3, Internal::num_dim<celltype> * i + 2) = 0;
      Bop(4, Internal::num_dim<celltype> * i + 0) = 0;
      Bop(4, Internal::num_dim<celltype> * i + 1) = jacobian_mapping.N_XYZ_(2, i);
      Bop(4, Internal::num_dim<celltype> * i + 2) = jacobian_mapping.N_XYZ_(1, i);
      Bop(5, Internal::num_dim<celltype> * i + 0) = jacobian_mapping.N_XYZ_(2, i);
      Bop(5, Internal::num_dim<celltype> * i + 1) = 0;
      Bop(5, Internal::num_dim<celltype> * i + 2) = jacobian_mapping.N_XYZ_(0, i);
    }

    return Bop;
  }

  template <Core::FE::CellType celltype>
  struct Stress
  {
    /// Second Piola-Kirchhoff stress tensor in stress-like voigt notation
    Core::LinAlg::Matrix<Internal::num_str<celltype>, 1> pk2_;

    /// Linearization of the 2. Piola Kirchhoff stress tensor w.r.t. Green-Lagrange strain tensor in
    /// mixed Voigt notation
    Core::LinAlg::Matrix<Internal::num_str<celltype>, Internal::num_str<celltype>> cmat_;
  };

  /*!
   * @brief Evaluates the material stress (2. Piola-Kirchhoff stress tensor and the linearization
   * w.r.t. Green-Lagrange strain)
   *
   * @tparam celltype : Cell type
   * @param material (in) : Reference to the material
   * @param spatial_material_mapping (in) : An object holding quantities of the spatial material
   * mapping (deformation_gradient, inverse_deformation_gradient,
   * determinant_deformation_gradient)
   * @param gl_strain (in) : Green-Lagrange strain
   * @param params (in) : List of additional parameter to pass quantities from the time integrator
   * to the material
   * @param gp (in) : Gauss point
   * @param eleGID (in) : Global element id
   * @return Stress<celltype> : Object holding the 2. Piola-Kirchhoff stress tensor and the
   * linearization w.r.t. Green-Lagrange strain tensor
   */
  template <Core::FE::CellType celltype>
  Stress<celltype> evaluate_material_stress(Mat::So3Material& material,
      const Core::LinAlg::Matrix<Internal::num_dim<celltype>, Internal::num_dim<celltype>>& defgrd,
      const Core::LinAlg::Matrix<Internal::num_str<celltype>, 1>& gl_strain,
      Teuchos::ParameterList& params, const int gp, const int eleGID)
  {
    Stress<celltype> stress;

    material.evaluate(&defgrd, &gl_strain, params, &stress.pk2_, &stress.cmat_, gp, eleGID);
    return stress;
  }

  /*!
   * @brief Adds the internal force vector contribution of one Gauss point
   *
   * @tparam celltype : Cell type
   * @param Bop (in) : Strain gradient (B-Operator)
   * @param stress (in) : Stress measures
   * @param integration_fac (in) : Integration factor (Gauss point weight times the determinant of
   * the jacobian)
   * @param force_vector (in/out) : Force vector where the local contribution is added to
   */
  template <Core::FE::CellType celltype>
  void add_internal_force_vector(
      const Core::LinAlg::Matrix<Internal::num_str<celltype>,
          Internal::num_dim<celltype> * Internal::num_nodes<celltype>>& Bop,
      const Stress<celltype>& stress, const double integration_fac,
      Core::LinAlg::Matrix<Internal::num_dim<celltype> * Internal::num_nodes<celltype>, 1>&
          force_vector)
  {
    force_vector.multiply_tn(integration_fac, Bop, stress.pk2_, 1.);
  }

  /*!
   * @brief Add elastic stiffness matrix contribution of one Gauss point
   *
   * @tparam celltype : Cell type
   * @param Bop (in) : Strain gradient (B-Operator)
   * @param stress (in) : Stress measures
   * @param integration_fac (in) : Integration factor (Gauss point weight times the determinant of
   * the jacobian)
   * @param stiffness_matrix (in/out) : stiffness matrix where the local contribution is added to
   */
  template <Core::FE::CellType celltype>
  void add_elastic_stiffness_matrix(
      const Core::LinAlg::Matrix<Internal::num_str<celltype>,
          Internal::num_dim<celltype> * Internal::num_nodes<celltype>>& Bop,
      const Stress<celltype>& stress, const double integration_fac,
      Core::LinAlg::Matrix<Internal::num_dim<celltype> * Internal::num_nodes<celltype>,
          Internal::num_dim<celltype> * Internal::num_nodes<celltype>>& stiffness_matrix)
  {
    Core::LinAlg::Matrix<Internal::num_str<celltype>,
        Internal::num_nodes<celltype> * Internal::num_dim<celltype>>
        cb;
    cb.multiply(stress.cmat_, Bop);
    stiffness_matrix.multiply_tn(integration_fac, Bop, cb, 1.0);
  }

  /*!
   * @brief Add geometric stiffness matrix contribution of one Gauss point
   *
   * @tparam celltype : Cell type
   * @param B_L (in) : B_L operator, i.e. derivatives of the shape functions w.r.t. XYZ
   *                   at the respective Gauss point
   * @param stress (in) : Stress measures
   * @param integration_fac (in) : Integration factor (Gauss point weight times the determinant of
   * the jacobian)
   * @param stiffness_matrix (in/out) : stiffness matrix where the local contribution is added to
   */
  template <Core::FE::CellType celltype>
  void add_geometric_stiffness_matrix(
      const Core::LinAlg::Matrix<Internal::num_dim<celltype>, Internal::num_nodes<celltype>>& B_L,
      const Stress<celltype>& stress, const double integration_fac,
      Core::LinAlg::Matrix<Internal::num_dim<celltype> * Internal::num_nodes<celltype>,
          Internal::num_dim<celltype> * Internal::num_nodes<celltype>>& stiffness_matrix)
  {
    std::array<double, 3> SmB_L;  // intermediate Sm.B_L
    // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
    for (int inod = 0; inod < Internal::num_nodes<celltype>; ++inod)
    {
      SmB_L[0] = stress.pk2_(0) * B_L(0, inod) + stress.pk2_(3) * B_L(1, inod) +
                 stress.pk2_(5) * B_L(2, inod);
      SmB_L[1] = stress.pk2_(3) * B_L(0, inod) + stress.pk2_(1) * B_L(1, inod) +
                 stress.pk2_(4) * B_L(2, inod);
      SmB_L[2] = stress.pk2_(5) * B_L(0, inod) + stress.pk2_(4) * B_L(1, inod) +
                 stress.pk2_(2) * B_L(2, inod);

      for (int jnod = 0; jnod < Internal::num_nodes<celltype>; ++jnod)
      {
        double bopstrbop = 0.0;  // intermediate value
        for (int idim = 0; idim < Internal::num_dim<celltype>; ++idim)
          bopstrbop += B_L(idim, jnod) * SmB_L[idim];

        for (int d = 0; d < Internal::num_dim<celltype>; ++d)
          stiffness_matrix(Internal::num_dim<celltype> * inod + d,
              Internal::num_dim<celltype> * jnod + d) += integration_fac * bopstrbop;
      }
    }
  }

  /*!
   * @brief Add mass matrix contribution of one Gauss point
   *
   * @tparam celltype : Cell type
   * @param shapefunctions (in) : Shape functions and derivatives evaluated at the respective point
   * in the parameter space
   * @param integration_factor (in) : Integration factor (Gauss point weight times the determinant
   * of the jacobian)
   * @param density (in) : density at the Gauss point
   * @param mass (in/out) : mass matrix where the local contribution is added to
   */
  template <Core::FE::CellType celltype>
  void add_mass_matrix(const ShapeFunctionsAndDerivatives<celltype>& shapefunctions,
      const double integration_factor, const double density,
      Core::LinAlg::Matrix<Internal::num_dim<celltype> * Internal::num_nodes<celltype>,
          Internal::num_dim<celltype> * Internal::num_nodes<celltype>>& mass)
  {
    for (int inod = 0; inod < Internal::num_nodes<celltype>; ++inod)
    {
      const double ifactor = shapefunctions.shapefunctions_(inod) * integration_factor * density;
      for (int jnod = 0; jnod < Internal::num_nodes<celltype>; ++jnod)
      {
        const double massfactor =
            shapefunctions.shapefunctions_(jnod) * ifactor;  // intermediate factor
        for (int d = 0; d < Internal::num_dim<celltype>; ++d)
          mass(Internal::num_dim<celltype> * inod + d, Internal::num_dim<celltype> * jnod + d) +=
              massfactor;
      }
    }
  }

  /*!
   * @brief Calls the @p gp_evaluator for each Gauss point with evaluated jacobian mapping using the
   * integration rule defined by @p integration.
   *
   * @tparam celltype : Cell type known at compile time
   * @tparam GaussPointEvaluator
   * @param nodal_coordinates (in) : The nodal coordinates of the element
   * @param integration (in) : The integration rule to be used.
   * @param gp_evaluator (in) : A callable object (e.g. lambda-function) with signature void(const
   * Core::LinAlg::Matrix<Internal::num_dim<celltype>, 1>& xi, const
   * ShapeFunctionsAndDerivatives<celltype>& shape_functions, const JacobianMapping<celltype>&
   * jacobian_mapping, double integration_factor, int gp) that will be called for each integration
   * point.
   */
  template <Core::FE::CellType celltype, typename GaussPointEvaluator>
  inline void for_each_gauss_point(const ElementNodes<celltype>& nodal_coordinates,
      const Core::FE::GaussIntegration& integration, GaussPointEvaluator gp_evaluator)
  {
    for (int gp = 0; gp < integration.num_points(); ++gp)
    {
      const Core::LinAlg::Matrix<Internal::num_dim<celltype>, 1> xi =
          evaluate_parameter_coordinate<celltype>(integration, gp);

      const ShapeFunctionsAndDerivatives<celltype> shape_functions =
          evaluate_shape_functions_and_derivs<celltype>(xi, nodal_coordinates);

      const JacobianMapping<celltype> jacobian_mapping =
          evaluate_jacobian_mapping(shape_functions, nodal_coordinates);

      const double integration_factor = jacobian_mapping.determinant_ * integration.weight(gp);

      gp_evaluator(xi, shape_functions, jacobian_mapping, integration_factor, gp);
    }
  }

  // create a struct to store the error computation components
  struct AnalyticalDisplacementErrorIntegrationResults
  {
    double integrated_squared_error = 0;
    double integrated_squared_displacements = 0;
    double integrated_volume = 0;
  };

  /*!
   * @brief compute displacement error and displacement integral
   *
   * @tparam celltype : Cell type
   * @param element  : element
   * @param discretization : discretization object
   * @param lm  : location vector
   * @param analytical_displacements_function  : reference to function describing the analytical
   * solution
   */

  template <Core::FE::CellType celltype>
  AnalyticalDisplacementErrorIntegrationResults compute_analytical_displacement_error_integration(
      Core::Elements::Element& element, const Core::FE::Discretization& discretization,
      const std::vector<int>& lm,
      const Core::Utils::FunctionOfSpaceTime& analytical_displacements_function)
  {
    constexpr auto numdim = Core::FE::dim<celltype>;
    const ElementNodes<celltype> nodal_coordinates =
        evaluate_element_nodes<celltype>(element, discretization, lm);
    Core::FE::GaussIntegration gauss_integration = Core::FE::create_gauss_integration<celltype>(
        Discret::Elements::get_gauss_rule_stiffness_matrix<celltype>());

    AnalyticalDisplacementErrorIntegrationResults error_result;
    Discret::Elements::for_each_gauss_point<celltype>(nodal_coordinates, gauss_integration,
        [&](const Core::LinAlg::Matrix<Internal::num_dim<celltype>, 1>& xi,
            const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
            const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp)
        {
          Core::LinAlg::Matrix<numdim, 1> gauss_point_reference_coordinates;
          gauss_point_reference_coordinates.multiply_nn(
              nodal_coordinates.reference_coordinates, shape_functions.shapefunctions_);

          Core::LinAlg::Matrix<numdim, 1> gauss_point_disp;
          gauss_point_disp.multiply(
              nodal_coordinates.displacements, shape_functions.shapefunctions_);

          Core::LinAlg::Matrix<3, 1> analytical_solution;
          // Calculate the analytical solution for each dimension
          for (int i_dim = 0; i_dim < 3; i_dim++)
          {
            analytical_solution(i_dim) = analytical_displacements_function.evaluate(
                gauss_point_reference_coordinates.data(), 0.0, i_dim);
          }

          // The data that will be added to the element_force_vector will be summed up for all
          // elements
          Core::LinAlg::Matrix<3, 1> error_pointwise = analytical_solution;
          error_pointwise -= gauss_point_disp;
          // Calculate the error squared integral
          error_result.integrated_squared_error +=
              error_pointwise.dot(error_pointwise) * integration_factor;
          // Calculate the L2 norm of displacements
          error_result.integrated_squared_displacements +=
              gauss_point_disp.dot(gauss_point_disp) * integration_factor;
          // Calculate the domain integral
          error_result.integrated_volume += integration_factor;
        });
    return error_result;
  }

  /*!
   * @brief Evaluates the Gauss point coordinates in reference configuration
   * and adds those to the parameter list
   *
   * @tparam celltype : Cell type
   * @param nodal_coordinates (in) : Reference and current coordinates of the nodes of the element
   * @param shape_functions_gp (in) : Shape functions evaluated at the Gauss point
   * @param params (in/out) : ParameterList the quantities are added to
   */
  template <Core::FE::CellType celltype>
  void evaluate_gp_coordinates_and_add_to_parameter_list(
      const ElementNodes<celltype>& nodal_coordinates,
      const ShapeFunctionsAndDerivatives<celltype>& shape_functions_gp,
      Teuchos::ParameterList& params)
  {
    auto gp_ref_coord = evaluate_reference_coordinate<celltype>(
        nodal_coordinates.reference_coordinates, shape_functions_gp.shapefunctions_);
    params.set("gp_coords_ref", gp_ref_coord);
  }

  /*!
   * @brief Evaluates the element centroid coordinates in reference configuration
   * and adds those to the parameter list
   *
   * @tparam celltype : Cell type
   * @param nodal_coordinates (in) : Reference and current coordinates of the nodes of the element
   * @param params (in/out) : ParameterList the quantities are added to
   */
  template <Core::FE::CellType celltype>
  void evaluate_centroid_coordinates_and_add_to_parameter_list(
      const ElementNodes<celltype>& nodal_coordinates, Teuchos::ParameterList& params)
  {
    auto element_center = evaluate_reference_coordinate_centroid<celltype>(nodal_coordinates);
    params.set("elecenter_coords_ref", element_center);
  }

  /*!
   * @brief For elements with fiber nodes, interpolate fibers to Gauss points and add to the
   * parameter list
   *
   * @tparam celltype : Cell type
   * @param stiffness_matrix_integration (in) : Container holding the integration points
   * @param ele (in) : Reference to the element, possibly having nodal fibers
   * @param params (in/out) : ParameterList the interpolated fibers are added to
   */
  template <Core::FE::CellType celltype>
  inline void interpolate_fibers_to_gauss_points_and_add_to_parameter_list(
      const Core::FE::GaussIntegration& stiffness_matrix_integration,
      const Core::Elements::Element& ele, Teuchos::ParameterList& params)
  {
    if (Core::Nodes::have_nodal_fibers<celltype>(ele.nodes()))
    {
      // This element has fiber nodes.
      // Interpolate fibers to the Gauss points and add them to the parameter list

      // Get shape functions
      const static std::vector<Core::LinAlg::Matrix<Internal::num_nodes<celltype>, 1>> shapefcts =
          std::invoke(
              [&]
              {
                std::vector<Core::LinAlg::Matrix<Internal::num_nodes<celltype>, 1>> shapefcns(
                    stiffness_matrix_integration.num_points());
                for (int gp = 0; gp < stiffness_matrix_integration.num_points(); ++gp)
                {
                  Core::LinAlg::Matrix<Internal::num_dim<celltype>, 1> xi(
                      stiffness_matrix_integration.point(gp), true);
                  Core::FE::shape_function<celltype>(xi, shapefcns[gp]);
                }
                return shapefcns;
              });

      // add fibers to the ParameterList
      Core::Nodes::NodalFiberHolder fiberHolder;

      // Do the interpolation
      Core::Nodes::project_fibers_to_gauss_points<celltype>(ele.nodes(), shapefcts, fiberHolder);

      params.set("fiberholder", fiberHolder);
    }
  }

  /*!
   * @brief Linearization of a solid formulation containing the derivatives of the deformation
   * gradient w.r.t. the nodal displacements, xi and the second derivative w.r.t. xi and
   * displacements
   *
   * @tparam celltype
   */
  template <Core::FE::CellType celltype>
  struct SolidFormulationLinearization
  {
    /// Derivative of the deformation gradient w.r.t. nodal displacements
    Core::LinAlg::Matrix<9, Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>> d_F_dd{};

    /// Derivative of the deformation gradient w.r.t. xi
    Core::LinAlg::Matrix<9, Core::FE::dim<celltype>> d_F_dxi{};

    /// 2. Derivative of the deformation gradient w.r.t. xi and nodal displacements
    Core::LinAlg::Matrix<9,
        Core::FE::num_nodes<celltype> * Core::FE::dim<celltype> * Core::FE::dim<celltype>>
        d2_F_dxi_dd{};
  };

  /*!
   * @brief Evaluates and returns the transformation matrix T^{-T} which maps voigt tensors from
   * parameter space to the material configuration
   *
   * For details, see Andelfinger et al., EAS-elements, 1993, doi: 10.1002/nme.1620360805.
   *
   * @tparam celltype : Cell type
   * @param jacobian(in) : Jacobian mapping evaluated
   * @return : transformation matrix
   */
  template <Core::FE::CellType celltype>
  Core::LinAlg::Matrix<Internal::num_str<celltype>, Internal::num_str<celltype>>
  evaluate_voigt_transformation_matrix(const Discret::Elements::JacobianMapping<celltype>& jacobian)
  {
    // build T^T (based on strain-like Voigt notation: xx,yy,zz,xy,yz,xz)
    // currently only works in 3D
    Core::LinAlg::Matrix<Internal::num_str<celltype>, Internal::num_str<celltype>> TinvT(false);
    TinvT(0, 0) = jacobian.jacobian_(0, 0) * jacobian.jacobian_(0, 0);
    TinvT(1, 0) = jacobian.jacobian_(1, 0) * jacobian.jacobian_(1, 0);
    TinvT(2, 0) = jacobian.jacobian_(2, 0) * jacobian.jacobian_(2, 0);
    TinvT(3, 0) = 2 * jacobian.jacobian_(0, 0) * jacobian.jacobian_(1, 0);
    TinvT(4, 0) = 2 * jacobian.jacobian_(1, 0) * jacobian.jacobian_(2, 0);
    TinvT(5, 0) = 2 * jacobian.jacobian_(0, 0) * jacobian.jacobian_(2, 0);

    TinvT(0, 1) = jacobian.jacobian_(0, 1) * jacobian.jacobian_(0, 1);
    TinvT(1, 1) = jacobian.jacobian_(1, 1) * jacobian.jacobian_(1, 1);
    TinvT(2, 1) = jacobian.jacobian_(2, 1) * jacobian.jacobian_(2, 1);
    TinvT(3, 1) = 2 * jacobian.jacobian_(0, 1) * jacobian.jacobian_(1, 1);
    TinvT(4, 1) = 2 * jacobian.jacobian_(1, 1) * jacobian.jacobian_(2, 1);
    TinvT(5, 1) = 2 * jacobian.jacobian_(0, 1) * jacobian.jacobian_(2, 1);

    TinvT(0, 2) = jacobian.jacobian_(0, 2) * jacobian.jacobian_(0, 2);
    TinvT(1, 2) = jacobian.jacobian_(1, 2) * jacobian.jacobian_(1, 2);
    TinvT(2, 2) = jacobian.jacobian_(2, 2) * jacobian.jacobian_(2, 2);
    TinvT(3, 2) = 2 * jacobian.jacobian_(0, 2) * jacobian.jacobian_(1, 2);
    TinvT(4, 2) = 2 * jacobian.jacobian_(1, 2) * jacobian.jacobian_(2, 2);
    TinvT(5, 2) = 2 * jacobian.jacobian_(0, 2) * jacobian.jacobian_(2, 2);

    TinvT(0, 3) = jacobian.jacobian_(0, 0) * jacobian.jacobian_(0, 1);
    TinvT(1, 3) = jacobian.jacobian_(1, 0) * jacobian.jacobian_(1, 1);
    TinvT(2, 3) = jacobian.jacobian_(2, 0) * jacobian.jacobian_(2, 1);
    TinvT(3, 3) = jacobian.jacobian_(0, 0) * jacobian.jacobian_(1, 1) +
                  jacobian.jacobian_(1, 0) * jacobian.jacobian_(0, 1);
    TinvT(4, 3) = jacobian.jacobian_(1, 0) * jacobian.jacobian_(2, 1) +
                  jacobian.jacobian_(2, 0) * jacobian.jacobian_(1, 1);
    TinvT(5, 3) = jacobian.jacobian_(0, 0) * jacobian.jacobian_(2, 1) +
                  jacobian.jacobian_(2, 0) * jacobian.jacobian_(0, 1);

    TinvT(0, 4) = jacobian.jacobian_(0, 1) * jacobian.jacobian_(0, 2);
    TinvT(1, 4) = jacobian.jacobian_(1, 1) * jacobian.jacobian_(1, 2);
    TinvT(2, 4) = jacobian.jacobian_(2, 1) * jacobian.jacobian_(2, 2);
    TinvT(3, 4) = jacobian.jacobian_(0, 1) * jacobian.jacobian_(1, 2) +
                  jacobian.jacobian_(1, 1) * jacobian.jacobian_(0, 2);
    TinvT(4, 4) = jacobian.jacobian_(1, 1) * jacobian.jacobian_(2, 2) +
                  jacobian.jacobian_(2, 1) * jacobian.jacobian_(1, 2);
    TinvT(5, 4) = jacobian.jacobian_(0, 1) * jacobian.jacobian_(2, 2) +
                  jacobian.jacobian_(2, 1) * jacobian.jacobian_(0, 2);

    TinvT(0, 5) = jacobian.jacobian_(0, 0) * jacobian.jacobian_(0, 2);
    TinvT(1, 5) = jacobian.jacobian_(1, 0) * jacobian.jacobian_(1, 2);
    TinvT(2, 5) = jacobian.jacobian_(2, 0) * jacobian.jacobian_(2, 2);
    TinvT(3, 5) = jacobian.jacobian_(0, 0) * jacobian.jacobian_(1, 2) +
                  jacobian.jacobian_(1, 0) * jacobian.jacobian_(0, 2);
    TinvT(4, 5) = jacobian.jacobian_(1, 0) * jacobian.jacobian_(2, 2) +
                  jacobian.jacobian_(2, 0) * jacobian.jacobian_(1, 2);
    TinvT(5, 5) = jacobian.jacobian_(0, 0) * jacobian.jacobian_(2, 2) +
                  jacobian.jacobian_(2, 0) * jacobian.jacobian_(0, 2);

    // evaluate the inverse T0^{-T} with solver
    Core::LinAlg::FixedSizeSerialDenseSolver<6, 6, 1> solve_for_inverse;
    solve_for_inverse.set_matrix(TinvT);

    int err_inv = solve_for_inverse.invert();
    FOUR_C_ASSERT_ALWAYS(!err_inv, "Inversion of matrix failed with LAPACK error code {}", err_inv);

    return TinvT;
  }



  /*!
   * @brief Compute the deformation gradient from the Green-Lagrange strain tensor

   * @param defgrd_disp(in) : displacement-based deformation gradient F^{u} (needed for the
   rotational part)
   * @param enhanced_gl_strain(in) : Green-Lagrange strains E^{enh} to compute the deformation
   gradient from
   * @return Core::LinAlg::Matrix<dim, dim> : deformation gradient F^{enh} computed from GL strains
   */
  static inline Core::LinAlg::Matrix<3, 3> compute_deformation_gradient_from_gl_strains(
      const Core::LinAlg::Matrix<3, 3>& defgrd_disp,
      const Core::LinAlg::Matrix<6, 1>& enhanced_gl_strain)
  {
    Core::LinAlg::Matrix<3, 3> R;       // rotation tensor
    Core::LinAlg::Matrix<3, 3> U_enh;   // enhanced right stretch tensor
    Core::LinAlg::Matrix<3, 3> U_disp;  // displacement-based right stretch tensor
    Core::LinAlg::Matrix<3, 3> EW;      // temporarily store eigenvalues
    Core::LinAlg::Matrix<3, 3> tmp;     // temporary matrix for matrix matrix matrix products
    Core::LinAlg::Matrix<3, 3> tmp2;    // temporary matrix for matrix matrix matrix products

    // calculate modified right stretch tensor
    for (unsigned i = 0; i < 3; i++) U_enh(i, i) = 2. * enhanced_gl_strain(i) + 1.;
    U_enh(0, 1) = enhanced_gl_strain(3);
    U_enh(1, 0) = enhanced_gl_strain(3);
    U_enh(1, 2) = enhanced_gl_strain(4);
    U_enh(2, 1) = enhanced_gl_strain(4);
    U_enh(0, 2) = enhanced_gl_strain(5);
    U_enh(2, 0) = enhanced_gl_strain(5);

    Core::LinAlg::syev(U_enh, EW, U_enh);
    for (unsigned i = 0; i < 3; ++i) EW(i, i) = sqrt(EW(i, i));
    tmp.multiply(U_enh, EW);
    tmp2.multiply_nt(tmp, U_enh);
    U_enh.update(tmp2);

    // calculate displacement-based right stretch tensor
    U_disp.multiply_tn(defgrd_disp, defgrd_disp);

    Core::LinAlg::syev(U_disp, EW, U_disp);
    for (unsigned i = 0; i < 3; ++i) EW(i, i) = sqrt(EW(i, i));
    tmp.multiply(U_disp, EW);
    tmp2.multiply_nt(tmp, U_disp);
    U_disp.update(tmp2);

    // compose consistent deformation gradient
    U_disp.invert();
    R.multiply(defgrd_disp, U_disp);

    Core::LinAlg::Matrix<3, 3> defgrd_enh;
    defgrd_enh.multiply(R, U_enh);
    return defgrd_enh;
  }

}  // namespace Discret::Elements

FOUR_C_NAMESPACE_CLOSE

#endif