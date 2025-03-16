// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_BEAM3_SPATIAL_DISCRETIZATION_UTILS_HPP
#define FOUR_C_BEAM3_SPATIAL_DISCRETIZATION_UTILS_HPP

#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_fem_general_utils_gausspoints.hpp"
#include "4C_fem_general_utils_integration.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret::Utils::Beam
{
  /** \brief evaluate shape functions at position \xi in element parameter space [-1,1]
   */
  template <unsigned int nnode, unsigned int vpernode, typename T>
  void evaluate_shape_functions_at_xi(const T& xi,
      Core::LinAlg::Matrix<1, vpernode * nnode, T>& I_i, const Core::FE::CellType& distype,
      double hermite_length_param = -1.0)
  {
    I_i.clear();

    switch (vpernode)
    {
      case 1:
      {
        // evaluate Lagrange shape functions at xi
        Core::FE::shape_function_1d(I_i, xi, distype);
        break;
      }
      case 2:
      {
        FOUR_C_ASSERT(hermite_length_param != -1.0,
            "you must provide a length parameter in case of "
            "Hermite interpolation!");

        // evaluate Hermite shape functions at xi: vpernode=2 means 3rd order, i.e. line2
        Core::FE::shape_function_hermite_1d(
            I_i, xi, hermite_length_param, Core::FE::CellType::line2);
        break;
      }
      default:
        FOUR_C_THROW("invalid value for vpernode (number of values per node) specified");
    }
  }

  /** \brief evaluate shape function derivatives at position \xi in element parameter space
   * [-1,1]
   */
  template <unsigned int nnode, unsigned int vpernode, typename T>
  void evaluate_shape_function_derivs_at_xi(const T& xi,
      Core::LinAlg::Matrix<1, vpernode * nnode, T>& I_i_xi, const Core::FE::CellType& distype,
      double hermite_length_param = -1.0)
  {
    I_i_xi.clear();

    switch (vpernode)
    {
      case 1:
      {
        // evaluate Lagrange shape function derivs at xi
        Core::FE::shape_function_1d_deriv1(I_i_xi, xi, distype);
        break;
      }
      case 2:
      {
        FOUR_C_ASSERT(hermite_length_param != -1.0,
            "you must provide a length parameter in case of "
            "Hermite interpolation!");

        // evaluate Hermite shape function derivs at xi: vpernode=2 means 3rd order, i.e. line2
        Core::FE::shape_function_hermite_1d_deriv1(
            I_i_xi, xi, hermite_length_param, Core::FE::CellType::line2);
        break;
      }
      default:
        FOUR_C_THROW("invalid value for vpernode (number of values per node) specified");
    }
  }

  /** \brief evaluate second derivatives of shape function at position \xi in element parameter
   *         space [-1,1]
   */
  template <unsigned int nnode, unsigned int vpernode, typename T>
  void evaluate_shape_function2nd_derivs_at_xi(const T& xi,
      Core::LinAlg::Matrix<1, vpernode * nnode, T>& I_i_xixi, const Core::FE::CellType& distype,
      double hermite_length_param = -1.0)
  {
    I_i_xixi.clear();

    switch (vpernode)
    {
      case 1:
      {
        // evaluate Lagrange shape function derivs at xi
        Core::FE::shape_function_1d_deriv2(I_i_xixi, xi, distype);
        break;
      }
      case 2:
      {
        FOUR_C_ASSERT(hermite_length_param != -1.0,
            "you must provide a length parameter in case of "
            "Hermite interpolation!");

        // evaluate Hermite shape function derivs at xi: vpernode=2 means 3rd order, i.e. line2
        Core::FE::shape_function_hermite_1d_deriv2(
            I_i_xixi, xi, hermite_length_param, Core::FE::CellType::line2);
        break;
      }
      default:
        FOUR_C_THROW("invalid value for vpernode (number of values per node) specified");
    }
  }

  /** \brief evaluate third derivatives of shape function at position \xi in element parameter
   *         space [-1,1]
   */
  template <unsigned int nnode, unsigned int vpernode, typename T>
  void evaluate_shape_function3rd_derivs_at_xi(const T& xi,
      Core::LinAlg::Matrix<1, vpernode * nnode, T>& I_i_xixixi, const Core::FE::CellType& distype,
      double hermite_length_param = -1.0)
  {
    I_i_xixixi.clear();

    switch (vpernode)
    {
      case 1:
      {
        // evaluate Lagrange shape function derivs at xi
        if (nnode >= 4)
          FOUR_C_THROW("Please implement 3rd derivatives of Lagrange shape functions!");

        break;
      }
      case 2:
      {
        FOUR_C_ASSERT(hermite_length_param != -1.0,
            "you must provide a length parameter in case of "
            "Hermite interpolation!");

        // evaluate Hermite shape function derivs at xi: vpernode=2 means 3rd order, i.e. line2
        Core::FE::shape_function_hermite_1d_deriv3(
            I_i_xixixi, xi, hermite_length_param, Core::FE::CellType::line2);
        break;
      }
      default:
        FOUR_C_THROW("invalid value for vpernode (number of values per node) specified");
    }
  }

  /** \brief evaluate shape functions and its first derivatives at position \xi in element
   *         parameter space [-1,1]
   */
  template <unsigned int nnode, unsigned int vpernode, typename T>
  void evaluate_shape_functions_and_derivs_at_xi(const T& xi,
      Core::LinAlg::Matrix<1, vpernode * nnode, T>& I_i,
      Core::LinAlg::Matrix<1, vpernode * nnode, T>& I_i_xi, const Core::FE::CellType& distype,
      double hermite_length_param = -1.0)
  {
    evaluate_shape_functions_at_xi<nnode, vpernode>(xi, I_i, distype, hermite_length_param);
    evaluate_shape_function_derivs_at_xi<nnode, vpernode>(
        xi, I_i_xi, distype, hermite_length_param);
  }

  /** \brief evaluate shape functions and its first and second derivatives at position \xi in
   *         element parameter space [-1,1]
   *
   */
  template <unsigned int nnode, unsigned int vpernode, typename T>
  void evaluate_shape_functions_and_derivs_and2nd_derivs_at_xi(const T& xi,
      Core::LinAlg::Matrix<1, vpernode * nnode, T>& I_i,
      Core::LinAlg::Matrix<1, vpernode * nnode, T>& I_i_xi,
      Core::LinAlg::Matrix<1, vpernode * nnode, T>& I_i_xixi, const Core::FE::CellType& distype,
      double hermite_length_param = -1.0)
  {
    evaluate_shape_functions_at_xi<nnode, vpernode>(xi, I_i, distype, hermite_length_param);
    evaluate_shape_function_derivs_at_xi<nnode, vpernode>(
        xi, I_i_xi, distype, hermite_length_param);
    evaluate_shape_function2nd_derivs_at_xi<nnode, vpernode>(
        xi, I_i_xixi, distype, hermite_length_param);
  }

  /** \brief evaluate shape functions at all specified Gauss points at once
   *
   */
  template <unsigned int nnode, unsigned int vpernode, typename T>
  void evaluate_shape_functions_all_gps(const Core::FE::IntegrationPoints1D& gausspoints,
      std::vector<Core::LinAlg::Matrix<1, vpernode * nnode, T>>& I_i,
      const Core::FE::CellType& distype, double hermite_length_param = -1.0,
      double integration_interval_lower_limit = -1.0, double integration_interval_upper_limit = 1.0)
  {
    if (I_i.size() != (unsigned int)gausspoints.nquad)
      FOUR_C_THROW(
          "vector for individual shape functions to be evaluated at {} GPs has wrong size: {}",
          gausspoints.nquad, I_i.size());

    for (int numgp = 0; numgp < gausspoints.nquad; numgp++)
    {
      // Get location of GP in element parameter space xi \in [-1;1]
      const T xi_tilde = gausspoints.qxg[numgp][0];

      /* do a mapping into integration interval, i.e. coordinate transformation
       * note: this has no effect if integration interval is [-1;1] */
      const T xi = 0.5 * ((1.0 - xi_tilde) * integration_interval_lower_limit +
                             (1.0 + xi_tilde) * integration_interval_upper_limit);

      evaluate_shape_functions_at_xi<nnode, vpernode>(
          xi, I_i[numgp], distype, hermite_length_param);
    }
  }

  /** \brief evaluate shape function derivatives at all specified Gauss points at once
   *
   */
  template <unsigned int nnode, unsigned int vpernode, typename T>
  void evaluate_shape_function_derivs_all_gps(const Core::FE::IntegrationPoints1D& gausspoints,
      std::vector<Core::LinAlg::Matrix<1, vpernode * nnode, T>>& I_i_xi,
      const Core::FE::CellType& distype, double hermite_length_param = -1.0,
      double integration_interval_lower_limit = -1.0, double integration_interval_upper_limit = 1.0)
  {
    if (I_i_xi.size() != (unsigned int)gausspoints.nquad)
      FOUR_C_THROW(
          "vector for individual shape function derivatives to be evaluated at {} GPs has "
          "wrong size: {}",
          gausspoints.nquad, I_i_xi.size());

    for (int numgp = 0; numgp < gausspoints.nquad; ++numgp)
    {
      // Get location of GP in element parameter space xi \in [-1;1]
      const T xi_tilde = gausspoints.qxg[numgp][0];

      /* do a mapping into integration interval, i.e. coordinate transformation
       * note: this has no effect if integration interval is [-1;1] */
      const T xi = 0.5 * ((1.0 - xi_tilde) * integration_interval_lower_limit +
                             (1.0 + xi_tilde) * integration_interval_upper_limit);

      evaluate_shape_function_derivs_at_xi<nnode, vpernode>(
          xi, I_i_xi[numgp], distype, hermite_length_param);
    }
  }

  /** \brief evaluate shape functions and derivatives at all specified Gauss points at once
   *
   */
  template <unsigned int nnode, unsigned int vpernode, typename T>
  void evaluate_shape_functions_and_derivs_all_gps(const Core::FE::IntegrationPoints1D& gausspoints,
      std::vector<Core::LinAlg::Matrix<1, vpernode * nnode, T>>& I_i,
      std::vector<Core::LinAlg::Matrix<1, vpernode * nnode, T>>& I_i_xi,
      const Core::FE::CellType& distype, double hermite_length_param = -1.0,
      double integration_interval_lower_limit = -1.0, double integration_interval_upper_limit = 1.0)
  {
    evaluate_shape_functions_all_gps<nnode, vpernode>(gausspoints, I_i, distype,
        hermite_length_param, integration_interval_lower_limit, integration_interval_upper_limit);
    evaluate_shape_function_derivs_all_gps<nnode, vpernode>(gausspoints, I_i_xi, distype,
        hermite_length_param, integration_interval_lower_limit, integration_interval_upper_limit);
  }

  /** \brief assemble one shape function matrix, such that: r=N*d
   */
  template <unsigned int numnodes, unsigned int numnodalvalues, typename T>
  void assemble_shape_functions(const Core::LinAlg::Matrix<1, numnodes * numnodalvalues, T>& N_i,
      Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, T>& N)
  {
    // assembly_N is just an array to help assemble the matrices of the shape functions
    // it determines, which shape function is used in which column of N
    std::array<std::array<unsigned int, 3 * numnodes * numnodalvalues>, 3> assembly_N{};

    /*
    Set number of shape functions for each 3*3 block:
    e.g. second order Lagrange polynomials (numnodes=3, numnodalvalues=1)
    int assembly_N[3][9]=  { {1,0,0,2,0,0,3,0,0},
                             {0,1,0,0,2,0,0,3,0},
                             {0,0,1,0,0,2,0,0,3}};

    e.g. cubic Hermite polynomials (numnodes=2, numnodalvalues=2)
    int assembly_N[3][12]=  {{1,0,0,2,0,0,3,0,0,4,0,0},
                             {0,1,0,0,2,0,0,3,0,0,4,0},
                             {0,0,1,0,0,2,0,0,3,0,0,4}};
    */

    for (unsigned int i = 0; i < numnodes * numnodalvalues; ++i)
    {
      assembly_N[0][3 * i] = i + 1;
      assembly_N[1][3 * i + 1] = i + 1;
      assembly_N[2][3 * i + 2] = i + 1;
    }

    // Assemble the matrices of the shape functions
    for (unsigned int i = 0; i < 3 * numnodes * numnodalvalues; ++i)
      for (unsigned int j = 0; j < 3; ++j)
      {
        if (assembly_N[j][i] == 0)
          N(j, i) = 0.0;
        else
          N(j, i) = N_i(assembly_N[j][i] - 1);
      }
  }

  /** \brief assemble shape function matrices, such that: r=N*d, r_xi=N_xi*d, r_xixi=N_xixi*d
   */
  template <unsigned int numnodes, unsigned int numnodalvalues, typename T>
  void assemble_shape_functions_and_derivs_and2nd_derivs(
      const Core::LinAlg::Matrix<1, numnodes * numnodalvalues, T>& N_i,
      const Core::LinAlg::Matrix<1, numnodes * numnodalvalues, T>& N_i_xi,
      const Core::LinAlg::Matrix<1, numnodes * numnodalvalues, T>& N_i_xixi,
      Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, T>& N,
      Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, T>& N_xi,
      Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, T>& N_xixi)
  {
    assemble_shape_functions<numnodes, numnodalvalues, T>(N_i, N);
    assemble_shape_functions<numnodes, numnodalvalues, T>(N_i_xi, N_xi);
    assemble_shape_functions<numnodes, numnodalvalues, T>(N_i_xixi, N_xixi);
  }

  /** \brief interpolation of nodal DoFs based on given shape function values
   */
  template <unsigned int nnode, unsigned int vpernode, unsigned int ndim, typename T, typename T2>
  void calc_interpolation(const Core::LinAlg::Matrix<ndim * vpernode * nnode, 1, T>& dof_vals,
      const Core::LinAlg::Matrix<1, vpernode * nnode, T2>& shapefcn_vals,
      Core::LinAlg::Matrix<ndim, 1, T>& result)
  {
    result.clear();

    for (unsigned int idim = 0; idim < ndim; ++idim)
      for (unsigned int idofperdim = 0; idofperdim < vpernode * nnode; ++idofperdim)
        result(idim) += shapefcn_vals(idofperdim) * dof_vals(ndim * idofperdim + idim);
  }

  /** \brief evaluate length and derivative at all specified Gauss points at once
   */
  template <unsigned int nnode, unsigned int vpernode>
  std::tuple<double, double> integrate_centerline_arc_length_and_arc_length_derivative(
      const Core::FE::IntegrationPoints1D& gausspoints,
      const Core::LinAlg::Matrix<3 * vpernode * nnode, 1, double>& disp_centerline,
      const Core::FE::CellType& distype, const double& reflength)
  {
    std::vector<Core::LinAlg::Matrix<1, nnode * vpernode, double>> H_i_xi(gausspoints.nquad);
    Core::LinAlg::Matrix<3, 1> r_xi;

    evaluate_shape_function_derivs_all_gps<nnode, vpernode>(
        gausspoints, H_i_xi, distype, reflength);

    double int_length = 0.0, deriv_length = 0.0;

    for (int numgp = 0; numgp < gausspoints.nquad; numgp++)
    {
      double deriv_int = 0.0;

      calc_interpolation<nnode, vpernode, 3, double>(disp_centerline, H_i_xi[numgp], r_xi);
      int_length += gausspoints.qwgt[numgp] * r_xi.norm2();

      for (int dim = 0; dim < 3; dim++)
        deriv_int +=
            (disp_centerline(3 + dim) * H_i_xi[numgp](1) / reflength +
                disp_centerline(3 * vpernode * 1 + 3 + dim) * H_i_xi[numgp](3) / reflength) *
            r_xi(dim);

      deriv_length += gausspoints.qwgt[numgp] * deriv_int / r_xi.norm2();
    }

    return {reflength - int_length, 1 - deriv_length};
  }
}  // namespace Discret::Utils::Beam

FOUR_C_NAMESPACE_CLOSE

#endif
