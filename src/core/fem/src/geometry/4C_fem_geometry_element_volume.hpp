// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FEM_GEOMETRY_ELEMENT_VOLUME_HPP
#define FOUR_C_FEM_GEOMETRY_ELEMENT_VOLUME_HPP


#include "4C_config.hpp"

#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_fem_general_utils_gausspoints.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::Geo
{
  //! calculates the length of an element in given configuration
  template <Core::FE::CellType distype, class Matrixtype>
  double element_length_t(const Matrixtype& xyze)  ///> xyze nsd = 3 coords, number of nodes
  {
    // gaussian points
    static constexpr int numnode = Core::FE::num_nodes<distype>;
    const Core::FE::GaussIntegration intpoints(distype);

    Core::LinAlg::Matrix<1, 1> eleCoord;
    Core::LinAlg::Matrix<1, numnode> deriv;
    Core::LinAlg::Matrix<1, 3> xjm;

    double length = 0.0;

    // integration loop
    for (int iquad = 0; iquad < intpoints.num_points(); ++iquad)
    {
      // coordinates of the current integration point in element coordinates \xi
      eleCoord(0) = intpoints.point(iquad)[0];

      // shape functions and their first derivatives
      Core::FE::shape_function_1d_deriv1(deriv, eleCoord(0), distype);

      // get transposed of the jacobian matrix d x / d \xi
      xjm = 0;

      for (int inode = 0; inode < numnode; ++inode)
        for (int i = 0; i < 1; ++i)
          for (int j = 0; j < 3; ++j) xjm(i, j) += deriv(i, inode) * xyze(j, inode);

      const double fac = intpoints.weight(iquad) * xjm.norm2();

      length += fac;
    }  // end loop over gauss points

    return length;
  }

  //! calculates the area of an element in given configuration
  template <Core::FE::CellType distype, class Matrixtype>
  double element_area_t(const Matrixtype& xyze)  ///> xyze nsd = 3 coords, number of nodes
  {
    // gaussian points
    static constexpr int numnode = Core::FE::num_nodes<distype>;
    const Core::FE::GaussIntegration intpoints(distype);

    Core::LinAlg::Matrix<2, 1> eleCoord;
    Core::LinAlg::Matrix<2, numnode> deriv;
    Core::LinAlg::Matrix<2, 3> xjm;
    Core::LinAlg::Matrix<2, 2> xjm_xjmt;

    double area = 0.0;

    // integration loop
    for (int iquad = 0; iquad < intpoints.num_points(); ++iquad)
    {
      // coordinates of the current integration point in element coordinates \xi
      eleCoord(0) = intpoints.point(iquad)[0];
      eleCoord(1) = intpoints.point(iquad)[1];

      // shape functions and their first derivatives
      Core::FE::shape_function_2d_deriv1(deriv, eleCoord(0), eleCoord(1), distype);

      // get transposed of the jacobian matrix d x / d \xi
      xjm = 0;

      for (int inode = 0; inode < numnode; ++inode)
        for (int i = 0; i < 2; ++i)
          for (int j = 0; j < 3; ++j) xjm(i, j) += deriv(i, inode) * xyze(j, inode);

      xjm_xjmt.multiply_nt<3>(xjm, xjm);

      const double det = xjm_xjmt.determinant();
      const double fac = intpoints.weight(iquad) * std::sqrt(det);

      area += fac;
    }  // end loop over gauss points

    return area;
  }

  //! calculates the volume of an element in given configuration          u.may
  template <Core::FE::CellType distype, class Matrixtype>
  double element_volume_t(const Matrixtype& xyze  ///> xyze nsd = 3 coords, number of nodes)
  )
  {
    // number of nodes for element
    const int numnode = Core::FE::num_nodes<distype>;
    // gaussian points
    const Core::FE::GaussIntegration intpoints(distype);

    Core::LinAlg::Matrix<3, 1> eleCoord;
    Core::LinAlg::Matrix<3, numnode> deriv;
    Core::LinAlg::Matrix<3, 3> xjm;

    double vol = 0.0;

    // integration loop
    for (int iquad = 0; iquad < intpoints.num_points(); ++iquad)
    {
      // coordinates of the current integration point in element coordinates \xi
      eleCoord(0) = intpoints.point(iquad)[0];
      eleCoord(1) = intpoints.point(iquad)[1];
      eleCoord(2) = intpoints.point(iquad)[2];

      // shape functions and their first derivatives
      Core::FE::shape_function_3d_deriv1(deriv, eleCoord(0), eleCoord(1), eleCoord(2), distype);

      // get transposed of the jacobian matrix d x / d \xi
      xjm = 0;

      for (int inode = 0; inode < numnode; ++inode)
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j) xjm(i, j) += deriv(i, inode) * xyze(j, inode);

      const double det = xjm.determinant();
      const double fac = intpoints.weight(iquad) * det;

      if (det <= 0.0) FOUR_C_THROW("NEGATIVE JACOBIAN DETERMINANT: {}", det);

      vol += fac;
    }  // end loop over gauss points
    return vol;
  }

  /** \brief calculates the length of a edge element in given configuration
   *
   *  \param distype (in) : discretization type of the given element
   *  \param xyze    (in) : spatial coordinates of the elememnt nodes
   *                         (row = dim, col = number of nodes)
   */
  template <class Matrixtype>
  double element_length(const Core::FE::CellType& distype, const Matrixtype& xyze)
  {
    if (distype != Core::FE::CellType::line2 or xyze.numCols() != 2)
      FOUR_C_THROW("Currently only line2 elements are supported!");

    // calculate the distance between the two given nodes and return
    // the value
    Core::LinAlg::Matrix<3, 1> d(&xyze(0, 0));
    const Core::LinAlg::Matrix<3, 1> x1(&xyze(0, 1));

    d.update(1.0, x1, -1.0);
    return d.norm2();
  }


  /** \brief calculates the area of a surface element in given configuration
   */
  template <class Matrixtype>
  double element_area(const Core::FE::CellType distype, const Matrixtype& xyze)
  {
    switch (distype)
    {
      // --- 2-D boundary elements
      case Core::FE::CellType::line2:
        return element_length(distype, xyze);
      case Core::FE::CellType::line3:
        return element_length_t<Core::FE::CellType::line3>(xyze);
      // --- 3-D boundary elements
      case Core::FE::CellType::tri3:
        return element_area_t<Core::FE::CellType::tri3>(xyze);
      case Core::FE::CellType::tri6:
        return element_area_t<Core::FE::CellType::tri6>(xyze);
      case Core::FE::CellType::quad4:
        return element_area_t<Core::FE::CellType::quad4>(xyze);
      case Core::FE::CellType::quad8:
        return element_area_t<Core::FE::CellType::quad8>(xyze);
      case Core::FE::CellType::quad9:
        return element_area_t<Core::FE::CellType::quad9>(xyze);
      default:
        FOUR_C_THROW("Unsupported surface element type!");
        exit(EXIT_FAILURE);
    }

    return -1.0;
  }

  //! calculates the volume of an element in given configuration          u.may
  template <class Matrixtype>
  double element_volume(const Core::FE::CellType distype,
      const Matrixtype& xyze  ///> xyze nsd = 3 coords, number of nodes
  )
  {
    switch (distype)
    {
      // --- 1-D elements -----------------------------------------------------
      case Core::FE::CellType::line2:
        return element_length(distype, xyze);
      // --- 2-D elements -----------------------------------------------------
      case Core::FE::CellType::quad4:
      case Core::FE::CellType::tri3:
        return element_area(distype, xyze);
      // --- 3-D elements -----------------------------------------------------
      case Core::FE::CellType::hex8:
        return element_volume_t<Core::FE::CellType::hex8>(xyze);
      case Core::FE::CellType::hex20:
        return element_volume_t<Core::FE::CellType::hex20>(xyze);
      case Core::FE::CellType::hex27:
        return element_volume_t<Core::FE::CellType::hex27>(xyze);
      case Core::FE::CellType::tet4:
        return element_volume_t<Core::FE::CellType::tet4>(xyze);
      case Core::FE::CellType::tet10:
        return element_volume_t<Core::FE::CellType::tet10>(xyze);
      case Core::FE::CellType::wedge6:
        return element_volume_t<Core::FE::CellType::wedge6>(xyze);
      case Core::FE::CellType::wedge15:
        return element_volume_t<Core::FE::CellType::wedge15>(xyze);
      case Core::FE::CellType::pyramid5:
        return element_volume_t<Core::FE::CellType::pyramid5>(xyze);
      default:
        FOUR_C_THROW("add you distype here");
    }
    return -1.0;
  }

}  // namespace Core::Geo


FOUR_C_NAMESPACE_CLOSE

#endif
