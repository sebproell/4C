// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_fem_general_element_integration.hpp"

#include "4C_fem_general_cell_type.hpp"
#include "4C_fem_general_element_integration_select.hpp"
#include "4C_fem_general_utils_gausspoints.hpp"
#include "4C_utils_exceptions.hpp"


namespace
{
  using namespace FourC;

  template <Core::FE::CellType celltype, unsigned dim>
  std::tuple<Core::Elements::ElementNodes<celltype, dim>, double>
  generate_element_nodes_and_reference_solution()
  {
    if constexpr (celltype == Core::FE::CellType::quad4)
    {
      if constexpr (dim == 2)
      {
        Core::Elements::ElementNodes<celltype, dim> element_nodes;
        element_nodes.coordinates(0, 0) = 1.2;
        element_nodes.coordinates(0, 1) = 1.5;
        element_nodes.coordinates(1, 0) = 4.3;
        element_nodes.coordinates(1, 1) = 1.8;
        element_nodes.coordinates(2, 0) = 4.0;
        element_nodes.coordinates(2, 1) = 3.6;
        element_nodes.coordinates(3, 0) = 1.5;
        element_nodes.coordinates(3, 1) = 3.2;

        return {element_nodes, 4.9};
      }

      if constexpr (dim == 3)
      {
        Core::Elements::ElementNodes<celltype, dim> element_nodes;
        element_nodes.coordinates(0, 0) = 1.2;
        element_nodes.coordinates(0, 1) = 1.5;
        element_nodes.coordinates(0, 2) = 1.5;
        element_nodes.coordinates(1, 0) = 4.3;
        element_nodes.coordinates(1, 1) = 1.8;
        element_nodes.coordinates(1, 2) = 1.5;
        element_nodes.coordinates(2, 0) = 4.0;
        element_nodes.coordinates(2, 1) = 3.6;
        element_nodes.coordinates(2, 2) = 1.5;
        element_nodes.coordinates(3, 0) = 1.5;
        element_nodes.coordinates(3, 1) = 3.2;
        element_nodes.coordinates(3, 2) = 1.5;

        return {element_nodes, 4.9};
      }
    }


    if constexpr (celltype == Core::FE::CellType::tri3)
    {
      if constexpr (dim == 2)
      {
        Core::Elements::ElementNodes<celltype, dim> element_nodes;
        element_nodes.coordinates(0, 0) = 1.5;
        element_nodes.coordinates(0, 1) = 2.3;
        element_nodes.coordinates(1, 0) = 3.1;
        element_nodes.coordinates(1, 1) = 1.8;
        element_nodes.coordinates(2, 0) = 4.7;
        element_nodes.coordinates(2, 1) = 5.2;

        return {element_nodes, 3.12};
      }

      if constexpr (dim == 3)
      {
        Core::Elements::ElementNodes<celltype, dim> element_nodes;
        element_nodes.coordinates(0, 0) = 1.5;
        element_nodes.coordinates(0, 1) = 2.3;
        element_nodes.coordinates(0, 2) = 2.1;
        element_nodes.coordinates(1, 0) = 3.1;
        element_nodes.coordinates(1, 1) = 1.8;
        element_nodes.coordinates(1, 2) = 2.2;
        element_nodes.coordinates(2, 0) = 4.7;
        element_nodes.coordinates(2, 1) = 5.2;
        element_nodes.coordinates(2, 2) = 2.4;

        return {element_nodes, 3.128769726266218};
      }
    }

    if constexpr (celltype == Core::FE::CellType::line2)
    {
      if constexpr (dim == 1)
      {
        Core::Elements::ElementNodes<celltype, dim> element_nodes;
        element_nodes.coordinates(0, 0) = 1.5;
        element_nodes.coordinates(1, 0) = 5.6;

        return {element_nodes, 4.1};
      }

      if constexpr (dim == 2)
      {
        Core::Elements::ElementNodes<celltype, dim> element_nodes;
        element_nodes.coordinates(0, 0) = 1.5;
        element_nodes.coordinates(0, 1) = 2.3;
        element_nodes.coordinates(1, 0) = 3.1;
        element_nodes.coordinates(1, 1) = 1.8;

        return {element_nodes, std::sqrt(2.81)};
      }

      if constexpr (dim == 3)
      {
        Core::Elements::ElementNodes<celltype, dim> element_nodes;
        element_nodes.coordinates(0, 0) = 1.5;
        element_nodes.coordinates(0, 1) = 2.3;
        element_nodes.coordinates(0, 2) = 2.1;
        element_nodes.coordinates(1, 0) = 3.1;
        element_nodes.coordinates(1, 1) = 1.8;
        element_nodes.coordinates(1, 2) = 2.2;

        return {element_nodes, std::sqrt(2.82)};
      }
    }

    FOUR_C_THROW("Unsupported celltype or dim");
  }


  template <Core::FE::CellType celltype, unsigned dim>
  void test_element_integration()
  {
    const auto [element_nodes, ref_solution] =
        generate_element_nodes_and_reference_solution<celltype, dim>();

    Core::FE::GaussIntegration gauss_integration = Core::FE::create_gauss_integration<celltype>(
        Discret::Elements::DisTypeToOptGaussRule<celltype>::rule);

    double unit_integral = 0.0;
    Core::Elements::for_each_gauss_point(element_nodes, gauss_integration,
        [&](auto xi, auto shape_functions, auto jacobian_mapping, auto integration_factor, auto gp)
        { unit_integral += integration_factor; });

    EXPECT_NEAR(unit_integral, ref_solution, 1.0e-10);
  }

  TEST(ElementIntegrationTest, TestIntegrationQuad4Dim2)
  {
    test_element_integration<Core::FE::CellType::quad4, 2>();
  }

  TEST(ElementIntegrationTest, TestIntegrationQuad4Dim3)
  {
    test_element_integration<Core::FE::CellType::quad4, 3>();
  }

  TEST(ElementIntegrationTest, TestIntegrationTri3Dim2)
  {
    test_element_integration<Core::FE::CellType::tri3, 2>();
  }

  TEST(ElementIntegrationTest, TestIntegrationTri3Dim3)
  {
    test_element_integration<Core::FE::CellType::tri3, 3>();
  }

  TEST(ElementIntegrationTest, TestIntegrationLine2Dim1)
  {
    test_element_integration<Core::FE::CellType::line2, 1>();
  }

  TEST(ElementIntegrationTest, TestIntegrationLine2Dim2)
  {
    test_element_integration<Core::FE::CellType::line2, 2>();
  }

  TEST(ElementIntegrationTest, TestIntegrationLine2Dim3)
  {
    test_element_integration<Core::FE::CellType::line2, 3>();
  }

}  // namespace
