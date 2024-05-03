/*---------------------------------------------------------------------------*/
/*! \file

\brief Unittests to check utility functions for the reference configuration

\level 3
*/
/*----------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include "4C_contact_element.hpp"
#include "4C_lib_utils_reference_configuration.hpp"
#include "4C_so3_hex8.hpp"
#include "4C_so3_tet4.hpp"
#include "4C_unittest_utils_assertions_test.hpp"

#include <Epetra_SerialComm.h>

namespace
{
  using namespace FourC;

  class UtilsRefConfigTest : public testing::Test
  {
   public:
    Teuchos::RCP<DRT::Discretization> testdis_;

    UtilsRefConfigTest()
    {
      // create a discretization, to store the created elements and nodes
      testdis_ =
          Teuchos::rcp(new DRT::Discretization("dummy", Teuchos::rcp(new Epetra_SerialComm)));

      // create hex8 element and store it in the test discretization
      const std::array<int, 8> nodeidshex8 = {0, 1, 2, 3, 4, 5, 6, 7};
      const std::vector<std::vector<double>> coordshex8 = {{-0.10, -0.20, -0.50},
          {1.25, 0.23, 0.66}, {1.20, 0.99, 0.50}, {-0.11, 1.20, 0.66}, {-0.10, -0.20, 1.90},
          {1.00, 0.00, 1.90}, {1.20, 0.99, 1.50}, {-0.11, -0.20, 1.66}};
      for (int i = 0; i < 8; ++i)
      {
        testdis_->AddNode(Teuchos::rcp(new DRT::Node(nodeidshex8[i], coordshex8[i], 0)));
      }
      Teuchos::RCP<DRT::ELEMENTS::SoHex8> testhex8ele =
          Teuchos::rcp(new DRT::ELEMENTS::SoHex8(0, 0));
      testhex8ele->SetNodeIds(8, nodeidshex8.data());
      testdis_->AddElement(testhex8ele);

      // create corresponding quad4 surface contact element and store it
      Teuchos::RCP<CONTACT::Element> testcontactquad4ele =
          Teuchos::rcp(new CONTACT::Element(testhex8ele->Id() + 1, testhex8ele->Owner(),
              testhex8ele->Shape(), testhex8ele->NumNode(), testhex8ele->NodeIds(), false, false));
      testdis_->AddElement(testcontactquad4ele);

      // create tet4 element and store it in the test discretization
      const std::array<int, 4> nodeidstet4 = {8, 9, 10, 11};
      const std::vector<std::vector<double>> coordstet4 = {
          {2.5, -0.5, 0.0}, {1.0, -1.1, 0.1}, {1.1, 0.11, 0.15}, {1.5, -0.5, 2.0}};
      for (int j = 0; j < 4; ++j)
      {
        testdis_->AddNode(Teuchos::rcp(new DRT::Node(nodeidstet4[j], coordstet4[j], 0)));
      }
      Teuchos::RCP<DRT::ELEMENTS::SoTet4> testtet4ele =
          Teuchos::rcp(new DRT::ELEMENTS::SoTet4(2, 0));
      testtet4ele->SetNodeIds(4, nodeidstet4.data());
      testdis_->AddElement(testtet4ele);

      // create corresponding tri3 surface contact element and store it
      Teuchos::RCP<CONTACT::Element> testcontacttri3ele =
          Teuchos::rcp(new CONTACT::Element(testtet4ele->Id() + 1, testtet4ele->Owner(),
              testtet4ele->Shape(), testtet4ele->NumNode(), testtet4ele->NodeIds(), false, false));
      testdis_->AddElement(testcontacttri3ele);
      testdis_->FillComplete(false, false, false);
    }
  };

  TEST_F(UtilsRefConfigTest, LocalToGlobalPositionAtXiRefConfig)
  {
    // get hex8 element and test it
    const DRT::Element* hex8ele = testdis_->gElement(0);
    CORE::LINALG::Matrix<3, 1> xicenterhex8ele(true);
    CORE::LINALG::Matrix<3, 1> hex8elecoords(true);
    CORE::LINALG::Matrix<3, 1> hex8refsolution(true);
    hex8refsolution(0, 0) = 423.0 / 800.0;
    hex8refsolution(1, 0) = 281.0 / 800.0;
    hex8refsolution(2, 0) = 207.0 / 200.0;
    DRT::UTILS::LocalToGlobalPositionAtXiRefConfig<3, CORE::FE::CellType::hex8>(
        hex8ele, xicenterhex8ele, hex8elecoords);

    FOUR_C_EXPECT_NEAR(hex8elecoords, hex8refsolution, 1e-14);

    // get quad4 element and test it
    const DRT::Element* quad4ele = testdis_->gElement(1);
    CORE::LINALG::Matrix<2, 1> xicenterquad4ele(true);
    CORE::LINALG::Matrix<3, 1> quad4elecoords(true);
    CORE::LINALG::Matrix<3, 1> quad4refsolution(true);
    quad4refsolution(0, 0) = 14.0 / 25.0;
    quad4refsolution(1, 0) = 111.0 / 200.0;
    quad4refsolution(2, 0) = 33.0 / 100.0;
    DRT::UTILS::LocalToGlobalPositionAtXiRefConfig<3, CORE::FE::CellType::quad4>(
        quad4ele, xicenterquad4ele, quad4elecoords);

    FOUR_C_EXPECT_NEAR(quad4elecoords, quad4refsolution, 1e-14);

    // get tet4 element stuff and test it
    const DRT::Element* tet4ele = testdis_->gElement(2);
    CORE::LINALG::Matrix<3, 1> xicentertet4ele(true);
    CORE::LINALG::Matrix<3, 1> tet4elecoords(true);
    CORE::LINALG::Matrix<3, 1> tet4refsolution(true);
    tet4refsolution(0, 0) = 61.0 / 40.0;
    tet4refsolution(1, 0) = -199.0 / 400.0;
    tet4refsolution(2, 0) = 9.0 / 16.0;
    xicentertet4ele.PutScalar(1.0 / 4.0);
    DRT::UTILS::LocalToGlobalPositionAtXiRefConfig<3, CORE::FE::CellType::tet4>(
        tet4ele, xicentertet4ele, tet4elecoords);

    FOUR_C_EXPECT_NEAR(tet4elecoords, tet4refsolution, 1e-14);

    // get tri3 element and test it
    const DRT::Element* tri3ele = testdis_->gElement(3);
    CORE::LINALG::Matrix<2, 1> xicentertri3ele(true);
    CORE::LINALG::Matrix<3, 1> tri3elecoords(true);
    CORE::LINALG::Matrix<3, 1> tri3refsolution(true);
    tri3refsolution(0, 0) = 23.0 / 15.0;
    tri3refsolution(1, 0) = -149.0 / 300.0;
    tri3refsolution(2, 0) = 1.0 / 12.0;
    xicentertri3ele.PutScalar(1.0 / 3.0);
    DRT::UTILS::LocalToGlobalPositionAtXiRefConfig<3, CORE::FE::CellType::tri3>(
        tri3ele, xicentertri3ele, tri3elecoords);

    FOUR_C_EXPECT_NEAR(tri3elecoords, tri3refsolution, 1e-14);
  }

  TEST_F(UtilsRefConfigTest, ComputeUnitNormalAtXiRefConfig)
  {
    // get quad4 element and test it
    const DRT::Element* quad4ele = testdis_->gElement(1);
    CORE::LINALG::Matrix<2, 1> xicenterquad4ele(true);
    CORE::LINALG::Matrix<3, 1> quad4elecoords(true);
    CORE::LINALG::Matrix<3, 1> quad4refsolution(true);
    quad4refsolution(0, 0) = -0.29138926578643;
    quad4refsolution(1, 0) = -0.40854577471087;
    quad4refsolution(2, 0) = 0.86497551742829;
    DRT::UTILS::ComputeUnitNormalAtXiRefConfig<CORE::FE::CellType::quad4>(
        quad4ele, xicenterquad4ele, quad4elecoords);

    FOUR_C_EXPECT_NEAR(quad4elecoords, quad4refsolution, 1e-14);

    // get tri3 element and test it
    const DRT::Element* tri3ele = testdis_->gElement(3);
    CORE::LINALG::Matrix<2, 1> xicentertri3ele(true);
    CORE::LINALG::Matrix<3, 1> tri3elecoords(true);
    CORE::LINALG::Matrix<3, 1> tri3refsolution(true);
    tri3refsolution(0, 0) = -0.085623542490578;
    tri3refsolution(1, 0) = 0.048198682858935;
    tri3refsolution(2, 0) = -0.995161040205065;
    xicentertri3ele.PutScalar(1.0 / 3.0);
    DRT::UTILS::ComputeUnitNormalAtXiRefConfig<CORE::FE::CellType::tri3>(
        tri3ele, xicentertri3ele, tri3elecoords);

    FOUR_C_EXPECT_NEAR(tri3elecoords, tri3refsolution, 1e-14);
  }
}  // namespace