// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_pre_exodus_validate.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_fem_general_utils_integration.hpp"
#include "4C_global_data.hpp"
#include "4C_global_data_read.hpp"
#include "4C_global_legacy_module.hpp"
#include "4C_io_input_file.hpp"
#include "4C_linalg_utils_densematrix_multiply.hpp"
#include "4C_pre_exodus_reader.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::validate_input_file(const MPI_Comm comm, const std::string datfile)
{
  using namespace FourC;

  // access our problem instance
  Global::Problem* problem = Global::Problem::instance();

  // create a InputFile
  Core::IO::InputFile input_file = Global::set_up_input_file(comm);
  input_file.read(datfile);

  // read and validate dynamic and solver sections
  std::cout << "...Read parameters" << std::endl;
  Global::read_parameter(*problem, input_file);

  // read and validate all material definitions
  std::cout << "...Read materials" << std::endl;
  Global::read_materials(*problem, input_file);

  // do NOT allocate the different fields (discretizations) here,
  // since RAM might be a problem for huge problems!
  // But, we have to perform at least the problem-specific setup since
  // some reading procedures depend on the number of fields (e.g., ReadKnots())
  std::cout << "...Read field setup" << std::endl;
  Global::read_fields(*problem, input_file, false);  // option false is important here!

  std::cout << "...";
  {
    Core::Utils::FunctionManager function_manager;
    global_legacy_module_callbacks().AttachFunctionDefinitions(function_manager);
    function_manager.read_input(input_file);
  }

  Global::read_result(*problem, input_file);
  Global::read_conditions(*problem, input_file);

  // input of materials of cloned fields (if needed)
  Global::read_cloning_material_map(*problem, input_file);

  // read all knot information for isogeometric analysis
  // and add it to the (derived) nurbs discretization
  Global::read_knots(*problem, input_file);

  // we wait till all procs are here. Otherwise a hang up might occur where
  // one proc ended with FOUR_C_THROW but other procs were not finished and waited...
  // we also want to have the printing above being finished.
  Core::Communication::barrier(comm);
  // the input file seems to be valid
  std::cout << "...OK\n\n";
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::validate_mesh_element_jacobians(Mesh& mymesh)
{
  if (mymesh.get_num_dim() != 3) FOUR_C_THROW("Element Validation only for 3 Dimensions");

  std::map<int, std::shared_ptr<ElementBlock>> myebs = mymesh.get_element_blocks();
  std::map<int, std::shared_ptr<ElementBlock>>::iterator i_eb;

  for (i_eb = myebs.begin(); i_eb != myebs.end(); ++i_eb)
  {
    std::shared_ptr<ElementBlock> eb = i_eb->second;
    const Core::FE::CellType distype = pre_shape_to_drt(eb->get_shape());
    // check and rewind if necessary
    validate_element_jacobian(mymesh, distype, *eb);
    // full check at all gausspoints
    int invalid_dets = validate_element_jacobian_fullgp(mymesh, distype, *eb);
    if (invalid_dets > 0)
      std::cout << invalid_dets << " negative Jacobian determinants in EB of shape "
                << shape_to_string(eb->get_shape()) << std::endl;
  }
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::validate_element_jacobian(
    Mesh& mymesh, const Core::FE::CellType distype, ElementBlock& eb)
{
  using namespace FourC;

  // use one point gauss rule to calculate jacobian at element center
  Core::FE::GaussRule3D integrationrule_1point = Core::FE::GaussRule3D::undefined;
  switch (distype)
  {
    case Core::FE::CellType::hex8:
    case Core::FE::CellType::hex20:
      integrationrule_1point = Core::FE::GaussRule3D::hex_1point;
      break;
    case Core::FE::CellType::hex27:
      integrationrule_1point =
          Core::FE::GaussRule3D::hex_27point;  // one point is not enough for hex27!!
      break;
    case Core::FE::CellType::tet4:
    case Core::FE::CellType::tet10:
      integrationrule_1point = Core::FE::GaussRule3D::tet_1point;
      break;
    case Core::FE::CellType::wedge6:
    case Core::FE::CellType::wedge15:
      integrationrule_1point = Core::FE::GaussRule3D::wedge_1point;
      break;
    case Core::FE::CellType::pyramid5:
      integrationrule_1point = Core::FE::GaussRule3D::pyramid_1point;
      break;
    // do nothing for 2D, 1D and 0D elements
    case Core::FE::CellType::quad4:
    case Core::FE::CellType::quad8:
    case Core::FE::CellType::quad9:
    case Core::FE::CellType::tri3:
    case Core::FE::CellType::tri6:
    case Core::FE::CellType::line2:
    case Core::FE::CellType::line3:
    case Core::FE::CellType::point1:
      return;
    default:
      FOUR_C_THROW("Unknown element type, validation failed!");
      break;
  }
  const Core::FE::IntegrationPoints3D intpoints(integrationrule_1point);
  const int iel = eb.get_ele_nodes(0).size();
  // shape functions derivatives
  const int NSD = 3;
  Core::LinAlg::SerialDenseMatrix deriv(NSD, iel);

  // go through all elements
  std::shared_ptr<std::map<int, std::vector<int>>> eleconn = eb.get_ele_conn();
  std::map<int, std::vector<int>>::iterator i_ele;
  int numrewindedeles = 0;
  for (i_ele = eleconn->begin(); i_ele != eleconn->end(); ++i_ele)
  {
    int rewcount = 0;
    for (int igp = 0; igp < intpoints.nquad; ++igp)
    {
      Core::FE::shape_function_3d_deriv1(
          deriv, intpoints.qxg[igp][0], intpoints.qxg[igp][1], intpoints.qxg[igp][2], distype);
      if (!positive_ele(i_ele->first, i_ele->second, mymesh, deriv))
      {
        // rewind the element nodes
        if (rewcount == 0)
        {
          i_ele->second = rewind_ele(i_ele->second, distype);
          numrewindedeles++;
        }
        // double check
        if (!positive_ele(i_ele->first, i_ele->second, mymesh, deriv))
          FOUR_C_THROW(
              "No proper rewinding for element id {} at gauss point {}", i_ele->first, igp);
        rewcount++;
      }
    }
  }
  if (numrewindedeles > 0)
    std::cout << "...Successfully rewinded " << numrewindedeles
              << " elements. For details see *.err file" << std::endl;

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int EXODUS::validate_element_jacobian_fullgp(
    Mesh& mymesh, const Core::FE::CellType distype, ElementBlock& eb)
{
  using namespace FourC;

  Core::FE::GaussRule3D integrationrule = Core::FE::GaussRule3D::undefined;
  switch (distype)
  {
    case Core::FE::CellType::hex8:
      integrationrule = Core::FE::GaussRule3D::hex_8point;
      break;
    case Core::FE::CellType::hex20:
      integrationrule = Core::FE::GaussRule3D::hex_27point;
      break;
    case Core::FE::CellType::hex27:
      integrationrule = Core::FE::GaussRule3D::hex_27point;
      break;
    case Core::FE::CellType::tet4:
      integrationrule = Core::FE::GaussRule3D::tet_4point;
      break;
    case Core::FE::CellType::tet10:
      integrationrule = Core::FE::GaussRule3D::tet_10point;
      break;
    case Core::FE::CellType::wedge6:
    case Core::FE::CellType::wedge15:
      integrationrule = Core::FE::GaussRule3D::wedge_6point;
      break;
    case Core::FE::CellType::pyramid5:
      integrationrule = Core::FE::GaussRule3D::pyramid_8point;
      break;
    // do nothing for 2D, 1D and 0D elements
    case Core::FE::CellType::quad4:
    case Core::FE::CellType::quad8:
    case Core::FE::CellType::quad9:
    case Core::FE::CellType::tri3:
    case Core::FE::CellType::tri6:
    case Core::FE::CellType::line2:
    case Core::FE::CellType::line3:
    case Core::FE::CellType::point1:
      return 0;
    default:
      FOUR_C_THROW("Unknown element type, validation failed!");
      break;
  }
  const Core::FE::IntegrationPoints3D intpoints(integrationrule);
  const int iel = eb.get_ele_nodes(0).size();
  // shape functions derivatives
  const int NSD = 3;
  Core::LinAlg::SerialDenseMatrix deriv(NSD, iel);

  // go through all elements
  int invalids = 0;
  std::shared_ptr<std::map<int, std::vector<int>>> eleconn = eb.get_ele_conn();
  std::map<int, std::vector<int>>::iterator i_ele;
  for (i_ele = eleconn->begin(); i_ele != eleconn->end(); ++i_ele)
  {
    for (int igp = 0; igp < intpoints.nquad; ++igp)
    {
      Core::FE::shape_function_3d_deriv1(
          deriv, intpoints.qxg[igp][0], intpoints.qxg[igp][1], intpoints.qxg[igp][2], distype);
      if (positive_ele(i_ele->first, i_ele->second, mymesh, deriv) == false)
      {
        invalids++;
      }
    }
  }

  return invalids;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool EXODUS::positive_ele(const int& eleid, const std::vector<int>& nodes, const Mesh& mymesh,
    const Core::LinAlg::SerialDenseMatrix& deriv)
{
  using namespace FourC;

  const int iel = deriv.numCols();
  const int NSD = deriv.numRows();
  Core::LinAlg::SerialDenseMatrix xyze(deriv.numRows(), iel);
  for (int inode = 0; inode < iel; inode++)
  {
    const std::vector<double> x = mymesh.get_node(nodes.at(inode));
    xyze(0, inode) = x[0];
    xyze(1, inode) = x[1];
    xyze(2, inode) = x[2];
  }
  // get Jacobian matrix and determinant
  // actually compute its transpose....
  if (NSD == 3)
  {
    Core::LinAlg::SerialDenseMatrix xjm(NSD, NSD);
    Core::LinAlg::multiply_nt(xjm, deriv, xyze);
    Core::LinAlg::Matrix<3, 3> jac(xjm.values(), true);
    const double det = jac.determinant();

    if (abs(det) < 1E-16)
      FOUR_C_THROW("ZERO JACOBIAN DETERMINANT FOR ELEMENT {}: DET = {}", eleid, det);

    if (det < 0.0)
    {
      return false;
    }
  }
  return true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int EXODUS::ele_sane_sign(
    const std::vector<int>& nodes, const std::map<int, std::vector<double>>& nodecoords)
{
  using namespace FourC;

  const int iel = nodes.size();
  // to be even stricter we test the Jacobian at every Node, not just at the gausspoints
  Core::LinAlg::SerialDenseMatrix local_nodecoords(iel, 3);
  Core::FE::CellType distype;
  switch (iel)
  {
    case 8:  // hex8
      local_nodecoords(0, 0) = -1.;
      local_nodecoords(0, 1) = -1.;
      local_nodecoords(0, 2) = -1.;
      local_nodecoords(1, 0) = 1.;
      local_nodecoords(1, 1) = -1.;
      local_nodecoords(1, 2) = -1.;
      local_nodecoords(2, 0) = 1.;
      local_nodecoords(2, 1) = 1.;
      local_nodecoords(2, 2) = -1.;
      local_nodecoords(3, 0) = -1.;
      local_nodecoords(3, 1) = 1.;
      local_nodecoords(3, 2) = -1.;
      local_nodecoords(4, 0) = -1.;
      local_nodecoords(4, 1) = -1.;
      local_nodecoords(4, 2) = 1.;
      local_nodecoords(5, 0) = 1.;
      local_nodecoords(5, 1) = -1.;
      local_nodecoords(5, 2) = 1.;
      local_nodecoords(6, 0) = 1.;
      local_nodecoords(6, 1) = 1.;
      local_nodecoords(6, 2) = 1.;
      local_nodecoords(7, 0) = -1.;
      local_nodecoords(7, 1) = 1.;
      local_nodecoords(7, 2) = 1.;
      distype = Core::FE::CellType::hex8;
      break;
    case 6:  // wedge6
      local_nodecoords(0, 0) = 0.;
      local_nodecoords(0, 1) = 0.;
      local_nodecoords(0, 2) = -1.;
      local_nodecoords(1, 0) = 1.;
      local_nodecoords(1, 1) = 0.;
      local_nodecoords(1, 2) = -1.;
      local_nodecoords(2, 0) = 0.;
      local_nodecoords(2, 1) = 1.;
      local_nodecoords(2, 2) = -1.;
      local_nodecoords(3, 0) = 0.;
      local_nodecoords(3, 1) = 0.;
      local_nodecoords(3, 2) = 1.;
      local_nodecoords(4, 0) = 1.;
      local_nodecoords(4, 1) = 0.;
      local_nodecoords(4, 2) = 1.;
      local_nodecoords(5, 0) = 0.;
      local_nodecoords(5, 1) = 1.;
      local_nodecoords(5, 2) = 1.;
      distype = Core::FE::CellType::wedge6;
      break;
    default:
      FOUR_C_THROW("No Element Sanity Check for this distype");
      break;
  }
  // shape functions derivatives
  const int NSD = 3;
  Core::LinAlg::SerialDenseMatrix deriv(NSD, iel);

  Core::LinAlg::SerialDenseMatrix xyze(deriv.numRows(), iel);
  for (int inode = 0; inode < iel; inode++)
  {
    const std::vector<double> x = nodecoords.find(nodes[inode])->second;
    xyze(0, inode) = x[0];
    xyze(1, inode) = x[1];
    xyze(2, inode) = x[2];
  }
  // get Jacobian matrix and determinant
  // actually compute its transpose....
  Core::LinAlg::SerialDenseMatrix xjm(NSD, NSD);
  int n_posdet = 0;
  int n_negdet = 0;

  for (int i = 0; i < iel; ++i)
  {
    Core::FE::shape_function_3d_deriv1(
        deriv, local_nodecoords(i, 0), local_nodecoords(i, 1), local_nodecoords(i, 2), distype);
    Core::LinAlg::multiply_nt(xjm, deriv, xyze);
    const double det = xjm(0, 0) * xjm(1, 1) * xjm(2, 2) + xjm(0, 1) * xjm(1, 2) * xjm(2, 0) +
                       xjm(0, 2) * xjm(1, 0) * xjm(2, 1) - xjm(0, 2) * xjm(1, 1) * xjm(2, 0) -
                       xjm(0, 0) * xjm(1, 2) * xjm(2, 1) - xjm(0, 1) * xjm(1, 0) * xjm(2, 2);
    if (abs(det) < 1E-16) FOUR_C_THROW("ZERO JACOBIAN DETERMINANT");
    if (det < 0)
    {
      ++n_negdet;
    }
    else
      ++n_posdet;
  }

  if (n_posdet == iel && n_negdet == 0)
    return 1;
  else if (n_posdet == 0 && n_negdet == iel)
    return -1;
  else
    return 0;

  return 0;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<int> EXODUS::rewind_ele(std::vector<int> old_nodeids, const Core::FE::CellType distype)
{
  using namespace FourC;

  std::vector<int> new_nodeids(old_nodeids.size());
  // rewinding of nodes to arrive at mathematically positive element
  switch (distype)
  {
    case Core::FE::CellType::tet4:
    {
      new_nodeids[0] = old_nodeids[0];
      new_nodeids[1] = old_nodeids[2];
      new_nodeids[2] = old_nodeids[1];
      new_nodeids[3] = old_nodeids[3];
      break;
    }
    case Core::FE::CellType::tet10:
    {
      new_nodeids[0] = old_nodeids[0];
      new_nodeids[1] = old_nodeids[2];
      new_nodeids[2] = old_nodeids[1];
      new_nodeids[3] = old_nodeids[3];
      new_nodeids[4] = old_nodeids[6];
      new_nodeids[5] = old_nodeids[5];
      new_nodeids[6] = old_nodeids[4];
      new_nodeids[7] = old_nodeids[7];
      new_nodeids[8] = old_nodeids[8];
      new_nodeids[9] = old_nodeids[9];
      break;
    }
    case Core::FE::CellType::hex8:
    {
      new_nodeids[0] = old_nodeids[4];
      new_nodeids[1] = old_nodeids[5];
      new_nodeids[2] = old_nodeids[6];
      new_nodeids[3] = old_nodeids[7];
      new_nodeids[4] = old_nodeids[0];
      new_nodeids[5] = old_nodeids[1];
      new_nodeids[6] = old_nodeids[2];
      new_nodeids[7] = old_nodeids[3];
      break;
    }
    case Core::FE::CellType::wedge6:
    {
      new_nodeids[0] = old_nodeids[3];
      new_nodeids[1] = old_nodeids[4];
      new_nodeids[2] = old_nodeids[5];
      new_nodeids[3] = old_nodeids[0];
      new_nodeids[4] = old_nodeids[1];
      new_nodeids[5] = old_nodeids[2];
      break;
    }
    case Core::FE::CellType::pyramid5:
    {
      new_nodeids[1] = old_nodeids[3];
      new_nodeids[3] = old_nodeids[1];
      // the other nodes can stay the same
      new_nodeids[0] = old_nodeids[0];
      new_nodeids[2] = old_nodeids[2];
      new_nodeids[4] = old_nodeids[4];
      break;
    }
    case Core::FE::CellType::hex27:
    {
      // nodes 1 - 20 can stay the same (no rewinding for hex20)
      for (int i = 0; i < 20; i++)
      {
        new_nodeids[i] = old_nodeids[i];
      }
      // rewind the nodes on the center of the 6 sides
      // and the center node of the actual hex27 element
      new_nodeids[20] = old_nodeids[21];
      new_nodeids[21] = old_nodeids[25];
      new_nodeids[22] = old_nodeids[24];
      new_nodeids[23] = old_nodeids[26];
      new_nodeids[24] = old_nodeids[23];
      new_nodeids[25] = old_nodeids[22];
      new_nodeids[26] = old_nodeids[20];
      break;
    }
    default:
      FOUR_C_THROW("no rewinding scheme for this type of element");
      break;
  }
  return new_nodeids;
}

FOUR_C_NAMESPACE_CLOSE
