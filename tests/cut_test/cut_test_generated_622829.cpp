// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_cut_combintersection.hpp"
#include "4C_cut_levelsetintersection.hpp"
#include "4C_cut_meshintersection.hpp"
#include "4C_cut_options.hpp"
#include "4C_cut_side.hpp"
#include "4C_cut_tetmeshintersection.hpp"
#include "4C_cut_volumecell.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "cut_test_utils.hpp"

void test_generated_622829()
{
  Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;
  std::vector<double> lsvs(8);
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.4891304347826088694;
    tri3_xyze(1, 0) = -0.17291304347826086385;
    tri3_xyze(2, 0) = -0.020000000000000000416;
    nids.push_back(68559);
    tri3_xyze(0, 1) = 1.4891304347826086474;
    tri3_xyze(1, 1) = -0.1620434782608695945;
    tri3_xyze(2, 1) = -0.020000000000000000416;
    nids.push_back(68575);
    tri3_xyze(0, 2) = 1.4945652173913044347;
    tri3_xyze(1, 2) = -0.16747826086956524305;
    tri3_xyze(2, 2) = -0.020000000000000000416;
    nids.push_back(-11259);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.17291304347826089161;
    tri3_xyze(2, 0) = -0.0066666666666666618904;
    nids.push_back(69306);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.16204347826086953899;
    tri3_xyze(2, 1) = -0.0066666666666666636251;
    nids.push_back(69322);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.1674782608695652153;
    tri3_xyze(2, 2) = -0.0053333333333333314205;
    nids.push_back(-11465);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.16204347826086953899;
    tri3_xyze(2, 0) = -0.0066666666666666636251;
    nids.push_back(69322);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.16204347826086953899;
    tri3_xyze(2, 1) = -0.0040000000000000000833;
    nids.push_back(69321);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.1674782608695652153;
    tri3_xyze(2, 2) = -0.0053333333333333314205;
    nids.push_back(-11465);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.17291304347826083609;
    tri3_xyze(2, 0) = -0.0093333333333333358406;
    nids.push_back(69307);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.16204347826086953899;
    tri3_xyze(2, 1) = -0.0093333333333333306364;
    nids.push_back(69323);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.16747826086956518754;
    tri3_xyze(2, 2) = -0.0079999999999999966971;
    nids.push_back(-11466);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.16204347826086953899;
    tri3_xyze(2, 0) = -0.0093333333333333306364;
    nids.push_back(69323);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.16204347826086953899;
    tri3_xyze(2, 1) = -0.0066666666666666636251;
    nids.push_back(69322);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.16747826086956518754;
    tri3_xyze(2, 2) = -0.0079999999999999966971;
    nids.push_back(-11466);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.16204347826086953899;
    tri3_xyze(2, 0) = -0.0066666666666666636251;
    nids.push_back(69322);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.17291304347826089161;
    tri3_xyze(2, 1) = -0.0066666666666666618904;
    nids.push_back(69306);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.16747826086956518754;
    tri3_xyze(2, 2) = -0.0079999999999999966971;
    nids.push_back(-11466);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.17291304347826086385;
    tri3_xyze(2, 0) = -0.01199999999999999678;
    nids.push_back(69308);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.16204347826086956674;
    tri3_xyze(2, 1) = -0.011999999999999995046;
    nids.push_back(69324);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.16747826086956518754;
    tri3_xyze(2, 2) = -0.010666666666666664576;
    nids.push_back(-11467);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.16204347826086956674;
    tri3_xyze(2, 0) = -0.011999999999999995046;
    nids.push_back(69324);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.16204347826086953899;
    tri3_xyze(2, 1) = -0.0093333333333333306364;
    nids.push_back(69323);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.16747826086956518754;
    tri3_xyze(2, 2) = -0.010666666666666664576;
    nids.push_back(-11467);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.16204347826086953899;
    tri3_xyze(2, 0) = -0.0093333333333333306364;
    nids.push_back(69323);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.17291304347826083609;
    tri3_xyze(2, 1) = -0.0093333333333333358406;
    nids.push_back(69307);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.16747826086956518754;
    tri3_xyze(2, 2) = -0.010666666666666664576;
    nids.push_back(-11467);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.17291304347826083609;
    tri3_xyze(2, 0) = -0.014666666666666668128;
    nids.push_back(69309);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.1620434782608695945;
    tri3_xyze(2, 1) = -0.014666666666666669863;
    nids.push_back(69325);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.16747826086956524305;
    tri3_xyze(2, 2) = -0.013333333333333332454;
    nids.push_back(-11468);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.1620434782608695945;
    tri3_xyze(2, 0) = -0.014666666666666669863;
    nids.push_back(69325);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.16204347826086956674;
    tri3_xyze(2, 1) = -0.011999999999999995046;
    nids.push_back(69324);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.16747826086956524305;
    tri3_xyze(2, 2) = -0.013333333333333332454;
    nids.push_back(-11468);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.16204347826086956674;
    tri3_xyze(2, 0) = -0.011999999999999995046;
    nids.push_back(69324);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.17291304347826086385;
    tri3_xyze(2, 1) = -0.01199999999999999678;
    nids.push_back(69308);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.16747826086956524305;
    tri3_xyze(2, 2) = -0.013333333333333332454;
    nids.push_back(-11468);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.17291304347826086385;
    tri3_xyze(2, 0) = -0.017333333333333332538;
    nids.push_back(69310);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.1620434782608695945;
    tri3_xyze(2, 1) = -0.017333333333333325599;
    nids.push_back(69326);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.16747826086956524305;
    tri3_xyze(2, 2) = -0.016000000000000000333;
    nids.push_back(-11469);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.1620434782608695945;
    tri3_xyze(2, 0) = -0.017333333333333325599;
    nids.push_back(69326);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.1620434782608695945;
    tri3_xyze(2, 1) = -0.014666666666666669863;
    nids.push_back(69325);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.16747826086956524305;
    tri3_xyze(2, 2) = -0.016000000000000000333;
    nids.push_back(-11469);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.1620434782608695945;
    tri3_xyze(2, 0) = -0.014666666666666669863;
    nids.push_back(69325);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.17291304347826083609;
    tri3_xyze(2, 1) = -0.014666666666666668128;
    nids.push_back(69309);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.16747826086956524305;
    tri3_xyze(2, 2) = -0.016000000000000000333;
    nids.push_back(-11469);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.17291304347826086385;
    tri3_xyze(2, 0) = -0.020000000000000000416;
    nids.push_back(69311);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.16204347826086956674;
    tri3_xyze(2, 1) = -0.020000000000000000416;
    nids.push_back(69327);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.1674782608695652153;
    tri3_xyze(2, 2) = -0.018666666666666664742;
    nids.push_back(-11470);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.16204347826086956674;
    tri3_xyze(2, 0) = -0.020000000000000000416;
    nids.push_back(69327);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.1620434782608695945;
    tri3_xyze(2, 1) = -0.017333333333333325599;
    nids.push_back(69326);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.1674782608695652153;
    tri3_xyze(2, 2) = -0.018666666666666664742;
    nids.push_back(-11470);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.1620434782608695945;
    tri3_xyze(2, 0) = -0.017333333333333325599;
    nids.push_back(69326);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.17291304347826086385;
    tri3_xyze(2, 1) = -0.017333333333333332538;
    nids.push_back(69310);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.1674782608695652153;
    tri3_xyze(2, 2) = -0.018666666666666664742;
    nids.push_back(-11470);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.15117391304347826964;
    tri3_xyze(2, 0) = -0.004000000000000001818;
    nids.push_back(69337);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.15117391304347826964;
    tri3_xyze(2, 1) = -0.0066666666666666662272;
    nids.push_back(69338);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.14573913043478259333;
    tri3_xyze(2, 2) = -0.0053333333333333331552;
    nids.push_back(-11495);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.15117391304347826964;
    tri3_xyze(2, 0) = -0.0066666666666666662272;
    nids.push_back(69338);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.14030434782608697253;
    tri3_xyze(2, 1) = -0.0066666666666666644925;
    nids.push_back(69354);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.14573913043478259333;
    tri3_xyze(2, 2) = -0.0053333333333333331552;
    nids.push_back(-11495);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.14030434782608697253;
    tri3_xyze(2, 0) = -0.0066666666666666644925;
    nids.push_back(69354);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.14030434782608691702;
    tri3_xyze(2, 1) = -0.0040000000000000009506;
    nids.push_back(69353);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.14573913043478259333;
    tri3_xyze(2, 2) = -0.0053333333333333331552;
    nids.push_back(-11495);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.15117391304347826964;
    tri3_xyze(2, 0) = -0.0066666666666666662272;
    nids.push_back(69338);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.15117391304347824188;
    tri3_xyze(2, 1) = -0.0093333333333333341059;
    nids.push_back(69339);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.14573913043478262108;
    tri3_xyze(2, 2) = -0.0080000000000000001665;
    nids.push_back(-11496);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.15117391304347824188;
    tri3_xyze(2, 0) = -0.0093333333333333341059;
    nids.push_back(69339);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.14030434782608697253;
    tri3_xyze(2, 1) = -0.0093333333333333341059;
    nids.push_back(69355);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.14573913043478262108;
    tri3_xyze(2, 2) = -0.0080000000000000001665;
    nids.push_back(-11496);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.14030434782608697253;
    tri3_xyze(2, 0) = -0.0093333333333333341059;
    nids.push_back(69355);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.14030434782608697253;
    tri3_xyze(2, 1) = -0.0066666666666666644925;
    nids.push_back(69354);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.14573913043478262108;
    tri3_xyze(2, 2) = -0.0080000000000000001665;
    nids.push_back(-11496);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.14030434782608697253;
    tri3_xyze(2, 0) = -0.0066666666666666644925;
    nids.push_back(69354);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.15117391304347826964;
    tri3_xyze(2, 1) = -0.0066666666666666662272;
    nids.push_back(69338);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.14573913043478262108;
    tri3_xyze(2, 2) = -0.0080000000000000001665;
    nids.push_back(-11496);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.15117391304347824188;
    tri3_xyze(2, 0) = -0.0093333333333333341059;
    nids.push_back(69339);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.15117391304347826964;
    tri3_xyze(2, 1) = -0.01199999999999999678;
    nids.push_back(69340);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.14573913043478259333;
    tri3_xyze(2, 2) = -0.010666666666666664576;
    nids.push_back(-11497);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.15117391304347826964;
    tri3_xyze(2, 0) = -0.01199999999999999678;
    nids.push_back(69340);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.14030434782608691702;
    tri3_xyze(2, 1) = -0.01200000000000000025;
    nids.push_back(69356);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.14573913043478259333;
    tri3_xyze(2, 2) = -0.010666666666666664576;
    nids.push_back(-11497);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.14030434782608691702;
    tri3_xyze(2, 0) = -0.01200000000000000025;
    nids.push_back(69356);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.14030434782608697253;
    tri3_xyze(2, 1) = -0.0093333333333333341059;
    nids.push_back(69355);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.14573913043478259333;
    tri3_xyze(2, 2) = -0.010666666666666664576;
    nids.push_back(-11497);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.14030434782608697253;
    tri3_xyze(2, 0) = -0.0093333333333333341059;
    nids.push_back(69355);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.15117391304347824188;
    tri3_xyze(2, 1) = -0.0093333333333333341059;
    nids.push_back(69339);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.14573913043478259333;
    tri3_xyze(2, 2) = -0.010666666666666664576;
    nids.push_back(-11497);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.15117391304347826964;
    tri3_xyze(2, 0) = -0.01199999999999999678;
    nids.push_back(69340);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.15117391304347826964;
    tri3_xyze(2, 1) = -0.014666666666666664659;
    nids.push_back(69341);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.14573913043478259333;
    tri3_xyze(2, 2) = -0.013333333333333332454;
    nids.push_back(-11498);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.15117391304347826964;
    tri3_xyze(2, 0) = -0.014666666666666664659;
    nids.push_back(69341);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.14030434782608694477;
    tri3_xyze(2, 1) = -0.014666666666666668128;
    nids.push_back(69357);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.14573913043478259333;
    tri3_xyze(2, 2) = -0.013333333333333332454;
    nids.push_back(-11498);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.14030434782608694477;
    tri3_xyze(2, 0) = -0.014666666666666668128;
    nids.push_back(69357);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.14030434782608691702;
    tri3_xyze(2, 1) = -0.01200000000000000025;
    nids.push_back(69356);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.14573913043478259333;
    tri3_xyze(2, 2) = -0.013333333333333332454;
    nids.push_back(-11498);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.14030434782608691702;
    tri3_xyze(2, 0) = -0.01200000000000000025;
    nids.push_back(69356);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.15117391304347826964;
    tri3_xyze(2, 1) = -0.01199999999999999678;
    nids.push_back(69340);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.14573913043478259333;
    tri3_xyze(2, 2) = -0.013333333333333332454;
    nids.push_back(-11498);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.15117391304347826964;
    tri3_xyze(2, 0) = -0.014666666666666664659;
    nids.push_back(69341);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.15117391304347824188;
    tri3_xyze(2, 1) = -0.017333333333333325599;
    nids.push_back(69342);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.14573913043478259333;
    tri3_xyze(2, 2) = -0.015999999999999993394;
    nids.push_back(-11499);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.15117391304347824188;
    tri3_xyze(2, 0) = -0.017333333333333325599;
    nids.push_back(69342);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.14030434782608688926;
    tri3_xyze(2, 1) = -0.017333333333333325599;
    nids.push_back(69358);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.14573913043478259333;
    tri3_xyze(2, 2) = -0.015999999999999993394;
    nids.push_back(-11499);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.14030434782608688926;
    tri3_xyze(2, 0) = -0.017333333333333325599;
    nids.push_back(69358);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.14030434782608694477;
    tri3_xyze(2, 1) = -0.014666666666666668128;
    nids.push_back(69357);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.14573913043478259333;
    tri3_xyze(2, 2) = -0.015999999999999993394;
    nids.push_back(-11499);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.14030434782608694477;
    tri3_xyze(2, 0) = -0.014666666666666668128;
    nids.push_back(69357);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.15117391304347826964;
    tri3_xyze(2, 1) = -0.014666666666666664659;
    nids.push_back(69341);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.14573913043478259333;
    tri3_xyze(2, 2) = -0.015999999999999993394;
    nids.push_back(-11499);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.16204347826086953899;
    tri3_xyze(2, 0) = -0.0040000000000000000833;
    nids.push_back(69321);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.16204347826086953899;
    tri3_xyze(2, 1) = -0.0066666666666666636251;
    nids.push_back(69322);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.15660869565217389043;
    tri3_xyze(2, 2) = -0.0053333333333333340226;
    nids.push_back(-11480);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.16204347826086953899;
    tri3_xyze(2, 0) = -0.0066666666666666636251;
    nids.push_back(69322);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.15117391304347826964;
    tri3_xyze(2, 1) = -0.0066666666666666662272;
    nids.push_back(69338);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.15660869565217389043;
    tri3_xyze(2, 2) = -0.0053333333333333340226;
    nids.push_back(-11480);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.15117391304347826964;
    tri3_xyze(2, 0) = -0.0066666666666666662272;
    nids.push_back(69338);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.15117391304347826964;
    tri3_xyze(2, 1) = -0.004000000000000001818;
    nids.push_back(69337);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.15660869565217389043;
    tri3_xyze(2, 2) = -0.0053333333333333340226;
    nids.push_back(-11480);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.16204347826086953899;
    tri3_xyze(2, 0) = -0.0066666666666666636251;
    nids.push_back(69322);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.16204347826086953899;
    tri3_xyze(2, 1) = -0.0093333333333333306364;
    nids.push_back(69323);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.15660869565217389043;
    tri3_xyze(2, 2) = -0.0079999999999999984318;
    nids.push_back(-11481);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.16204347826086953899;
    tri3_xyze(2, 0) = -0.0093333333333333306364;
    nids.push_back(69323);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.15117391304347824188;
    tri3_xyze(2, 1) = -0.0093333333333333341059;
    nids.push_back(69339);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.15660869565217389043;
    tri3_xyze(2, 2) = -0.0079999999999999984318;
    nids.push_back(-11481);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.15117391304347824188;
    tri3_xyze(2, 0) = -0.0093333333333333341059;
    nids.push_back(69339);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.15117391304347826964;
    tri3_xyze(2, 1) = -0.0066666666666666662272;
    nids.push_back(69338);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.15660869565217389043;
    tri3_xyze(2, 2) = -0.0079999999999999984318;
    nids.push_back(-11481);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.15117391304347826964;
    tri3_xyze(2, 0) = -0.0066666666666666662272;
    nids.push_back(69338);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.16204347826086953899;
    tri3_xyze(2, 1) = -0.0066666666666666636251;
    nids.push_back(69322);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.15660869565217389043;
    tri3_xyze(2, 2) = -0.0079999999999999984318;
    nids.push_back(-11481);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.16204347826086953899;
    tri3_xyze(2, 0) = -0.0093333333333333306364;
    nids.push_back(69323);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.16204347826086956674;
    tri3_xyze(2, 1) = -0.011999999999999995046;
    nids.push_back(69324);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.15660869565217391819;
    tri3_xyze(2, 2) = -0.010666666666666664576;
    nids.push_back(-11482);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.16204347826086956674;
    tri3_xyze(2, 0) = -0.011999999999999995046;
    nids.push_back(69324);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.15117391304347826964;
    tri3_xyze(2, 1) = -0.01199999999999999678;
    nids.push_back(69340);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.15660869565217391819;
    tri3_xyze(2, 2) = -0.010666666666666664576;
    nids.push_back(-11482);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.15117391304347826964;
    tri3_xyze(2, 0) = -0.01199999999999999678;
    nids.push_back(69340);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.15117391304347824188;
    tri3_xyze(2, 1) = -0.0093333333333333341059;
    nids.push_back(69339);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.15660869565217391819;
    tri3_xyze(2, 2) = -0.010666666666666664576;
    nids.push_back(-11482);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.15117391304347824188;
    tri3_xyze(2, 0) = -0.0093333333333333341059;
    nids.push_back(69339);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.16204347826086953899;
    tri3_xyze(2, 1) = -0.0093333333333333306364;
    nids.push_back(69323);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.15660869565217391819;
    tri3_xyze(2, 2) = -0.010666666666666664576;
    nids.push_back(-11482);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.16204347826086956674;
    tri3_xyze(2, 0) = -0.011999999999999995046;
    nids.push_back(69324);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.1620434782608695945;
    tri3_xyze(2, 1) = -0.014666666666666669863;
    nids.push_back(69325);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.15660869565217391819;
    tri3_xyze(2, 2) = -0.013333333333333332454;
    nids.push_back(-11483);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.1620434782608695945;
    tri3_xyze(2, 0) = -0.014666666666666669863;
    nids.push_back(69325);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.15117391304347826964;
    tri3_xyze(2, 1) = -0.014666666666666664659;
    nids.push_back(69341);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.15660869565217391819;
    tri3_xyze(2, 2) = -0.013333333333333332454;
    nids.push_back(-11483);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.15117391304347826964;
    tri3_xyze(2, 0) = -0.014666666666666664659;
    nids.push_back(69341);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.15117391304347826964;
    tri3_xyze(2, 1) = -0.01199999999999999678;
    nids.push_back(69340);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.15660869565217391819;
    tri3_xyze(2, 2) = -0.013333333333333332454;
    nids.push_back(-11483);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.15117391304347826964;
    tri3_xyze(2, 0) = -0.01199999999999999678;
    nids.push_back(69340);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.16204347826086956674;
    tri3_xyze(2, 1) = -0.011999999999999995046;
    nids.push_back(69324);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.15660869565217391819;
    tri3_xyze(2, 2) = -0.013333333333333332454;
    nids.push_back(-11483);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.1620434782608695945;
    tri3_xyze(2, 0) = -0.014666666666666669863;
    nids.push_back(69325);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.1620434782608695945;
    tri3_xyze(2, 1) = -0.017333333333333325599;
    nids.push_back(69326);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.15660869565217391819;
    tri3_xyze(2, 2) = -0.015999999999999996864;
    nids.push_back(-11484);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.1620434782608695945;
    tri3_xyze(2, 0) = -0.017333333333333325599;
    nids.push_back(69326);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.15117391304347824188;
    tri3_xyze(2, 1) = -0.017333333333333325599;
    nids.push_back(69342);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.15660869565217391819;
    tri3_xyze(2, 2) = -0.015999999999999996864;
    nids.push_back(-11484);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.15117391304347824188;
    tri3_xyze(2, 0) = -0.017333333333333325599;
    nids.push_back(69342);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.15117391304347826964;
    tri3_xyze(2, 1) = -0.014666666666666664659;
    nids.push_back(69341);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.15660869565217391819;
    tri3_xyze(2, 2) = -0.015999999999999996864;
    nids.push_back(-11484);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.15117391304347826964;
    tri3_xyze(2, 0) = -0.014666666666666664659;
    nids.push_back(69341);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.1620434782608695945;
    tri3_xyze(2, 1) = -0.014666666666666669863;
    nids.push_back(69325);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.15660869565217391819;
    tri3_xyze(2, 2) = -0.015999999999999996864;
    nids.push_back(-11484);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.1620434782608695945;
    tri3_xyze(2, 0) = -0.017333333333333325599;
    nids.push_back(69326);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.16204347826086956674;
    tri3_xyze(2, 1) = -0.020000000000000000416;
    nids.push_back(69327);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.15660869565217391819;
    tri3_xyze(2, 2) = -0.018666666666666664742;
    nids.push_back(-11485);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.16204347826086956674;
    tri3_xyze(2, 0) = -0.020000000000000000416;
    nids.push_back(69327);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.15117391304347826964;
    tri3_xyze(2, 1) = -0.020000000000000000416;
    nids.push_back(69343);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.15660869565217391819;
    tri3_xyze(2, 2) = -0.018666666666666664742;
    nids.push_back(-11485);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.15117391304347826964;
    tri3_xyze(2, 0) = -0.020000000000000000416;
    nids.push_back(69343);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.15117391304347824188;
    tri3_xyze(2, 1) = -0.017333333333333325599;
    nids.push_back(69342);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.15660869565217391819;
    tri3_xyze(2, 2) = -0.018666666666666664742;
    nids.push_back(-11485);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.15117391304347824188;
    tri3_xyze(2, 0) = -0.017333333333333325599;
    nids.push_back(69342);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.1620434782608695945;
    tri3_xyze(2, 1) = -0.017333333333333325599;
    nids.push_back(69326);
    tri3_xyze(0, 2) = 1.5;
    tri3_xyze(1, 2) = -0.15660869565217391819;
    tri3_xyze(2, 2) = -0.018666666666666664742;
    nids.push_back(-11485);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.4782608695652172948;
    tri3_xyze(1, 0) = -0.1620434782608695945;
    tri3_xyze(2, 0) = -0.020000000000000000416;
    nids.push_back(67823);
    tri3_xyze(0, 1) = 1.4891304347826086474;
    tri3_xyze(1, 1) = -0.1620434782608695945;
    tri3_xyze(2, 1) = -0.020000000000000000416;
    nids.push_back(68575);
    tri3_xyze(0, 2) = 1.4836956521739130821;
    tri3_xyze(1, 2) = -0.16747826086956524305;
    tri3_xyze(2, 2) = -0.020000000000000000416;
    nids.push_back(-11137);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.4891304347826086474;
    tri3_xyze(1, 0) = -0.1620434782608695945;
    tri3_xyze(2, 0) = -0.020000000000000000416;
    nids.push_back(68575);
    tri3_xyze(0, 1) = 1.4891304347826088694;
    tri3_xyze(1, 1) = -0.17291304347826086385;
    tri3_xyze(2, 1) = -0.020000000000000000416;
    nids.push_back(68559);
    tri3_xyze(0, 2) = 1.4836956521739130821;
    tri3_xyze(1, 2) = -0.16747826086956524305;
    tri3_xyze(2, 2) = -0.020000000000000000416;
    nids.push_back(-11137);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.4782608695652172948;
    tri3_xyze(1, 0) = -0.15117391304347826964;
    tri3_xyze(2, 0) = -0.020000000000000000416;
    nids.push_back(67839);
    tri3_xyze(0, 1) = 1.4891304347826084253;
    tri3_xyze(1, 1) = -0.15117391304347826964;
    tri3_xyze(2, 1) = -0.020000000000000000416;
    nids.push_back(68591);
    tri3_xyze(0, 2) = 1.4836956521739128601;
    tri3_xyze(1, 2) = -0.15660869565217394594;
    tri3_xyze(2, 2) = -0.020000000000000000416;
    nids.push_back(-11139);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.4891304347826084253;
    tri3_xyze(1, 0) = -0.15117391304347826964;
    tri3_xyze(2, 0) = -0.020000000000000000416;
    nids.push_back(68591);
    tri3_xyze(0, 1) = 1.4891304347826086474;
    tri3_xyze(1, 1) = -0.1620434782608695945;
    tri3_xyze(2, 1) = -0.020000000000000000416;
    nids.push_back(68575);
    tri3_xyze(0, 2) = 1.4836956521739128601;
    tri3_xyze(1, 2) = -0.15660869565217394594;
    tri3_xyze(2, 2) = -0.020000000000000000416;
    nids.push_back(-11139);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.4891304347826086474;
    tri3_xyze(1, 0) = -0.1620434782608695945;
    tri3_xyze(2, 0) = -0.020000000000000000416;
    nids.push_back(68575);
    tri3_xyze(0, 1) = 1.4782608695652172948;
    tri3_xyze(1, 1) = -0.1620434782608695945;
    tri3_xyze(2, 1) = -0.020000000000000000416;
    nids.push_back(67823);
    tri3_xyze(0, 2) = 1.4836956521739128601;
    tri3_xyze(1, 2) = -0.15660869565217394594;
    tri3_xyze(2, 2) = -0.020000000000000000416;
    nids.push_back(-11139);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.4891304347826086474;
    tri3_xyze(1, 0) = -0.1620434782608695945;
    tri3_xyze(2, 0) = -0.020000000000000000416;
    nids.push_back(68575);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.16204347826086956674;
    tri3_xyze(2, 1) = -0.020000000000000000416;
    nids.push_back(69327);
    tri3_xyze(0, 2) = 1.4945652173913044347;
    tri3_xyze(1, 2) = -0.16747826086956524305;
    tri3_xyze(2, 2) = -0.020000000000000000416;
    nids.push_back(-11259);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.16204347826086956674;
    tri3_xyze(2, 0) = -0.020000000000000000416;
    nids.push_back(69327);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.17291304347826086385;
    tri3_xyze(2, 1) = -0.020000000000000000416;
    nids.push_back(69311);
    tri3_xyze(0, 2) = 1.4945652173913044347;
    tri3_xyze(1, 2) = -0.16747826086956524305;
    tri3_xyze(2, 2) = -0.020000000000000000416;
    nids.push_back(-11259);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.4891304347826086474;
    tri3_xyze(1, 0) = -0.1620434782608695945;
    tri3_xyze(2, 0) = -0.020000000000000000416;
    nids.push_back(68575);
    tri3_xyze(0, 1) = 1.4891304347826084253;
    tri3_xyze(1, 1) = -0.15117391304347826964;
    tri3_xyze(2, 1) = -0.020000000000000000416;
    nids.push_back(68591);
    tri3_xyze(0, 2) = 1.4945652173913042127;
    tri3_xyze(1, 2) = -0.15660869565217394594;
    tri3_xyze(2, 2) = -0.020000000000000000416;
    nids.push_back(-11261);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.4891304347826084253;
    tri3_xyze(1, 0) = -0.15117391304347826964;
    tri3_xyze(2, 0) = -0.020000000000000000416;
    nids.push_back(68591);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.15117391304347826964;
    tri3_xyze(2, 1) = -0.020000000000000000416;
    nids.push_back(69343);
    tri3_xyze(0, 2) = 1.4945652173913042127;
    tri3_xyze(1, 2) = -0.15660869565217394594;
    tri3_xyze(2, 2) = -0.020000000000000000416;
    nids.push_back(-11261);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.15117391304347826964;
    tri3_xyze(2, 0) = -0.020000000000000000416;
    nids.push_back(69343);
    tri3_xyze(0, 1) = 1.5;
    tri3_xyze(1, 1) = -0.16204347826086956674;
    tri3_xyze(2, 1) = -0.020000000000000000416;
    nids.push_back(69327);
    tri3_xyze(0, 2) = 1.4945652173913042127;
    tri3_xyze(1, 2) = -0.15660869565217394594;
    tri3_xyze(2, 2) = -0.020000000000000000416;
    nids.push_back(-11261);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.5;
    tri3_xyze(1, 0) = -0.16204347826086956674;
    tri3_xyze(2, 0) = -0.020000000000000000416;
    nids.push_back(69327);
    tri3_xyze(0, 1) = 1.4891304347826086474;
    tri3_xyze(1, 1) = -0.1620434782608695945;
    tri3_xyze(2, 1) = -0.020000000000000000416;
    nids.push_back(68575);
    tri3_xyze(0, 2) = 1.4945652173913042127;
    tri3_xyze(1, 2) = -0.15660869565217394594;
    tri3_xyze(2, 2) = -0.020000000000000000416;
    nids.push_back(-11261);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = 1.4888837989634799985;
    hex8_xyze(1, 0) = -0.17777693745167399975;
    hex8_xyze(2, 0) = -0.0055550085050327501629;
    nids.push_back(1653793);
    hex8_xyze(0, 1) = 1.4888838267851700614;
    hex8_xyze(1, 1) = -0.17777695280934899258;
    hex8_xyze(2, 1) = -0.016665064572966199058;
    nids.push_back(1653794);
    hex8_xyze(0, 2) = 1.4888841561094099397;
    hex8_xyze(1, 2) = -0.16666619296931300953;
    hex8_xyze(2, 2) = -0.016665065045985600484;
    nids.push_back(1653797);
    hex8_xyze(0, 3) = 1.4888841103540999544;
    hex8_xyze(1, 3) = -0.16666617600415301048;
    hex8_xyze(2, 3) = -0.005555010668345490045;
    nids.push_back(1653796);
    hex8_xyze(0, 4) = 1.4999993064943799581;
    hex8_xyze(1, 4) = -0.17777652881863101331;
    hex8_xyze(2, 4) = -0.0055550171282643398887;
    nids.push_back(1653802);
    hex8_xyze(0, 5) = 1.4999994621857399846;
    hex8_xyze(1, 5) = -0.17777654050145599851;
    hex8_xyze(2, 5) = -0.016665059052471998396;
    nids.push_back(1653803);
    hex8_xyze(0, 6) = 1.5000000240384498973;
    hex8_xyze(1, 6) = -0.16666591955666298919;
    hex8_xyze(2, 6) = -0.016665059299914599528;
    nids.push_back(1653806);
    hex8_xyze(0, 7) = 1.4999998115573600632;
    hex8_xyze(1, 7) = -0.16666590305784700909;
    hex8_xyze(2, 7) = -0.0055550208030329603637;
    nids.push_back(1653805);

    intersection.add_element(622295, nids, hex8_xyze, Core::FE::CellType::hex8);
  }

  {
    Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = 1.4777663419814799362;
    hex8_xyze(1, 0) = -0.16666637043361701154;
    hex8_xyze(2, 0) = -0.0055550106717242699397;
    nids.push_back(1653787);
    hex8_xyze(0, 1) = 1.4777662917684499799;
    hex8_xyze(1, 1) = -0.16666638549249698786;
    hex8_xyze(2, 1) = -0.016665084007053600906;
    nids.push_back(1653788);
    hex8_xyze(0, 2) = 1.4777664314381200317;
    hex8_xyze(1, 2) = -0.15555542022984999995;
    hex8_xyze(2, 2) = -0.016665084382067499313;
    nids.push_back(1654340);
    hex8_xyze(0, 3) = 1.4777664823350400436;
    hex8_xyze(1, 3) = -0.15555540311010399024;
    hex8_xyze(2, 3) = -0.0055550086982320704201;
    nids.push_back(1654339);
    hex8_xyze(0, 4) = 1.4888841103540999544;
    hex8_xyze(1, 4) = -0.16666617600415301048;
    hex8_xyze(2, 4) = -0.005555010668345490045;
    nids.push_back(1653796);
    hex8_xyze(0, 5) = 1.4888841561094099397;
    hex8_xyze(1, 5) = -0.16666619296931300953;
    hex8_xyze(2, 5) = -0.016665065045985600484;
    nids.push_back(1653797);
    hex8_xyze(0, 6) = 1.4888842858838400307;
    hex8_xyze(1, 6) = -0.15555532048654499566;
    hex8_xyze(2, 6) = -0.016665065474201001122;
    nids.push_back(1654349);
    hex8_xyze(0, 7) = 1.4888842572649398921;
    hex8_xyze(1, 7) = -0.15555530184635499302;
    hex8_xyze(2, 7) = -0.0055550089610411400656;
    nids.push_back(1654348);

    intersection.add_element(622820, nids, hex8_xyze, Core::FE::CellType::hex8);
  }

  {
    Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = 1.4888841103953200928;
    hex8_xyze(1, 0) = -0.16666615692926001202;
    hex8_xyze(2, 0) = 0.0055550113112029396098;
    nids.push_back(1653795);
    hex8_xyze(0, 1) = 1.4888841103540999544;
    hex8_xyze(1, 1) = -0.16666617600415301048;
    hex8_xyze(2, 1) = -0.005555010668345490045;
    nids.push_back(1653796);
    hex8_xyze(0, 2) = 1.4888842572649398921;
    hex8_xyze(1, 2) = -0.15555530184635499302;
    hex8_xyze(2, 2) = -0.0055550089610411400656;
    nids.push_back(1654348);
    hex8_xyze(0, 3) = 1.4888842572979799073;
    hex8_xyze(1, 3) = -0.15555528272000199164;
    hex8_xyze(2, 3) = 0.0055550092788237302385;
    nids.push_back(1654347);
    hex8_xyze(0, 4) = 1.499999811793909954;
    hex8_xyze(1, 4) = -0.16666588367014800731;
    hex8_xyze(2, 4) = 0.005555022038289450341;
    nids.push_back(1653804);
    hex8_xyze(0, 5) = 1.4999998115573600632;
    hex8_xyze(1, 5) = -0.16666590305784700909;
    hex8_xyze(2, 5) = -0.0055550208030329603637;
    nids.push_back(1653805);
    hex8_xyze(0, 6) = 1.4999999278943199066;
    hex8_xyze(1, 6) = -0.15555513430187101198;
    hex8_xyze(2, 6) = -0.0055550176470937602871;
    nids.push_back(1654357);
    hex8_xyze(0, 7) = 1.4999999280317299899;
    hex8_xyze(1, 7) = -0.15555511477299099887;
    hex8_xyze(2, 7) = 0.0055550181922827601352;
    nids.push_back(1654356);

    intersection.add_element(622828, nids, hex8_xyze, Core::FE::CellType::hex8);
  }

  {
    Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = 1.4888841103540999544;
    hex8_xyze(1, 0) = -0.16666617600415301048;
    hex8_xyze(2, 0) = -0.005555010668345490045;
    nids.push_back(1653796);
    hex8_xyze(0, 1) = 1.4888841561094099397;
    hex8_xyze(1, 1) = -0.16666619296931300953;
    hex8_xyze(2, 1) = -0.016665065045985600484;
    nids.push_back(1653797);
    hex8_xyze(0, 2) = 1.4888842858838400307;
    hex8_xyze(1, 2) = -0.15555532048654499566;
    hex8_xyze(2, 2) = -0.016665065474201001122;
    nids.push_back(1654349);
    hex8_xyze(0, 3) = 1.4888842572649398921;
    hex8_xyze(1, 3) = -0.15555530184635499302;
    hex8_xyze(2, 3) = -0.0055550089610411400656;
    nids.push_back(1654348);
    hex8_xyze(0, 4) = 1.4999998115573600632;
    hex8_xyze(1, 4) = -0.16666590305784700909;
    hex8_xyze(2, 4) = -0.0055550208030329603637;
    nids.push_back(1653805);
    hex8_xyze(0, 5) = 1.5000000240384498973;
    hex8_xyze(1, 5) = -0.16666591955666298919;
    hex8_xyze(2, 5) = -0.016665059299914599528;
    nids.push_back(1653806);
    hex8_xyze(0, 6) = 1.5000000843443499488;
    hex8_xyze(1, 6) = -0.15555515548912299262;
    hex8_xyze(2, 6) = -0.016665059622743800399;
    nids.push_back(1654358);
    hex8_xyze(0, 7) = 1.4999999278943199066;
    hex8_xyze(1, 7) = -0.15555513430187101198;
    hex8_xyze(2, 7) = -0.0055550176470937602871;
    nids.push_back(1654357);

    intersection.add_element(622829, nids, hex8_xyze, Core::FE::CellType::hex8);
  }

  {
    Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = 1.4888842572649398921;
    hex8_xyze(1, 0) = -0.15555530184635499302;
    hex8_xyze(2, 0) = -0.0055550089610411400656;
    nids.push_back(1654348);
    hex8_xyze(0, 1) = 1.4888842858838400307;
    hex8_xyze(1, 1) = -0.15555532048654499566;
    hex8_xyze(2, 1) = -0.016665065474201001122;
    nids.push_back(1654349);
    hex8_xyze(0, 2) = 1.4888843921406000881;
    hex8_xyze(1, 2) = -0.14444430714100200963;
    hex8_xyze(2, 2) = -0.016665065715478599717;
    nids.push_back(1654352);
    hex8_xyze(0, 3) = 1.4888843634408799321;
    hex8_xyze(1, 3) = -0.14444429163695798879;
    hex8_xyze(2, 3) = -0.0055550091920389997949;
    nids.push_back(1654351);
    hex8_xyze(0, 4) = 1.4999999278943199066;
    hex8_xyze(1, 4) = -0.15555513430187101198;
    hex8_xyze(2, 4) = -0.0055550176470937602871;
    nids.push_back(1654357);
    hex8_xyze(0, 5) = 1.5000000843443499488;
    hex8_xyze(1, 5) = -0.15555515548912299262;
    hex8_xyze(2, 5) = -0.016665059622743800399;
    nids.push_back(1654358);
    hex8_xyze(0, 6) = 1.5000002376542600491;
    hex8_xyze(1, 6) = -0.14444420306060198889;
    hex8_xyze(2, 6) = -0.016665060125803200092;
    nids.push_back(1654361);
    hex8_xyze(0, 7) = 1.500000081045760103;
    hex8_xyze(1, 7) = -0.14444419147942599846;
    hex8_xyze(2, 7) = -0.0055550180418534201437;
    nids.push_back(1654360);

    intersection.add_element(622832, nids, hex8_xyze, Core::FE::CellType::hex8);
  }

  {
    Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = 1.4888841561094099397;
    hex8_xyze(1, 0) = -0.16666619296931300953;
    hex8_xyze(2, 0) = -0.016665065045985600484;
    nids.push_back(1653797);
    hex8_xyze(0, 1) = 1.4888842535221900043;
    hex8_xyze(1, 1) = -0.16666620405004300975;
    hex8_xyze(2, 1) = -0.027775221752330700453;
    nids.push_back(1653822);
    hex8_xyze(0, 2) = 1.488884401433290039;
    hex8_xyze(1, 2) = -0.15555533015350000992;
    hex8_xyze(2, 2) = -0.027775224281203998722;
    nids.push_back(1654374);
    hex8_xyze(0, 3) = 1.4888842858838400307;
    hex8_xyze(1, 3) = -0.15555532048654499566;
    hex8_xyze(2, 3) = -0.016665065474201001122;
    nids.push_back(1654349);
    hex8_xyze(0, 4) = 1.5000000240384498973;
    hex8_xyze(1, 4) = -0.16666591955666298919;
    hex8_xyze(2, 4) = -0.016665059299914599528;
    nids.push_back(1653806);
    hex8_xyze(0, 5) = 1.5000003147648799384;
    hex8_xyze(1, 5) = -0.16666592794941101352;
    hex8_xyze(2, 5) = -0.027775153874800999343;
    nids.push_back(1653831);
    hex8_xyze(0, 6) = 1.5000004316765200851;
    hex8_xyze(1, 6) = -0.1555551590795289929;
    hex8_xyze(2, 6) = -0.027775157712656900477;
    nids.push_back(1654383);
    hex8_xyze(0, 7) = 1.5000000843443499488;
    hex8_xyze(1, 7) = -0.15555515548912299262;
    hex8_xyze(2, 7) = -0.016665059622743800399;
    nids.push_back(1654358);

    intersection.add_element(622854, nids, hex8_xyze, Core::FE::CellType::hex8);
  }

  {
    Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = 1.4999998115573600632;
    hex8_xyze(1, 0) = -0.16666590305784700909;
    hex8_xyze(2, 0) = -0.0055550208030329603637;
    nids.push_back(1653805);
    hex8_xyze(0, 1) = 1.5000000240384498973;
    hex8_xyze(1, 1) = -0.16666591955666298919;
    hex8_xyze(2, 1) = -0.016665059299914599528;
    nids.push_back(1653806);
    hex8_xyze(0, 2) = 1.5000000843443499488;
    hex8_xyze(1, 2) = -0.15555515548912299262;
    hex8_xyze(2, 2) = -0.016665059622743800399;
    nids.push_back(1654358);
    hex8_xyze(0, 3) = 1.4999999278943199066;
    hex8_xyze(1, 3) = -0.15555513430187101198;
    hex8_xyze(2, 3) = -0.0055550176470937602871;
    nids.push_back(1654357);
    hex8_xyze(0, 4) = 1.5111168475913099307;
    hex8_xyze(1, 4) = -0.16666554260500698881;
    hex8_xyze(2, 4) = -0.0055549949031295899407;
    nids.push_back(1666345);
    hex8_xyze(0, 5) = 1.5111171109621499564;
    hex8_xyze(1, 5) = -0.16666555546697500723;
    hex8_xyze(2, 5) = -0.016665101355255299009;
    nids.push_back(1666346);
    hex8_xyze(0, 6) = 1.5111174408790200019;
    hex8_xyze(1, 6) = -0.15555481939484200327;
    hex8_xyze(2, 6) = -0.016665101458739898371;
    nids.push_back(1666898);
    hex8_xyze(0, 7) = 1.5111170460271099447;
    hex8_xyze(1, 7) = -0.15555485502832899769;
    hex8_xyze(2, 7) = -0.0055550495644879202203;
    nids.push_back(1666897);

    intersection.add_element(634961, nids, hex8_xyze, Core::FE::CellType::hex8);
  }

  intersection.cut_test_cut(
      true, Cut::VCellGaussPts_DirectDivergence, Cut::BCellGaussPts_Tessellation);
  intersection.cut_finalize(
      true, Cut::VCellGaussPts_DirectDivergence, Cut::BCellGaussPts_Tessellation, false, true);
}
