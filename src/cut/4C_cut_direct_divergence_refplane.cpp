// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_cut_direct_divergence_refplane.hpp"

#include "4C_cut_kernel.hpp"
#include "4C_cut_options.hpp"
#include "4C_cut_output.hpp"
#include "4C_cut_side.hpp"
#include "4C_cut_volumecell.hpp"

FOUR_C_NAMESPACE_OPEN

/*-----------------------------------------------------------------------------------*
 * Perform all the operations related to computing reference plane           sudhakar 06/15
 *-----------------------------------------------------------------------------------*/
std::vector<double> Cut::DirectDivergenceGlobalRefplane::get_reference_plane()
{
  if (elem1_->shape() != Core::FE::CellType::hex8)
  {
    FOUR_C_THROW("Currently can handle only hexagonal family, not {}\n", elem1_->shape());
  }

  std::vector<double> RefPlaneEqn(4, 0.0);

  bool comp_ref_plane = false;

  // In the first round check if the projection of element corner points is inside
  std::vector<Point*> points = elem1_->points();

  static const int max_attempts = 9;
  double tol = 1e-8;

  for (int refplaneattempt = 0; refplaneattempt < max_attempts; ++refplaneattempt)
  {
    //---
    // First estimate -- Compute reference plane based on the facets of the volumecell information
    //---
    comp_ref_plane = facet_based_ref(RefPlaneEqn, points, tol);
    if (comp_ref_plane)
    {
      return RefPlaneEqn;
    }
    ref_pts_gmsh_.clear();

    //---
    // Second estimate -- Compute reference plane based on the diagonal information
    //---
    comp_ref_plane = diagonal_based_ref(RefPlaneEqn, points, tol);
    if (comp_ref_plane)
    {
      return RefPlaneEqn;
    }
    ref_pts_gmsh_.clear();

    //---
    // Third estimate -- Compute reference plane based on Side information
    //---
    comp_ref_plane = side_based_ref(RefPlaneEqn, points, tol);
    if (comp_ref_plane)
    {
      return RefPlaneEqn;
    }

    if (not comp_ref_plane)
    {
      std::stringstream str;
      str << ".no_refplane_" << refplaneattempt << "_CUTFAIL.pos";
      std::string filename(Cut::Output::generate_gmsh_output_filename(str.str()));
      std::ofstream file(filename.c_str());

      Cut::Output::gmsh_complete_cut_element(file, elem1_);
      Cut::Output::gmsh_new_section(file, "VolumeCell");
      Cut::Output::gmsh_volumecell_dump(file, volcell_);
      Cut::Output::gmsh_end_section(file, true);
    }

    if (!refplaneattempt)
    {
      // We couldn't find a reference plane where gauss points from the whole element are inside.
      // But anyway this is not required, since we are just integrating a vc here --> therefore
      // project all points from the vc to check if it is a good reference plane --> this is more
      // expensive depending on the number of points describing the vc but will just be necessary in
      // some cases ...
      points.clear();
      for (plain_facet_set::const_iterator fit = volcell_->facets().begin();
          fit != volcell_->facets().end(); ++fit)
      {
        for (std::size_t p = 0; p < (*fit)->points().size(); ++p)
        {
          bool insert = true;
          for (std::size_t ap = 0; ap < points.size(); ++ap)
            if (points[ap]->id() == (*fit)->points()[p]->id())
            {
              insert = false;
              break;
            }
          if (insert) points.push_back((*fit)->points()[p]);
        }
      }
    }
    else
    {
      tol *= 10;
      std::cout << "Warning: Increasing my tolerance to find a proper reference plane to tol = "
                << tol << std::endl;
    }
  }
  FOUR_C_THROW("Proper reference plane not found");

  return RefPlaneEqn;
}

/*-------------------------------------------------------------------------------------------------------*
 * Computation of reference plane based on the diagonal of background element sudhakar 06/15 This
 *considers all 6 diagonals of background Hex element, and choose the one that has maximum normal
 *component in x-direction
 *-------------------------------------------------------------------------------------------------------*/
bool Cut::DirectDivergenceGlobalRefplane::diagonal_based_ref(
    std::vector<double>& RefPlaneEqn, std::vector<Point*> points, double tol)
{
  if (options_.direct_divergence_refplane() != DirDiv_refplane_all &&
      options_.direct_divergence_refplane() != DirDiv_refplane_diagonal &&
      options_.direct_divergence_refplane() != DirDiv_refplane_diagonal_side)
    return false;

  // Any method should involve the following two steps

  //---
  // STEP 1: Estimate equation of reference plane and store the corresponding points for gmsh
  // output. Here we take all possible 6 diagonals of the Hex element, and choose the diagonal which
  // has maximum normal component in x-direction as the reference plane
  //---
  std::vector<Point*> ptslist = elem1_->points();

  std::vector<Point*> diag;
  std::vector<std::vector<Point*>> diagonals;
  for (unsigned i = 0; i < 24; ++i)  // 24 is the number of considered diagonals in
                                     // tri_diags_[24][3]
  {
    diag.clear();
    for (unsigned plid = 0; plid < 3; ++plid) diag.push_back(ptslist[tri_diags_[i][plid]]);
    diagonals.push_back(diag);
  }

  double xnormal = 0.0;
  bool found_refplane = false;
  for (std::vector<std::vector<Point*>>::iterator itd = diagonals.begin(); itd != diagonals.end();
      itd++)
  {
    std::vector<Point*> ptl = *itd;
    std::vector<double> RefPlaneTemp = Kernel::eqn_plane_of_polygon(ptl);

    scale_equation_of_plane(RefPlaneTemp);

    if (fabs(RefPlaneTemp[0]) < REF_PLANE_DIRDIV) continue;

    //---
    // STEP 2: Project all the corner points of the Hex element onto the reference plane.
    // If all these projected points are within the background element, then we consider this as a
    // possible reference plane!
    if (is_all_projected_corners_inside_ele(RefPlaneTemp, points, tol))
    {
      double fac =
          sqrt(pow(RefPlaneTemp[0], 2) + pow(RefPlaneTemp[1], 2) + pow(RefPlaneTemp[2], 2));
      //---
      // STEP 3: Take the reference plane with biggest component in x-direction of the normal
      // vector!
      double xn = fabs(RefPlaneTemp[0]) / fac;

      if (xn > xnormal)
      {
        xnormal = xn;
        RefPlaneEqn = RefPlaneTemp;
        ref_pts_gmsh_ = ptl;
        found_refplane = true;
      }
    }
  }

  return found_refplane;
}

/*-------------------------------------------------------------------------------------------------------*
 * Computation of reference plane based on the facets of the volumecell                         ager
 *02/16 Sort all the sides based on n_x (the one has more n_x gets on the top) Iterate through all
 *the sides to get the correct reference plane
 *-------------------------------------------------------------------------------------------------------*/
bool Cut::DirectDivergenceGlobalRefplane::facet_based_ref(
    std::vector<double>& RefPlaneEqn, std::vector<Point*> points, double tol)
{
  if (options_.direct_divergence_refplane() != DirDiv_refplane_all &&
      options_.direct_divergence_refplane() != DirDiv_refplane_facet)
    return false;

  const plain_facet_set& allfacets = volcell_->facets();

  //---
  // STEP 1: Estimate equation of reference plane and store the corresponding points for gmsh
  // output. First get all the facets of the volumecell and compute the equation of plane for each
  // facet Store them in a data structure which stores the sides based on n_x in descending order
  //---
  std::multimap<double, std::pair<std::vector<double>, std::vector<Point*>>, CompareClass>
      facet_data;
  for (plain_facet_set::const_iterator it = allfacets.begin(); it != allfacets.end(); it++)
  {
    std::vector<double> RefPlaneTemp = Kernel::eqn_plane_of_polygon((*it)->points());
    scale_equation_of_plane(RefPlaneTemp);
    if (fabs(RefPlaneTemp[0]) < REF_PLANE_DIRDIV) continue;

    facet_data.insert(
        std::make_pair(RefPlaneTemp[0], std::make_pair(RefPlaneTemp, (*it)->points())));
  }

  double xnormal = 0.0;
  bool found_refplane = false;
  for (auto it = facet_data.begin(); it != facet_data.end(); it++)
  {
    std::vector<double> RefPlaneTemp = it->second.first;

    scale_equation_of_plane(RefPlaneTemp);

    if (fabs(RefPlaneTemp[0]) < REF_PLANE_DIRDIV) continue;

    //---
    // STEP 2: Project all the corner points of the Hex element onto the reference plane.
    // If all these projected points are within the background element, then we consider this as a
    // possible reference plane!

    if (is_all_projected_corners_inside_ele(RefPlaneTemp, points, tol))
    {
      double fac =
          sqrt(pow(RefPlaneTemp[0], 2) + pow(RefPlaneTemp[1], 2) + pow(RefPlaneTemp[2], 2));
      //---
      // STEP 3: Take the reference plane with biggest component in x-direction of the normal
      // vector!
      double xn = fabs(RefPlaneTemp[0]) / fac;
      if (xn > xnormal)
      {
        xnormal = xn;
        RefPlaneEqn = RefPlaneTemp;
        ref_pts_gmsh_ = it->second.second;
        found_refplane = true;
      }
    }
  }

  //---
  // Basically using a diagonal reference plane should be enough, otherwise we have to look into
  // that again!
  return found_refplane;
}

/*-------------------------------------------------------------------------------------------------------*
 * Computation of reference plane based on the sides of background element sudhakar 06/15 Sort all
 *the sides based on n_x (the one has more n_x gets on the top) Iterate through all the sides to get
 *the correct reference plane
 *-------------------------------------------------------------------------------------------------------*/
bool Cut::DirectDivergenceGlobalRefplane::side_based_ref(
    std::vector<double>& RefPlaneEqn, std::vector<Point*> points, double tol)
{
  if (options_.direct_divergence_refplane() != DirDiv_refplane_all &&
      options_.direct_divergence_refplane() != DirDiv_refplane_side &&
      options_.direct_divergence_refplane() != DirDiv_refplane_diagonal_side)
    return false;

  const std::vector<Side*>& allsides = elem1_->sides();

  //---
  // STEP 1: Estimate equation of reference plane and store the corresponding points for gmsh
  // output. First get all the sides of the element and compute the equation of plane for each side
  // Store them in a data structure which stores the sides based on n_x in descending order
  //---
  std::multimap<double, std::pair<std::vector<double>, std::vector<Point*>>, CompareClass>
      side_data;
  for (std::vector<Side*>::const_iterator it = allsides.begin(); it != allsides.end(); it++)
  {
    const Side* s = *it;
    const std::vector<Node*> nds = s->nodes();


    for (std::size_t split_quadidx = 0; split_quadidx < 3 * nds.size() - 8; ++split_quadidx)
    {
      std::vector<Point*> ptside;
      if (nds.size() == 4)
      {
        ptside.push_back(nds[side_split_[split_quadidx][0]]->point());
        ptside.push_back(nds[side_split_[split_quadidx][1]]->point());
        ptside.push_back(nds[side_split_[split_quadidx][2]]->point());
      }
      else if (nds.size() == 3)
      {
        for (std::vector<Node*>::const_iterator itn = nds.begin(); itn != nds.end(); itn++)
          ptside.push_back((*itn)->point());
      }
      else
        FOUR_C_THROW("Side with another number of nodes than 3 or 4?");

      std::vector<double> RefPlaneTemp = Kernel::eqn_plane_of_polygon(ptside);
      scale_equation_of_plane(RefPlaneTemp);
      if (fabs(RefPlaneTemp[0]) < REF_PLANE_DIRDIV) continue;

      side_data.insert(std::make_pair(RefPlaneTemp[0], std::make_pair(RefPlaneTemp, ptside)));
    }
  }

  double xnormal = 0.0;
  bool found_refplane = false;
  for (std::multimap<double, std::pair<std::vector<double>, std::vector<Point*>>>::iterator it =
           side_data.begin();
      it != side_data.end(); it++)
  {
    std::vector<double> RefPlaneTemp = it->second.first;

    scale_equation_of_plane(RefPlaneTemp);

    if (fabs(RefPlaneTemp[0]) < REF_PLANE_DIRDIV) continue;

    //---
    // STEP 2: Project all the corner points of the Hex element onto the reference plane.
    // If all these projected points are within the background element, then we consider this as a
    // possible reference plane!
    if (is_all_projected_corners_inside_ele(RefPlaneTemp, points, tol))
    {
      double fac =
          sqrt(pow(RefPlaneTemp[0], 2) + pow(RefPlaneTemp[1], 2) + pow(RefPlaneTemp[2], 2));
      //---
      // STEP 3: Take the reference plane with biggest component in x-direction of the normal
      // vector!
      double xn = fabs(RefPlaneTemp[0]) / fac;
      if (xn > xnormal)
      {
        xnormal = xn;
        RefPlaneEqn = RefPlaneTemp;
        ref_pts_gmsh_ = it->second.second;
        found_refplane = true;
      }
    }
  }

  //---
  // Basically using a diagonal reference plane should be enough, otherwise we have to look into
  // that again!
  return found_refplane;
}

/*---------------------------------------------------------------------------------------------------------------*
 * In order to check whether the chosen reference plane is a correct choice,
 * we project all the corner points of the element onto this reference plane sudhakar 06/15 If all
 *these projected points are within the element, then we got the correct ref plane
 *---------------------------------------------------------------------------------------------------------------*/
bool Cut::DirectDivergenceGlobalRefplane::is_all_projected_corners_inside_ele(
    std::vector<double>& RefPlaneEqn, std::vector<Point*> points, double tol)
{
  for (std::vector<Point*>::iterator it = points.begin(); it != points.end(); it++)
  {
    Point* pt = *it;

    Core::LinAlg::Matrix<3, 1> coo;
    pt->coordinates(coo.data());

    Core::LinAlg::Matrix<3, 1> xyz_proj(coo), rst_proj;

    // x-coordinate of corner pt projected in the reference plane
    // y- and z-coordinate remain the same
    xyz_proj(0, 0) =
        (RefPlaneEqn[3] - RefPlaneEqn[1] * coo(1, 0) - RefPlaneEqn[2] * coo(2, 0)) / RefPlaneEqn[0];

    // get the local coordinates of the projected point
    elem1_->local_coordinates(xyz_proj, rst_proj);

    // Check whether the local coordinate of the projected point is within the specified limits
    if (std::abs(rst_proj(0, 0)) > 1.0 + tol or std::abs(rst_proj(1, 0)) > 1.0 + tol or
        std::abs(rst_proj(2, 0)) > 1.0 + tol or std::isnan(rst_proj(0, 0)) or
        std::isnan(rst_proj(1, 0)) or std::isnan(rst_proj(2, 0)))
    {
      return false;
    }
  }
  return true;
}

/*----------------------------------------------------------------------------------------------------------*
 * Scale the equation of plane to enable comparison of normals sudhakar 07/15
 *----------------------------------------------------------------------------------------------------------*/
void Cut::DirectDivergenceGlobalRefplane::scale_equation_of_plane(std::vector<double>& RefPlaneEqn)
{
  double scale =
      sqrt(pow(RefPlaneEqn[0], 2.0) + pow(RefPlaneEqn[1], 2.0) + pow(RefPlaneEqn[2], 2.0));
  for (unsigned i = 0; i < 4; i++) RefPlaneEqn[i] /= scale;
}

const unsigned Cut::DirectDivergenceGlobalRefplane::tri_diags_[24][3] = {{1, 6, 7}, {0, 6, 7},
    {0, 1, 7}, {0, 1, 6}, {3, 4, 5}, {2, 4, 5}, {2, 3, 5}, {2, 3, 4}, {6, 3, 0}, {5, 3, 0},
    {5, 6, 0}, {5, 6, 3}, {7, 2, 1}, {4, 2, 1}, {4, 7, 1}, {4, 7, 2}, {4, 6, 2}, {0, 6, 2},
    {0, 4, 2}, {0, 4, 6}, {1, 3, 7}, {5, 3, 7}, {5, 1, 7}, {5, 1, 3}};

const unsigned Cut::DirectDivergenceGlobalRefplane::side_split_[4][3] = {
    {1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}};

FOUR_C_NAMESPACE_CLOSE
