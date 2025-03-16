// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_cut_pointpool.hpp"

#include "4C_cut_output.hpp"
#include "4C_global_data.hpp"

FOUR_C_NAMESPACE_OPEN


/*-----------------------------------------------------------------------------------------*
 * If a point with the coordinates "x" does not exists, it creates a new point correspondingly
 *-----------------------------------------------------------------------------------------*/
Cut::Point* Cut::OctTreeNode::new_point(const double* x, Edge* cut_edge, Side* cut_side,
    double tolerance, PointpoolMergeStrategy merge_strategy)
{
  // check if the point already exists
#if CUT_CREATION_INFO
  bool new_point = false;
#endif
  Core::LinAlg::Matrix<3, 1> px(x);

  Point* p = get_point(x, cut_edge, cut_side, tolerance, merge_strategy);

  if (p == nullptr)
  {
    p = &*create_point(points_.size(), x, cut_edge, cut_side, tolerance);  // create the point
#if CUT_CREATION_INFO
    new_point = true;
#endif

    if (points_.size() % 1000 == 0)  // split the node starting from level 0
    {
      split(0);
    }
  }


#if CUT_CREATION_INFO
  // if it was merged
  if (not new_point)
  {
    std::stringstream info;
    info << "// Another point was merged in this one with\n";
    info << "// Initial coordinates" << std::setprecision(15) << px << std::endl;
    // merged_to point coordinates
    Core::LinAlg::Matrix<3, 1> nx;
    p->coordinates(nx.data());
    px.update(-1, nx, 1);
    info << "// Merged with tolerance of " << std::setprecision(15) << px.norm2() << std::endl;
    p->AddAdditionalCreationInfo(info.str());
  }

  else
  {
    p->AddAdditionalCreationInfo(
        "// This point was really created with the exact coordinates given");
  }
#endif

  return p;
}

/*-----------------------------------------------------------------------------------------*
 * Get the point with the specified coordinates "x" from the pointpool
 *-----------------------------------------------------------------------------------------*/
Cut::Point* Cut::OctTreeNode::get_point(const double* x, Edge* cut_edge, Side* cut_side,
    double tolerance, PointpoolMergeStrategy merge_strategy)
{
  // try to find the point in all of the 8 children nodes
  if (not is_leaf())
  {
    // stop finding the point when point is not included in the current bounding box
    if (!bb_->within(1.0, x)) return nullptr;

    for (int i = 0; i < 8; ++i)
    {
      Point* p = nodes_[i]->get_point(x, cut_edge, cut_side, tolerance, merge_strategy);
      if (p != nullptr)
      {
        return p;
      }
    }
  }
  else
  {
    Core::LinAlg::Matrix<3, 1> px(x);
    Core::LinAlg::Matrix<3, 1> nx;

    double tol = TOPOLOGICAL_TOLERANCE * norm_;


    switch (merge_strategy)
    {
      case PointpoolMergeStrategy::SelfCutLoad:
      {
        tol = TOPOLOGICAL_TOLERANCE * norm_ * NODAL_POINT_TOLERANCE_SELFCUT_SCALE;
        break;
      }
      case PointpoolMergeStrategy::InitialLoad:
      {
        // when we are loading the geometry into the cut  we want the distance to be large than
        // topological tolerance. Otherwise we might experience huge problems
        tol = TOPOLOGICAL_TOLERANCE * norm_ * NODAL_POINT_TOLERANCE_SCALE;
        // safety check
#ifdef NODAL_POINT_TOLERANCE_SCALE
        // this should happen with both cut_side and cut_edge equal to nullptr
        // note: be careful with other cases, apart from mesh loading
        if (cut_side and cut_edge)
        {
          FOUR_C_THROW(
              "Scaling is {} for non-nullptr cut_side and cut_edge. This should not be possible!",
              NODAL_POINT_TOLERANCE_SCALE);
        }
#endif
        break;
      }
      case PointpoolMergeStrategy::NormalCutLoad:
      {
        // no additional scale
        break;
      }
      default:
        FOUR_C_THROW("Unknown merge strategy is equal to {}", merge_strategy);
    }

    // linear search for the node in the current leaf
    std::vector<Point*> merge_candidates;  // candidates for merge

    for (RCPPointSet::iterator i = points_.begin(); i != points_.end(); ++i)
    {
      Point* p = &**i;

      p->coordinates(nx.data());
      nx.update(-1, px, 1);
      if (nx.norm2() <= tol)
      {
        merge_candidates.push_back(p);
      }
    }


    // if there are merge candidates
    if (merge_candidates.size() > 0)
    {
      // first consider topologically connected points
      // and try to find any candidate there
      std::vector<Point*> topological_candidates;
      Point* mymerge = nullptr;
      int match_number = 0;
      // base tolerance, all the candidates should fit into that one
      for (unsigned int i = 0; i < merge_candidates.size(); ++i)
      {
        // if this point was cut but same side and edge.
        // can happen if this edge is shared between multiple sides, we want the intersection point
        // to be consistent and not to merged somewhere else, as then it will cause problem
        if (merge_candidates[i]->is_cut(cut_side, cut_edge))
        {
          topological_candidates.push_back(merge_candidates[i]);
          mymerge = merge_candidates[i];
          match_number++;
        }
      }

      if (match_number >= 2)
      {
        std::stringstream err_msg;
        err_msg << "More then one merge candidate in the pointpool, that both are intersected by "
                   "the same side"
                   "and edge and fit into tolerance. This should not happen. Some points are "
                   "probably created twice or not merged";
        // do gmsh side and edge dump for analysis
        std::ofstream file("pointpool_conflict_topologically_connected.pos");
        for (std::vector<Point*>::iterator it = topological_candidates.begin();
            it != topological_candidates.end(); ++it)
        {
#if CUT_CREATION_INFO
          std::cout << "Id" << (*it)->Id()
                    << ((*it)->IsReallyCut(std::pair<Side*, Edge*>(cut_side, cut_edge))
                               ? " is really cut"
                               : "connection is added later")
                    << std::endl;
#endif
          (*it)->dump_connectivity_info();
        }

        if (cut_side)
        {
          Cut::Output::gmsh_new_section(file, "CurrentIntersectionSide");
          Cut::Output::gmsh_side_dump(file, cut_side, false, nullptr);
          Cut::Output::gmsh_end_section(file, false);
        }
        if (cut_edge)
        {
          Cut::Output::gmsh_new_section(file, "CurrentIntersectionEdge");
          Cut::Output::gmsh_edge_dump(file, cut_edge, false, nullptr);
          Cut::Output::gmsh_end_section(file, false);
        }

        file.close();
        FOUR_C_THROW("{}", err_msg.str());
      }

      // if there are no points intersected by same side and edge but there  are still other points
      // that fit into the tolerance
      if (mymerge == nullptr)
      {
        // there only unknown non-topologically connected case
        if (merge_candidates.empty()) return nullptr;

        // if nothing was intersected before return first matching point. else do nothing
        mymerge = merge_candidates[0];

#if EXTENDED_CUT_DEBUG_OUTPUT
        // If there are more matching candidates we basically don't really know where to merge.
        // For now we just merge into the closest one
        // Otherwise One possible idea is is would be to
        // merge into the point that is "the most close". To do this we just need to find minimum
        // merge_candidates[i] ->coordinates( nx.data() )
        // min (nx  - px).norm2()  and merge into that one.
        if (merge_candidates.size() > 1)
        {
          std::cout << "NOTE: More then one merge candidate in the pointpool that fit into the "
                       "tolerance. There are possibly"
                    << merge_candidates.size() << "merging candidates. Merging it into "
                    << mymerge->Id() << std::endl;

          // do gmsh side and edge dump for analysis
          std::ofstream file("pointpool_conflict_side_no_shared_cut.pos", std::ios_base::app);

          if (cut_side)
          {
            Cut::Output::gmsh_new_section(file, "InterSides");
            Cut::Output::gmsh_side_dump(file, cut_side, false, nullptr);
            Cut::Output::gmsh_end_section(file, false);
          }
          if (cut_edge)
          {
            Cut::Output::gmsh_new_section(file, "InterEdges");
            Cut::Output::GmshEdgeDump(file, cut_edge, false, nullptr);
            Cut::Output::gmsh_end_section(file, false);
          }

          file.close();
          std::cout << "Actual point coordinates, (the real one, before merging) ar "
                    << std::setprecision(15) << px << std::endl;
        }
#endif
      }

      // after this we found proper megin candidates
      if (mymerge != nullptr)
      {  // just a safety check

        mymerge->coordinates(nx.data());
        nx.update(-1, px, 1);
        if (not mymerge->is_cut(cut_side, cut_edge))
        {
          if (cut_edge != nullptr)
          {
            // Here we need to set this point local coordinates on the edge, based on the unmerged
            // point, that it why "t" is called explicitly before AddEdge
            mymerge->t(cut_edge, px);
            mymerge->add_edge(cut_edge);
          }
          if (cut_side != nullptr)
          {
            mymerge->add_side(cut_side);
          }
          if ((cut_side != nullptr) && (cut_edge != nullptr))
          {
            mymerge->add_pair(cut_side, cut_edge);
#if CUT_CREATION_INFO
            std::pair<Side*, Edge*> int_pair = std::make_pair(cut_side, cut_edge);
            mymerge->AddMergedPair(int_pair, &px);
#endif
          }
        }

        // This is done because we want to merge cut_mesh into the normal mesh first
        if (merge_strategy == PointpoolMergeStrategy::InitialLoad)
        {
          mymerge->move_point(x);
        }
        return mymerge;
      }
    }
  }
  return nullptr;
}


/*-----------------------------------------------------------------------------------------*
 * Get the point with the specified coordinates "x" from the pointpool
 *-----------------------------------------------------------------------------------------*/
std::shared_ptr<Cut::Point> Cut::OctTreeNode::create_point(
    unsigned newid, const double* x, Edge* cut_edge, Side* cut_side, double tolerance)
{
  if (not is_leaf())
  {
    // call recursively create_point for the child where the Point shall lie in
    std::shared_ptr<Point> p = leaf(x)->create_point(newid, x, cut_edge, cut_side, tolerance);
    // add the pointer not only in the leaf but also on the current level
    add_point(x, p);
    return p;
  }
  else
  {
    // create a new point and add the point at the lowest level
    std::shared_ptr<Point> p = Cut::create_point(newid, x, cut_edge, cut_side, tolerance);
    add_point(x, p);
    return p;
  }
}


/*-----------------------------------------------------------------------------------------*
 * Simply insert p into the pointpool and correspondingly modify the boundingbox size
 *-----------------------------------------------------------------------------------------*/
void Cut::OctTreeNode::add_point(const double* x, std::shared_ptr<Point> p)
{
  points_.insert(p);  // insert the point in the pointpool
  bb_->add_point(x);  // modify the boundingbox size
}


/*-----------------------------------------------------------------------------------------*
 * get the leaf where the point with the given coordinates lies in
 *-----------------------------------------------------------------------------------------*/
Cut::OctTreeNode* Cut::OctTreeNode::leaf(const double* x)
{
  // navigate to the right one of the 8 children nodes
  //
  //    z <0            1  |  3         z > 0        5  |  7
  //        ____ y     ____|____                    ____|____
  //       |               |                            |
  //       |            2  |  4                      6  |  8
  //       x
  //

  int idx = 0;
  if (x[0] > splitpoint_(0))  // add an index of one to move in x direction
  {
    idx += 1;
  }
  if (x[1] > splitpoint_(1))  // add an index of two to move in y direction
  {
    idx += 2;
  }
  if (x[2] > splitpoint_(2))  // add an index of four to move in z direction
  {
    idx += 4;
  }
  return &*nodes_[idx];
}


/*-----------------------------------------------------------------------------------------*
 * split the current boounding box (tree-node)
 *-----------------------------------------------------------------------------------------*/
void Cut::OctTreeNode::split(int level)
{
  // We must not end up with a OctTreeNode that holds just nodes from the
  // cutter mesh. However, there is no real way to test this right now.

  if (points_.size() > 125)  /// 125 = 1/8 *1000 -> see NewPoint
  {
    Core::LinAlg::Matrix<3, 1> x;
    bool first = true;

    for (RCPPointSet::iterator i = points_.begin(); i != points_.end(); ++i)
    {
      Point* p = &**i;
      if (first)
      {
        first = false;
        p->coordinates(splitpoint_.data());
      }
      else
      {
        p->coordinates(x.data());
        splitpoint_.update(1, x, 1);
      }
    }

    splitpoint_.scale(1. / points_.size());

    for (int i = 0; i < 8; ++i)
    {
      nodes_[i] = std::make_shared<OctTreeNode>(norm_);
    }

    // avoid empty room (room not covered by boundary boxes)
    for (int i = 0; i < 8; ++i)
    {
      // always have the split point in all boxes
      nodes_[i]->bb_->add_point(splitpoint_);

      // always have the outmost point in each box
      double x[3];
      bb_->corner_point(i, x);
      leaf(x)->bb_->add_point(x);
    }

    for (RCPPointSet::iterator i = points_.begin(); i != points_.end(); ++i)
    {
      std::shared_ptr<Point> p = *i;
      double x[3];
      p->coordinates(x);
      leaf(x)->add_point(x, p);
    }

    for (int i = 0; i < 8; ++i)
    {
      nodes_[i]->split(level + 1);
    }
  }
}


/*-----------------------------------------------------------------------------------------*
 * collect all edges
 *-----------------------------------------------------------------------------------------*/
void Cut::OctTreeNode::collect_edges(const BoundingBox& edgebox, plain_edge_set& edges)
{
  if (not is_leaf())
  {
    if (edgebox.within(norm_, *bb_))
    {
      for (int i = 0; i < 8; ++i)
      {
        nodes_[i]->collect_edges(edgebox, edges);
      }
    }
  }
  else
  {
    std::shared_ptr<BoundingBox> sbox(BoundingBox::create());
    for (RCPPointSet::iterator i = points_.begin(); i != points_.end(); ++i)
    {
      Point* p = &**i;
      const plain_edge_set& sds = p->cut_edges();
      for (plain_edge_set::const_iterator i = sds.begin(); i != sds.end(); ++i)
      {
        Edge* s = *i;
        if (edges.count(s) == 0)
        {
          sbox->assign(*s);
          if (sbox->within(norm_, edgebox))
          {
            edges.insert(s);
          }
        }
      }
    }
  }
}


/*-----------------------------------------------------------------------------------------*
 * collect all sides
 *-----------------------------------------------------------------------------------------*/
void Cut::OctTreeNode::collect_sides(const BoundingBox& sidebox, plain_side_set& sides)
{
  if (not is_leaf())
  {
    if (sidebox.within(norm_, *bb_))
    {
      for (int i = 0; i < 8; ++i)
      {
        nodes_[i]->collect_sides(sidebox, sides);
      }
    }
  }
  else
  {
    std::shared_ptr<BoundingBox> sbox(BoundingBox::create());
    for (RCPPointSet::iterator i = points_.begin(); i != points_.end(); ++i)
    {
      Point* p = &**i;
      const plain_side_set& sds = p->cut_sides();
      for (plain_side_set::const_iterator i = sds.begin(); i != sds.end(); ++i)
      {
        Side* s = *i;
        if (sides.count(s) == 0)
        {
          sbox->assign(*s);
          if (sbox->within(norm_, sidebox))
          {
            sides.insert(s);
          }
        }
      }
    }
  }
}


/*-----------------------------------------------------------------------------------------*
 * collect all elements near the sidebox
 * (unused, does not work properly when there is no point adjacent to elements in a tree's leaf,
 * e.g. when side lies within an element)
 *-----------------------------------------------------------------------------------------*/
void Cut::OctTreeNode::collect_elements(const BoundingBox& sidebox, plain_element_set& elements)
{
  // see REMARK in cut_mesh.cpp
  FOUR_C_THROW(
      "collecting elements via the OctTreeNode does not find all possible element-side "
      "intersections");

  if (not is_leaf())
  {
    if (sidebox.within(
            norm_, *bb_))  // within check is a check of overlap between the 2 bounding boxes
    {
      for (int i = 0; i < 8; ++i)
      {
        nodes_[i]->collect_elements(sidebox, elements);
      }
    }
  }
  else
  {
    std::shared_ptr<BoundingBox> elementbox(BoundingBox::create());
    for (RCPPointSet::iterator i = points_.begin(); i != points_.end(); ++i)
    {
      Point* p = &**i;
      const plain_element_set& els = p->elements();

      // add all elements adjacent to the current point
      // REMARK: this does not find all elements that have an overlap with the sidebox!!!
      for (plain_element_set::const_iterator i = els.begin(); i != els.end(); ++i)
      {
        Element* e = *i;
        if (elements.count(e) == 0)
        {
          elementbox->assign(*e);
          if (elementbox->within(norm_, sidebox))
          {
            elements.insert(e);
          }
        }
      }
    }
  }
}


/*-----------------------------------------------------------------------------------------*
 * reset the Point::Position of outside points
 *-----------------------------------------------------------------------------------------*/
void Cut::OctTreeNode::reset_outside_points()
{
  for (RCPPointSet::iterator i = points_.begin(); i != points_.end(); ++i)
  {
    Point* p = &**i;
    if (p->position() == Point::outside)
    {
      p->position(Point::undecided);
    }
  }
}


/*-----------------------------------------------------------------------------------------*
 * print the tree at a given level
 *-----------------------------------------------------------------------------------------*/
void Cut::OctTreeNode::print(int level, std::ostream& stream)
{
  if (not is_leaf())
  {
    for (int i = 0; i < 8; ++i)
    {
      nodes_[i]->print(level + 1, stream);
    }
  }
  else
  {
    for (RCPPointSet::iterator i = points_.begin(); i != points_.end(); ++i)
    {
      Point* p = &**i;
      p->plot(stream);
    }
    stream << "\n";
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Cut::PointPool::PointPool(double norm) : tree_(norm), probdim_(Global::Problem::instance()->n_dim())
{
}


// find node of the octree where current point resides, search start from "this" node and goes
// downwards
Cut::OctTreeNode* Cut::OctTreeNode::find_node(const double* coord, Point* p)
{
  if (p)
  {
    if (not is_leaf())
    {
      // stop finding the point when point is not included in the current bounding box
      if (!bb_->within(1.0, coord)) return nullptr;
      for (int i = 0; i < 8; ++i)
      {
        Cut::OctTreeNode* node = (nodes_[i]->find_node(coord, p));
        if (node != nullptr)
        {
          return node;
        }
      }
    }
    // if it is leaf
    else
    {
      RCPPointSet::iterator i = points_.begin();
      for (; i != points_.end(); ++i)
      {
        Point* b = &**i;
        if (b == p)
        {
          // return current node
          return this;
        }
      }
      if (i == points_.end())
      {
        return nullptr;
      }
    }
  }
  else
  {
    FOUR_C_THROW("Invalid point");
    return nullptr;
  }
  return nullptr;
}

FOUR_C_NAMESPACE_CLOSE
