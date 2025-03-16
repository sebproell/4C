// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_cut_tetmeshintersection.hpp"

#include "4C_cut_boundarycell.hpp"
#include "4C_cut_integrationcell.hpp"
#include "4C_cut_node.hpp"
#include "4C_cut_options.hpp"
#include "4C_cut_output.hpp"
#include "4C_cut_pointpool.hpp"
#include "4C_cut_side.hpp"
#include "4C_cut_volumecell.hpp"

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyTraits.hpp>

FOUR_C_NAMESPACE_OPEN

/* Initialize a Mesh within an element. This is done if a TetMesh is created but can't be filled in
   fill_facet_mesh(). Basically this constructor converts the data structure used for TetMesh into
   the data structures used for Mesh.

   The TetMesh produced from the in the "parent element" is decomposed into a "small" Mesh object
   with its TETs. When Cut is Called these TETs are cut as a normal Mesh would with cut-sides etc...
   Then every TET of this Mesh, is is tested to see if a TetMesh can be produced out of its new cut
   configuration. If Yes -> We are done! The cut element has been tesselated If No  -> Another call
   to the TetMeshIntersection class is necessary PROBLEM: For problematic cut situations we run into
   machine precision error with this type of algorithm. It is NOT very stable... An alternative
   should be pursued. It is likely that ALE will not work for LevelSet with this algorithm,
   eventhough it works (sort of) all right for a cartesian constellation.

   For a LevelSet-Cut the triangulated facet is decomposed into TRI-surfaces here.

   o mesh_      is a normal Mesh object composed of the the accepted TETs
   o cut_mesh_  is created from the facets and triangulated facets????
   o pp_        a new local point pool for the decomposition of this tetmesh of this element into a
   mesh object.

 */
Cut::TetMeshIntersection::TetMeshIntersection(Options& options, Element* element,
    const std::vector<std::vector<int>>& tets, const std::vector<int>& accept_tets,
    const std::vector<Point*>& points, const plain_side_set& cut_sides)
    : pp_(std::make_shared<PointPool>()), mesh_(options, 1, pp_), cut_mesh_(options, 1, pp_, true)
{
  // Create the nodes and make the connectivity to the parent_mesh.
  for (std::vector<Point*>::const_iterator i = points.begin(); i != points.end(); ++i)
  {
    Point* p = *i;
    Node* n = mesh_.get_node(std::distance(points.begin(), i), p->x());
    Point* np = n->point();
    np->position(p->position());
    do_register(p, np);
  }

  // Create the tets and register if the nodes of the tets are on a cut surface (i.e. register
  // edges)
  for (std::vector<std::vector<int>>::const_iterator i = tets.begin(); i != tets.end(); ++i)
  {
    const std::vector<int>& tet = *i;
    unsigned id = std::distance(tets.begin(), i);
    if (accept_tets[id])
    {
      Element* e = mesh_.get_element(
          id, tet, *shards::getCellTopologyData<shards::Tetrahedron<4>>(), accept_tets[id]);
      const std::vector<Node*>& nodes = e->nodes();
      for (std::vector<Node*>::const_iterator i = nodes.begin(); i != nodes.end(); ++i)
      {
        Node* n = *i;
        n->register_cuts();
      }
    }
  }

  const plain_facet_set& element_facets = element->facets();

  // Triangulated cut facets need to be converted to tri cut sides. This is
  // done after the normal cut sides.

  std::vector<Facet*> triangulated;
  std::map<Point*, Node*> nodemap;

  for (plain_side_set::const_iterator i = cut_sides.begin(); i != cut_sides.end(); ++i)
  {
    Side* s = *i;

    plain_facet_set facets;
    const std::vector<Facet*>& side_facets = s->facets();
    // Extract the facets of the cut-side
    for (std::vector<Facet*>::const_iterator i = side_facets.begin(); i != side_facets.end(); ++i)
    {
      Facet* f = *i;
      if (element_facets.count(f) > 0)
      {
        facets.insert(f);
      }
    }

    for (plain_facet_set::iterator i = facets.begin(); i != facets.end(); ++i)
    {
      Facet* f = *i;

      // Is this facet a LevelSetSide or is it Triangulated
      // add this to triangulated and the new nodes to nodemap
      if (f->belongs_to_level_set_side() or f->is_triangulated())
      {
        triangulated.push_back(f);
        PointSet points;
        f->all_points(points);
        for (PointSet::iterator i = points.begin(); i != points.end(); ++i)
        {
          Point* p = *i;
          nodemap[to_child(p)];
        }
      }
      else
      {
        copy_cut_side(s, f);
      }
    }
  }

  // create nodes

  cut_mesh_.new_nodes_from_points(nodemap);

  // do triangulated facets (create extra cut sides)
  // and initialize cut_mesh_

  for (std::vector<Facet*>::iterator i = triangulated.begin(); i != triangulated.end(); ++i)
  {
    Facet* f = *i;
    Side* s = f->parent_side();

    if (f->is_triangulated())
    {
      const std::vector<std::vector<Point*>>& triangulation = f->triangulation();
      for (std::vector<std::vector<Point*>>::const_iterator i = triangulation.begin();
          i != triangulation.end(); ++i)
      {
        const std::vector<Point*>& tri = *i;
        if (tri.size() != 3) FOUR_C_THROW("tri3 expected");
        std::vector<Node*> nodes;
        nodes.reserve(3);
        for (std::vector<Point*>::const_iterator i = tri.begin(); i != tri.end(); ++i)
        {
          Point* p = *i;
          nodes.push_back(nodemap[to_child(p)]);
        }
        Side* cs =
            cut_mesh_.get_side(s->id(), nodes, shards::getCellTopologyData<shards::Triangle<3>>());
        side_parent_to_child_[s].push_back(cs);
        for (std::vector<Node*>::iterator i = nodes.begin(); i != nodes.end(); ++i)
        {
          Node* n = *i;
          n->register_cuts();
        }
      }
    }
    else
    {
      const std::vector<Point*>& points = f->corner_points();

      std::vector<Node*> nodes;
      nodes.reserve(points.size());
      for (std::vector<Point*>::const_iterator i = points.begin(); i != points.end(); ++i)
      {
        Point* p = *i;
        nodes.push_back(nodemap[to_child(p)]);
      }

      switch (points.size())
      {
        case 2:
          // Degenerated nonsense. Why does that happen?
          FOUR_C_THROW(
              "In TetMeshIntersection() the corner points of a facet is points.size()=2. Should't "
              "happen?");

          // make sure the mapping entry exists
          side_parent_to_child_[s];
          break;
        case 3:
        {
          Side* cs = cut_mesh_.get_side(
              s->id(), nodes, shards::getCellTopologyData<shards::Triangle<3>>());
          side_parent_to_child_[s].push_back(cs);
          for (std::vector<Node*>::iterator i = nodes.begin(); i != nodes.end(); ++i)
          {
            Node* n = *i;
            n->register_cuts();
          }
          break;
        }
        default:
          FOUR_C_THROW("facet with more that three points");
      }
    }
  }
}

void Cut::TetMeshIntersection::find_edge_cuts()
{
  plain_edge_set cut_edges;
  const std::map<plain_int_set, std::shared_ptr<Edge>>& c_edges = mesh_.edges();
  for (std::map<plain_int_set, std::shared_ptr<Edge>>::const_iterator i = c_edges.begin();
      i != c_edges.end(); ++i)
  {
    Edge* e = &*i->second;
    cut_edges.insert(e);
  }

  for (plain_edge_set::iterator i = cut_edges.begin(); i != cut_edges.end(); ++i)
  {
    Edge* ce = *i;
    std::shared_ptr<BoundingBox> edgebox(BoundingBox::create(*ce));
    plain_edge_set edges;
    pp_->collect_edges(*edgebox, edges);
    // edges.erase( ce );

    PointSet cp;
    cp.insert(ce->begin_node()->point());
    cp.insert(ce->end_node()->point());

    for (plain_edge_set::iterator i = edges.begin(); i != edges.end(); ++i)
    {
      Edge* e = *i;
      if (cut_edges.count(e) == 0)
      {
        // Find cut points between edges. Some might be new.

        if (cp.count(e->begin_node()->point()) == 0 and cp.count(e->end_node()->point()) == 0)
        {
          double tolerance;
          try
          {
            e->compute_cut(&mesh_, ce, nullptr, nullptr, tolerance);
          }
          catch (Core::Exception& err)
          {
            std::cout << "\n-------------------\n";
            std::cout << "\nCut Edge\n";
            ce->print();
            std::cout << "\nOther Edge\n";
            e->print();
            std::cout << "\n-------------------\n";

            plain_edge_set edges;
            edges.reserve(2);
            edges.insert(e);
            edges.insert(ce);
            Cut::Output::gmsh_edges_only(edges);
            throw;
          }
        }
      }
    }
  }
}

/* Cut the created mesh in a similar fashion as is done for ParentIntersection.
 * However it is not done exactly the same way... It seems to be quite complex...
 * One would have to spend quite some time to find out the thinking behind it and
 * maybe connect it to the existing algo. */
void Cut::TetMeshIntersection::cut(Mesh& parent_mesh, Element* element,
    const plain_volumecell_set& parent_cells, int count, bool tetcellsonly)
{
  find_edge_cuts();

  plain_element_set elements_done;
  // increase counter
  ++count;
  // note: purpose of count is not really known, guess it was related to recursion
  cut_mesh_.cut(mesh_, elements_done /*, count */);
  // cut_mesh_.Cut( mesh_, elements_done, count );

  // New way of cut? Might need to join the the TetMeshIntersection with the new Cut algorithm?
  // THIS IS NOT WORKING AT THE MOMENT!!!
  //  mesh_.SearchCollisions(cut_mesh_);
  //  mesh_.find_cut_points(count+1);

  cut_mesh_.rectify_cut_numerics();
  mesh_.rectify_cut_numerics();

  mesh_.make_cut_lines();
  mesh_.make_facets();
  mesh_.make_volume_cells();

  std::map<VolumeCell*, ChildCell> cellmap;

  map_volume_cells(parent_mesh, element, parent_cells, cellmap);

  mesh_.create_integration_cells(count, tetcellsonly);

  fill(parent_mesh, element, parent_cells, cellmap);
}

/// Make sure connectivity between the new VolumeCell (of the children of the VC) are correctly
/// connected to its parents. However not sure exactly how this is done....
void Cut::TetMeshIntersection::map_volume_cells(Mesh& parent_mesh, Element* element,
    const plain_volumecell_set& parent_cells, std::map<VolumeCell*, ChildCell>& cellmap)
{
  plain_volumecell_set done_child_cells;

  seed_cells(parent_mesh, parent_cells, cellmap, done_child_cells);

  int nonnodecells = 0;

  for (plain_volumecell_set::const_iterator i = parent_cells.begin(); i != parent_cells.end(); ++i)
  {
    VolumeCell* vc = *i;
    ChildCell& cc = cellmap[vc];
    plain_volumecell_set& childset = cc.cells_;

    if (childset.size() > 0)
    {
      fill(vc, cc);
      std::copy(childset.begin(), childset.end(),
          std::inserter(done_child_cells, done_child_cells.begin()));
    }
    else
    {
      nonnodecells += 1;
    }
  }

  // emergency seed cell filling

  while (nonnodecells > 0)
  {
    int backup = nonnodecells;

    for (std::map<VolumeCell*, ChildCell>::iterator i = cellmap.begin(); i != cellmap.end(); ++i)
    {
      // VolumeCell * vc = i->first;
      ChildCell& cc = i->second;
      // plain_volumecell_set & childset = cc.cells_;
      std::map<Side*, std::vector<Facet*>>& facetsonsurface = cc.facetsonsurface_;

      // match parent and child volumes at cut surface

      for (std::map<Side*, std::vector<Facet*>>::iterator i = facetsonsurface.begin();
          i != facetsonsurface.end(); ++i)
      {
        Side* child_side = i->first;
        std::vector<Facet*>& parent_facets = i->second;

        // match only possible if the side has only one facet on this volume
        // there should be no other case?

        if (parent_facets.size() == 1)
        {
          Facet* facet = parent_facets[0];
          const plain_volumecell_set& parent_cells = facet->cells();
          std::vector<ChildCell*> parent_cell_info;
          parent_cell_info.reserve(2);
          for (plain_volumecell_set::const_iterator i = parent_cells.begin();
              i != parent_cells.end(); ++i)
          {
            VolumeCell* vc = *i;

            std::map<VolumeCell*, ChildCell>::iterator j = cellmap.find(vc);
            if (j != cellmap.end())
            {
              parent_cell_info.push_back(&j->second);
            }
          }

          // If there are less than two volumes at this cut, we have a touch
          // at a boundary.

          if (parent_cell_info.size() == 1)
          {
            if (not parent_cell_info[0]->done_)
            {
              const std::vector<Facet*>& child_facets = child_side->facets();
              for (std::vector<Facet*>::const_iterator i = child_facets.begin();
                  i != child_facets.end(); ++i)
              {
                Facet* f = *i;
                const plain_volumecell_set& child_cells = f->cells();
                if (child_cells.size() == 1)
                {
                  VolumeCell* c = *child_cells.begin();
                  if (done_child_cells.count(c) == 0)
                  {
                    if (not parent_cell_info[0]->contains_child(c))
                    {
                      parent_cell_info[0]->cells_.insert(c);
                      done_child_cells.insert(c);
                    }
                  }
                }
                else if (child_cells.size() == 2)
                {
                  // odd.
                  FOUR_C_THROW(
                      "illegal number of neighbouring volume cells: child_cells.size()==2");
                }
                else
                {
                  std::stringstream str;
                  str << "illegal number of neighbouring volume cells: child_cells.size()=="
                      << child_cells.size();
                  FOUR_C_THROW("{}", str.str());
                }
              }

              ChildCell& cc = *parent_cell_info[0];
              fill(cc.parent_, cc);
              std::copy(cc.cells_.begin(), cc.cells_.end(),
                  std::inserter(done_child_cells, done_child_cells.begin()));
              nonnodecells -= 1;
            }
          }
          else if (parent_cell_info.size() == 2)
          {
            // Only useful if one volume is already done and the other one is not.

            int doneindex = -1;
            int otherindex = -1;
            if (parent_cell_info[0]->done_ and not parent_cell_info[1]->done_)
            {
              doneindex = 0;
              otherindex = 1;
            }
            else if (parent_cell_info[1]->done_ and not parent_cell_info[0]->done_)
            {
              doneindex = 1;
              otherindex = 0;
            }

            if (doneindex > -1)
            {
              // There are many child facets on the child side. But those facets
              // match the one parent facet. There must be no other facet
              // outside that region.

              bool found = false;

              const std::vector<Facet*>& child_facets = child_side->facets();
              //               if ( child_facets.size()==0 )
              //               {
              //                 FOUR_C_THROW( "Wo sind die facets?" );
              //               }

              for (std::vector<Facet*>::const_iterator i = child_facets.begin();
                  i != child_facets.end(); ++i)
              {
                Facet* f = *i;
                const plain_volumecell_set& child_cells = f->cells();
                if (child_cells.size() == 2)
                {
                  std::vector<VolumeCell*> child_cell_vector;
                  child_cell_vector.reserve(2);
                  child_cell_vector.assign(child_cells.begin(), child_cells.end());
                  if (parent_cell_info[doneindex]->contains_child(child_cell_vector[0]))
                  {
                    VolumeCell* c = child_cell_vector[1];
                    if (done_child_cells.count(c) == 0)
                    {
                      parent_cell_info[otherindex]->cells_.insert(c);
                      done_child_cells.insert(c);
                      found = true;
                    }
                  }
                  else if (parent_cell_info[doneindex]->contains_child(child_cell_vector[1]))
                  {
                    VolumeCell* c = child_cell_vector[0];
                    if (done_child_cells.count(c) == 0)
                    {
                      parent_cell_info[otherindex]->cells_.insert(c);
                      done_child_cells.insert(c);
                      found = true;
                    }
                  }
                  else
                  {
                    FOUR_C_THROW("child must be part of done parent cell");
                  }
                }
                else if (child_cells.size() == 1)
                {
                  VolumeCell* c = *child_cells.begin();
                  if (not parent_cell_info[doneindex]->contains_child(c))
                  {
                    if (done_child_cells.count(c) == 0)
                    {
                      parent_cell_info[otherindex]->cells_.insert(c);
                      done_child_cells.insert(c);
                      found = true;
                    }
                  }
                }
                else
                {
                  std::stringstream str;
                  str << "illegal number of neighbouring volume cells: child_cells.size() == "
                      << child_cells.size();
                  FOUR_C_THROW("{}", str.str());
                }
              }

              if (found)
              {
                ChildCell& cc = *parent_cell_info[otherindex];
                fill(cc.parent_, cc);
                std::copy(cc.cells_.begin(), cc.cells_.end(),
                    std::inserter(done_child_cells, done_child_cells.begin()));
                nonnodecells -= 1;
              }
              else
              {
                // std::cout << "not found after " << child_facets.size() << " facets from " <<
                // child_side << "\n";
              }
            }
            else
            {
              // std::cout << "parent_cell_info[0]->done_=" << parent_cell_info[0]->done_ << "  "
              //          << "parent_cell_info[1]->done_=" << parent_cell_info[1]->done_ << "\n";
            }
          }
          else
          {
            // std::cout << "parent_cell_info.size()=" << parent_cell_info.size() << "\n";
          }
        }
        else
        {
          // std::cout << "parent_facets.size()=" << parent_facets.size() << "\n";
        }
      }
    }

    if (nonnodecells == 1)
    {
      VolumeCell* parent_vc = nullptr;
      ChildCell* child_cells = nullptr;
      for (std::map<VolumeCell*, ChildCell>::iterator i = cellmap.begin(); i != cellmap.end(); ++i)
      {
        VolumeCell* vc = i->first;
        ChildCell& cc = i->second;
        if (not cc.done_)
        {
          if (parent_vc == nullptr)
          {
            parent_vc = vc;
            child_cells = &cc;
          }
          else
          {
            FOUR_C_THROW("more than one open parent cells");
          }
        }
      }
      if (parent_vc == nullptr)
      {
        FOUR_C_THROW("no open parent cell");
      }

      ChildCell& cc = *child_cells;
      plain_volumecell_set& childset = cc.cells_;

      plain_volumecell_set done_child_cells;

      for (std::map<VolumeCell*, ChildCell>::iterator i = cellmap.begin(); i != cellmap.end(); ++i)
      {
        // VolumeCell * vc = i->first;
        ChildCell& cc = i->second;
        std::copy(cc.cells_.begin(), cc.cells_.end(),
            std::inserter(done_child_cells, done_child_cells.begin()));
      }

      const std::list<std::shared_ptr<VolumeCell>>& all_child_cells = mesh_.volume_cells();
      for (std::list<std::shared_ptr<VolumeCell>>::const_iterator i = all_child_cells.begin();
          i != all_child_cells.end(); ++i)
      {
        VolumeCell* child_vc = &**i;
        if (done_child_cells.count(child_vc) == 0)
        {
          childset.insert(child_vc);
        }
      }
      if (childset.size() > 0)
      {
        fill(parent_vc, cc);
        nonnodecells -= 1;
      }
      else
      {
        // FOUR_C_THROW( "no child cell for open parent cell" );

        // Empty parent cell. We did not get any children. The cell is most
        // probably too small.
        cc.done_ = true;
        nonnodecells -= 1;
      }
    }

    if (nonnodecells > 0)
    {
      // test if there are any volume cells left

      plain_volumecell_set done_child_cells;
      for (std::map<VolumeCell*, ChildCell>::iterator i = cellmap.begin(); i != cellmap.end(); ++i)
      {
        ChildCell& cc = i->second;
        std::copy(cc.cells_.begin(), cc.cells_.end(),
            std::inserter(done_child_cells, done_child_cells.begin()));
      }

      bool found = false;
      const std::list<std::shared_ptr<VolumeCell>>& all_child_cells = mesh_.volume_cells();
      for (std::list<std::shared_ptr<VolumeCell>>::const_iterator i = all_child_cells.begin();
          i != all_child_cells.end(); ++i)
      {
        VolumeCell* child_vc = &**i;
        if (done_child_cells.count(child_vc) == 0)
        {
          found = true;
          break;
        }
      }

      if (not found)
      {
        // Done. There are a few empty parent cells. We do not mind.
        nonnodecells = 0;
        for (std::map<VolumeCell*, ChildCell>::iterator i = cellmap.begin(); i != cellmap.end();
            ++i)
        {
          ChildCell& cc = i->second;
          cc.done_ = true;
        }
      }
    }

    if (backup == nonnodecells) FOUR_C_THROW("no progress in child cell--parent cell mapping");
  }

  for (std::map<VolumeCell*, ChildCell>::iterator i = cellmap.begin(); i != cellmap.end(); ++i)
  {
    VolumeCell* vc = i->first;
    ChildCell& cc = i->second;
    plain_volumecell_set& childset = cc.cells_;

    if (not cc.done_)
    {
      // finish partly filled volume cells
      fill(vc, cc);
    }

    register_new_points(parent_mesh, childset);
  }

  // copy volume cell position

  if (mesh_.create_options().find_positions())
  {
    for (std::map<VolumeCell*, ChildCell>::iterator i = cellmap.begin(); i != cellmap.end(); ++i)
    {
      VolumeCell* vc = i->first;
      ChildCell& cc = i->second;
      plain_volumecell_set& childset = cc.cells_;

      Point::PointPosition pos = vc->position();
      //       if ( pos==Point::undecided )
      //         FOUR_C_THROW( "undecided volume cell" );
      if (pos != Point::undecided)
        for (plain_volumecell_set::iterator i = childset.begin(); i != childset.end(); ++i)
        {
          VolumeCell* c = *i;
          c->position(pos);
        }
    }
  }
}

void Cut::TetMeshIntersection::seed_cells(Mesh& parent_mesh,
    const plain_volumecell_set& parent_cells, std::map<VolumeCell*, ChildCell>& cellmap,
    plain_volumecell_set& done_child_cells)
{
  std::map<Point*, std::vector<VolumeCell*>> parent_point_cells;

  for (plain_volumecell_set::const_iterator i = parent_cells.begin(); i != parent_cells.end(); ++i)
  {
    VolumeCell* vc = *i;
    ChildCell& cc = cellmap[vc];
    cc.parent_ = vc;
    plain_volumecell_set& childset = cc.cells_;

    PointSet volume_points;
    vc->get_all_points(parent_mesh, volume_points);

    // seed cells at parent element nodes (if unique)

    for (PointSet::iterator i = volume_points.begin(); i != volume_points.end(); ++i)
    {
      Point* p = *i;
      if (p->position() != Point::oncutsurface)
      {
        Point* np = to_child(p);
        find_volume_cell(np, childset);
      }
    }

    for (PointSet::iterator i = volume_points.begin(); i != volume_points.end(); ++i)
    {
      Point* p = *i;
      parent_point_cells[p].push_back(vc);
    }
  }

  // seed cells with unique point

  for (std::map<Point*, std::vector<VolumeCell*>>::iterator i = parent_point_cells.begin();
      i != parent_point_cells.end(); ++i)
  {
    Point* p = i->first;
    std::vector<VolumeCell*>& vcs = i->second;
    if (vcs.size() == 1)
    {
      VolumeCell* vc = vcs[0];
      Point* np = to_child(p);
      ChildCell& cc = cellmap[vc];
      plain_volumecell_set& childset = cc.cells_;
      find_volume_cell(np, childset);
    }
  }

  // collect done cells

  for (std::map<VolumeCell*, ChildCell>::iterator i = cellmap.begin(); i != cellmap.end(); ++i)
  {
    ChildCell& cc = i->second;
    plain_volumecell_set& childset = cc.cells_;
    std::copy(childset.begin(), childset.end(),
        std::inserter(done_child_cells, done_child_cells.begin()));
  }

  // look at all points of each free child volume cell and see if there is a
  // unique parent volume cell to these points

  const std::list<std::shared_ptr<VolumeCell>>& all_child_cells = mesh_.volume_cells();
  for (std::list<std::shared_ptr<VolumeCell>>::const_iterator i = all_child_cells.begin();
      i != all_child_cells.end(); ++i)
  {
    VolumeCell* child_vc = &**i;
    if (done_child_cells.count(child_vc) == 0)
    {
      PointSet child_cut_points;
      child_vc->get_all_points(mesh_, child_cut_points);

      // Remove all points that are new in the child mesh. Those do not at
      // all help to find the parent cell.
      for (PointSet::iterator i = child_cut_points.begin(); i != child_cut_points.end();)
      {
        Point* p = *i;
        if (child_to_parent_.count(p) == 0)
        {
          set_erase(child_cut_points, i);
        }
        else
        {
          ++i;
        }
      }

      if (child_cut_points.size() > 0)
      {
        plain_volumecell_set used_parent_cells;

        PointSet::iterator j = child_cut_points.begin();
        Point* p = *j;
        find_volume_cell(to_parent(p), used_parent_cells);

        for (++j; j != child_cut_points.end(); ++j)
        {
          Point* p = *j;
          plain_volumecell_set upc;
          find_volume_cell(to_parent(p), upc);

          plain_volumecell_set intersection;
          std::set_intersection(used_parent_cells.begin(), used_parent_cells.end(), upc.begin(),
              upc.end(), std::inserter(intersection, intersection.begin()));

          std::swap(used_parent_cells, intersection);

          if (used_parent_cells.size() == 0) FOUR_C_THROW("no possible parent cell");
        }

        if (used_parent_cells.size() == 1)
        {
          VolumeCell* parent_vc = *used_parent_cells.begin();
          ChildCell& cc = cellmap[parent_vc];
          if (cc.done_)
          {
            // FOUR_C_THROW( "free child cell to done parent cell?" );
          }
          else
          {
            plain_volumecell_set& childset = cc.cells_;
            childset.insert(child_vc);
            done_child_cells.insert(child_vc);
          }
        }
      }
      else
      {
        // FOUR_C_THROW( "child cell with all new points?" );
      }
    }
  }
}

void Cut::TetMeshIntersection::build_surface_cell_map(VolumeCell* vc, ChildCell& cc)
{
  // find parent facets on cut surface

  std::map<Side*, std::vector<Facet*>>& facetsonsurface = cc.facetsonsurface_;

  const plain_facet_set& parent_facets = vc->facets();
  for (plain_facet_set::const_iterator i = parent_facets.begin(); i != parent_facets.end(); ++i)
  {
    Facet* f = *i;
    if (f->on_boundary_cell_side())
    {
      Side* s = f->parent_side();
      std::map<Side*, std::vector<Side*>>::iterator j = side_parent_to_child_.find(s);
      if (j == side_parent_to_child_.end()) FOUR_C_THROW("unknown parent cut facet");
      std::vector<Side*>& side_vector = j->second;
      for (std::vector<Side*>::iterator i = side_vector.begin(); i != side_vector.end(); ++i)
      {
        Side* cs = *i;
        facetsonsurface[cs].push_back(f);
      }
    }
  }
}

void Cut::TetMeshIntersection::fill(Mesh& parent_mesh, Element* element,
    const plain_volumecell_set& parent_cells, std::map<VolumeCell*, ChildCell>& cellmap)
{
  for (std::map<VolumeCell*, ChildCell>::iterator i = cellmap.begin(); i != cellmap.end(); ++i)
  {
    VolumeCell* parent_cell = i->first;
    ChildCell& cc = i->second;
    plain_volumecell_set& childset = cc.cells_;
    std::map<Side*, std::vector<Facet*>>& facetsonsurface = cc.facetsonsurface_;

    for (plain_volumecell_set::iterator i = childset.begin(); i != childset.end(); ++i)
    {
      VolumeCell* vc = *i;
      const plain_integrationcell_set& cells = vc->integration_cells();
      const plain_boundarycell_set& bcells = vc->boundary_cells();

      for (plain_integrationcell_set::const_iterator i = cells.begin(); i != cells.end(); ++i)
      {
        IntegrationCell* ic = *i;
        const std::vector<Point*>& points = ic->points();

        std::vector<Point*> parent_points(points);
        to_parent(parent_mesh, parent_points);

        parent_cell->new_integration_cell(parent_mesh, ic->shape(), parent_points);
      }
      for (plain_boundarycell_set::const_iterator i = bcells.begin(); i != bcells.end(); ++i)
      {
        BoundaryCell* bc = *i;
        const std::vector<Point*>& points = bc->points();

        std::vector<Point*> parent_points(points);
        to_parent(parent_mesh, parent_points);

        Facet* parent_facet = nullptr;
        Facet* child_facet = bc->get_facet();

        if (not child_facet->on_boundary_cell_side())
          FOUR_C_THROW("boundary cell not on cut surface");

        std::vector<Facet*> facets;
        std::stringstream str;

        std::map<Side*, std::vector<Facet*>>::iterator j =
            facetsonsurface.find(child_facet->parent_side());
        if (j == facetsonsurface.end())
        {
          // Parent side not included in facetsonsurface. This might be a
          // numerical problem. This might be a real bug.

          str << (*child_facet) << "\non boundary cell " << (*child_facet->parent_side())
              << "\non unknown cut surface: ";

          // Try to recover.
          if (not child_facet->has_holes() and not child_facet->is_triangulated())
          {
            for (std::map<Side*, std::vector<Side*>>::iterator i = side_parent_to_child_.begin();
                i != side_parent_to_child_.end(); ++i)
            {
              Side* parent_side = i->first;
              std::vector<Side*>& children = i->second;
              if (std::find(children.begin(), children.end(), child_facet->parent_side()) !=
                  children.end())
              {
                const std::vector<Facet*>& fs = parent_side->facets();
                for (std::vector<Facet*>::const_iterator i = fs.begin(); i != fs.end(); ++i)
                {
                  Facet* f = *i;

                  if (parent_cell->facets().count(f) > 0)
                  {
                    facets.push_back(f);
                  }
                }
              }
            }
          }

          //           continue;
        }
        else
        {
          facets = j->second;
        }

        if (facets.size() == 1)
        {
          parent_facet = facets[0];
        }
        else if (facets.size() > 1)
        {
          // this can happen with levelset

          // get facet points in case those differ from boundary cell points
          std::vector<Point*> facet_points = child_facet->points();
          to_parent(parent_mesh, facet_points);

          // search for matching parent facet
          for (std::vector<Facet*>::iterator i = facets.begin(); i != facets.end(); ++i)
          {
            Facet* f = *i;
            if (f->contains_some(facet_points))
            {
              if (parent_facet == nullptr)
              {
                parent_facet = f;
              }
              else
              {
                str << "parent facet not unique";
                FOUR_C_THROW("{}", str.str());
              }
            }
          }
          if (parent_facet == nullptr)
          {
            str << "no parent facet found";
            FOUR_C_THROW("{}", str.str());
          }
        }
        else
        {
          str << "empty list bug";
          FOUR_C_THROW("{}", str.str());
        }

        parent_cell->new_boundary_cell(parent_mesh, bc->shape(), parent_facet, parent_points);
      }
    }
  }
}

void Cut::TetMeshIntersection::fill(VolumeCell* parent_cell, ChildCell& childcell)
{
  plain_volumecell_set& child_cells = childcell.cells_;

  if (child_cells.size() == 0)
  {
    FOUR_C_THROW("failed to find seed cells");
  }

  plain_volumecell_set done_child_cells;

  while (child_cells.size() > 0)
  {
    plain_facet_set open_facets;

    while (child_cells.size() > 0)
    {
      for (plain_volumecell_set::iterator i = child_cells.begin(); i != child_cells.end(); ++i)
      {
        VolumeCell* vc = *i;
        if (done_child_cells.count(vc) == 0)
        {
          child_cells.erase(vc);
          done_child_cells.insert(vc);
          const plain_facet_set& facets = vc->facets();
          for (plain_facet_set::const_iterator i = facets.begin(); i != facets.end(); ++i)
          {
            Facet* f = *i;
            if (not f->on_boundary_cell_side())
            {
              VolumeCell* nc = f->neighbor(vc);
              if (nc != nullptr)
              {
                if (done_child_cells.count(nc) == 0)
                {
                  child_cells.insert(nc);
                }
              }
              else
              {
                if (not f->has_holes() and not f->is_triangulated() and f->points().size() == 3)
                {
                  open_facets.insert(f);
                }
              }
            }
          }
          break;
        }
      }
    }

    // search for four point parent "facets" that are split by a flat tet

    if (open_facets.size() > 0)
    {
      FacetMesh facet_mesh;

      for (plain_facet_set::iterator i = open_facets.begin(); i != open_facets.end(); ++i)
      {
        Facet* f = *i;
        facet_mesh.add(f);
      }

      for (std::map<std::pair<Point*, Point*>, std::vector<Facet*>>::iterator i =
               facet_mesh.facet_mesh_.begin();
          i != facet_mesh.facet_mesh_.end(); ++i)
      {
        const std::pair<Point*, Point*>& line = i->first;
        const std::vector<Facet*>& facets = i->second;

        if (facets.size() == 2)
        {
          Facet* f1 = facets[0];
          Facet* f2 = facets[1];

          Point* p1 = f1->other_point(line.first, line.second);
          Point* p2 = f2->other_point(line.first, line.second);

          plain_facet_set facets1;
          plain_facet_set facets2;

          find_common_facets(p1, p2, line.first, facets1);
          find_common_facets(p1, p2, line.second, facets2);

          if (facets1.size() == 1 and facets2.size() == 1)
          {
            Facet* f3 = *facets1.begin();
            Facet* f4 = *facets2.begin();

            if (not f3->on_boundary_cell_side() and not f4->on_boundary_cell_side())
            {
              const plain_volumecell_set& cells3 = f3->cells();
              const plain_volumecell_set& cells4 = f4->cells();

              if (cells3.size() == 1 and cells4.size() == 1)
              {
                VolumeCell* vc3 = *cells3.begin();
                VolumeCell* vc4 = *cells4.begin();

                if (done_child_cells.count(vc3) == 0 and done_child_cells.count(vc4) == 0)
                {
                  child_cells.insert(vc3);
                  child_cells.insert(vc4);

                  facet_mesh.erase(f1);
                  facet_mesh.erase(f2);
                }
              }
            }
          }
        }
      }
    }
  }

  child_cells.swap(done_child_cells);
  childcell.done_ = true;

  build_surface_cell_map(parent_cell, childcell);
}

void Cut::TetMeshIntersection::register_new_points(
    Mesh& parent_mesh, const plain_volumecell_set& childset)
{
  for (plain_volumecell_set::const_iterator i = childset.begin(); i != childset.end(); ++i)
  {
    VolumeCell* vc = *i;
    const plain_facet_set& facets = vc->facets();
    for (plain_facet_set::const_iterator i = facets.begin(); i != facets.end(); ++i)
    {
      Facet* f = *i;
      if (f->on_boundary_cell_side())
      {
        PointSet points;
        f->all_points(points);
        for (PointSet::iterator i = points.begin(); i != points.end(); ++i)
        {
          Point* p = *i;
          if (child_to_parent_.count(p) == 0)
          {
            Point* pp = parent_mesh.new_point(p->x(), nullptr, nullptr, p->tolerance());
            do_register(pp, p);
          }
        }
      }
    }
  }
}

void Cut::TetMeshIntersection::find_volume_cell(Point* p, plain_volumecell_set& childset)
{
  const plain_facet_set& facets = p->facets();
  for (plain_facet_set::const_iterator i = facets.begin(); i != facets.end(); ++i)
  {
    Facet* f = *i;
    const plain_volumecell_set& facet_cells = f->cells();
    std::copy(facet_cells.begin(), facet_cells.end(), std::inserter(childset, childset.begin()));
  }
}

void Cut::TetMeshIntersection::swap_points(
    Mesh& mesh, const std::map<Point*, Point*>& pointmap, std::vector<Point*>& points)
{
  std::vector<Point*> new_points;
  new_points.reserve(points.size());
  for (std::vector<Point*>::iterator i = points.begin(); i != points.end(); ++i)
  {
    Point* p = *i;
    std::map<Point*, Point*>::const_iterator j = pointmap.find(p);
    if (j == pointmap.end())
    {
      Point* np = mesh.new_point(p->x(), nullptr, nullptr, p->tolerance());
      new_points.push_back(np);
    }
    else
    {
      new_points.push_back(j->second);
    }
  }
  std::swap(new_points, points);
}

void Cut::TetMeshIntersection::swap_points(
    const std::map<Point*, Point*>& pointmap, std::vector<Point*>& points)
{
  std::vector<Point*> new_points;
  new_points.reserve(points.size());
  for (std::vector<Point*>::iterator i = points.begin(); i != points.end(); ++i)
  {
    Point* p = *i;
    std::map<Point*, Point*>::const_iterator j = pointmap.find(p);
    if (j == pointmap.end())
    {
      FOUR_C_THROW("no such point");
    }
    new_points.push_back(j->second);
  }
  std::swap(new_points, points);
}

void Cut::TetMeshIntersection::swap_points(
    const std::map<Point*, Point*>& pointmap, PointSet& points)
{
  PointSet new_points;
  for (PointSet::iterator i = points.begin(); i != points.end(); ++i)
  {
    Point* p = *i;
    std::map<Point*, Point*>::const_iterator j = pointmap.find(p);
    if (j == pointmap.end())
    {
      FOUR_C_THROW("no such point");
    }
    new_points.insert(j->second);
  }
  std::swap(new_points, points);
}

Cut::Point* Cut::TetMeshIntersection::swap_point(
    const std::map<Point*, Point*>& pointmap, Point* point)
{
  std::map<Point*, Point*>::const_iterator j = pointmap.find(point);
  if (j == pointmap.end())
  {
    // FOUR_C_THROW( "no such point" );
    return nullptr;
  }
  return j->second;
}

void Cut::TetMeshIntersection::do_register(Point* parent_point, Point* child_point)
{
  child_to_parent_[child_point] = parent_point;
  parent_to_child_[parent_point] = child_point;
}

void Cut::TetMeshIntersection::copy_cut_side(Side* s, Facet* f)
{
  const std::vector<Node*>& nodes = s->nodes();
  std::vector<int> nids;
  nids.reserve(nodes.size());
  for (std::vector<Node*>::const_iterator i = nodes.begin(); i != nodes.end(); ++i)
  {
    Node* n = *i;
    nids.push_back(n->id());
    Point* p = n->point();

    Node* new_node = cut_mesh_.get_node(n->id(), p->x());
    Point* np = to_child(p);
    if (np != nullptr)
    {
      if (new_node->point() != np)
      {
        FOUR_C_THROW("did not catch known cut point");
      }
    }
    else
    {
      do_register(p, new_node->point());
    }
  }

  Side* cs = cut_mesh_.get_side(s->id(), nids, s->topology());

  side_parent_to_child_[s].push_back(cs);

  // Copy cut point to cut surfaces, since a second cut search could result
  // in different cut points.

  const std::vector<Edge*>& old_edges = s->edges();
  const std::vector<Edge*>& new_edges = cs->edges();

  for (std::vector<Edge*>::const_iterator ei = old_edges.begin(); ei != old_edges.end(); ++ei)
  {
    Edge* e = *ei;
    Edge* ne = new_edges[std::distance(old_edges.begin(), ei)];
    const PointPositionSet& cutpoints = e->cut_points();
    for (PointPositionSet::const_iterator i = cutpoints.begin(); i != cutpoints.end(); ++i)
    {
      Point* p = *i;
      Point* np = to_child(p);

      if (np != nullptr)
      {
        np->add_edge(ne);
        np->position(Point::oncutsurface);
      }
      else
      {
        np = Point::new_point(mesh_, p->x(), p->t(e), ne, nullptr, p->tolerance());
        np->position(Point::oncutsurface);
        do_register(p, np);
      }
    }
  }

  // Copy cut points from facets. If the facets is triangulated, there is a
  // middle point that needs to be introduces as a cut point.

  PointSet points;
  f->all_points(points);
  for (PointSet::iterator i = points.begin(); i != points.end(); ++i)
  {
    Point* p = *i;
    Point* np = to_child(p);
    if (np != nullptr)
    {
      np->add_side(cs);
      np->position(Point::oncutsurface);
    }
  }

  // Copy cut points from cut lines to cut surfaces. We cannot copy cut lines
  // here, since there might be additional cut points.

  const std::vector<Line*>& cutlines = s->cut_lines();

  for (std::vector<Line*>::const_iterator i = cutlines.begin(); i != cutlines.end(); ++i)
  {
    Line* l = *i;
    Point* p1 = to_child(l->begin_point());
    Point* p2 = to_child(l->end_point());

    if (p1 != nullptr)
    {
      p1->add_side(cs);
    }

    if (p2 != nullptr)
    {
      p2->add_side(cs);
    }
  }
}

FOUR_C_NAMESPACE_CLOSE
