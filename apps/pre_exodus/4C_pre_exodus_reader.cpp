// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_pre_exodus_reader.hpp"

#include "4C_fem_general_utils_local_connectivity_matrices.hpp"

#include <exodusII.h>
#include <Teuchos_Time.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <fstream>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 12/07|
 *----------------------------------------------------------------------*/
EXODUS::Mesh::Mesh(const std::string exofilename)
{
  int error;
  int CPU_word_size, IO_word_size;
  float exoversion;               /* version of exodus */
  CPU_word_size = sizeof(double); /* size of a double */
  IO_word_size = 0;               /* use what is stored in file */

  const char* exofilenamechar = exofilename.c_str();

  // open EXODUS II file
  int exo_handle = ex_open(exofilenamechar, EX_READ, &CPU_word_size, &IO_word_size, &exoversion);
  if (exo_handle <= 0) FOUR_C_THROW("Error while opening EXODUS II file {}", exofilenamechar);

  // print version
  std::cout << "File " << exofilename << " was created with EXODUS II library version "
            << exoversion << std::endl;

  // read database parameters
  int num_elem_blk, num_node_sets, num_side_sets, num_nodes;
  char title[MAX_LINE_LENGTH + 1];
  error = ex_get_init(exo_handle, title, &four_c_dim_, &num_nodes, &num_elem_, &num_elem_blk,
      &num_node_sets, &num_side_sets);
  title_ = std::string(title);

  num_dim_ = 3;

  // get nodal coordinates
  {
    std::vector<double> x(num_nodes);
    std::vector<double> y(num_nodes);
    std::vector<double> z(num_nodes);
    error = ex_get_coord(exo_handle, x.data(), y.data(), z.data());
    if (error != 0) FOUR_C_THROW("exo error returned");

    // store nodes in map
    nodes_ = std::make_shared<std::map<int, std::vector<double>>>();
    for (int i = 0; i < num_nodes; ++i)
    {
      std::vector<double> coords;
      coords.push_back(x[i]);
      coords.push_back(y[i]);
      coords.push_back(z[i]);
      nodes_->insert(std::pair<int, std::vector<double>>(
          i + 1, coords));  // to store the EXO-ID starting with 1
    }
  }  // free coordinate vectors x, y ,z

  // Get all ElementBlocks
  {
    std::vector<int> epropID(num_elem_blk);
    std::vector<int> ebids(num_elem_blk);
    error = ex_get_ids(exo_handle, EX_ELEM_BLOCK, ebids.data());
    if (error != 0) FOUR_C_THROW("exo error returned");
    error = ex_get_prop_array(exo_handle, EX_ELEM_BLOCK, "ID", epropID.data());
    if (error != 0) FOUR_C_THROW("exo error returned");
    for (int i = 0; i < num_elem_blk; ++i)
    {
      // Read Element Blocks into Map
      char mychar[MAX_STR_LENGTH + 1];
      int num_el_in_blk, num_nod_per_elem, num_attr;
      // error = ex_get_elem_block (exo_handle, epropID[i], mychar, &num_el_in_blk,
      // &num_nod_per_elem, &num_attr);
      error = ex_get_block(exo_handle, EX_ELEM_BLOCK, ebids[i], mychar, &num_el_in_blk,
          &num_nod_per_elem, nullptr, nullptr, &num_attr);
      if (error != 0) FOUR_C_THROW("exo error returned");
      // prefer std::string to store element type
      std::string ele_type(mychar);

      // get ElementBlock name
      error = ex_get_name(exo_handle, EX_ELEM_BLOCK, ebids[i], mychar);
      if (error != 0) FOUR_C_THROW("exo error returned");
      // prefer std::string to store name
      std::string blockname(mychar);

      // get element connectivity
      std::vector<int> allconn(num_nod_per_elem * num_el_in_blk);
      error = ex_get_conn(exo_handle, EX_ELEM_BLOCK, ebids[i], allconn.data(), nullptr, nullptr);
      if (error != 0) FOUR_C_THROW("exo error returned");
      std::shared_ptr<std::map<int, std::vector<int>>> eleconn =
          std::make_shared<std::map<int, std::vector<int>>>();
      for (int j = 0; j < num_el_in_blk; ++j)
      {
        std::vector<int> actconn;
        actconn.reserve(num_nod_per_elem);
        for (int k = 0; k < num_nod_per_elem; ++k)
        {
          actconn.push_back(allconn[k + j * num_nod_per_elem]);
        }
        eleconn->insert(std::pair<int, std::vector<int>>(j, actconn));
      }
      std::shared_ptr<ElementBlock> actEleBlock =
          std::make_shared<ElementBlock>(string_to_shape(ele_type), eleconn, blockname);

      // Add this ElementBlock into Mesh map
      element_blocks_.insert(std::pair<int, std::shared_ptr<ElementBlock>>(ebids[i], actEleBlock));
    }
  }  // end of element section

  // get all NodeSets
  {
    std::map<int, NodeSet> prelimNodeSets;  // prelim due to possible prop names
    std::vector<int> npropID(num_node_sets);
    error = ex_get_prop_array(exo_handle, EX_NODE_SET, "ID", npropID.data());
    for (int i = 0; i < num_node_sets; ++i)
    {
      // Read NodeSet params
      int num_nodes_in_set, num_df_in_set;
      error =
          ex_get_set_param(exo_handle, EX_NODE_SET, npropID[i], &num_nodes_in_set, &num_df_in_set);

      // get NodeSet name
      char mychar[MAX_STR_LENGTH + 1];
      error = ex_get_name(exo_handle, EX_NODE_SET, npropID[i], mychar);
      // prefer std::string to store name
      std::string nodesetname(mychar);

      // get nodes in node set
      std::vector<int> node_set_node_list(num_nodes_in_set);
      error = ex_get_set(exo_handle, EX_NODE_SET, npropID[i], node_set_node_list.data(), nullptr);
      if (error > 0)
        std::cout << "'ex_get_set' for EX_NODE_SET returned warning while reading node set "
                  << npropID[i] << std::endl;
      else if (error < 0)
        FOUR_C_THROW("error reading node set");
      std::set<int> nodes_in_set;
      for (int j = 0; j < num_nodes_in_set; ++j) nodes_in_set.insert(node_set_node_list[j]);
      NodeSet actNodeSet(nodes_in_set, nodesetname, "none");

      // Add this NodeSet into Mesh map (here prelim due to pro names)
      prelimNodeSets.insert(std::pair<int, NodeSet>(npropID[i], actNodeSet));
    }

    /* Read NodeSet property names ***********************************************
     * They are assigned by ICEM and provide recognition */
    int num_props;
    float fdum;
    char cdum;  // dummy argument
    error = ex_inquire(exo_handle, EX_INQ_NS_PROP, &num_props, &fdum, &cdum);
    // allocate memory for NodeSet property names
    char** prop_names = new char*[num_props];
    for (int i = 0; i < num_props; ++i)
    {
      prop_names[i] = new char[MAX_STR_LENGTH + 1];
    }
    // get prop names of node sets
    error = ex_get_prop_names(exo_handle, EX_NODE_SET, prop_names);

    // Add prop names to final Mesh NodeSet if available
    std::map<int, NodeSet>::const_iterator i_ns;
    if ((num_props - 1) == num_node_sets)
    {
      int i = 1;  // id of propname, starts with 1 because 0 is "ID"
      for (i_ns = prelimNodeSets.begin(); i_ns != prelimNodeSets.end(); ++i_ns)
      {
        std::string propname(prop_names[i]);
        const NodeSet actNodeSet = i_ns->second;
        std::string ns_name = actNodeSet.get_name();
        if (ns_name.size() == 0) ns_name = propname;
        const NodeSet newNodeSet(actNodeSet.get_node_set(), ns_name, propname);
        node_sets_.insert(std::pair<int, NodeSet>(i_ns->first, newNodeSet));
        ++i;  // next propname refers to next NodeSet
      }
    }
    else
    {
      // this is the standard case without prop names
      node_sets_ = prelimNodeSets;
    }

    // clean up node set names
    for (int i = 0; i < num_props; i++)
    {
      delete[] prop_names[i];
    }
    delete[] prop_names;
  }  // end of nodeset section
  // ***************************************************************************

  // get all SideSets
  if (num_side_sets > 0)
  {
    std::vector<int> spropID(num_side_sets);
    error = ex_get_prop_array(exo_handle, EX_SIDE_SET, "ID", spropID.data());
    for (int i = 0; i < num_side_sets; ++i)
    {
      // get SideSet name
      char mychar[MAX_STR_LENGTH + 1];
      error = ex_get_name(exo_handle, EX_SIDE_SET, spropID[i], mychar);
      // prefer std::string to store name
      std::string sidesetname(mychar);

      // Read SideSet params
      int num_side_in_set, num_dist_fact_in_set;
      error = ex_get_set_param(
          exo_handle, EX_SIDE_SET, spropID[i], &num_side_in_set, &num_dist_fact_in_set);

      // get SideSet
      std::vector<int> side_set_elem_list(num_side_in_set);
      std::vector<int> side_set_side_list(num_side_in_set);
      error = ex_get_set(exo_handle, EX_SIDE_SET, spropID[i], side_set_elem_list.data(),
          side_set_side_list.data());
      if (error != 0) FOUR_C_THROW("error reading side set");
      std::map<int, std::vector<int>> sides_in_set;
      for (int j = 0; j < num_side_in_set; ++j)
      {
        std::vector<int> side(2);  // first entry is element, second side
        side[0] = side_set_elem_list[j];
        side[1] = side_set_side_list[j];
        sides_in_set.insert(std::pair<int, std::vector<int>>(j, side));
      }

      SideSet actSideSet(sides_in_set, sidesetname);

      // Add this SideSet into Mesh map
      side_sets_.insert(std::pair<int, SideSet>(spropID[i], actSideSet));
    }
  }  // end of sideset section

  error = ex_close(exo_handle);
  if (error < 0) FOUR_C_THROW("error while closing exodus II file");
}


/*----------------------------------------------------------------------*
 |  Print method (public)                                      maf 12/07|
 *----------------------------------------------------------------------*/
void EXODUS::Mesh::print(std::ostream& os, bool verbose) const
{
  os << "Mesh consists of ";
  os << get_num_nodes() << " Nodes, ";
  os << num_elem_ << " Elements, organized in " << std::endl;
  os << get_num_element_blocks() << " ElementBlocks, ";
  os << get_num_node_sets() << " NodeSets, ";
  os << get_num_side_sets() << " SideSets ";
  os << std::endl << std::endl;
  if (verbose)
  {
    os << "ElementBlocks" << std::endl;
    std::map<int, std::shared_ptr<ElementBlock>>::const_iterator it;
    std::map<int, std::shared_ptr<ElementBlock>> eleBlocks = get_element_blocks();
    for (it = eleBlocks.begin(); it != eleBlocks.end(); it++)
    {
      os << it->first << ": ";
      it->second->print(os);
    }
    os << std::endl << "NodeSets" << std::endl;
    std::map<int, NodeSet>::const_iterator it2;
    std::map<int, NodeSet> nodeSets = get_node_sets();
    for (it2 = nodeSets.begin(); it2 != nodeSets.end(); it2++)
    {
      os << "NodeSet " << it2->first << ": ";
      it2->second.print(os);
    }
    os << std::endl << "SideSets" << std::endl;
    os << "Warning: SideSets are not yet fully supported by PreExodus!" << std::endl;
    std::map<int, SideSet>::const_iterator it3;
    std::map<int, SideSet> sideSets = get_side_sets();
    for (it3 = sideSets.begin(); it3 != sideSets.end(); it3++)
    {
      os << "SideSet " << it3->first << ": ";
      it3->second.print(os);
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<EXODUS::ElementBlock> EXODUS::Mesh::get_element_block(const int id) const
{
  if (element_blocks_.find(id) == element_blocks_.end())
    FOUR_C_THROW("ElementBlock {} not found.", id);

  return (element_blocks_.find(id))->second;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
EXODUS::NodeSet EXODUS::Mesh::get_node_set(const int id) const
{
  if (node_sets_.find(id) == node_sets_.end()) FOUR_C_THROW("NodeSet {} not found.", id);

  return (node_sets_.find(id))->second;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
EXODUS::SideSet EXODUS::Mesh::get_side_set(const int id) const
{
  if (side_sets_.find(id) == side_sets_.end()) FOUR_C_THROW("SideSet {} not found.", id);

  return (side_sets_.find(id))->second;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::Mesh::print_nodes(std::ostream& os, bool storeid) const
{
  std::map<int, std::vector<double>>::const_iterator it;
  for (it = nodes_->begin(); it != nodes_->end(); it++)
  {
    if (storeid) os << "MapID: " << it->first;
    int exoid = it->first + 1;
    os << " ExoID: " << exoid << " : ";
    const std::vector<double> mycoords = it->second;
    for (int i = 0; i < signed(mycoords.size()); i++)
    {
      os << mycoords[i] << ",";
    }
    os << std::endl;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> EXODUS::Mesh::get_node(const int NodeID) const
{
  std::map<int, std::vector<double>>::const_iterator it = nodes_->find(NodeID);
  return it->second;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::Mesh::set_node(const int NodeID, const std::vector<double> coord)
{
  // if entry exits already , delete it first and the insert the new value
  // other wise nothing is inserted
  if (nodes_->find(NodeID) != nodes_->end()) nodes_->erase(NodeID);

  nodes_->insert(std::make_pair(NodeID, coord));
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::Mesh::set_nsd(const int nsd)
{
  if (nsd != 2 && nsd != 3) return;

  four_c_dim_ = nsd;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> EXODUS::Mesh::normal(const int head1, const int origin, const int head2) const
{
  std::vector<double> normal(3);
  std::vector<double> h1 = get_node(head1);
  std::vector<double> h2 = get_node(head2);
  std::vector<double> o = get_node(origin);

  normal[0] = ((h1[1] - o[1]) * (h2[2] - o[2]) - (h1[2] - o[2]) * (h2[1] - o[1]));
  normal[1] = -((h1[0] - o[0]) * (h2[2] - o[2]) - (h1[2] - o[2]) * (h2[0] - o[0]));
  normal[2] = ((h1[0] - o[0]) * (h2[1] - o[1]) - (h1[1] - o[1]) * (h2[0] - o[0]));

  double length = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
  normal[0] = normal[0] / length;
  normal[1] = normal[1] / length;
  normal[2] = normal[2] / length;

  return normal;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> EXODUS::Mesh::node_vec(const int tail, const int head) const
{
  std::vector<double> nv(3);
  std::vector<double> t = get_node(tail);
  std::vector<double> h = get_node(head);
  nv[0] = h[0] - t[0];
  nv[1] = h[1] - t[1];
  nv[2] = h[2] - t[2];
  double length = sqrt(nv[0] * nv[0] + nv[1] * nv[1] + nv[2] * nv[2]);
  nv[0] = nv[0] / length;
  nv[1] = nv[1] / length;
  nv[2] = nv[2] / length;
  return nv;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
EXODUS::ElementBlock::ElementBlock(ElementBlock::Shape Distype,
    std::shared_ptr<std::map<int, std::vector<int>>>& eleconn, std::string name)
    : distype_(Distype), eleconn_(eleconn), name_(name.c_str())
{
  // do a sanity check
  for (std::map<int, std::vector<int>>::const_iterator elem = eleconn->begin();
      elem != eleconn->end(); ++elem)
  {
    if (Core::FE::get_number_of_element_nodes(pre_shape_to_drt(Distype)) !=
        (int)elem->second.size())
    {
      FOUR_C_THROW("number of read nodes does not fit the distype");
    }
  }
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<int> EXODUS::ElementBlock::get_ele_nodes(int i) const
{
  std::map<int, std::vector<int>>::const_iterator it = eleconn_->find(i);
  return it->second;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int EXODUS::ElementBlock::get_ele_node(int ele, int node) const
{
  std::map<int, std::vector<int>>::const_iterator it = eleconn_->find(ele);
  if (it == eleconn_->end()) FOUR_C_THROW("Element not found");
  std::vector<int> elenodes = get_ele_nodes(ele);
  return elenodes[node];
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::ElementBlock::print(std::ostream& os, bool verbose) const
{
  os << "Element Block, named: " << name_ << std::endl
     << "of Shape: " << shape_to_string(distype_) << std::endl
     << "has " << get_num_ele() << " Elements" << std::endl;
  if (verbose)
  {
    std::map<int, std::vector<int>>::const_iterator it;
    for (it = eleconn_->begin(); it != eleconn_->end(); it++)
    {
      os << "Ele " << it->first << ": ";
      const std::vector<int> myconn = it->second;  // GetEleNodes(int(it));
      for (int i = 0; i < signed(myconn.size()); i++)
      {
        os << myconn[i] << ",";
      }
      os << std::endl;
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
EXODUS::NodeSet::NodeSet(
    const std::set<int>& nodeids, const std::string& name, const std::string& propname)
    : nodeids_(nodeids), name_(name.c_str()), propname_(propname.c_str())
{
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::NodeSet::print(std::ostream& os, bool verbose) const
{
  os << "Node Set, named: " << name_ << std::endl
     << "Property Name: " << propname_ << std::endl
     << "has " << get_num_nodes() << " Nodes" << std::endl;
  if (verbose)
  {
    os << "Contains Nodes:" << std::endl;
    std::set<int>::iterator it;
    for (it = nodeids_.begin(); it != nodeids_.end(); it++) os << *it << ",";
    os << std::endl;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
EXODUS::SideSet::SideSet(const std::map<int, std::vector<int>>& sides, const std::string& name)
    : sides_(sides), name_(name.c_str())
{
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::SideSet::print(std::ostream& os, bool verbose) const
{
  os << "SideSet, named: " << name_ << std::endl
     << "has " << get_num_sides() << " Sides" << std::endl;
  if (verbose)
  {
    std::map<int, std::vector<int>>::const_iterator it;
    for (it = sides_.begin(); it != sides_.end(); it++)
    {
      os << "Side " << it->first << ": ";
      os << "Ele: " << it->second.at(0) << ", Side: " << it->second.at(1) << std::endl;
    }
  }
}


FOUR_C_NAMESPACE_CLOSE
