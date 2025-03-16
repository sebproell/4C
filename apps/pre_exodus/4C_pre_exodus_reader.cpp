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
  exoid_ = ex_open(exofilenamechar, EX_READ, &CPU_word_size, &IO_word_size, &exoversion);
  if (exoid_ <= 0) FOUR_C_THROW("Error while opening EXODUS II file {}", exofilenamechar);

  // print version
  std::cout << "File " << exofilename << " was created with EXODUS II library version "
            << exoversion << std::endl;

  // read database parameters
  int num_elem_blk, num_node_sets, num_side_sets, num_nodes;
  char title[MAX_LINE_LENGTH + 1];
  error = ex_get_init(exoid_, title, &four_c_dim_, &num_nodes, &num_elem_, &num_elem_blk,
      &num_node_sets, &num_side_sets);
  title_ = std::string(title);

  num_dim_ = 3;

  // get nodal coordinates
  {
    std::vector<double> x(num_nodes);
    std::vector<double> y(num_nodes);
    std::vector<double> z(num_nodes);
    error = ex_get_coord(exoid_, x.data(), y.data(), z.data());
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
    error = ex_get_ids(exoid_, EX_ELEM_BLOCK, ebids.data());
    if (error != 0) FOUR_C_THROW("exo error returned");
    error = ex_get_prop_array(exoid_, EX_ELEM_BLOCK, "ID", epropID.data());
    if (error != 0) FOUR_C_THROW("exo error returned");
    for (int i = 0; i < num_elem_blk; ++i)
    {
      // Read Element Blocks into Map
      char mychar[MAX_STR_LENGTH + 1];
      int num_el_in_blk, num_nod_per_elem, num_attr;
      // error = ex_get_elem_block (exoid_, epropID[i], mychar, &num_el_in_blk, &num_nod_per_elem,
      // &num_attr);
      error = ex_get_block(exoid_, EX_ELEM_BLOCK, ebids[i], mychar, &num_el_in_blk,
          &num_nod_per_elem, nullptr, nullptr, &num_attr);
      if (error != 0) FOUR_C_THROW("exo error returned");
      // prefer std::string to store element type
      std::string ele_type(mychar);

      // get ElementBlock name
      error = ex_get_name(exoid_, EX_ELEM_BLOCK, ebids[i], mychar);
      if (error != 0) FOUR_C_THROW("exo error returned");
      // prefer std::string to store name
      std::string blockname(mychar);

      // get element connectivity
      std::vector<int> allconn(num_nod_per_elem * num_el_in_blk);
      error = ex_get_conn(exoid_, EX_ELEM_BLOCK, ebids[i], allconn.data(), nullptr, nullptr);
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
    error = ex_get_prop_array(exoid_, EX_NODE_SET, "ID", npropID.data());
    for (int i = 0; i < num_node_sets; ++i)
    {
      // Read NodeSet params
      int num_nodes_in_set, num_df_in_set;
      error = ex_get_set_param(exoid_, EX_NODE_SET, npropID[i], &num_nodes_in_set, &num_df_in_set);

      // get NodeSet name
      char mychar[MAX_STR_LENGTH + 1];
      error = ex_get_name(exoid_, EX_NODE_SET, npropID[i], mychar);
      // prefer std::string to store name
      std::string nodesetname(mychar);

      // get nodes in node set
      std::vector<int> node_set_node_list(num_nodes_in_set);
      error = ex_get_set(exoid_, EX_NODE_SET, npropID[i], node_set_node_list.data(), nullptr);
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
    error = ex_inquire(exoid_, EX_INQ_NS_PROP, &num_props, &fdum, &cdum);
    // allocate memory for NodeSet property names
    char** prop_names = new char*[num_props];
    for (int i = 0; i < num_props; ++i)
    {
      prop_names[i] = new char[MAX_STR_LENGTH + 1];
    }
    // get prop names of node sets
    error = ex_get_prop_names(exoid_, EX_NODE_SET, prop_names);

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
    error = ex_get_prop_array(exoid_, EX_SIDE_SET, "ID", spropID.data());
    for (int i = 0; i < num_side_sets; ++i)
    {
      // get SideSet name
      char mychar[MAX_STR_LENGTH + 1];
      error = ex_get_name(exoid_, EX_SIDE_SET, spropID[i], mychar);
      // prefer std::string to store name
      std::string sidesetname(mychar);

      // Read SideSet params
      int num_side_in_set, num_dist_fact_in_set;
      error = ex_get_set_param(
          exoid_, EX_SIDE_SET, spropID[i], &num_side_in_set, &num_dist_fact_in_set);

      // get SideSet
      std::vector<int> side_set_elem_list(num_side_in_set);
      std::vector<int> side_set_side_list(num_side_in_set);
      error = ex_get_set(
          exoid_, EX_SIDE_SET, spropID[i], side_set_elem_list.data(), side_set_side_list.data());
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

  // close ExoFile
  close_exo();

  return;
}

EXODUS::Mesh::Mesh()
{
  nodes_ = std::make_shared<std::map<int, std::vector<double>>>();
  num_dim_ = 3;
  four_c_dim_ = 3;
  num_elem_ = 0;
  exoid_ = 0;
  title_ = "emptymesh";
}


/*----------------------------------------------------------------------*
 |  Extension constructor (public)                             maf 01/08|
 *----------------------------------------------------------------------*/
EXODUS::Mesh::Mesh(const EXODUS::Mesh& basemesh,
    const std::shared_ptr<std::map<int, std::vector<double>>> extNodes,
    const std::map<int, std::shared_ptr<ElementBlock>>& extBlocks,
    const std::map<int, NodeSet>& extNodesets, const std::map<int, SideSet>& extSidesets,
    const std::string newtitle)
    : title_(newtitle.c_str())
{
  // get all data from basemesh
  const int basedim = basemesh.get_num_dim();
  const int fourcdim = basemesh.get_four_c_dim();
  const int basenumele = basemesh.get_num_ele();
  std::shared_ptr<std::map<int, std::vector<double>>> baseNodes =
      std::make_shared<std::map<int, std::vector<double>>>();
  baseNodes = basemesh.get_nodes();
  // std::shared_ptr<std::map<int,std::vector<double> > > baseNodes = basemesh.GetNodes();
  std::map<int, std::shared_ptr<ElementBlock>> baseEblocks = basemesh.get_element_blocks();
  std::map<int, NodeSet> baseNodesets = basemesh.get_node_sets();
  std::map<int, SideSet> baseSidesets = basemesh.get_side_sets();

  //  // get infos from extension
  //  int extnumele = extBlocks.size();

  /********************* merge everything into new mesh ***********************/
  num_dim_ = basedim;
  four_c_dim_ = fourcdim;
  int total_num_elem = basenumele;
  //  num_elem_ = basenumele; // + extnumele;
  exoid_ = basemesh.get_exo_id();  // basefile still used for writing minor infos, e.g. qa record or
                                   // coordnames

  // merge nodes
  std::map<int, std::vector<double>>::const_iterator i_node;
  for (i_node = extNodes->begin(); i_node != extNodes->end(); ++i_node)
  {
    std::pair<std::map<int, std::vector<double>>::iterator, bool> check;
    check = baseNodes->insert(std::pair<int, std::vector<double>>(i_node->first, i_node->second));
    // happens when concatenating: if (check.second == false)  FOUR_C_THROW("Extension node already
    // exists!");
  }
  nodes_ = baseNodes;

  // merge ElementBlocks
  std::map<int, std::shared_ptr<ElementBlock>>::const_iterator i_block;
  for (i_block = baseEblocks.begin(); i_block != baseEblocks.end(); ++i_block)
  {
    element_blocks_.insert(
        std::pair<int, std::shared_ptr<ElementBlock>>(i_block->first, i_block->second));
  }
  for (i_block = extBlocks.begin(); i_block != extBlocks.end(); ++i_block)
  {
    std::pair<std::map<int, std::shared_ptr<ElementBlock>>::iterator, bool> check;
    check = element_blocks_.insert(
        std::pair<int, std::shared_ptr<ElementBlock>>(i_block->first, i_block->second));
    if (check.second == false)
      FOUR_C_THROW("Extension ElementBlock already exists!");
    else
      total_num_elem += i_block->second->get_num_ele();
  }
  num_elem_ = total_num_elem;

  // merge NodeSets
  std::map<int, NodeSet>::const_iterator i_ns;
  for (i_ns = baseNodesets.begin(); i_ns != baseNodesets.end(); ++i_ns)
  {
    node_sets_.insert(std::pair<int, NodeSet>(i_ns->first, i_ns->second));
  }
  for (i_ns = extNodesets.begin(); i_ns != extNodesets.end(); ++i_ns)
  {
    std::pair<std::map<int, NodeSet>::iterator, bool> check;
    check = node_sets_.insert(std::pair<int, NodeSet>(i_ns->first, i_ns->second));
    if (check.second == false) FOUR_C_THROW("Extension NodeSet already exists!");
  }

  // merge SideSets
  std::map<int, SideSet>::const_iterator i_ss;
  for (i_ss = baseSidesets.begin(); i_ss != baseSidesets.end(); ++i_ss)
  {
    side_sets_.insert(std::pair<int, SideSet>(i_ss->first, i_ss->second));
  }
  for (i_ss = extSidesets.begin(); i_ss != extSidesets.end(); ++i_ss)
  {
    std::pair<std::map<int, SideSet>::iterator, bool> check;
    check = side_sets_.insert(std::pair<int, SideSet>(i_ss->first, i_ss->second));
    if (check.second == false) FOUR_C_THROW("Extension SideSet already exists!");
  }
}


/*----------------------------------------------------------------------*
 |  Close corresponding Exofile(public)                        maf 12/07|
 *----------------------------------------------------------------------*/
void EXODUS::Mesh::close_exo() const
{
  // close exodus II file
  int exoid = get_exo_id();
  int error = ex_close(exoid);
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
std::string EXODUS::Mesh::get_title() const
{
  std::string title(title_, int(MAX_LINE_LENGTH + 1));
  return title;
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
std::map<int, std::vector<int>> EXODUS::Mesh::get_side_set_conn(const SideSet sideset) const
{
  std::cout << "Creating SideSet Connectivity... " << std::endl;
  fflush(stdout);

  Teuchos::Time time("", true);
  auto timetot = Teuchos::TimeMonitor::getNewTimer("Side Set Connect total");
  auto time1 = Teuchos::TimeMonitor::getNewTimer("One Side Set");
  auto time2 = Teuchos::TimeMonitor::getNewTimer("Get one Element Block");
  auto time3 = Teuchos::TimeMonitor::getNewTimer("Get Ele Conn");
  auto time4 = Teuchos::TimeMonitor::getNewTimer("Get one Ele");
  auto time5 = Teuchos::TimeMonitor::getNewTimer("Build one Side Conn");
  auto time6 = Teuchos::TimeMonitor::getNewTimer("Get all Eblocks and Econns");
  std::shared_ptr<Teuchos::TimeMonitor> tm_total = std::make_shared<Teuchos::TimeMonitor>(*timetot);

  std::map<int, std::vector<int>> conn;
  std::map<int, std::vector<int>> mysides = sideset.get_side_set();
  std::map<int, std::vector<int>>::iterator i_side;

  std::map<int, std::shared_ptr<EXODUS::ElementBlock>> ebs = get_element_blocks();
  std::map<int, std::shared_ptr<EXODUS::ElementBlock>>::const_iterator i_ebs;

  // Range Vector for global eleID identification in SideSet
  std::vector<int> glob_eb_erange(1, 0);
  int rangebreak = 0;
  // Also we once get all EBlocks and EConns to enable quick access
  std::vector<EXODUS::ElementBlock> eblocks;
  std::vector<std::map<int, std::vector<int>>> econns;
  std::shared_ptr<Teuchos::TimeMonitor> tm6 = std::make_shared<Teuchos::TimeMonitor>(*time6);
  for (i_ebs = ebs.begin(); i_ebs != ebs.end(); ++i_ebs)
  {
    rangebreak += i_ebs->second->get_num_ele();
    glob_eb_erange.push_back(rangebreak);
    eblocks.push_back(*i_ebs->second);
    econns.push_back(*(i_ebs->second->get_ele_conn()));
  }
  tm6 = nullptr;

  // fill SideSet Connectivity
  for (i_side = mysides.begin(); i_side != mysides.end(); ++i_side)
  {
    Teuchos::TimeMonitor tm1(*time1);
    int actele = i_side->second.at(0) -
                 1;  // ExoIds start from 1, but we from 0 //ToDo: remove -1 idconfusion
    int actface = i_side->second.at(1) - 1;  // ExoIds start from 1, but we from 0
    // find actual EBlock where actele lies in
    int actebid = -1;
    for (unsigned int i = 0; i < glob_eb_erange.size(); ++i)
      if (actele < glob_eb_erange[i])
      {
        actebid = i - 1;
        break;
      }
    std::shared_ptr<Teuchos::TimeMonitor> tm2 = std::make_shared<Teuchos::TimeMonitor>(*time2);
    // EXODUS::ElementBlock acteb = ebs.find(actebid)->second;
    tm2 = nullptr;
    // EXODUS::ElementBlock::Shape actshape = acteb.GetShape();
    if (actebid < 0) FOUR_C_THROW("invalid element block id");
    EXODUS::ElementBlock::Shape actshape = eblocks[actebid].get_shape();
    std::shared_ptr<Teuchos::TimeMonitor> tm3 = std::make_shared<Teuchos::TimeMonitor>(*time3);
    // map<int,std::vector<int> > acteconn = acteb.GetEleConn();
    tm3 = nullptr;
    // get act parent ele from actual Side
    int parent_ele_id = actele - glob_eb_erange[actebid];
    std::shared_ptr<Teuchos::TimeMonitor> tm4 = std::make_shared<Teuchos::TimeMonitor>(*time4);
    // vector<int> parent_ele = acteconn.find(parent_ele_id)->second;
    std::vector<int> parent_ele = econns[actebid].find(parent_ele_id)->second;
    tm4 = nullptr;
    // Face to ElementNode Map
    //// **** temporary hex map due to conflicts between side numbering exo<->4C
    std::shared_ptr<Teuchos::TimeMonitor> tm5 = std::make_shared<Teuchos::TimeMonitor>(*time5);
    switch (actshape)
    {
      case ElementBlock::tet4:
      {
        break;
      }
      case ElementBlock::hex8:
      {
        actface = hex_side_number_exo_to_four_c(actface);
        break;
      }
      case ElementBlock::pyramid5:
      {
        actface = pyr_side_number_exo_to_four_c(actface);
        break;
      }
      case ElementBlock::wedge6:
      {
        break;
      }
      default:
      {
        std::cout << shape_to_string(actshape) << ":" << std::endl;
        FOUR_C_THROW("Parent Element Type not supported");
        break;
      }
    }
    std::vector<int> childmap =
        Core::FE::get_ele_node_numbering_surfaces(pre_shape_to_drt(actshape))[actface];
    // child gets its node ids
    std::vector<int> child;
    for (unsigned int j = 0; j < childmap.size(); ++j) child.push_back(parent_ele[childmap[j]]);
    //    PrintVec(std::cout,childmap);
    //    PrintVec(std::cout,child);
    // some checking
    if ((child.size() != 3) && (child.size() != 4))
    {
      print_vec(std::cout, child);
      print_vec(std::cout, childmap);
      print_vec(std::cout, parent_ele);
      std::cout << shape_to_string(actshape) << ",Face: " << actface
                << ",childsize:" << child.size() << std::endl;
      FOUR_C_THROW("Child Ele error");
    }
    // insert child into SideSet Connectivity
    conn.insert(std::pair<int, std::vector<int>>(i_side->first, child));
    tm5 = nullptr;

    //    // progress output
    //    tm1 = nullptr;
    //    if (signed(i_side->first) == perc * signed(mysides.size())/100){
    //      Teuchos::TimeMonitor::summarize();
    //      std::cout << perc << " % of " << mysides.size() << " Sides done" << std::endl;
    //      perc ++;
    //      fflush(stdout);
    //    }
  }
  tm_total = nullptr;
  //  Teuchos::TimeMonitor::summarize();
  std::cout << "...done" << std::endl;
  fflush(stdout);

  return conn;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<int, std::vector<int>> EXODUS::Mesh::get_side_set_conn(
    const SideSet sideset, bool checkoutside) const
{
  std::cout << "Creating SideSet Connectivity with outside-check... " << std::endl;
  fflush(stdout);

  Teuchos::Time time("", true);
  auto timetot = Teuchos::TimeMonitor::getNewTimer("Side Set Connect total");
  auto time1 = Teuchos::TimeMonitor::getNewTimer("One Side Set");
  std::shared_ptr<Teuchos::TimeMonitor> tm_total = std::make_shared<Teuchos::TimeMonitor>(*timetot);

  std::map<int, std::vector<int>> conn;
  std::map<int, std::vector<int>> mysides = sideset.get_side_set();
  std::map<int, std::vector<int>>::iterator i_side;

  std::map<int, std::shared_ptr<EXODUS::ElementBlock>> ebs = get_element_blocks();
  std::map<int, std::shared_ptr<EXODUS::ElementBlock>>::const_iterator i_ebs;

  // Range Vector for global eleID identification in SideSet
  std::vector<int> glob_eb_erange(1, 0);
  int rangebreak = 0;
  // Also we once get all EBlocks and EConns to enable quick access
  std::vector<EXODUS::ElementBlock> eblocks;
  std::vector<std::map<int, std::vector<int>>> econns;
  for (i_ebs = ebs.begin(); i_ebs != ebs.end(); ++i_ebs)
  {
    rangebreak += i_ebs->second->get_num_ele();
    glob_eb_erange.push_back(rangebreak);
    eblocks.push_back(*i_ebs->second);
    econns.push_back(*(i_ebs->second->get_ele_conn()));
  }

  // fill SideSet Connectivity
  for (i_side = mysides.begin(); i_side != mysides.end(); ++i_side)
  {
    std::shared_ptr<Teuchos::TimeMonitor> tm1 = std::make_shared<Teuchos::TimeMonitor>(*time1);
    int actele = i_side->second.at(0) -
                 1;  // ExoIds start from 1, but we from 0 //ToDo: remove -1 idconfusion
    int actface = i_side->second.at(1) - 1;  // ExoIds start from 1, but we from 0
    // find actual EBlock where actele lies in
    int actebid = -1;
    for (unsigned int i = 0; i < glob_eb_erange.size(); ++i)
      if (actele < glob_eb_erange[i])
      {
        actebid = i - 1;
        break;
      }
    if (actebid < 0) FOUR_C_THROW("invalid element block id");
    EXODUS::ElementBlock::Shape actshape = eblocks[actebid].get_shape();

    // get act parent ele from actual Side
    int parent_ele_id = actele - glob_eb_erange[actebid];
    std::vector<int> parent_ele = econns[actebid].find(parent_ele_id)->second;

    // Face to ElementNode Map
    //// **** temporary hex map due to conflicts between side numbering exo<->4C
    switch (actshape)
    {
      case ElementBlock::tet4:
      {
        break;
      }
      case ElementBlock::hex8:
      {
        actface = hex_side_number_exo_to_four_c(actface);
        break;
      }
      case ElementBlock::pyramid5:
      {
        actface = pyr_side_number_exo_to_four_c(actface);
        break;
      }
      case ElementBlock::wedge6:
      {
        break;
      }
      default:
      {
        std::cout << shape_to_string(actshape) << ":" << std::endl;
        FOUR_C_THROW("Parent Element Type not supported");
        break;
      }
    }
    std::vector<int> childmap =
        Core::FE::get_ele_node_numbering_surfaces(pre_shape_to_drt(actshape))[actface];

    std::vector<int> child;
    if (checkoutside)
      child = outside_oriented_side(parent_ele, childmap);
    else
      for (unsigned int j = 0; j < childmap.size(); ++j) child.push_back(parent_ele[childmap[j]]);

    // some checking
    if ((child.size() != 3) && (child.size() != 4))
    {
      print_vec(std::cout, child);
      print_vec(std::cout, childmap);
      print_vec(std::cout, parent_ele);
      std::cout << shape_to_string(actshape) << ",Face: " << actface
                << ",childsize:" << child.size() << std::endl;
      FOUR_C_THROW("Child Ele error");
    }

    // insert child into SideSet Connectivity
    conn.insert(std::pair<int, std::vector<int>>(i_side->first, child));
    tm1 = nullptr;
  }
  tm_total = nullptr;
  std::cout << "...done" << std::endl;
  fflush(stdout);

  return conn;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<int> EXODUS::Mesh::outside_oriented_side(
    const std::vector<int> parentele, const std::vector<int> sidemap) const
{
  // first guess of child
  std::vector<int> child;
  // set simplifies later inverse search
  std::set<int> childset;
  for (unsigned int j = 0; j < sidemap.size(); ++j)
  {
    child.push_back(parentele[sidemap[j]]);
    childset.insert(parentele[sidemap[j]]);
  }

  // set of parentele for later inverse search
  std::set<int> parentset;
  for (unsigned int i = 0; i < parentele.size(); ++i) parentset.insert(parentele.at(i));

  // find parentele node not within side
  int insidenode = -1;
  std::set<int>::iterator it;
  for (it = parentset.begin(); it != parentset.end(); ++it)
  {
    if (childset.find(*it) == childset.end())
    {
      insidenode = *it;
      break;
    }
  }

  // build normal at first side node
  std::vector<double> sidenormal = normal(child.back(), child.front(), child.at(1));
  // build vector from first side node to inside element node
  std::vector<double> insidevec = node_vec(child.front(), insidenode);

  // scalar-product
  double scp =
      sidenormal[0] * insidevec[0] + sidenormal[1] * insidevec[1] + sidenormal[2] * insidevec[2];

  std::vector<int> out_side;
  if (scp < 0)
  {
    std::vector<int> reversechild;
    std::vector<int>::reverse_iterator rit;
    for (rit = child.rbegin(); rit < child.rend(); ++rit) out_side.push_back(*rit);
  }
  else
    out_side = child;

  return out_side;
}

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
std::vector<EXODUS::ElementBlock> EXODUS::Mesh::side_set_to_e_blocks(
    const EXODUS::SideSet& sideset, const std::map<int, std::vector<int>>& sidesetconn) const
{
  std::vector<ElementBlock> eblocks;
  std::map<int, std::vector<int>>::const_iterator i_ele;
  std::shared_ptr<std::map<int, std::vector<int>>> quadconn =
      std::make_shared<std::map<int, std::vector<int>>>();
  int quadcounter = 0;
  std::shared_ptr<std::map<int, std::vector<int>>> triconn =
      std::make_shared<std::map<int, std::vector<int>>>();
  int tricounter = 0;
  for (i_ele = sidesetconn.begin(); i_ele != sidesetconn.end(); ++i_ele)
  {
    int numnodes = i_ele->second.size();
    if (numnodes == 4)
    {
      quadconn->insert(std::pair<int, std::vector<int>>(quadcounter, i_ele->second));
      quadcounter++;
    }
    else if (numnodes == 3)
    {
      triconn->insert(std::pair<int, std::vector<int>>(tricounter, i_ele->second));
      tricounter++;
    }
    else
      FOUR_C_THROW("Number of basenodes for conversion from SideSet to EBlock not supported");
  }
  if (quadcounter > 0)
  {
    std::ostringstream quadblockname;
    quadblockname << sideset.get_name() << "quad";
    EXODUS::ElementBlock neweblock(ElementBlock::quad4, quadconn, quadblockname.str());
    eblocks.push_back(neweblock);
  }
  if (tricounter > 0)
  {
    std::ostringstream triblockname;
    triblockname << sideset.get_name() << "tri";
    EXODUS::ElementBlock neweblock(ElementBlock::tri3, triconn, triblockname.str());
    eblocks.push_back(neweblock);
  }

  return eblocks;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
EXODUS::NodeSet EXODUS::Mesh::side_set_to_node_set(
    const EXODUS::SideSet& sideset, const std::map<int, std::vector<int>>& sidesetconn) const
{
  std::map<int, std::vector<int>>::const_iterator i_side;
  std::vector<int>::const_iterator i_node;
  std::set<int> nodes;
  for (i_side = sidesetconn.begin(); i_side != sidesetconn.end(); ++i_side)
    for (i_node = i_side->second.begin(); i_node != i_side->second.end(); ++i_node)
      nodes.insert(*i_node);  // nodes.insert(i_side->second.at(i_node));
  std::ostringstream nodesetname;
  nodesetname << "nodes";  // sideset.GetName() << "nodes";
  std::string propname = "";
  EXODUS::NodeSet nodeset(nodes, nodesetname.str(), propname);

  return nodeset;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::set<int> EXODUS::Mesh::get_side_set_nodes(
    const EXODUS::SideSet& sideset, const std::map<int, std::vector<int>>& sidesetconn) const
{
  std::map<int, std::vector<int>>::const_iterator i_side;
  std::vector<int>::const_iterator i_node;
  std::set<int> nodes;
  for (i_side = sidesetconn.begin(); i_side != sidesetconn.end(); ++i_side)
    for (i_node = i_side->second.begin(); i_node != i_side->second.end(); ++i_node)
      nodes.insert(*i_node);  // nodes.insert(i_side->second.at(i_node));
  return nodes;
}

/*----------------------------------------------------------------------*
 |  Write Mesh into exodus file (public)                       maf 01/08|
 *----------------------------------------------------------------------*/
void EXODUS::Mesh::write_mesh(const std::string newexofilename) const
{
  std::cout << "Writing Mesh into: " << newexofilename << std::endl;
  // std::string newexofile(newexofilename);
  // newexofile += ".exo";
  const char* newexofilechar = newexofilename.c_str();

  int CPU_word_size, IO_word_size, exoid, error;
  CPU_word_size = sizeof(double); /* use double*/
  IO_word_size = 8;               /* store variables as doubles */

  /* create EXODUS II file */
  exoid = ex_create(newexofilechar, /* filename path */
      EX_CLOBBER,                   /* create mode */
      &CPU_word_size,               /* CPU double word size in bytes */
      &IO_word_size);               /* I/O double word size in bytes */

  int num_elem_blk = get_num_element_blocks();
  int num_node_sets = get_num_node_sets();
  int num_side_sets = get_num_side_sets();
  int num_nodes = get_num_nodes();
  /* initialize file with parameters */
  const char* title = title_.c_str();
  error = ex_put_init(
      exoid, title, num_dim_, num_nodes, num_elem_, num_elem_blk, num_node_sets, num_side_sets);
  if (error != 0) FOUR_C_THROW("error in exfile init");

  /* Write QA record based on original exofile */
  int num_qa_rec;
  char* qa_record[MAX_STR_LENGTH][4];  // should be MAX_QA_REC][4], but this is nowhere defined!;
  char cdum;                           // dummy variable
  float fdum;
  /* read QA records */
  ex_inquire(exoid_, EX_INQ_QA, &num_qa_rec, &fdum, &cdum); /* write QA records */
  for (int i = 0; i < num_qa_rec; i++)
    for (int j = 0; j < 4; j++) qa_record[i][j] = (char*)calloc((MAX_STR_LENGTH + 1), sizeof(char));
  // for (int j=0; j<4; j++) qa_record[i][j] = new (char)[MAX_STR_LENGTH+1];
  error = ex_get_qa(exoid_, qa_record);
  error = ex_put_qa(exoid, num_qa_rec, qa_record);
  for (int i = 0; i < num_qa_rec; i++)
    for (int j = 0; j < 4; j++)
      // delete [] qa_record[i][j];
      free(qa_record[i][j]);

  // Write coord names based on original exofile
  char* coord_names[3];
  for (int i = 0; i < num_dim_; i++)
    coord_names[i] = (char*)calloc((MAX_STR_LENGTH + 1), sizeof(char));
  error = ex_get_coord_names(exoid_, coord_names);
  error = ex_put_coord_names(exoid, coord_names);
  for (int i = 0; i < num_dim_; i++) free(coord_names[i]);

  // Write nodal coordinates
  std::vector<double> xc(num_nodes);
  std::vector<double> yc(num_nodes);
  std::vector<double> zc(num_nodes);
  std::map<int, std::vector<double>>::const_iterator it;
  std::shared_ptr<std::map<int, std::vector<double>>> nodes = get_nodes();
  for (it = nodes->begin(); it != nodes->end(); ++it)
  {
    xc[it->first - 1] = it->second[0];  // vector starts with 0
    yc[it->first - 1] = it->second[1];  // vector starts with 0
    zc[it->first - 1] = it->second[2];  // vector starts with 0
  }
  error = ex_put_coord(exoid, xc.data(), yc.data(), zc.data());

  // Write NodeSets ************************************************************
  std::map<int, NodeSet>::const_iterator ins;
  const std::map<int, NodeSet> nss = get_node_sets();
  for (ins = nss.begin(); ins != nss.end(); ++ins)
  {
    const int nsID = ins->first;
    const NodeSet ns = ins->second;
    const int num_nodes_in_set = ns.get_num_nodes();
    if (num_nodes_in_set > 0)  // do not bother if nodeset empty
    {
      const std::string name = ns.get_name();
      const char* nsname = name.c_str();
      const std::string propname = ns.get_prop_name();
      const std::set<int> nodes = ns.get_node_set();
      error = ex_put_set_param(exoid,  // of write file
          EX_NODE_SET,
          nsID,  // node set id
          num_nodes_in_set,
          0);  // yet no distribution factors
      if (error != 0) FOUR_C_THROW("error writing node set params");
      std::vector<int> nodelist(num_nodes_in_set);
      ns.fill_nodelist_array(nodelist.data());
      error = ex_put_set(exoid, EX_NODE_SET, nsID, nodelist.data(), nullptr);
      if (error != 0) FOUR_C_THROW("error writing node set \"{}\" ", nsname);
      error = ex_put_name(exoid, EX_NODE_SET, nsID, nsname);
      if (error != 0) FOUR_C_THROW("error writing node set name");
    }
  }

  // Write ElementBlocks  ******************************************************
  std::map<int, std::shared_ptr<ElementBlock>>::const_iterator iebs;
  const std::map<int, std::shared_ptr<ElementBlock>> ebs = get_element_blocks();
  for (iebs = ebs.begin(); iebs != ebs.end(); iebs++)
  {
    const int blockID = iebs->first;
    const ElementBlock eb = (*iebs->second);
    const ElementBlock::Shape shape = eb.get_shape();
    const std::string shapestring = shape_to_string(shape);
    const std::vector<int> exampleconn = eb.get_ele_nodes(0);  // iebs->first);
    const int num_nod_per_elem = exampleconn.size();
    const int numele = eb.get_num_ele();
    const char* elem_type = shapestring.c_str();
    error = ex_put_block(exoid,  // of write file
        EX_ELEM_BLOCK,
        blockID,           // element block id
        elem_type,         // its name
        numele,            // number of element in block
        num_nod_per_elem,  // num of nodes per ele
        0, 0,
        1);  // num of attributes, not supported yet ->1
    if (error != 0) FOUR_C_THROW("error writing element block");
    // Write Element Connectivity
    std::vector<int> conn(num_nod_per_elem * numele);
    eb.fill_econn_array(conn.data());
    error = ex_put_conn(exoid, EX_ELEM_BLOCK, blockID, conn.data(), nullptr, nullptr);
    if (error != 0) FOUR_C_THROW("error writing element block conns");
    // write block name
    const std::string bname = eb.get_name();
    const char* blockname = bname.c_str();
    error = ex_put_name(exoid, EX_ELEM_BLOCK, blockID, blockname);
    if (error != 0) FOUR_C_THROW("error writing element block name");
  }

  // Write SideSets ************************************************************
  std::map<int, SideSet>::const_iterator iss;
  const std::map<int, SideSet> sss = get_side_sets();
  for (iss = sss.begin(); iss != sss.end(); ++iss)
  {
    const int ssID = iss->first;
    const SideSet ss = iss->second;
    const int num_side_in_set = ss.get_num_sides();
    const std::string name = ss.get_name();
    const char* ssname = name.c_str();
    error = ex_put_set_param(exoid,  // of write file
        EX_SIDE_SET,
        ssID,  // side set id
        num_side_in_set,
        0);  // yet no distribution factors
    if (error != 0) FOUR_C_THROW("error writing side set params");
    std::vector<int> side_set_elem_list(num_side_in_set);
    std::vector<int> side_set_side_list(num_side_in_set);
    // in case the sideset is newly created we have to adjust element ids to global numbering
    std::map<int, std::vector<int>> globalsides;
    if (iss->second.get_first_side_set().size() == 3)
    {
      globalsides = globalify_s_seleids(ssID);
      ss.fill_side_lists(side_set_elem_list.data(), side_set_side_list.data(), globalsides);
    }
    else
      ss.fill_side_lists(side_set_elem_list.data(), side_set_side_list.data());

    error =
        ex_put_set(exoid, EX_SIDE_SET, ssID, side_set_elem_list.data(), side_set_side_list.data());
    if (error != 0) FOUR_C_THROW("error writing side set");
    error = ex_put_name(exoid, EX_SIDE_SET, ssID, ssname);
    if (error != 0) FOUR_C_THROW("error writing sideset name");
  }

  // ***************************************************************************

  // close file
  error = ex_close(exoid);
  if (error != 0) FOUR_C_THROW("error closing exodus file");
  std::cout << ".. finished" << std::endl;
}

/*----------------------------------------------------------------------*
 |  Add Element Block to mesh(public)                          maf 01/08|
 *----------------------------------------------------------------------*/
void EXODUS::Mesh::add_element_block(const std::shared_ptr<EXODUS::ElementBlock> eblock) const
{
  std::map<int, std::shared_ptr<ElementBlock>> eblocks = get_element_blocks();
  eblocks.insert(
      std::pair<int, std::shared_ptr<ElementBlock>>(get_num_element_blocks() + 1, eblock));
}

/*----------------------------------------------------------------------*
 |  Erase Element Block from mesh(public)                      maf 07/08|
 *----------------------------------------------------------------------*/
void EXODUS::Mesh::erase_element_block(const int id)
{
  int red_numele = get_element_block(id)->get_num_ele();
  element_blocks_.erase(id);
  num_elem_ = num_elem_ - red_numele;
}

/*----------------------------------------------------------------------*
 |  Erase SideSet from mesh(public)                            maf 07/08|
 *----------------------------------------------------------------------*/
void EXODUS::Mesh::erase_side_set(const int id) { side_sets_.erase(id); }

/*------------------------------------------------------------------------*
 | - calculates the midpoint of each element                               |
 | - returns map <midpoint-ID,std::pair<eblock-ID,element-ID> >         SP 06/08|
 *------------------------------------------------------------------------*/
std::map<int, std::pair<int, int>> EXODUS::Mesh::create_midpoints(
    std::map<int, std::vector<double>>& midpoints, const std::vector<int>& eb_ids) const
{
  // map that will be returned
  std::map<int, std::pair<int, int>> conn_mpID_elID;

  // auxiliary variables
  int counter_elements = 0;
  int nodes_per_element = 0;

  std::vector<double> sumVector(3, 0);
  std::vector<double> midPoint(3, 0);

  std::map<int, std::shared_ptr<ElementBlock>> EBlocks;
  std::map<int, std::shared_ptr<ElementBlock>>::const_iterator it;
  // work only on eblocks in eb_ids
  std::vector<int>::const_iterator id;
  std::map<int, std::shared_ptr<EXODUS::ElementBlock>> ebs;
  for (id = eb_ids.begin(); id != eb_ids.end(); ++id)
  {
    std::shared_ptr<EXODUS::ElementBlock> acteb = this->get_element_block(*id);
    EBlocks.insert(std::pair<int, std::shared_ptr<EXODUS::ElementBlock>>(*id, acteb));
  }

  std::map<int, std::vector<int>> EleConn;
  std::map<int, std::vector<int>>::const_iterator it_2;

  std::vector<int>::const_iterator it_3;

  // loop over ElementBlocks
  for (it = EBlocks.begin(); it != EBlocks.end(); ++it)
  {
    EleConn = *(it->second->get_ele_conn());

    // loop over element connectivity
    for (it_2 = EleConn.begin(); it_2 != EleConn.end(); ++it_2)
    {
      counter_elements++;
      nodes_per_element = 0;

      sumVector[0] = 0;
      sumVector[1] = 0;
      sumVector[2] = 0;

      // loop over each Node in one element connectivity
      for (it_3 = it_2->second.begin(); it_3 != it_2->second.end(); ++it_3)
      {
        nodes_per_element++;

        // sum of two vectors
        sumVector[0] += get_node(*it_3)[0];
        sumVector[1] += get_node(*it_3)[1];
        sumVector[2] += get_node(*it_3)[2];
      }

      // midpoint of element i
      midPoint[0] = sumVector[0] / nodes_per_element;
      midPoint[1] = sumVector[1] / nodes_per_element;
      midPoint[2] = sumVector[2] / nodes_per_element;

      // insert calculated midpoint in midpoints_
      midpoints.insert(std::pair<int, std::vector<double>>(counter_elements, midPoint));
      // conn_mpID_elID = (midpoint-ID, eblock-ID, element-ID)
      std::pair<int, int> eb_e = std::make_pair(it->first, it_2->first);
      conn_mpID_elID.insert(std::pair<int, std::pair<int, int>>(counter_elements, eb_e));
    }
  }
  // EXODUS::PrintMap(std::cout,conn_mpID_elID);
  return conn_mpID_elID;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<int, std::vector<int>> EXODUS::Mesh::globalify_s_seleids(const int ssid) const
{
  SideSet ss = get_side_set(ssid);

  std::map<int, std::shared_ptr<EXODUS::ElementBlock>> ebs = get_element_blocks();
  std::map<int, std::shared_ptr<EXODUS::ElementBlock>>::const_iterator i_ebs;

  // Range Vector for global eleID identification in SideSet
  std::vector<int> glob_eb_erange(1, 0);
  int rangebreak = 0;

  std::map<int, int> ebid_rangepos;
  int rangepos = 0;
  for (i_ebs = ebs.begin(); i_ebs != ebs.end(); ++i_ebs)
  {
    rangebreak += i_ebs->second->get_num_ele();
    glob_eb_erange.push_back(rangebreak);
    ebid_rangepos.insert(std::pair<int, int>(i_ebs->first, rangepos));
    ++rangepos;
  }

  std::map<int, std::vector<int>> sideset = ss.get_side_set();
  std::map<int, std::vector<int>>::iterator i_ss;

  for (i_ss = sideset.begin(); i_ss != sideset.end(); ++i_ss)
  {
    if (i_ss->second.size() != 3) FOUR_C_THROW("Problem in new SideSet!");  // double check
    int ebid = i_ss->second.at(2);
    int lowerbound = glob_eb_erange.at(ebid_rangepos.find(ebid)->second);
    i_ss->second.at(0) = i_ss->second.at(0) + lowerbound + 1;
  }

  return sideset;
}


/*------------------------------------------------------------------------*
 |creates gmsh-file to visualize mesh                             MF 07/08|
 *------------------------------------------------------------------------*/
void EXODUS::Mesh::plot_element_blocks_gmsh(
    const std::string fname, const EXODUS::Mesh& mymesh) const
{
  std::shared_ptr<std::map<int, std::vector<double>>> nodes = mymesh.get_nodes();
  std::ofstream f_system(fname.c_str());
  std::stringstream gmshfilecontent;
  gmshfilecontent << "View \" Mesh \" {" << std::endl;

  std::map<int, std::shared_ptr<EXODUS::ElementBlock>> ebs = mymesh.get_element_blocks();
  std::map<int, std::shared_ptr<EXODUS::ElementBlock>>::const_iterator eb_it;
  std::map<int, std::vector<int>> conn;
  std::map<int, std::vector<int>>::const_iterator it;

  for (eb_it = ebs.begin(); eb_it != ebs.end(); ++eb_it)
  {
    std::shared_ptr<std::map<int, std::vector<int>>> actconn = eb_it->second->get_ele_conn();
    for (it = actconn->begin(); it != actconn->end(); ++it)
    {
      int eleid = it->first;
      const std::vector<int> elenodes = it->second;
      int numnodes = elenodes.size();
      if (numnodes == 6)
        gmshfilecontent << "SI(";
      else if (numnodes == 8)
        gmshfilecontent << "SH(";
      for (unsigned int i = 0; i < elenodes.size(); ++i)
      {
        gmshfilecontent << nodes->find(elenodes.at(i))->second[0] << ",";
        gmshfilecontent << nodes->find(elenodes.at(i))->second[1] << ",";
        gmshfilecontent << nodes->find(elenodes.at(i))->second[2];
        if (i == (elenodes.size() - 1))
          gmshfilecontent << ")";
        else
          gmshfilecontent << ",";
      }
      gmshfilecontent << "{";
      for (unsigned int i = 0; i < (elenodes.size() - 1); ++i) gmshfilecontent << eleid << ",";
      gmshfilecontent << eleid << "};" << std::endl;
    }
  }
  gmshfilecontent << "};" << std::endl;
  f_system << gmshfilecontent.str();
  f_system.close();
  return;
}

/*------------------------------------------------------------------------*
 |creates gmsh-file to visualize mesh                             MF 07/08|
 *------------------------------------------------------------------------*/
void EXODUS::Mesh::plot_element_blocks_gmsh(
    const std::string fname, const EXODUS::Mesh& mymesh, const std::vector<int>& ebids) const
{
  std::shared_ptr<std::map<int, std::vector<double>>> nodes = mymesh.get_nodes();
  std::ofstream f_system(fname.c_str());
  std::stringstream gmshfilecontent;
  gmshfilecontent << "View \" Mesh \" {" << std::endl;

  std::vector<int>::const_iterator id;
  std::map<int, std::shared_ptr<EXODUS::ElementBlock>> ebs;
  for (id = ebids.begin(); id != ebids.end(); ++id)
  {
    std::shared_ptr<EXODUS::ElementBlock> acteb = mymesh.get_element_block(*id);
    ebs.insert(std::pair<int, std::shared_ptr<EXODUS::ElementBlock>>(*id, acteb));
  }

  std::map<int, std::shared_ptr<EXODUS::ElementBlock>>::const_iterator eb_it;
  std::map<int, std::vector<int>> conn;
  std::map<int, std::vector<int>>::const_iterator it;

  for (eb_it = ebs.begin(); eb_it != ebs.end(); ++eb_it)
  {
    std::shared_ptr<std::map<int, std::vector<int>>> actconn = eb_it->second->get_ele_conn();
    for (it = actconn->begin(); it != actconn->end(); ++it)
    {
      int eleid = it->first;
      const std::vector<int> elenodes = it->second;
      int numnodes = elenodes.size();
      if (numnodes == 6)
        gmshfilecontent << "SI(";
      else if (numnodes == 8)
        gmshfilecontent << "SH(";
      for (unsigned int i = 0; i < elenodes.size(); ++i)
      {
        gmshfilecontent << nodes->find(elenodes.at(i))->second[0] << ",";
        gmshfilecontent << nodes->find(elenodes.at(i))->second[1] << ",";
        gmshfilecontent << nodes->find(elenodes.at(i))->second[2];
        if (i == (elenodes.size() - 1))
          gmshfilecontent << ")";
        else
          gmshfilecontent << ",";
      }
      gmshfilecontent << "{";
      for (unsigned int i = 0; i < (elenodes.size() - 1); ++i) gmshfilecontent << eleid << ",";
      gmshfilecontent << eleid << "};" << std::endl;
    }
  }
  gmshfilecontent << "};" << std::endl;
  f_system << gmshfilecontent.str();
  f_system.close();
  return;
}

/*------------------------------------------------------------------------*
 |creates gmsh-file to visualize all nodes                        SP 06/08|
 *------------------------------------------------------------------------*/
void EXODUS::Mesh::plot_nodes_gmsh() const
{
  std::ofstream f_system("mesh_all_nodes.gmsh");
  std::stringstream gmshfilecontent;
  gmshfilecontent << "View \" Nodes \" {" << std::endl;

  std::map<int, std::vector<double>>::const_iterator it;
  // loop over all nodes
  for (it = nodes_->begin(); it != nodes_->end(); it++)
  {
    const std::vector<double> mycoords = it->second;

    // writing of coordinates of each node
    gmshfilecontent << "SP(";
    gmshfilecontent << mycoords[0] << ",";
    gmshfilecontent << mycoords[1] << ",";
    gmshfilecontent << mycoords[2];
    gmshfilecontent << ")";

    // writing of node-ID
    gmshfilecontent << "{";
    gmshfilecontent << it->first << "," << it->first << "," << it->first << "};" << std::endl;
  }
  gmshfilecontent << "};" << std::endl;
  f_system << gmshfilecontent.str();
  f_system.close();
}

/*------------------------------------------------------------------------*
 |creates gmsh-file to visualize connectivity                     MF 12/08|
 *------------------------------------------------------------------------*/
void EXODUS::Mesh::plot_conn_gmsh(const std::string fname, const EXODUS::Mesh& mymesh,
    const std::map<int, std::vector<int>>& conn) const
{
  std::shared_ptr<std::map<int, std::vector<double>>> nodes = mymesh.get_nodes();
  std::ofstream f_system(fname.c_str());
  std::stringstream gmshfilecontent;
  gmshfilecontent << "View \" Connectivity \" {" << std::endl;

  std::map<int, std::vector<int>>::const_iterator it;

  for (it = conn.begin(); it != conn.end(); ++it)
  {
    int eleid = it->first;
    const std::vector<int> elenodes = it->second;
    int numnodes = elenodes.size();
    if (numnodes == 6)
      gmshfilecontent << "SI(";
    else if (numnodes == 8)
      gmshfilecontent << "SH(";
    else if (numnodes == 3)
      gmshfilecontent << "ST(";
    else if (numnodes == 4)
      gmshfilecontent << "SQ(";
    for (unsigned int i = 0; i < elenodes.size(); ++i)
    {
      gmshfilecontent << nodes->find(elenodes.at(i))->second[0] << ",";
      gmshfilecontent << nodes->find(elenodes.at(i))->second[1] << ",";
      gmshfilecontent << nodes->find(elenodes.at(i))->second[2];
      if (i == (elenodes.size() - 1))
        gmshfilecontent << ")";
      else
        gmshfilecontent << ",";
    }
    gmshfilecontent << "{";
    for (unsigned int i = 0; i < (elenodes.size() - 1); ++i) gmshfilecontent << eleid << ",";
    gmshfilecontent << eleid << "};" << std::endl;
  }
  gmshfilecontent << "};" << std::endl;
  f_system << gmshfilecontent.str();
  f_system.close();
  return;
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
void EXODUS::ElementBlock::fill_econn_array(int* connarray) const
{
  const std::map<int, std::vector<int>>::const_iterator iele;
  int numele = eleconn_->size();
  for (int i = 0; i < numele; ++i)
  {
    std::vector<int> ele = get_ele_nodes(i);
    int num_nod_per_elem = ele.size();
    for (int j = 0; j < num_nod_per_elem; ++j)
    {
      connarray[i * num_nod_per_elem + j] = get_ele_node(i, j);
    }
  }
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
void EXODUS::NodeSet::fill_nodelist_array(int* nodelist) const
{
  std::set<int> nlist = get_node_set();
  std::set<int>::iterator it;
  int i = 0;
  for (it = nlist.begin(); it != nlist.end(); ++it)
  {
    nodelist[i] = (*it);
    ++i;
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
void EXODUS::SideSet::fill_side_lists(int* elemlist, int* sidelist) const
{
  std::map<int, std::vector<int>> sides = get_side_set();
  std::map<int, std::vector<int>>::iterator it;
  int i = 0;
  for (it = sides.begin(); it != sides.end(); ++it)
  {
    elemlist[i] = it->second[0];
    sidelist[i] = it->second[1];
    ++i;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::SideSet::fill_side_lists(
    int* elemlist, int* sidelist, const std::map<int, std::vector<int>>& sides) const
{
  std::map<int, std::vector<int>>::const_iterator it;
  int i = 0;
  for (it = sides.begin(); it != sides.end(); ++it)
  {
    elemlist[i] = it->second[0];
    sidelist[i] = it->second[1];
    ++i;
  }
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


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
EXODUS::Mesh EXODUS::quadto_tri(EXODUS::Mesh& basemesh)
{
  std::map<int, std::shared_ptr<EXODUS::ElementBlock>>
      neweblocks;  // here the new EBlocks are stored
  std::map<int, std::shared_ptr<EXODUS::ElementBlock>> ebs = basemesh.get_element_blocks();
  std::map<int, std::shared_ptr<EXODUS::ElementBlock>>::const_iterator i_ebs;

  for (i_ebs = ebs.begin(); i_ebs != ebs.end(); ++i_ebs)
  {
    std::shared_ptr<EXODUS::ElementBlock> quadblock = i_ebs->second;
    EXODUS::ElementBlock::Shape quadshape = quadblock->get_shape();
    if ((quadshape != EXODUS::ElementBlock::quad4) && (quadshape != EXODUS::ElementBlock::shell4))
    {
      // FOUR_C_THROW("Only quad4 or shell4 in quad->tri conversion");
      std::cout << "Warning! Only quad4 or shell4 in quad->tri conversion. Skipping EBlock"
                << std::endl;
    }
    else
    {
      std::shared_ptr<std::map<int, std::vector<int>>> quad_conn = quadblock->get_ele_conn();
      std::map<int, std::vector<int>>::const_iterator i_quad;
      std::shared_ptr<std::map<int, std::vector<int>>> triconn =
          std::make_shared<std::map<int, std::vector<int>>>();

      for (i_quad = quad_conn->begin(); i_quad != quad_conn->end(); ++i_quad)
      {
        std::vector<int> quad = i_quad->second;
        std::vector<int> tri1(3);
        std::vector<int> tri2(3);
        tri1[0] = quad[0];
        tri1[1] = quad[1];
        tri1[2] = quad[2];
        tri2[0] = quad[2];
        tri2[1] = quad[3];
        tri2[2] = quad[0];
        int tri1_id = 2 * i_quad->first;
        int tri2_id = 2 * i_quad->first + 1;
        triconn->insert(std::pair<int, std::vector<int>>(tri1_id, tri1));
        triconn->insert(std::pair<int, std::vector<int>>(tri2_id, tri2));
      }

      std::shared_ptr<EXODUS::ElementBlock> triblock = std::make_shared<EXODUS::ElementBlock>(
          EXODUS::ElementBlock::tri3, triconn, quadblock->get_name());
      neweblocks.insert(
          std::pair<int, std::shared_ptr<EXODUS::ElementBlock>>(i_ebs->first, triblock));
      basemesh.erase_element_block(i_ebs->first);
    }
  }

  std::string newtitle = "trimesh";
  std::map<int, EXODUS::NodeSet> emptynodeset;
  std::map<int, EXODUS::SideSet> emptysideset;
  if (basemesh.get_num_side_sets() > 0)
    std::cout << "Warning! SideSets will not be transferred by quad->tri!" << std::endl;
  std::shared_ptr<std::map<int, std::vector<double>>> emptynodes =
      std::make_shared<std::map<int, std::vector<double>>>();
  EXODUS::Mesh trimesh(basemesh, emptynodes, neweblocks, emptynodeset, emptysideset, newtitle);
  return trimesh;
}



/*----------------------------------------------------------------------*
  TINY HELPER FUNCTIONS
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::print_map(std::ostream& os, const std::map<int, std::vector<int>> mymap)
{
  std::map<int, std::vector<int>>::const_iterator iter;
  for (iter = mymap.begin(); iter != mymap.end(); ++iter)
  {
    os << iter->first << ": ";
    std::vector<int> actvec = iter->second;
    std::vector<int>::iterator i;
    for (i = actvec.begin(); i < actvec.end(); ++i)
    {
      os << *i << ",";
    }
    os << std::endl;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::print_map(std::ostream& os, const std::map<int, std::set<int>> mymap)
{
  std::map<int, std::set<int>>::const_iterator iter;
  for (iter = mymap.begin(); iter != mymap.end(); ++iter)
  {
    os << iter->first << ": ";
    std::set<int> actset = iter->second;
    std::set<int>::iterator i;
    for (i = actset.begin(); i != actset.end(); ++i)
    {
      os << *i << ",";
    }
    os << std::endl;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::print_map(std::ostream& os, const std::map<int, std::vector<double>> mymap)
{
  std::map<int, std::vector<double>>::const_iterator iter;
  for (iter = mymap.begin(); iter != mymap.end(); ++iter)
  {
    os << iter->first << ": ";
    std::vector<double> actvec = iter->second;
    std::vector<double>::iterator i;
    for (i = actvec.begin(); i < actvec.end(); ++i)
    {
      os << *i << ",";
    }
    os << std::endl;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::print_map(std::ostream& os, const std::map<double, int> mymap)
{
  std::map<double, int>::const_iterator iter;
  for (iter = mymap.begin(); iter != mymap.end(); ++iter)
  {
    os << iter->first << ": " << iter->second;
    os << std::endl;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::print_map(std::ostream& os, const std::map<int, std::map<int, int>> mymap)
{
  std::map<int, std::map<int, int>>::const_iterator iter;
  for (iter = mymap.begin(); iter != mymap.end(); ++iter)
  {
    os << iter->first << ": ";
    std::map<int, int> actmap = iter->second;
    std::map<int, int>::iterator i;
    for (i = actmap.begin(); i != actmap.end(); ++i)
    {
      os << i->first << ": " << i->second;
    }
    os << std::endl;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::print_map(std::ostream& os, const std::map<int, std::pair<int, int>> mymap)
{
  std::map<int, std::pair<int, int>>::const_iterator iter;
  for (iter = mymap.begin(); iter != mymap.end(); ++iter)
  {
    os << iter->first << ": ";
    std::pair<int, int> actpair = iter->second;
    os << actpair.first << " <=> " << actpair.second;
    os << std::endl;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::print_map(std::ostream& os, const std::map<int, int> mymap)
{
  std::map<int, int>::const_iterator iter;
  for (iter = mymap.begin(); iter != mymap.end(); ++iter)
  {
    os << iter->first << ": ";
    os << iter->second;
    os << std::endl;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::print_map(std::ostream& os, const std::map<int, double> mymap)
{
  std::map<int, double>::const_iterator iter;
  for (iter = mymap.begin(); iter != mymap.end(); ++iter)
  {
    os << iter->first << ": ";
    os << iter->second;
    os << std::endl;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::print_vec(std::ostream& os, const std::vector<int> actvec)
{
  std::vector<int>::const_iterator i;
  for (i = actvec.begin(); i < actvec.end(); ++i)
  {
    os << *i << ",";
  }
  os << std::endl;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::print_vec(std::ostream& os, const std::vector<double> actvec)
{
  std::vector<double>::const_iterator i;
  for (i = actvec.begin(); i < actvec.end(); ++i)
  {
    os << *i << ",";
  }
  os << std::endl;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::print_set(std::ostream& os, const std::set<int> actset)
{
  std::set<int>::iterator i;
  for (i = actset.begin(); i != actset.end(); ++i)
  {
    os << *i << ",";
  }
  os << std::endl;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int EXODUS::hex_side_number_exo_to_four_c(const int exoface)
{
  constexpr std::array<int, 6> map = {1, 2, 3, 4, 0, 5};
  return map[exoface];
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int EXODUS::pyr_side_number_exo_to_four_c(const int exoface)
{
  constexpr std::array<int, 5> map = {1, 2, 3, 4, 0};
  return map[exoface];
}

FOUR_C_NAMESPACE_CLOSE
