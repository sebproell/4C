// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_IO_EXODUS_HPP
#define FOUR_C_IO_EXODUS_HPP

#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"

#include <set>
#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::IO::Exodus
{
  // forward declaration
  class ElementBlock;
  class NodeSet;
  class SideSet;

  /**
   * A class that stores the mesh information read from an Exodus file.
   */
  class Mesh
  {
   public:
    //! constructor
    Mesh(std::string exofilename);

    //! Print mesh info
    void print(std::ostream& os, bool verbose = false) const;

    //! Get number of nodes in mesh
    int get_num_nodes() const { return nodes_.size(); }

    //! Get number of elements in mesh
    int get_num_ele() const { return num_elem_; }

    //! Get number of dimensions
    int get_num_dim() const { return num_dim_; }

    //! Get number of dimensions
    int get_four_c_dim() const { return four_c_dim_; }

    //! Get ElementBlock map
    const std::map<int, ElementBlock>& get_element_blocks() const { return element_blocks_; }

    //! Get Number of ElementBlocks
    int get_num_element_blocks() const { return element_blocks_.size(); }

    //! Get one ElementBlock
    const ElementBlock& get_element_block(const int id) const;

    //! Get NodeSet map
    std::map<int, NodeSet> get_node_sets() const { return node_sets_; }

    //! Get Number of NodeSets
    int get_num_node_sets() const { return node_sets_.size(); }

    //! Get one NodeSet
    NodeSet get_node_set(const int id) const;

    //! Get SideSet map
    std::map<int, SideSet> get_side_sets() const { return side_sets_; }

    //! Get Number of SideSets
    int get_num_side_sets() const { return side_sets_.size(); }

    //! Get one SideSet
    SideSet get_side_set(const int id) const;

    //! Get edge Normal at node
    std::vector<double> normal(const int head1, const int origin, const int head2) const;

    //! Get normalized Vector between 2 nodes
    std::vector<double> node_vec(const int tail, const int head) const;

    //! Get Node map
    const std::map<int, std::vector<double>>& get_nodes() const { return nodes_; }

    //! Get one nodal coords
    std::vector<double> get_node(const int NodeID) const;

    //! Set one nodal coords
    void set_node(const int NodeID, const std::vector<double> coord);

    //! Set number of space dimensions
    void set_nsd(const int nsd);

   private:
    std::map<int, std::vector<double>> nodes_;

    std::map<int, ElementBlock> element_blocks_;

    std::map<int, NodeSet> node_sets_;

    std::map<int, SideSet> side_sets_;

    //! number of dimensions
    int num_dim_;

    //! number of dimensions for 4C problem (wall and fluid2 elements require 2d, although we have
    //! spatial dimensions)
    int four_c_dim_;

    //! number of elements
    int num_elem_;

    //! title
    std::string title_;
  };


  /*!
  \brief ElementBlock is a set of Elements of same discretization Type

  A Element Block is a tiny class storing element-type, name, etc. of a ElementBlock
  It implements its printout.

  */
  class ElementBlock
  {
   public:
    ElementBlock(FE::CellType cell_type, std::map<int, std::vector<int>> eleconn, std::string name);

    FE::CellType get_shape() const { return cell_type_; }

    int get_num_ele() const { return eleconn_.size(); }

    const std::map<int, std::vector<int>>& get_ele_conn() const { return eleconn_; }

    const std::vector<int>& get_ele_nodes(int i) const;

    std::string get_name() const { return name_; }

    void print(std::ostream& os, bool verbose = false) const;

   private:
    FE::CellType cell_type_;

    //! Element Connectivity
    std::map<int, std::vector<int>> eleconn_;

    std::string name_;
  };

  class NodeSet
  {
   public:
    NodeSet(const std::set<int>& nodeids, const std::string& name);

    const std::set<int>& get_node_set() const { return nodeids_; };

    std::string get_name() const { return name_; };

    inline int get_num_nodes() const { return nodeids_.size(); }

    void print(std::ostream& os, bool verbose = false) const;

   private:
    std::set<int> nodeids_;  // nodids in NodeSet
    std::string name_;       // NodeSet name
  };

  class SideSet
  {
   public:
    SideSet(const std::map<int, std::vector<int>>& sides, const std::string& name);

    inline int get_num_sides() const { return sides_.size(); }

    std::string get_name() const { return name_; }

    const std::map<int, std::vector<int>>& get_side_set() const { return sides_; }

    void print(std::ostream& os, bool verbose = false) const;

   private:
    std::map<int, std::vector<int>> sides_;
    std::string name_;
  };
}  // namespace Core::IO::Exodus

FOUR_C_NAMESPACE_CLOSE

#endif
