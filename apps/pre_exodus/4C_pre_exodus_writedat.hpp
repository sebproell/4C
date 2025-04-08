// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_PRE_EXODUS_WRITEDAT_HPP
#define FOUR_C_PRE_EXODUS_WRITEDAT_HPP

#include "4C_config.hpp"

#include "4C_pre_exodus_readbc.hpp"

FOUR_C_NAMESPACE_OPEN

namespace EXODUS
{
  //! write outfile
  int write_dat_file(const std::string& outfile, const EXODUS::Mesh& mymesh,
      const std::string& headfile, const std::vector<EXODUS::ElemDef>& eledefs,
      const std::vector<EXODUS::CondDef>& condefs);

  //! conditions
  void write_dat_conditions(
      const std::vector<EXODUS::CondDef>& condefs, const EXODUS::Mesh& mymesh, std::ostream& dat);

  //! special condition SurfLocsys: calculates Normal and Tangent
  std::vector<double> calc_normal_surf_locsys(const int ns_id, const EXODUS::Mesh& m);

  //! DesignNode - Node Topology
  void write_dat_design_topology(
      const std::vector<EXODUS::CondDef>& condefs, const EXODUS::Mesh& mymesh, std::ostream& dat);

  //! return set of nodes in bc_entity
  std::set<int> get_ns_from_bc_entity(const EXODUS::CondDef& e, const EXODUS::Mesh& m);

  //! Nodes into outfile
  void write_dat_nodes(const EXODUS::Mesh& mymesh, std::ostream& dat);

  //! Elements into outfile
  void write_dat_eles(
      const std::vector<EXODUS::ElemDef>& eledefs, const EXODUS::Mesh& mymesh, std::ostream& dat);

  //! Elements from eblock
  void dat_eles(const EXODUS::ElementBlock& eb, const EXODUS::ElemDef& acte, int& startele,
      std::ostream& outfile, const int eb_id);

}  // namespace EXODUS

FOUR_C_NAMESPACE_CLOSE

#endif
