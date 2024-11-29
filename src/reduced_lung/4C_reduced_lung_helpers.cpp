// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_config.hpp"

#include "4C_reduced_lung_helpers.hpp"

#include "4C_fem_condition.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_red_airways_elementbase.hpp"


FOUR_C_NAMESPACE_OPEN

namespace ReducedLung
{

  Epetra_Map create_domain_map(const Epetra_Comm& comm, const std::vector<Airway>& airways,
      const std::vector<TerminalUnit>& terminal_units)
  {
    std::vector<int> locally_owned_dof_indices;
    for (const auto& airway : airways)
    {
      locally_owned_dof_indices.insert(locally_owned_dof_indices.end(),
          airway.global_dof_ids.begin(), airway.global_dof_ids.end());
    }
    for (const auto& terminal_unit : terminal_units)
    {
      locally_owned_dof_indices.insert(locally_owned_dof_indices.end(),
          terminal_unit.global_dof_ids.begin(), terminal_unit.global_dof_ids.end());
    }
    const Epetra_Map domain_map(
        -1, locally_owned_dof_indices.size(), locally_owned_dof_indices.data(), 0, comm);

    return domain_map;
  }

  Epetra_Map create_row_map(const Epetra_Comm& comm, const std::vector<Airway>& airways,
      const std::vector<TerminalUnit>& terminal_units, const std::vector<Connection>& connections,
      const std::vector<Bifurcation>& bifurcations,
      const std::vector<BoundaryCondition>& boundary_conditions)
  {
    int n_local_state_equations = 0;
    for (const auto& airway : airways)
    {
      n_local_state_equations += airway.n_state_equations;
    }
    for (const auto& terminal_unit : terminal_units)
    {
      n_local_state_equations += terminal_unit.n_state_equations;
    }
    // Intermediate maps for the different equation types
    const Epetra_Map state_equations(-1, n_local_state_equations, 0, comm);
    const Epetra_Map couplings(-1, connections.size() * 2 + bifurcations.size() * 3,
        state_equations.NumGlobalElements(), comm);
    const Epetra_Map boundaries(-1, boundary_conditions.size(),
        couplings.NumGlobalElements() + state_equations.NumGlobalElements(), comm);

    //  Merge all maps to the full local matrix row map
    std::vector<int> global_row_indices;
    int n_local_indices =
        state_equations.NumMyElements() + couplings.NumMyElements() + boundaries.NumMyElements();
    global_row_indices.reserve(n_local_indices);
    global_row_indices.insert(global_row_indices.end(), state_equations.MyGlobalElements(),
        state_equations.MyGlobalElements() + state_equations.NumMyElements());
    global_row_indices.insert(global_row_indices.end(), couplings.MyGlobalElements(),
        couplings.MyGlobalElements() + couplings.NumMyElements());
    global_row_indices.insert(global_row_indices.end(), boundaries.MyGlobalElements(),
        boundaries.MyGlobalElements() + boundaries.NumMyElements());

    const Epetra_Map row_map(-1, n_local_indices, global_row_indices.data(), 0, comm);

    return row_map;
  }

  Epetra_Map create_column_map(const Epetra_Comm& comm, const std::vector<Airway>& airways,
      const std::vector<TerminalUnit>& terminal_units, const std::map<int, int>& global_dof_per_ele,
      const std::map<int, int>& first_global_dof_of_ele, const std::vector<Connection>& connections,
      const std::vector<Bifurcation>& bifurcations,
      const std::vector<BoundaryCondition>& boundary_conditions)
  {
    // Vector for intermediate storage of necessary dof ids
    std::vector<int> locally_relevant_dof_indices;

    // Loop over all elements and add their global dof ids
    for (const auto& airway : airways)
    {
      locally_relevant_dof_indices.insert(locally_relevant_dof_indices.end(),
          airway.global_dof_ids.begin(), airway.global_dof_ids.end());
    }
    for (const auto& terminal_unit : terminal_units)
    {
      locally_relevant_dof_indices.insert(locally_relevant_dof_indices.end(),
          terminal_unit.global_dof_ids.begin(), terminal_unit.global_dof_ids.end());
    }

    // Loop over all connections of two elements and add relevant dof ids (p and q associated with
    // the end of the parent element and p and q associated with the start of the child element)
    for (const auto& conn : connections)
    {
      int n_dofs_parent = global_dof_per_ele.find(conn.global_parent_element_id)->second;
      // p2 always second dof of element
      locally_relevant_dof_indices.insert(locally_relevant_dof_indices.end(),
          first_global_dof_of_ele.find(conn.global_parent_element_id)->second + 1);
      // q_out always last dof of element
      locally_relevant_dof_indices.insert(locally_relevant_dof_indices.end(),
          first_global_dof_of_ele.find(conn.global_parent_element_id)->second + n_dofs_parent - 1);
      // p1 always first dof of element
      locally_relevant_dof_indices.insert(locally_relevant_dof_indices.end(),
          first_global_dof_of_ele.find(conn.global_child_element_id)->second);
      // q_in always third dof of element
      locally_relevant_dof_indices.insert(locally_relevant_dof_indices.end(),
          first_global_dof_of_ele.find(conn.global_child_element_id)->second + 2);
    }

    // Loop over all bifurcations and add relevant dof ids (p and q associated with the
    // end of the parent element and p and q associated with the start of the child elements)
    for (const auto& bif : bifurcations)
    {
      int n_dofs_parent = global_dof_per_ele.find(bif.global_parent_element_id)->second;
      // p2 always second dof of element
      locally_relevant_dof_indices.insert(locally_relevant_dof_indices.end(),
          first_global_dof_of_ele.find(bif.global_parent_element_id)->second + 1);
      // q_out always last dof of element
      locally_relevant_dof_indices.insert(locally_relevant_dof_indices.end(),
          first_global_dof_of_ele.find(bif.global_parent_element_id)->second + n_dofs_parent - 1);
      // p1 always first dof of element
      locally_relevant_dof_indices.insert(locally_relevant_dof_indices.end(),
          first_global_dof_of_ele.find(bif.global_child_1_element_id)->second);
      locally_relevant_dof_indices.insert(locally_relevant_dof_indices.end(),
          first_global_dof_of_ele.find(bif.global_child_2_element_id)->second);
      // q_in always third dof of element
      locally_relevant_dof_indices.insert(locally_relevant_dof_indices.end(),
          first_global_dof_of_ele.find(bif.global_child_1_element_id)->second + 2);
      locally_relevant_dof_indices.insert(locally_relevant_dof_indices.end(),
          first_global_dof_of_ele.find(bif.global_child_2_element_id)->second + 2);
    }
    // Loop over all boundary conditions and add relevant dof ids (dof where bc is applied)
    for (const auto& bc : boundary_conditions)
    {
      switch (bc.bc_type)
      {
        case BoundaryConditionType::pressure_in:
          locally_relevant_dof_indices.insert(locally_relevant_dof_indices.end(),
              first_global_dof_of_ele.find(bc.global_element_id)->second);
          break;
        case BoundaryConditionType::pressure_out:
          locally_relevant_dof_indices.insert(locally_relevant_dof_indices.end(),
              first_global_dof_of_ele.find(bc.global_element_id)->second + 1);
          break;
        case BoundaryConditionType::flow_in:
          locally_relevant_dof_indices.insert(locally_relevant_dof_indices.end(),
              first_global_dof_of_ele.find(bc.global_element_id)->second + 2);
          break;
        case BoundaryConditionType::flow_out:
          locally_relevant_dof_indices.insert(locally_relevant_dof_indices.end(),
              first_global_dof_of_ele.find(bc.local_bc_id)->second +
                  global_dof_per_ele.find(bc.global_element_id)->second - 1);
          break;
      }
    }

    // Erase duplicate dof indices and sort the remaining ids
    std::sort(locally_relevant_dof_indices.begin(), locally_relevant_dof_indices.end());
    locally_relevant_dof_indices.erase(
        std::unique(locally_relevant_dof_indices.begin(), locally_relevant_dof_indices.end()),
        locally_relevant_dof_indices.end());

    const Epetra_Map column_map(
        -1, locally_relevant_dof_indices.size(), locally_relevant_dof_indices.data(), 0, comm);

    return column_map;
  }
}  // namespace ReducedLung

FOUR_C_NAMESPACE_CLOSE