// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_DEAL_II_CONTEXT_IMPLEMENTATION_HPP
#define FOUR_C_DEAL_II_CONTEXT_IMPLEMENTATION_HPP

#include "4C_config.hpp"

#include "4C_deal_ii_context.hpp"
#include "4C_fem_discretization.hpp"

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/mapping_collection.h>

#include <unordered_map>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace DealiiWrappers::Internal
{
  template <int dim, int spacedim = dim>
  struct ContextImplementation
  {
    //! Store the local mapping between deal.II cells and 4C elements
    std::unordered_map<int, int> cell_index_to_element_lid;

    //! All dealii::FiniteElement objects that are required for the original 4C discretization.
    dealii::hp::FECollection<dim, spacedim> finite_elements;

    //! The names of the FiniteElements in #finite_elements.
    std::vector<std::string> finite_element_names;
  };
}  // namespace DealiiWrappers::Internal

FOUR_C_NAMESPACE_CLOSE

#endif
