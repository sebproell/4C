// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fem_general_shape_function_type.hpp"

#include "4C_utils_exceptions.hpp"

#include <algorithm>
#include <map>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  namespace
  {
    const std::map<std::string, ShapeFunctionType> string2shapefuntype{
        {"Polynomial", ShapeFunctionType::polynomial}, {"Nurbs", ShapeFunctionType::nurbs},
        {"HDG", ShapeFunctionType::hdg}};
  }  // namespace
  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  ShapeFunctionType string_to_shape_function_type(std::string name)
  {
    const auto it = string2shapefuntype.find(name);
    if (it != string2shapefuntype.end()) return it->second;
    FOUR_C_THROW(
        "'{}' does not name a shape function type. Check for typos or consider adding the shape "
        "function type to the map.",
        name.c_str());
  }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  std::string shape_function_type_to_string(ShapeFunctionType shapefunctiontype)
  {
    const auto it = std::find_if(string2shapefuntype.begin(), string2shapefuntype.end(),
        [&](const auto& kv) { return kv.second == shapefunctiontype; });


    if (it != string2shapefuntype.end()) return it->first;
    FOUR_C_THROW(
        "Could not find the name of the given shape function type or the shapefunction is "
        "undefined.");
  }

  const std::map<std::string, ShapeFunctionType>& string_to_shape_function_type_map()
  {
    return string2shapefuntype;
  }
}  // namespace Core::FE
FOUR_C_NAMESPACE_CLOSE
