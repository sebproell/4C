// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_inpar_boolifyparameters.hpp"

#include "4C_utils_exceptions.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// Auxiliary routine to boolify integral Yes/No data
void Input::boolify_valid_input_parameters(
    Teuchos::ParameterList& list  ///< the valid input parameter list
)
{
  // collect parameter names with Yes/No entries
  // parse sub-lists as well
  std::vector<std::string> boolnames;  // collector
  for (Teuchos::ParameterList::ConstIterator i = list.begin(); i != list.end(); ++i)
  {
    if (list.isSublist(list.name(i)))
      boolify_valid_input_parameters(list.sublist(list.name(i)));
    else
    {
      if (list.isType<std::string>(list.name(i)))
      {
        const std::string value = list.get<std::string>(list.name(i));
        if ((value == "Yes") or (value == "YES") or (value == "yes") or (value == "No") or
            (value == "NO") or (value == "no"))
        {
          boolnames.push_back(list.name(i));
        }
      }
    }
  }

  // remove integral Yes/No parameters and replace them by Boolean
  for (std::vector<std::string>::iterator name = boolnames.begin(); name != boolnames.end(); ++name)
  {
    const std::string value = list.get<std::string>(*name);
    list.remove(*name);
    if ((value == "Yes") or (value == "YES") or (value == "yes"))
      list.set<bool>(*name, true);
    else if ((value == "No") or (value == "NO") or (value == "no"))
      list.set<bool>(*name, false);
    else
      FOUR_C_THROW("Cannot deal with entry \"{}\"", value.c_str());
  }
}


/*----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_CLOSE
