// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_porofluidmultiphase_resulttest.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_global_data.hpp"
#include "4C_porofluidmultiphase_meshtying_strategy_base.hpp"
#include "4C_porofluidmultiphase_timint_implicit.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | ctor                                                     vuong 08/16 |
 *----------------------------------------------------------------------*/
POROFLUIDMULTIPHASE::ResultTest::ResultTest(TimIntImpl& porotimint)
    : Core::Utils::ResultTest("POROFLUIDMULTIPHASE"), porotimint_(porotimint)
{
  return;
}


/*----------------------------------------------------------------------*
 | test node                                                vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::ResultTest::test_node(
    const Core::IO::InputParameterContainer& container, int& nerr, int& test_count)
{
  // care for the case of multiple discretizations of the same field type
  std::string dis = container.get<std::string>("DIS");
  if (dis != porotimint_.discretization()->name()) return;

  int node = container.get<int>("NODE");
  node -= 1;

  int havenode(porotimint_.discretization()->have_global_node(node));
  int isnodeofanybody(0);
  Core::Communication::sum_all(
      &havenode, &isnodeofanybody, 1, porotimint_.discretization()->get_comm());

  if (isnodeofanybody == 0)
  {
    FOUR_C_THROW("Node {} does not belong to discretization {}", node + 1,
        porotimint_.discretization()->name().c_str());
  }
  else
  {
    if (porotimint_.discretization()->have_global_node(node))
    {
      Core::Nodes::Node* actnode = porotimint_.discretization()->g_node(node);

      // Here we are just interested in the nodes that we own (i.e. a row node)!
      if (actnode->owner() !=
          Core::Communication::my_mpi_rank(porotimint_.discretization()->get_comm()))
        return;

      // extract name of quantity to be tested
      std::string quantity = container.get<std::string>("QUANTITY");

      // get result to be tested
      const double result = result_node(quantity, actnode);

      nerr += compare_values(result, "NODE", container);
      test_count++;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | test element                                        kremheller 10/19 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::ResultTest::test_element(
    const Core::IO::InputParameterContainer& container, int& nerr, int& test_count)
{
  // care for the case of multiple discretizations of the same field type
  std::string dis = container.get<std::string>("DIS");

  if (dis != porotimint_.discretization()->name()) return;

  int element = container.get<int>("ELEMENT");
  element -= 1;

  int haveelement(porotimint_.discretization()->have_global_element(element));
  int iselementofanybody(0);
  Core::Communication::sum_all(
      &haveelement, &iselementofanybody, 1, porotimint_.discretization()->get_comm());

  if (iselementofanybody == 0)
  {
    FOUR_C_THROW("Element {} does not belong to discretization {}", element + 1,
        porotimint_.discretization()->name().c_str());
  }
  else
  {
    if (porotimint_.discretization()->have_global_element(element))
    {
      const Core::Elements::Element* actelement = porotimint_.discretization()->g_element(element);

      // Here we are just interested in the elements that we own (i.e. a row element)!
      if (actelement->owner() !=
          Core::Communication::my_mpi_rank(porotimint_.discretization()->get_comm()))
        return;

      // extract name of quantity to be tested
      std::string quantity = container.get<std::string>("QUANTITY");

      // get result to be tested
      const double result = result_element(quantity, actelement);

      nerr += compare_values(result, "ELEMENT", container);
      test_count++;
    }
  }
}


/*----------------------------------------------------------------------*
 | get nodal result to be tested                            vuong 08/16 |
 *----------------------------------------------------------------------*/
double POROFLUIDMULTIPHASE::ResultTest::result_node(
    const std::string quantity, Core::Nodes::Node* node) const
{
  // initialize variable for result
  double result(0.);

  // extract row map from solution vector
  const Epetra_BlockMap& phinpmap = porotimint_.phinp()->get_map();

  // test result value of phi field
  if (quantity == "phi")
    result = (*porotimint_.phinp())[phinpmap.LID(porotimint_.discretization()->dof(0, node, 0))];

  // test result value for a system of scalars
  else if (!quantity.compare(0, 3, "phi"))
  {
    // read species ID
    std::string k_string = quantity.substr(3);
    char* locator(nullptr);
    int k = strtol(k_string.c_str(), &locator, 10) - 1;
    if (locator == k_string.c_str()) FOUR_C_THROW("Couldn't read species ID!");

    if (porotimint_.discretization()->num_dof(0, node) <= k)
      FOUR_C_THROW("Species ID is larger than number of DOFs of node!");

    // extract result
    result = (*porotimint_.phinp())[phinpmap.LID(porotimint_.discretization()->dof(0, node, k))];
  }

  // test result value of phi field
  else if (quantity == "pressure")
    result = (*porotimint_.pressure())[phinpmap.LID(porotimint_.discretization()->dof(0, node, 0))];

  // test result value for a system of scalars
  else if (!quantity.compare(0, 8, "pressure"))
  {
    // read species ID
    std::string k_string = quantity.substr(8);
    char* locator(nullptr);
    int k = strtol(k_string.c_str(), &locator, 10) - 1;
    if (locator == k_string.c_str()) FOUR_C_THROW("Couldn't read pressure ID!");

    if (porotimint_.discretization()->num_dof(0, node) <= k)
      FOUR_C_THROW("Pressure ID is larger than number of DOFs of node!");

    // extract result
    result = (*porotimint_.pressure())[phinpmap.LID(porotimint_.discretization()->dof(0, node, k))];
  }

  // test result value of phi field
  else if (quantity == "saturation")
    result =
        (*porotimint_.saturation())[phinpmap.LID(porotimint_.discretization()->dof(0, node, 0))];

  // test result value for a system of scalars
  else if (!quantity.compare(0, 10, "saturation"))
  {
    // read species ID
    std::string k_string = quantity.substr(10);
    char* locator(nullptr);
    int k = strtol(k_string.c_str(), &locator, 10) - 1;
    if (locator == k_string.c_str()) FOUR_C_THROW("Couldn't read saturation ID!");

    if (porotimint_.discretization()->num_dof(0, node) <= k)
      FOUR_C_THROW("Saturation ID is larger than number of DOFs of node!");

    // extract result
    result =
        (*porotimint_.saturation())[phinpmap.LID(porotimint_.discretization()->dof(0, node, k))];
  }

  // catch unknown quantity strings
  else
    FOUR_C_THROW("Quantity '{}' not supported in result test!", quantity.c_str());

  return result;
}  // POROFLUIDMULTIPHASE::ResultTest::ResultNode

/*----------------------------------------------------------------------*
 | get element result to be tested                     kremheller 10/19 |
 *----------------------------------------------------------------------*/
double POROFLUIDMULTIPHASE::ResultTest::result_element(
    const std::string quantity, const Core::Elements::Element* element) const
{
  // initialize variable for result
  double result(0.);

  if (quantity == "bloodvesselvolfrac")
  {
    result = (*porotimint_.mesh_tying_strategy()
            ->blood_vessel_volume_fraction())[porotimint_.discretization()->element_row_map()->LID(
        element->id())];
  }
  else if (!quantity.compare(0, 13, "phasevelocity"))
  {
    const int num_dim = Global::Problem::instance()->n_dim();
    // get phase ID
    // example: "phasevelocity3x" -> k = 2 (phase IDs start at index 0)
    std::string k_string = quantity.substr(13);
    char* locator(nullptr);
    auto idx_poro_dof = int(strtol(k_string.c_str(), &locator, 13) - 1);
    if (locator == k_string.c_str()) FOUR_C_THROW("Could not read phase ID in result test!");

    // get spatial dimension
    int idx_dim(-1);
    if (!quantity.compare(14, 15, "x"))
      idx_dim = 0;
    else if (!quantity.compare(14, 15, "y"))
      idx_dim = 1;
    else if (!quantity.compare(14, 15, "z"))
      idx_dim = 2;

    result = ((*porotimint_.phase_velocity())(
        idx_poro_dof * num_dim +
        idx_dim))[porotimint_.discretization()->element_row_map()->LID(element->id())];
  }
  // catch unknown quantity strings
  else
    FOUR_C_THROW("Quantity '{}' not supported in result test!", quantity.c_str());

  return result;
}

/*-------------------------------------------------------------------------------------*
 | test special quantity not associated with a particular element or node  vuong 08/16 |
 *-------------------------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::ResultTest::test_special(
    const Core::IO::InputParameterContainer& container, int& nerr, int& test_count)
{
  // make sure that quantity is tested only once
  if (Core::Communication::my_mpi_rank(porotimint_.discretization()->get_comm()) == 0)
  {
    // extract name of quantity to be tested
    std::string quantity = container.get<std::string>("QUANTITY");

    // get result to be tested
    const double result = result_special(quantity);

    // compare values
    const int err = compare_values(result, "SPECIAL", container);
    nerr += err;
    test_count++;
  }
}


/*----------------------------------------------------------------------*
 | get special result to be tested                          vuong 08/16 |
 *----------------------------------------------------------------------*/
double POROFLUIDMULTIPHASE::ResultTest::result_special(
    const std::string quantity  //! name of quantity to be tested
) const
{
  // initialize variable for result
  double result(0.);

  if (quantity == "numiterlastnewton") result = (double)porotimint_.iter_num();
  // result test of domain integrals
  else if (!quantity.compare(0, 22, "domain_integral_value_"))
  {
    // get the index of the value which should be checked
    std::string suffix = quantity.substr(22);
    int idx = -1;
    try
    {
      idx = std::stoi(suffix);
    }
    catch (const std::invalid_argument& e)
    {
      FOUR_C_THROW(
          "You provided the wrong format for output of domain_integral_values. The integer number "
          "must be at the very last position of the name, separated by an underscore.\n"
          "The correct format is: domain_integral_value_<number>");
    }

    // index should be in range [0, number_functions - 1]
    if (idx < 0 || idx >= porotimint_.num_domain_int_functions())
      FOUR_C_THROW("detected wrong index {}, index should be in range [0,{}]", idx,
          porotimint_.num_domain_int_functions() - 1);

    // return the result
    result = (*porotimint_.domain_int_values())[idx];
  }
  // catch unknown quantity strings
  else
    FOUR_C_THROW("Quantity '{}' not supported in result test!", quantity.c_str());

  return result;
}  // POROFLUIDMULTIPHASE::ResultTest::result_special

FOUR_C_NAMESPACE_CLOSE
