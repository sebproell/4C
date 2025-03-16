// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_solver_nonlin_nox_statustest_factory.hpp"  // class definition

#include "4C_solver_nonlin_nox_enum_lists.hpp"
#include "4C_solver_nonlin_nox_statustest_activeset.hpp"
#include "4C_solver_nonlin_nox_statustest_combo.hpp"
#include "4C_solver_nonlin_nox_statustest_normf.hpp"
#include "4C_solver_nonlin_nox_statustest_normupdate.hpp"
#include "4C_solver_nonlin_nox_statustest_normwrms.hpp"

#include <NOX_Abstract_Vector.H>
#include <NOX_StatusTest_Factory.H>
#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::StatusTest::Factory::Factory()
    : noxfactory_(Teuchos::make_rcp<::NOX::StatusTest::Factory>())
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::StatusTest::Generic> NOX::Nln::StatusTest::Factory::build_outer_status_tests(
    Teuchos::ParameterList& p, const ::NOX::Utils& u,
    std::map<std::string, Teuchos::RCP<::NOX::StatusTest::Generic>>* tagged_tests) const
{
  Teuchos::RCP<::NOX::StatusTest::Generic> status_test;

  std::string test_type = "???";

  if (Teuchos::isParameterType<std::string>(p, "Test Type"))
    test_type = Teuchos::get<std::string>(p, "Test Type");
  else
  {
    throw_error("build_outer_status_tests()", "The \"Test Type\" is a required parameter!");
  }

  if (test_type == "Combo")
    status_test = this->build_combo_test(p, u, tagged_tests);
  else if (test_type == "NormF")
    status_test = this->build_norm_f_test(p, u);
  else if (test_type == "RelativeNormF")
    status_test = this->build_norm_f_test(p, u, true);
  else if (test_type == "NormWRMS")
    status_test = this->build_norm_wrms_test(p, u);
  else if (test_type == "NormUpdate" or test_type == "NormUpdateSkipFirstIter")
    status_test = this->build_norm_update_test(p, u);
  else if (test_type == "ActiveSet")
    status_test = this->build_active_set_test(p, u);
  else
  {
    status_test = noxfactory_->buildStatusTests(p, u, tagged_tests);
  }

  this->check_and_tag_test(p, status_test, tagged_tests);

  return status_test;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::StatusTest::Generic> NOX::Nln::StatusTest::Factory::build_norm_f_test(
    Teuchos::ParameterList& p, const ::NOX::Utils& u) const
{
  return build_norm_f_test(p, u, false);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::StatusTest::Generic> NOX::Nln::StatusTest::Factory::build_norm_f_test(
    Teuchos::ParameterList& p, const ::NOX::Utils& u, const bool& relativeNormF) const
{
  // initialized to size zero
  std::vector<double> tolerances = std::vector<double>(0);
  std::vector<::NOX::StatusTest::NormF::ToleranceType> toltypes(0);
  std::vector<QuantityType> quantitytypes(0);
  std::vector<::NOX::Abstract::Vector::NormType> normtypes(0);
  std::vector<::NOX::StatusTest::NormF::ScaleType> scaletypes(0);

  int i = 0;
  std::ostringstream quantity_string;
  quantity_string << "Quantity " << i;
  while ((i == 0) or (p.isSublist(quantity_string.str())))
  {
    // If we use only one quantity per test, there is no need to define a
    // extra sub-list.
    Teuchos::RCP<Teuchos::ParameterList> ql = Teuchos::rcpFromRef(p);
    if (p.isSublist(quantity_string.str()))
      ql = Teuchos::rcpFromRef(p.sublist(quantity_string.str(), true));

    // ------------------------------------------
    // get the quantity type
    // ------------------------------------------
    std::string quantity_type_string = ql->get("Quantity Type", "???");
    QuantityType qType = string_to_quantity_type(quantity_type_string);
    if (qType == quantity_unknown)
    {
      std::ostringstream msg;
      msg << "The \"Quantity Type\" is a required parameter \n"
          << "and the chosen option \"" << quantity_type_string << "\" is invalid!";
      throw_error("build_norm_f_test()", msg.str());
    }
    else
      quantitytypes.push_back(qType);

    // ------------------------------------------
    // get the tolerance type
    // ------------------------------------------
    if (relativeNormF)
      toltypes.push_back(::NOX::StatusTest::NormF::Relative);
    else
    {
      std::string tol_type_string = ql->get("Tolerance Type", "Absolute");
      if (tol_type_string == "Absolute")
        toltypes.push_back(::NOX::StatusTest::NormF::Absolute);
      else if (tol_type_string == "Relative")
        toltypes.push_back(::NOX::StatusTest::NormF::Relative);
      else
      {
        std::string msg = "\"Tolerance Type\" must be either \"Absolute\" or \"Relative\"!";
        throw_error("build_norm_f_test()", msg);
      }
    }

    // ------------------------------------------
    // get the tolerance
    // ------------------------------------------
    if (toltypes.at(i) == ::NOX::StatusTest::NormF::Absolute)
      tolerances.push_back(ql->get("Tolerance", 1.0e-8));
    else
      tolerances.push_back(ql->get("Tolerance", 1.0e-5));

    // ------------------------------------------
    // get the norm type
    // ------------------------------------------
    std::string norm_type_string = ql->get("Norm Type", "Two Norm");
    if (norm_type_string == "Two Norm")
      normtypes.push_back(::NOX::Abstract::Vector::TwoNorm);
    else if (norm_type_string == "One Norm")
      normtypes.push_back(::NOX::Abstract::Vector::OneNorm);
    else if (norm_type_string == "Max Norm")
      normtypes.push_back(::NOX::Abstract::Vector::MaxNorm);
    else
    {
      std::string msg = "\"Norm Type\" must be either \"Two Norm\", \"One Norm\", or \"Max Norm\"!";
      throw_error("build_norm_f_test()", msg);
    }

    // ------------------------------------------
    // get the scale type
    // ------------------------------------------
    std::string scale_type_string = ql->get("Scale Type", "Unscaled");
    if (scale_type_string == "Unscaled")
      scaletypes.push_back(::NOX::StatusTest::NormF::Unscaled);
    else if (scale_type_string == "Scaled")
      scaletypes.push_back(::NOX::StatusTest::NormF::Scaled);
    else
    {
      std::string msg = "\"Scale Type\" must be either \"Unscaled\" or \"Scaled\"!";
      throw_error("build_norm_f_test()", msg);
    }
    // ------------------------------------------
    // increase iterator
    // ------------------------------------------
    ++i;
    quantity_string.str("");
    quantity_string << "Quantity " << i;
  }  // loop over all quantity types

  // build the normF status test
  Teuchos::RCP<NOX::Nln::StatusTest::NormF> status_test =
      Teuchos::make_rcp<NOX::Nln::StatusTest::NormF>(
          quantitytypes, toltypes, tolerances, normtypes, scaletypes, &u);

  return status_test;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::StatusTest::Generic> NOX::Nln::StatusTest::Factory::build_norm_update_test(
    Teuchos::ParameterList& p, const ::NOX::Utils& u) const
{
  // initialized to size zero
  std::vector<double> tolerances(0.0);
  std::vector<NOX::Nln::StatusTest::NormUpdate::ToleranceType> toltypes(0);
  std::vector<QuantityType> quantitytypes(0);
  std::vector<::NOX::Abstract::Vector::NormType> normtypes(0);
  std::vector<NOX::Nln::StatusTest::NormUpdate::ScaleType> scaletypes(0);

  // ------------------------------------------
  // Get Alpha
  // Minimum step size allowed during a
  // line search for incr norm to be flagged
  // as converged.
  // ------------------------------------------
  double alpha = p.get("Alpha", 1.0);

  // ------------------------------------------
  // Get Beta
  // Maximum linear solve tolerance allowed
  // for incr norm to be flagged as converged.
  // ------------------------------------------
  double beta = p.get("Beta", 0.5);

  int i = 0;
  std::ostringstream quantity_string;
  quantity_string << "Quantity " << i;
  while ((i == 0) or (p.isSublist(quantity_string.str())))
  {
    // If we use only one quantity per test, there is no need to define a
    // extra sub-list.
    Teuchos::RCP<Teuchos::ParameterList> ql = Teuchos::rcpFromRef(p);
    if (p.isSublist(quantity_string.str()))
      ql = Teuchos::rcpFromRef(p.sublist(quantity_string.str(), true));

    // ------------------------------------------
    // get the quantity type
    // ------------------------------------------
    std::string quantity_type_string = ql->get("Quantity Type", "???");
    QuantityType qType = string_to_quantity_type(quantity_type_string);
    if (qType == quantity_unknown)
    {
      std::ostringstream msg;
      msg << "The \"Quantity Type\" is a required parameter \n"
          << "and the chosen option \"" << quantity_type_string << "\" is invalid!";
      throw_error("build_norm_update_test()", msg.str());
    }
    else
      quantitytypes.push_back(qType);

    // ------------------------------------------
    // get the tolerance type
    // ------------------------------------------
    std::string tol_type_string = ql->get("Tolerance Type", "Absolute");
    if (tol_type_string == "Absolute")
      toltypes.push_back(NOX::Nln::StatusTest::NormUpdate::Absolute);
    else if (tol_type_string == "Relative")
      toltypes.push_back(NOX::Nln::StatusTest::NormUpdate::Relative);
    else
    {
      std::string msg = "\"Tolerance Type\" must be either \"Absolute\" or \"Relative\"!";
      throw_error("build_norm_update_test()", msg);
    }

    // ------------------------------------------
    // get the tolerance
    // ------------------------------------------
    if (toltypes.at(i) == NOX::Nln::StatusTest::NormUpdate::Absolute)
      tolerances.push_back(ql->get("Tolerance", 1.0e-8));
    else
      tolerances.push_back(ql->get("Tolerance", 1.0e-5));

    // ------------------------------------------
    // get the norm type
    // ------------------------------------------
    std::string norm_type_string = ql->get("Norm Type", "Two Norm");
    if (norm_type_string == "Two Norm")
      normtypes.push_back(::NOX::Abstract::Vector::TwoNorm);
    else if (norm_type_string == "One Norm")
      normtypes.push_back(::NOX::Abstract::Vector::OneNorm);
    else if (norm_type_string == "Max Norm")
      normtypes.push_back(::NOX::Abstract::Vector::MaxNorm);
    else
    {
      std::string msg = "\"Norm Type\" must be either \"Two Norm\", \"One Norm\", or \"Max Norm\"!";
      throw_error("build_norm_update_test()", msg);
    }

    // ------------------------------------------
    // get the scale type
    // ------------------------------------------
    std::string scale_type_string = ql->get("Scale Type", "Unscaled");
    if (scale_type_string == "Unscaled")
      scaletypes.push_back(NOX::Nln::StatusTest::NormUpdate::Unscaled);
    else if (scale_type_string == "Scaled")
      scaletypes.push_back(NOX::Nln::StatusTest::NormUpdate::Scaled);
    else
    {
      std::string msg = "\"Scale Type\" must be either \"Unscaled\" or \"Scaled\"!";
      throw_error("build_norm_update_test()", msg);
    }
    // ------------------------------------------
    // increase iterator
    // ------------------------------------------
    ++i;
    quantity_string.str("");
    quantity_string << "Quantity " << i;
  }  // loop over all quantity types

  // build the normUpdate status test
  if (Teuchos::get<std::string>(p, "Test Type") == "NormUpdateSkipFirstIter")
    return Teuchos::make_rcp<NOX::Nln::StatusTest::NormUpdateSkipFirstIter>(
        quantitytypes, toltypes, tolerances, normtypes, scaletypes, alpha, beta, &u);
  else
    return Teuchos::make_rcp<NOX::Nln::StatusTest::NormUpdate>(
        quantitytypes, toltypes, tolerances, normtypes, scaletypes, alpha, beta, &u);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::StatusTest::Generic> NOX::Nln::StatusTest::Factory::build_norm_wrms_test(
    Teuchos::ParameterList& p, const ::NOX::Utils& u) const
{
  // initialized to size zero
  std::vector<QuantityType> quantitytypes = std::vector<QuantityType>(0);
  std::vector<double> atol = std::vector<double>(0);
  std::vector<double> rtol = std::vector<double>(0);
  std::vector<double> tol = std::vector<double>(0);
  std::vector<double> BDFMultiplier = std::vector<double>(0);
  std::vector<bool> disable_implicit_weighting = std::vector<bool>(0);

  // ------------------------------------------
  // Get Alpha
  // Minimum step size allowed during a
  // line search for WRMS norm to be flagged
  // as converged.
  // ------------------------------------------
  double alpha = p.get("Alpha", 1.0);

  // ------------------------------------------
  // Get Beta
  // Maximum linear solve tolerance allowed
  // for WRMS norm to be flagged as converged.
  // ------------------------------------------
  double beta = p.get("Beta", 0.5);

  int i = 0;
  std::ostringstream quantity_string;
  quantity_string << "Quantity " << i;
  while ((i == 0) or (p.isSublist(quantity_string.str())))
  {
    // If we use only one quantity per test, there is no need to define a
    // extra sub-list.
    Teuchos::RCP<Teuchos::ParameterList> ql = Teuchos::rcpFromRef(p);
    if (p.isSublist(quantity_string.str()))
      ql = Teuchos::rcpFromRef(p.sublist(quantity_string.str(), true));

    // ------------------------------------------
    // get the Quantity Type
    // ------------------------------------------
    std::string quantity_type_string = ql->get("Quantity Type", "???");
    QuantityType qType = string_to_quantity_type(quantity_type_string);
    if (qType == quantity_unknown)
    {
      std::ostringstream msg;
      msg << "The \"Quantity Type\" is a required parameter \n"
          << "and the chosen option \"" << quantity_type_string << "\" is invalid!";
      throw_error("build_norm_wrms_test()", msg.str());
    }
    else
      quantitytypes.push_back(qType);

    // ------------------------------------------
    // get Absolute Tolerance
    // ------------------------------------------
    atol.push_back(ql->get("Absolute Tolerance", 1.0e-8));

    // ------------------------------------------
    // get Relative Tolerance
    // ------------------------------------------
    rtol.push_back(ql->get("Relative Tolerance", 1.0e-5));

    // ------------------------------------------
    // get Tolerance
    // Required tolerance for the NormWRMS to be
    // declared converged.
    // ------------------------------------------
    tol.push_back(ql->get("Tolerance", 1.0e-5));

    // ------------------------------------------
    // get BDF multiplier
    // ------------------------------------------
    BDFMultiplier.push_back(ql->get("BDF Multiplier", 1.0));

    // ------------------------------------------
    // get Disable Implicit Weighting
    // ------------------------------------------
    if ((ql->isParameter("Disable Implicit Weighting") and
            ql->isType<bool>("Disable Implicit Weighting")) or
        not ql->isParameter("Disable Implicit Weighting"))
    {
      const bool isdiw = ql->get<bool>("Disable Implicit Weighting", true);
      disable_implicit_weighting.push_back(isdiw);
    }
    else
    {
      std::string msg =
          "\"Disable Implicit Weighting\" must be either "
          "\"Yes\", \"yes\", \"No\" or \"no\"!";
      throw_error("build_norm_wrms_test()", msg);
    }

    // ------------------------------------------
    // increase iterator
    // ------------------------------------------
    ++i;
    quantity_string.str("");
    quantity_string << "Quantity " << i;
  }  // loop over all quantity types

  Teuchos::RCP<NOX::Nln::StatusTest::NormWRMS> status_test =
      Teuchos::make_rcp<NOX::Nln::StatusTest::NormWRMS>(
          quantitytypes, rtol, atol, BDFMultiplier, tol, alpha, beta, disable_implicit_weighting);

  return status_test;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::StatusTest::Generic> NOX::Nln::StatusTest::Factory::build_active_set_test(
    Teuchos::ParameterList& p, const ::NOX::Utils& u) const
{
  // ------------------------------------------
  // get the Quantity Type
  // ------------------------------------------
  std::string quantity_type_string = p.get("Quantity Type", "???");
  QuantityType qtype = string_to_quantity_type(quantity_type_string);
  if (qtype == quantity_unknown)
  {
    std::ostringstream msg;
    msg << "The \"Quantity Type\" is a required parameter \n"
        << "and the chosen option \"" << quantity_type_string << "\" is invalid!";
    throw_error("build_active_set_test()", msg.str());
  }

  // ------------------------------------------
  // get the maximal cycle size of the active
  // set, which is tested
  // ------------------------------------------
  int max_cycle_size = p.get<int>("Max Cycle Size", 0);


  Teuchos::RCP<NOX::Nln::StatusTest::ActiveSet> status_test =
      Teuchos::make_rcp<NOX::Nln::StatusTest::ActiveSet>(qtype, max_cycle_size);

  return status_test;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::StatusTest::Generic> NOX::Nln::StatusTest::Factory::build_combo_test(
    Teuchos::ParameterList& p, const ::NOX::Utils& u,
    std::map<std::string, Teuchos::RCP<::NOX::StatusTest::Generic>>* tagged_tests) const
{
  std::string combo_type_string = Teuchos::get<std::string>(p, "Combo Type");
  ::NOX::StatusTest::Combo::ComboType combo_type = ::NOX::StatusTest::Combo::AND;
  if (combo_type_string == "AND")
    combo_type = ::NOX::StatusTest::Combo::AND;
  else if (combo_type_string == "OR")
    combo_type = ::NOX::StatusTest::Combo::OR;
  else
  {
    throw_error("build_combo_test()", "The \"Combo Type\" must be \"AND\" or \"OR\"!");
  }

  Teuchos::RCP<NOX::Nln::StatusTest::Combo> combo_test =
      Teuchos::make_rcp<NOX::Nln::StatusTest::Combo>(combo_type, &u);

  int i = 0;
  std::ostringstream subtest_name;
  subtest_name << "Test " << i;
  while (p.isSublist(subtest_name.str()))
  {
    Teuchos::ParameterList& subtest_list = p.sublist(subtest_name.str(), true);

    Teuchos::RCP<::NOX::StatusTest::Generic> subtest =
        this->build_outer_status_tests(subtest_list, u, tagged_tests);

    combo_test->addStatusTest(subtest);

    // increase iterator
    ++i;
    subtest_name.str("");
    subtest_name << "Test " << i;
  }

  return combo_test;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::Nln::StatusTest::Factory::check_and_tag_test(const Teuchos::ParameterList& p,
    const Teuchos::RCP<::NOX::StatusTest::Generic>& test,
    std::map<std::string, Teuchos::RCP<::NOX::StatusTest::Generic>>* tagged_tests) const
{
  if ((Teuchos::isParameterType<std::string>(p, "Tag")) && (tagged_tests != nullptr))
  {
    (*tagged_tests)[Teuchos::getParameter<std::string>(p, "Tag")] = test;
    return true;
  }

  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::StatusTest::Factory::throw_error(
    const std::string& functionName, const std::string& errorMsg) const
{
  std::ostringstream msg;
  msg << "ERROR - NOX::Nln::StatusTest::Factory::" << functionName << " - " << errorMsg
      << std::endl;
  FOUR_C_THROW("{}", msg.str());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::StatusTest::Generic> NOX::Nln::StatusTest::build_outer_status_tests(
    Teuchos::ParameterList& p, const ::NOX::Utils& u,
    std::map<std::string, Teuchos::RCP<::NOX::StatusTest::Generic>>* tagged_tests)
{
  Factory factory;
  return factory.build_outer_status_tests(p, u, tagged_tests);
}

FOUR_C_NAMESPACE_CLOSE
