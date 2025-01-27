// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_inpar_elch.hpp"

#include "4C_fem_condition_definition.hpp"
#include "4C_inpar_scatra.hpp"
#include "4C_io_input_spec_builders.hpp"
#include "4C_linalg_sparseoperator.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

void Inpar::ElCh::set_valid_parameters(Teuchos::ParameterList& list)
{
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& elchcontrol =
      list.sublist("ELCH CONTROL", false, "control parameters for electrochemistry problems\n");

  Core::Utils::int_parameter("MOVBOUNDARYITEMAX", 10,
      "Maximum number of outer iterations in electrode shape change computations", &elchcontrol);
  Core::Utils::double_parameter("MOVBOUNDARYCONVTOL", 1e-6,
      "Convergence check tolerance for outer loop in electrode shape change computations",
      &elchcontrol);
  Core::Utils::double_parameter(
      "TEMPERATURE", 298.0, "Constant temperature (Kelvin)", &elchcontrol);
  Core::Utils::int_parameter("TEMPERATURE_FROM_FUNCT", -1,
      "Homogeneous temperature within electrochemistry field that can be time dependent according "
      "to function definition",
      &elchcontrol);
  Core::Utils::double_parameter("FARADAY_CONSTANT", 9.64853399e4,
      "Faraday constant (in unit system as chosen in input file)", &elchcontrol);
  Core::Utils::double_parameter("GAS_CONSTANT", 8.314472,
      "(universal) gas constant (in unit system as chosen in input file)", &elchcontrol);
  // parameter for possible types of ELCH algorithms for deforming meshes
  setStringToIntegralParameter<Inpar::ElCh::ElchMovingBoundary>("MOVINGBOUNDARY", "No",
      "ELCH algorithm for deforming meshes",
      tuple<std::string>("No", "pseudo-transient", "fully-transient"),
      tuple<std::string>("no moving boundary algorithm",
          "pseudo-transient moving boundary algorithm",
          "full moving boundary algorithm including fluid solve"),
      tuple<Inpar::ElCh::ElchMovingBoundary>(
          elch_mov_bndry_no, elch_mov_bndry_pseudo_transient, elch_mov_bndry_fully_transient),
      &elchcontrol);
  Core::Utils::double_parameter(
      "MOLARVOLUME", 0.0, "Molar volume for electrode shape change computations", &elchcontrol);
  Core::Utils::double_parameter("MOVBOUNDARYTHETA", 0.0,
      "One-step-theta factor in electrode shape change computations", &elchcontrol);
  Core::Utils::bool_parameter("GALVANOSTATIC", "No", "flag for galvanostatic mode", &elchcontrol);
  setStringToIntegralParameter<Inpar::ElCh::ApproxElectResist>("GSTAT_APPROX_ELECT_RESIST",
      "relation_pot_cur", "relation of potential and current flow",
      tuple<std::string>("relation_pot_cur", "effective_length_with_initial_cond",
          "effective_length_with_integrated_cond"),
      tuple<Inpar::ElCh::ApproxElectResist>(approxelctresist_relpotcur,
          approxelctresist_effleninitcond, approxelctresist_efflenintegcond),
      &elchcontrol);
  Core::Utils::int_parameter(
      "GSTATCONDID_CATHODE", 0, "condition id of electrode kinetics for cathode", &elchcontrol);
  Core::Utils::int_parameter(
      "GSTATCONDID_ANODE", 1, "condition id of electrode kinetics for anode", &elchcontrol);
  Core::Utils::double_parameter(
      "GSTATCONVTOL", 1.e-5, "Convergence check tolerance for galvanostatic mode", &elchcontrol);
  Core::Utils::double_parameter("GSTATCURTOL", 1.e-15, "Current Tolerance", &elchcontrol);
  Core::Utils::int_parameter(
      "GSTATFUNCTNO", -1, "function number defining the imposed current curve", &elchcontrol);
  Core::Utils::int_parameter(
      "GSTATITEMAX", 10, "maximum number of iterations for galvanostatic mode", &elchcontrol);
  Core::Utils::double_parameter(
      "GSTAT_LENGTH_CURRENTPATH", 0.0, "average length of the current path", &elchcontrol);

  setStringToIntegralParameter<Inpar::ElCh::EquPot>("EQUPOT", "Undefined",
      "type of closing equation for electric potential",
      tuple<std::string>(
          "Undefined", "ENC", "ENC_PDE", "ENC_PDE_ELIM", "Poisson", "Laplace", "divi"),
      tuple<Inpar::ElCh::EquPot>(equpot_undefined, equpot_enc, equpot_enc_pde, equpot_enc_pde_elim,
          equpot_poisson, equpot_laplace, equpot_divi),
      &elchcontrol);
  Core::Utils::bool_parameter(
      "DIFFCOND_FORMULATION", "No", "Activation of diffusion-conduction formulation", &elchcontrol);
  Core::Utils::bool_parameter("INITPOTCALC", "No",
      "Automatically calculate initial field for electric potential", &elchcontrol);
  Core::Utils::bool_parameter("ONLYPOTENTIAL", "no",
      "Coupling of general ion transport equation with Laplace equation", &elchcontrol);
  Core::Utils::bool_parameter("COUPLE_BOUNDARY_FLUXES", "Yes",
      "Coupling of lithium-ion flux density and electric current density at Dirichlet and Neumann "
      "boundaries",
      &elchcontrol);
  Core::Utils::double_parameter(
      "CYCLING_TIMESTEP", -1., "modified time step size for CCCV cell cycling", &elchcontrol);
  Core::Utils::bool_parameter("ELECTRODE_INFO_EVERY_STEP", "No",
      "the cell voltage, SOC, and C-Rate will be written to the csv file every step, even if "
      "RESULTSEVERY is not 1",
      &elchcontrol);

  /*----------------------------------------------------------------------*/
  // attention: this list is a sublist of elchcontrol
  Teuchos::ParameterList& elchdiffcondcontrol = elchcontrol.sublist(
      "DIFFCOND", false, "control parameters for electrochemical diffusion conduction problems\n");

  Core::Utils::bool_parameter(
      "CURRENT_SOLUTION_VAR", "No", "Current as a solution variable", &elchdiffcondcontrol);
  Core::Utils::bool_parameter("MAT_DIFFCOND_DIFFBASED", "Yes",
      "Coupling terms of chemical diffusion for current equation are based on t and kappa",
      &elchdiffcondcontrol);

  /// dilute solution theory (diffusion potential in current equation):
  ///    A          B
  ///   |--|  |----------|
  ///   z_1 + (z_2 - z_1) t_1
  /// ------------------------ (RT/F kappa (1+f+-) 1/c_k grad c_k)
  ///      z_1 z_2
  ///     |________|
  ///         C
  //
  // default: concentrated solution theory according to Newman
  Core::Utils::double_parameter("MAT_NEWMAN_CONST_A", 2.0,
      "Constant A for the Newman model(term for the concentration overpotential)",
      &elchdiffcondcontrol);
  Core::Utils::double_parameter("MAT_NEWMAN_CONST_B", -2.0,
      "Constant B for the Newman model(term for the concentration overpotential)",
      &elchdiffcondcontrol);
  Core::Utils::double_parameter("MAT_NEWMAN_CONST_C", -1.0,
      "Constant C for the Newman model(term for the concentration overpotential)",
      &elchdiffcondcontrol);
  Core::Utils::double_parameter(
      "PERMITTIVITY_VACUUM", 8.8541878128e-12, "Vacuum permittivity", &elchdiffcondcontrol);

  /*----------------------------------------------------------------------*/
  // sublist for space-charge layers
  auto& sclcontrol = elchcontrol.sublist(
      "SCL", false, "control parameters for coupled problems with space-charge layer formation\n");

  Core::Utils::bool_parameter(
      "ADD_MICRO_MACRO_COUPLING", "No", "flag for micro macro coupling with scls", &sclcontrol);
  Core::Utils::bool_parameter("COUPLING_OUTPUT", "No",
      "write coupled node gids and node coordinates to csv file", &sclcontrol);
  Core::Utils::bool_parameter(
      "INITPOTCALC", "No", "calculate initial potential field?", &sclcontrol);
  Core::Utils::int_parameter("SOLVER", -1, "solver for coupled SCL problem", &sclcontrol);
  setStringToIntegralParameter<Core::LinAlg::MatrixType>("MATRIXTYPE", "undefined",
      "type of global system matrix in global system of equations",
      tuple<std::string>("undefined", "block", "sparse"),
      tuple<Core::LinAlg::MatrixType>(Core::LinAlg::MatrixType::undefined,
          Core::LinAlg::MatrixType::block_field, Core::LinAlg::MatrixType::sparse),
      &sclcontrol);
  Core::Utils::int_parameter("ADAPT_TIME_STEP", -1,
      "time step when time step size should be updated to 'ADAPTED_TIME_STEP_SIZE'.", &sclcontrol);
  Core::Utils::double_parameter("ADAPTED_TIME_STEP_SIZE", -1.0, "new time step size.", &sclcontrol);

  setStringToIntegralParameter<ScaTra::InitialField>("INITIALFIELD", "zero_field",
      "Initial Field for scalar transport problem",
      tuple<std::string>("zero_field", "field_by_function", "field_by_condition"),
      tuple<ScaTra::InitialField>(ScaTra::initfield_zero_field, ScaTra::initfield_field_by_function,
          ScaTra::initfield_field_by_condition),
      &sclcontrol);

  Core::Utils::int_parameter(
      "INITFUNCNO", -1, "function number for scalar transport initial field", &sclcontrol);
}


void Inpar::ElCh::set_valid_conditions(std::vector<Core::Conditions::ConditionDefinition>& condlist)
{
  using namespace Core::IO::InputSpecBuilders;

  /*--------------------------------------------------------------------*/
  // electrode state of charge
  {
    // definition of electrode state of charge surface and volume conditions
    Core::Conditions::ConditionDefinition electrodesocline(
        "DESIGN ELECTRODE STATE OF CHARGE LINE CONDITIONS", "ElectrodeSOC",
        "electrode state of charge line condition", Core::Conditions::ElectrodeSOC, true,
        Core::Conditions::geometry_type_line);

    Core::Conditions::ConditionDefinition electrodesocsurf(
        "DESIGN ELECTRODE STATE OF CHARGE SURF CONDITIONS", "ElectrodeSOC",
        "electrode state of charge surface condition", Core::Conditions::ElectrodeSOC, true,
        Core::Conditions::geometry_type_surface);
    Core::Conditions::ConditionDefinition electrodesocvol(
        "DESIGN ELECTRODE STATE OF CHARGE VOL CONDITIONS", "ElectrodeSOC",
        "electrode state of charge volume condition", Core::Conditions::ElectrodeSOC, true,
        Core::Conditions::geometry_type_volume);

    const auto make_electrodesoc = [&condlist](Core::Conditions::ConditionDefinition& cond)
    {
      // insert input file line components into condition definitions
      cond.add_component(entry<int>("ConditionID"));
      cond.add_component(entry<double>("C_0%"));
      cond.add_component(entry<double>("C_100%"));
      cond.add_component(entry<double>("ONE_HOUR"));

      // insert condition definitions into global list of valid condition definitions
      condlist.emplace_back(cond);
    };

    make_electrodesoc(electrodesocline);
    make_electrodesoc(electrodesocsurf);
    make_electrodesoc(electrodesocvol);
  }

  /*--------------------------------------------------------------------*/
  // cell voltage

  {
    // definition of cell voltage point, line, and surface conditions
    Core::Conditions::ConditionDefinition cellvoltagepoint("DESIGN CELL VOLTAGE POINT CONDITIONS",
        "CellVoltagePoint", "cell voltage point condition", Core::Conditions::CellVoltage, false,
        Core::Conditions::geometry_type_point);

    Core::Conditions::ConditionDefinition cellvoltageline("DESIGN CELL VOLTAGE LINE CONDITIONS",
        "CellVoltage", "cell voltage line condition", Core::Conditions::CellVoltage, true,
        Core::Conditions::geometry_type_line);

    Core::Conditions::ConditionDefinition cellvoltagesurf("DESIGN CELL VOLTAGE SURF CONDITIONS",
        "CellVoltage", "cell voltage surface condition", Core::Conditions::CellVoltage, true,
        Core::Conditions::geometry_type_surface);

    const auto make_cellvoltage = [&condlist](Core::Conditions::ConditionDefinition& cond)
    {
      // insert input file line components into condition definitions
      cond.add_component(entry<int>("ConditionID"));

      // insert condition definitions into global list of valid condition definitions
      condlist.emplace_back(cond);
    };

    make_cellvoltage(cellvoltagepoint);
    make_cellvoltage(cellvoltageline);
    make_cellvoltage(cellvoltagesurf);
  }


  using namespace Core::IO::InputSpecBuilders;

  /*--------------------------------------------------------------------*/
  // electrode kinetics as boundary condition on electrolyte
  {
    auto reaction_model_choices = one_of({
        all_of({
            selection<int>("KINETIC_MODEL",
                {
                    {"Butler-Volmer", Inpar::ElCh::ElectrodeKinetics::butler_volmer},
                    {"Butler-Volmer-Yang1997",
                        Inpar::ElCh::ElectrodeKinetics::butler_volmer_yang1997},
                }),
            entry<double>("ALPHA_A"),
            entry<double>("ALPHA_C"),
            entry<double>("I0"),
            entry<double>("GAMMA"),
            entry<double>("REFCON"),
            entry<double>("DL_SPEC_CAP"),
        }),
        all_of({
            selection<int>("KINETIC_MODEL", {{"Tafel", Inpar::ElCh::ElectrodeKinetics::tafel}}),
            entry<double>("ALPHA"),
            entry<double>("I0"),
            entry<double>("GAMMA"),
            entry<double>("REFCON"),
            entry<double>("DL_SPEC_CAP"),
        }),
        all_of({
            selection<int>("KINETIC_MODEL", {{"linear", Inpar::ElCh::ElectrodeKinetics::linear}}),
            entry<double>("ALPHA"),
            entry<double>("I0"),
            entry<double>("GAMMA"),
            entry<double>("REFCON"),
            entry<double>("DL_SPEC_CAP"),
        }),
        all_of({
            selection<int>("KINETIC_MODEL",
                {{"Butler-Volmer-Newman", Inpar::ElCh::ElectrodeKinetics::butler_volmer_newman}}),
            entry<double>("K_A"),
            entry<double>("K_C"),
            entry<double>("BETA"),
            entry<double>("DL_SPEC_CAP"),
        }),
        all_of({
            selection<int>("KINETIC_MODEL",
                {{"Butler-Volmer-Bard", Inpar::ElCh::ElectrodeKinetics::butler_volmer_bard}}),
            entry<double>("E0"),
            entry<double>("K0"),
            entry<double>("BETA"),
            entry<double>("C_C0"),
            entry<double>("C_A0"),
            entry<double>("DL_SPEC_CAP"),
        }),
        all_of({
            selection<int>("KINETIC_MODEL", {{"Nernst", Inpar::ElCh::ElectrodeKinetics::nernst}}),
            entry<double>("E0"),
            entry<double>("C0"),
            entry<double>("DL_SPEC_CAP"),
        }),
    });


    Core::Conditions::ConditionDefinition electrodeboundarykineticspoint(
        "ELECTRODE BOUNDARY KINETICS POINT CONDITIONS", "ElchBoundaryKineticsPoint",
        "point electrode boundary kinetics", Core::Conditions::ElchBoundaryKinetics, false,
        Core::Conditions::geometry_type_point);

    Core::Conditions::ConditionDefinition electrodeboundarykineticsline(
        "ELECTRODE BOUNDARY KINETICS LINE CONDITIONS", "ElchBoundaryKinetics",
        "line electrode boundary kinetics", Core::Conditions::ElchBoundaryKinetics, true,
        Core::Conditions::geometry_type_line);

    Core::Conditions::ConditionDefinition electrodeboundarykineticssurf(
        "ELECTRODE BOUNDARY KINETICS SURF CONDITIONS", "ElchBoundaryKinetics",
        "surface electrode boundary kinetics", Core::Conditions::ElchBoundaryKinetics, true,
        Core::Conditions::geometry_type_surface);

    const auto make_electrodeboundarykinetics = [&condlist, &reaction_model_choices](
                                                    Core::Conditions::ConditionDefinition& cond)
    {
      cond.add_component(entry<int>("ConditionID"));
      cond.add_component(entry<double>("POT"));
      cond.add_component(entry<Noneable<int>>("FUNCT", {.description = ""}));
      cond.add_component(entry<int>("NUMSCAL"));
      cond.add_component(entry<std::vector<int>>(
          "STOICH", {.description = "", .size = from_parameter<int>("NUMSCAL")}));
      cond.add_component(entry<int>("E-"));
      cond.add_component(
          entry<double>("EPSILON", {.description = "porosity of electrode boundary, set to -1 if "
                                                   "equal to porosity of electrolyte domain"}));
      cond.add_component(entry<int>("ZERO_CUR"));
      cond.add_component(reaction_model_choices);
      condlist.emplace_back(cond);
    };

    make_electrodeboundarykinetics(electrodeboundarykineticspoint);
    make_electrodeboundarykinetics(electrodeboundarykineticsline);
    make_electrodeboundarykinetics(electrodeboundarykineticssurf);
  }

  /*--------------------------------------------------------------------*/
  // electrode kinetics as domain condition within electrolyte
  {
    // definition of line, surface, and volume conditions for electrode domain kinetics
    Core::Conditions::ConditionDefinition electrodedomainkineticsline(
        "ELECTRODE DOMAIN KINETICS LINE CONDITIONS", "ElchDomainKinetics",
        "line electrode domain kinetics", Core::Conditions::ElchDomainKinetics, true,
        Core::Conditions::geometry_type_line);

    Core::Conditions::ConditionDefinition electrodedomainkineticssurf(
        "ELECTRODE DOMAIN KINETICS SURF CONDITIONS", "ElchDomainKinetics",
        "surface electrode domain kinetics", Core::Conditions::ElchDomainKinetics, true,
        Core::Conditions::geometry_type_surface);

    Core::Conditions::ConditionDefinition electrodedomainkineticsvol(
        "ELECTRODE DOMAIN KINETICS VOL CONDITIONS", "ElchDomainKinetics",
        "volume electrode domain kinetics", Core::Conditions::ElchDomainKinetics, true,
        Core::Conditions::geometry_type_volume);

    // equip condition definition with input file line components
    auto electrodedomainkineticscomponents = all_of({
        entry<int>("ConditionID"),
        entry<double>("POT"),
        entry<Core::IO::Noneable<int>>("FUNCT"),
        entry<int>("NUMSCAL"),
        entry<std::vector<int>>("STOICH", {.size = from_parameter<int>("NUMSCAL")}),
        entry<int>("E-"),
        entry<int>("ZERO_CUR"),
        selection<int>("KINETIC_MODEL",
            {
                {"Butler-Volmer", Inpar::ElCh::ElectrodeKinetics::butler_volmer},
            }),
        entry<double>("A_S"),
        entry<double>("ALPHA_A"),
        entry<double>("ALPHA_C"),
        entry<double>("I0"),
        entry<double>("GAMMA"),
        entry<double>("REFCON"),
        entry<double>("DL_SPEC_CAP"),
    });

    {
      electrodedomainkineticsline.add_component(electrodedomainkineticscomponents);
      electrodedomainkineticssurf.add_component(electrodedomainkineticscomponents);
      electrodedomainkineticsvol.add_component(electrodedomainkineticscomponents);
    }

    // insert condition definitions into global list of valid condition definitions
    condlist.emplace_back(electrodedomainkineticsline);
    condlist.emplace_back(electrodedomainkineticssurf);
    condlist.emplace_back(electrodedomainkineticsvol);
  }

  /*--------------------------------------------------------------------*/
  // boundary condition for constant-current constant-voltage (CCCV) cell cycling

  // definition of point, line and surface conditions for CCCV cell cycling
  Core::Conditions::ConditionDefinition cccvcyclingpoint(
      "DESIGN CCCV CELL CYCLING POINT CONDITIONS", "CCCVCycling",
      "line boundary condition for constant-current constant-voltage (CCCV) cell cycling",
      Core::Conditions::CCCVCycling, true, Core::Conditions::geometry_type_point);

  Core::Conditions::ConditionDefinition cccvcyclingline("DESIGN CCCV CELL CYCLING LINE CONDITIONS",
      "CCCVCycling",
      "line boundary condition for constant-current constant-voltage (CCCV) cell cycling",
      Core::Conditions::CCCVCycling, true, Core::Conditions::geometry_type_line);

  Core::Conditions::ConditionDefinition cccvcyclingsurf("DESIGN CCCV CELL CYCLING SURF CONDITIONS",
      "CCCVCycling",
      "surface boundary condition for constant-current constant-voltage (CCCV) cell cycling",
      Core::Conditions::CCCVCycling, true, Core::Conditions::geometry_type_surface);

  const auto make_cccvcycling = [&condlist](Core::Conditions::ConditionDefinition& cond)
  {
    // insert input file line components into condition definitions
    {
      cond.add_component(entry<int>("NUMBER_OF_HALF_CYCLES"));
      cond.add_component(
          entry<int>("BEGIN_WITH_CHARGING"));  // Boolean parameter represented by integer parameter
      cond.add_component(entry<Noneable<int>>("CONDITION_ID_FOR_CHARGE", {.description = ""}));
      cond.add_component(entry<Noneable<int>>("CONDITION_ID_FOR_DISCHARGE", {.description = ""}));
      cond.add_component(entry<double>("INIT_RELAX_TIME"));
      cond.add_component(entry<int>("ADAPTIVE_TIME_STEPPING_INIT_RELAX"));
      cond.add_component(entry<Noneable<int>>("NUM_ADD_ADAPT_TIME_STEPS", {.description = ""}));
      cond.add_component(
          entry<Noneable<int>>("MIN_TIME_STEPS_DURING_INIT_RELAX", {.description = ""}));

      // insert condition definitions into global list of valid condition definitions
      condlist.emplace_back(cond);
    }
  };

  make_cccvcycling(cccvcyclingpoint);
  make_cccvcycling(cccvcyclingline);
  make_cccvcycling(cccvcyclingsurf);

  /*--------------------------------------------------------------------*/
  // boundary condition for constant-current constant-voltage (CCCV) half-cycle

  // definition of point, line and surface conditions for CCCV half-cycle
  Core::Conditions::ConditionDefinition cccvhalfcyclepoint(
      "DESIGN CCCV HALF-CYCLE POINT CONDITIONS", "CCCVHalfCycle",
      "line boundary condition for constant-current constant-voltage (CCCV) half-cycle",
      Core::Conditions::CCCVHalfCycle, true, Core::Conditions::geometry_type_point);

  Core::Conditions::ConditionDefinition cccvhalfcycleline("DESIGN CCCV HALF-CYCLE LINE CONDITIONS",
      "CCCVHalfCycle",
      "line boundary condition for constant-current constant-voltage (CCCV) half-cycle",
      Core::Conditions::CCCVHalfCycle, true, Core::Conditions::geometry_type_line);

  Core::Conditions::ConditionDefinition cccvhalfcyclesurf("DESIGN CCCV HALF-CYCLE SURF CONDITIONS",
      "CCCVHalfCycle",
      "surface boundary condition for constant-current constant-voltage (CCCV) half-cycle",
      Core::Conditions::CCCVHalfCycle, true, Core::Conditions::geometry_type_surface);

  const auto make_cccvhalfcycle = [&condlist](Core::Conditions::ConditionDefinition& cond)
  {
    // insert input file line components into condition definitions
    cond.add_component(entry<int>("ConditionID"));
    cond.add_component(entry<double>("CURRENT"));
    cond.add_component(entry<double>("CUT_OFF_VOLTAGE"));
    cond.add_component(entry<double>("CUT_OFF_C_RATE"));
    cond.add_component(entry<double>("RELAX_TIME"));
    // switch adaptive time stepping on for different phases of half cycle: 1st: end of constant
    // current, 2nd: end of constant voltage, 3rd: end of relaxation
    cond.add_component(entry<std::vector<int>>(
        "ADAPTIVE_TIME_STEPPING_PHASE_ON_OFF", {.description = "", .size = 3}));

    // insert condition definitions into global list of valid condition definitions
    condlist.emplace_back(cond);
  };

  make_cccvhalfcycle(cccvhalfcyclepoint);
  make_cccvhalfcycle(cccvhalfcycleline);
  make_cccvhalfcycle(cccvhalfcyclesurf);
}

FOUR_C_NAMESPACE_CLOSE
