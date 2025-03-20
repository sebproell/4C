// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_inpar_fluid.hpp"

#include "4C_fem_condition_definition.hpp"
#include "4C_io_geometry_type.hpp"
#include "4C_io_input_spec_builders.hpp"
#include "4C_legacy_enum_definitions_conditions.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN



void Inpar::FLUID::set_valid_parameters(std::map<std::string, Core::IO::InputSpec>& list)
{
  using Teuchos::tuple;
  using namespace Core::IO::InputSpecBuilders;

  Core::Utils::SectionSpecs fdyn{"FLUID DYNAMIC"};

  // physical type of fluid flow (incompressible, varying density, loma, Boussinesq approximation,
  // temperature-dependent water)
  fdyn.specs.emplace_back(deprecated_selection<Inpar::FLUID::PhysicalType>("PHYSICAL_TYPE",
      {
          {"Incompressible", incompressible},
          {"Weakly_compressible", weakly_compressible},
          {"Weakly_compressible_stokes", weakly_compressible_stokes},
          {"Weakly_compressible_dens_mom", weakly_compressible_dens_mom},
          {"Weakly_compressible_stokes_dens_mom", weakly_compressible_stokes_dens_mom},
          {"Artificial_compressibility", artcomp},
          {"Varying_density", varying_density},
          {"Loma", loma},
          {"Temp_dep_water", tempdepwater},
          {"Boussinesq", boussinesq},
          {"Stokes", stokes},
          {"Oseen", oseen},
      },
      {.description = "Physical Type", .default_value = incompressible}));

  // number of linear solver used for fluid problem
  fdyn.specs.emplace_back(parameter<int>("LINEAR_SOLVER",
      {.description = "number of linear solver used for fluid dynamics", .default_value = -1}));

  // number of linear solver used for fluid problem (former fluid pressure solver for SIMPLER
  // preconditioning with fluid)
  fdyn.specs.emplace_back(parameter<int>(
      "SIMPLER_SOLVER", {.description = "number of linear solver used for fluid dynamic (ONLY "
                                        "NECESSARY FOR BlockGaussSeidel solver block within fluid "
                                        "mehstying case any more!!!!)",
                            .default_value = -1}));

  // Flag to define the way of calculating stresses and wss
  fdyn.specs.emplace_back(deprecated_selection<Inpar::FLUID::WSSType>("WSS_TYPE",
      {
          {"Standard", wss_standard},
          {"Aggregation", wss_aggregation},
          {"Mean", wss_mean},
      },
      {.description = "which type of stresses and wall shear stress",
          .default_value = wss_standard}));

  // Set ML-solver number for smoothing of residual-based calculated wallshearstress via plain
  // aggregation.
  fdyn.specs.emplace_back(parameter<int>("WSS_ML_AGR_SOLVER",
      {.description = "Set ML-solver number for smoothing of residual-based calculated "
                      "wallshearstress via plain aggregation.",
          .default_value = -1}));

  fdyn.specs.emplace_back(deprecated_selection<Inpar::FLUID::TimeIntegrationScheme>("TIMEINTEGR",
      {
          {"Stationary", timeint_stationary},
          {"Np_Gen_Alpha", timeint_npgenalpha},
          {"Af_Gen_Alpha", timeint_afgenalpha},
          {"One_Step_Theta", timeint_one_step_theta},
          {"BDF2", timeint_bdf2},
      },
      {.description = "Time Integration Scheme", .default_value = timeint_one_step_theta}));

  fdyn.specs.emplace_back(parameter<Inpar::FLUID::OstContAndPress>("OST_CONT_PRESS",
      {.description =
              "One step theta option for time discretization of continuity eq. and pressure",
          .default_value = Cont_normal_Press_normal}));

  fdyn.specs.emplace_back(deprecated_selection<Core::IO::GeometryType>("GEOMETRY",
      {
          {"full", Core::IO::geometry_full},
          {"box", Core::IO::geometry_box},
          {"file", Core::IO::geometry_file},
      },
      {.description = "How the geometry is specified", .default_value = Core::IO::geometry_full}));

  fdyn.specs.emplace_back(parameter<Inpar::FLUID::LinearisationAction>("NONLINITER",
      {.description = "Nonlinear iteration scheme", .default_value = fixed_point_like}));

  std::vector<std::string> predictor_valid_input = {"steady_state", "zero_acceleration",
      "constant_acceleration", "constant_increment", "explicit_second_order_midpoint", "TangVel"};
  fdyn.specs.emplace_back(deprecated_selection<std::string>("PREDICTOR", predictor_valid_input,
      {.description = "Predictor for first guess in nonlinear iteration",
          .default_value = "steady_state"}));


  fdyn.specs.emplace_back(deprecated_selection<Inpar::FLUID::ItNorm>("CONVCHECK",
      {
          {"L_2_norm", fncc_L2},
      },
      {.description = "Norm for convergence check (relative increments and absolute residuals)",
          .default_value = fncc_L2}));

  fdyn.specs.emplace_back(parameter<bool>("INCONSISTENT_RESIDUAL",
      {.description = "do not evaluate residual after solution has converged (->faster)",
          .default_value = false}));

  {
    std::map<std::string, InitialField> map{
        {"zero_field", initfield_zero_field},
        {"field_by_function", initfield_field_by_function},
        {"disturbed_field_from_function", initfield_disturbed_field_from_function},
        {"FLAME_VORTEX_INTERACTION", initfield_flame_vortex_interaction},
        {"BELTRAMI-FLOW", initfield_beltrami_flow},
        {"KIM-MOIN-FLOW", initfield_kim_moin_flow},
        {"hit_comte_bellot_corrsin_initial_field", initfield_hit_comte_bellot_corrsin},
        {"forced_hit_simple_algebraic_spectrum", initfield_forced_hit_simple_algebraic_spectrum},
        {"forced_hit_numeric_spectrum", initfield_forced_hit_numeric_spectrum},
        {"forced_hit_passive", initfield_passive_hit_const_input},
        {"channel_weakly_compressible", initfield_channel_weakly_compressible},
    };

    fdyn.specs.emplace_back(deprecated_selection<Inpar::FLUID::InitialField>("INITIALFIELD", map,
        {.description = "Initial Field for transport problem",
            .default_value = initfield_zero_field}));
  }

  fdyn.specs.emplace_back(parameter<int>("OSEENFIELDFUNCNO",
      {.description = "function number of Oseen advective field", .default_value = -1}));

  fdyn.specs.emplace_back(parameter<bool>(
      "LIFTDRAG", {.description = "Calculate lift and drag forces along specified boundary",
                      .default_value = false}));

  std::vector<std::string> convform_valid_input = {"convective", "conservative"};
  fdyn.specs.emplace_back(deprecated_selection<std::string>("CONVFORM", convform_valid_input,
      {.description = "form of convective term", .default_value = "convective"}));

  fdyn.specs.emplace_back(parameter<bool>("NONLINEARBC",
      {.description = "Flag to activate check for potential nonlinear boundary conditions",
          .default_value = false}));

  fdyn.specs.emplace_back(deprecated_selection<Inpar::FLUID::MeshTying>("MESHTYING",
      {
          {"no", no_meshtying},
          {"Condensed_Smat", condensed_smat},
          {"Condensed_Bmat", condensed_bmat},
          {"Condensed_Bmat_merged", condensed_bmat_merged},
      },
      {.description = "Flag to (de)activate mesh tying algorithm", .default_value = no_meshtying}));

  fdyn.specs.emplace_back(parameter<Inpar::FLUID::Gridvel>(
      "GRIDVEL", {.description = "scheme for determination of gridvelocity from displacements",
                     .default_value = BE}));

  fdyn.specs.emplace_back(parameter<bool>("ALLDOFCOUPLED",
      {.description = "all dof (incl. pressure) are coupled", .default_value = true}));

  fdyn.specs.emplace_back(parameter<CalcError>("CALCERROR",
      {.description = "Flag to (de)activate error calculations", .default_value = no}));
  fdyn.specs.emplace_back(parameter<int>(
      "CALCERRORFUNCNO", {.description = "Function for Error Calculation", .default_value = -1}));

  fdyn.specs.emplace_back(parameter<int>("CORRTERMFUNCNO",
      {.description =
              "Function for calculation of the correction term for the weakly compressible problem",
          .default_value = -1}));

  fdyn.specs.emplace_back(parameter<int>("BODYFORCEFUNCNO",
      {.description =
              "Function for calculation of the body force for the weakly compressible problem",
          .default_value = -1}));

  fdyn.specs.emplace_back(parameter<double>(
      "STAB_DEN_REF", {.description = "Reference stabilization parameter for the density for the "
                                      "HDG weakly compressible formulation",
                          .default_value = 0.0}));

  fdyn.specs.emplace_back(parameter<double>(
      "STAB_MOM_REF", {.description = "Reference stabilization parameter for the momentum for the "
                                      "HDG weakly compressible formulation",
                          .default_value = 0.0}));

  fdyn.specs.emplace_back(parameter<int>(
      "VARVISCFUNCNO", {.description = "Function for calculation of a variable viscosity for the "
                                       "weakly compressible problem",
                           .default_value = -1}));

  fdyn.specs.emplace_back(parameter<bool>(
      "PRESSAVGBC", {.description = "Flag to (de)activate imposition of boundary condition for the "
                                    "considered element average pressure",
                        .default_value = false}));

  fdyn.specs.emplace_back(
      parameter<double>("REFMACH", {.description = "Reference Mach number", .default_value = 1.0}));

  fdyn.specs.emplace_back(parameter<bool>("BLOCKMATRIX",
      {.description =
              "Indicates if system matrix should be assembled into a sparse block matrix type.",
          .default_value = false}));

  fdyn.specs.emplace_back(parameter<bool>("ADAPTCONV",
      {.description =
              "Switch on adaptive control of linear solver tolerance for nonlinear solution",
          .default_value = false}));
  fdyn.specs.emplace_back(parameter<double>("ADAPTCONV_BETTER",
      {.description = "The linear solver shall be this much better than the current nonlinear "
                      "residual in the nonlinear convergence limit",
          .default_value = 0.1}));

  fdyn.specs.emplace_back(parameter<bool>("INFNORMSCALING",
      {.description = "Scale blocks of matrix with row infnorm?", .default_value = false}));

  fdyn.specs.emplace_back(parameter<bool>(
      "GMSH_OUTPUT", {.description = "write output to gmsh files", .default_value = false}));
  fdyn.specs.emplace_back(parameter<bool>(
      "COMPUTE_DIVU", {.description = "Compute divergence of velocity field at the element center",
                          .default_value = false}));
  fdyn.specs.emplace_back(parameter<bool>("COMPUTE_EKIN",
      {.description = "Compute kinetic energy at the end of each time step and write it to file.",
          .default_value = false}));
  fdyn.specs.emplace_back(parameter<bool>("NEW_OST",
      {.description = "Solve the Navier-Stokes equation with the new One Step Theta algorithm",
          .default_value = false}));  // TODO: To be removed.
  fdyn.specs.emplace_back(parameter<int>(
      "RESULTSEVERY", {.description = "Increment for writing solution", .default_value = 1}));
  fdyn.specs.emplace_back(parameter<int>(
      "RESTARTEVERY", {.description = "Increment for writing restart", .default_value = 20}));
  fdyn.specs.emplace_back(
      parameter<int>("NUMSTEP", {.description = "Total number of Timesteps", .default_value = 1}));
  fdyn.specs.emplace_back(parameter<int>(
      "STEADYSTEP", {.description = "steady state check every step", .default_value = -1}));
  fdyn.specs.emplace_back(parameter<int>(
      "NUMSTASTEPS", {.description = "Number of Steps for Starting Scheme", .default_value = 0}));
  fdyn.specs.emplace_back(parameter<int>(
      "STARTFUNCNO", {.description = "Function for Initial Starting Field", .default_value = -1}));
  fdyn.specs.emplace_back(parameter<int>(
      "ITEMAX", {.description = "max. number of nonlin. iterations", .default_value = 10}));
  fdyn.specs.emplace_back(parameter<int>("INITSTATITEMAX",
      {.description = "max number of nonlinear iterations for initial stationary solution",
          .default_value = 5}));
  fdyn.specs.emplace_back(
      parameter<double>("TIMESTEP", {.description = "Time increment dt", .default_value = 0.01}));
  fdyn.specs.emplace_back(parameter<double>(
      "MAXTIME", {.description = "Total simulation time", .default_value = 1000.0}));
  fdyn.specs.emplace_back(parameter<double>(
      "ALPHA_M", {.description = "Time integration factor", .default_value = 1.0}));
  fdyn.specs.emplace_back(parameter<double>(
      "ALPHA_F", {.description = "Time integration factor", .default_value = 1.0}));
  fdyn.specs.emplace_back(
      parameter<double>("GAMMA", {.description = "Time integration factor", .default_value = 1.0}));
  fdyn.specs.emplace_back(parameter<double>(
      "THETA", {.description = "Time integration factor", .default_value = 0.66}));

  fdyn.specs.emplace_back(parameter<double>("START_THETA",
      {.description = "Time integration factor for starting scheme", .default_value = 1.0}));

  fdyn.specs.emplace_back(parameter<bool>("STRONG_REDD_3D_COUPLING_TYPE",
      {.description = "Flag to (de)activate potential Strong 3D redD coupling",
          .default_value = false}));

  fdyn.specs.emplace_back(parameter<int>("VELGRAD_PROJ_SOLVER",
      {.description = "Number of linear solver used for L2 projection", .default_value = -1}));

  fdyn.specs.emplace_back(
      deprecated_selection<Inpar::FLUID::GradientReconstructionMethod>("VELGRAD_PROJ_METHOD",
          {
              {"none", gradreco_none},
              {"superconvergent_patch_recovery", gradreco_spr},
              {"L2_projection", gradreco_l2},
          },
          {.description = "Flag to (de)activate gradient reconstruction.",
              .default_value = gradreco_none}));

  fdyn.specs.emplace_back(parameter<bool>(
      "OFF_PROC_ASSEMBLY", {.description = "Do not evaluate ghosted elements but communicate them "
                                           "--> faster if element call is expensive",
                               .default_value = false}));

  fdyn.move_into_collection(list);

  /*----------------------------------------------------------------------*/
  Core::Utils::SectionSpecs fdyn_nln{fdyn, "NONLINEAR SOLVER TOLERANCES"};

  fdyn_nln.specs.emplace_back(parameter<double>(
      "TOL_VEL_RES", {.description = "Tolerance for convergence check of velocity residual",
                         .default_value = 1e-6}));

  fdyn_nln.specs.emplace_back(parameter<double>(
      "TOL_VEL_INC", {.description = "Tolerance for convergence check of velocity increment",
                         .default_value = 1e-6}));

  fdyn_nln.specs.emplace_back(parameter<double>(
      "TOL_PRES_RES", {.description = "Tolerance for convergence check of pressure residual",
                          .default_value = 1e-6}));

  fdyn_nln.specs.emplace_back(parameter<double>(
      "TOL_PRES_INC", {.description = "Tolerance for convergence check of pressure increment",
                          .default_value = 1e-6}));

  fdyn_nln.move_into_collection(list);

  /*----------------------------------------------------------------------*/
  Core::Utils::SectionSpecs fdyn_stab{fdyn, "RESIDUAL-BASED STABILIZATION"};
  // this parameter defines various stabilized methods
  fdyn_stab.specs.emplace_back(deprecated_selection<Inpar::FLUID::StabType>("STABTYPE",
      {
          {"no_stabilization", stabtype_nostab},
          {"residual_based", stabtype_residualbased},
          {"edge_based", stabtype_edgebased},
          {"pressure_projection", stabtype_pressureprojection},
      },
      {.description = "Apply (un)stabilized fluid formulation. No stabilization is only possible "
                      "for inf-sup stable elements.Use a residual-based stabilization or, more "
                      "generally, a stabilization \nbased on the concept of the residual-based "
                      "variational multiscale method...\nExpecting additional input. Use an "
                      "edge-based stabilization, especially for XFEM. Alternative: Element/cell "
                      "based polynomial pressure projection, see Dohrmann/Bochev 2004, IJNMF",
          .default_value = stabtype_residualbased}));

  fdyn_stab.specs.emplace_back(parameter<bool>(
      "INCONSISTENT", {.description = "residual based without second derivatives (i.e. only "
                                      "consistent for tau->0, but faster)",
                          .default_value = false}));

  fdyn_stab.specs.emplace_back(parameter<bool>("Reconstruct_Sec_Der",
      {.description = "residual computed with a reconstruction of the second derivatives via "
                      "projection or superconvergent patch recovery",
          .default_value = false}));

  // the following parameters are necessary only if a residual based stabilized method is applied
  fdyn_stab.specs.emplace_back(deprecated_selection<SubscalesTD>("TDS",
      {
          {"quasistatic", subscales_quasistatic},
          {"time_dependent", subscales_time_dependent},
      },
      {.description =
              "Flag to allow time dependency of subscales for residual-based stabilization.",
          .default_value = subscales_quasistatic}));

  fdyn_stab.specs.emplace_back(deprecated_selection<Transient>("TRANSIENT",
      {
          {"no_transient", inertia_stab_drop},
          {"yes_transient", inertia_stab_keep},
          {"transient_complete", inertia_stab_keep_complete},
      },
      {.description =
              "Specify how to treat the transient term. Use transient term (recommended for time "
              "dependent subscales) or use transient term including a linearisation of 1/tau",
          .default_value = inertia_stab_drop}));

  fdyn_stab.specs.emplace_back(parameter<bool>(
      "PSPG", {.description = "Flag to (de)activate PSPG stabilization.", .default_value = true}));
  fdyn_stab.specs.emplace_back(parameter<bool>(
      "SUPG", {.description = "Flag to (de)activate SUPG stabilization.", .default_value = true}));
  fdyn_stab.specs.emplace_back(parameter<bool>(
      "GRAD_DIV", {.description = "Flag to (de)activate grad-div term.", .default_value = true}));

  fdyn_stab.specs.emplace_back(deprecated_selection<VStab>("VSTAB",
      {
          {"no_vstab", viscous_stab_none},
          {"vstab_gls", viscous_stab_gls},
          {"vstab_gls_rhs", viscous_stab_gls_only_rhs},
          {"vstab_usfem", viscous_stab_usfem},
          {"vstab_usfem_rhs", viscous_stab_usfem_only_rhs},
      },
      {.description = "Flag to (de)activate viscous term in residual-based stabilization. Options: "
                      "No viscous term in stabilization, or, Viscous stabilization of GLS type, "
                      "or, Viscous stabilization of GLS type, included only on the right hand "
                      "side, or, Viscous stabilization of USFEM type, or, Viscous stabilization of "
                      "USFEM type, included only on the right hand side",
          .default_value = viscous_stab_none}));

  fdyn_stab.specs.emplace_back(deprecated_selection<RStab>("RSTAB",
      {
          {"no_rstab", reactive_stab_none},
          {"rstab_gls", reactive_stab_gls},
          {"rstab_usfem", reactive_stab_usfem},
      },
      {.description = "Flag to (de)activate reactive term in residual-based stabilization.",
          .default_value = reactive_stab_none}));

  fdyn_stab.specs.emplace_back(deprecated_selection<CrossStress>("CROSS-STRESS",
      {
          {"no_cross", cross_stress_stab_none},
          {"yes_cross", cross_stress_stab},
          {"cross_rhs", cross_stress_stab_only_rhs},
      },
      {.description = "Flag to (de)activate cross-stress term -> residual-based VMM. Options:No "
                      "cross-stress term, or,Include the cross-stress term with a linearization of "
                      "the convective part, or, Include cross-stress term, but only explicitly on "
                      "right hand side",
          .default_value = cross_stress_stab_none}));

  fdyn_stab.specs.emplace_back(deprecated_selection<ReynoldsStress>("REYNOLDS-STRESS",
      {
          {"no_reynolds", reynolds_stress_stab_none},
          {"yes_reynolds", reynolds_stress_stab},
          {"reynolds_rhs", reynolds_stress_stab_only_rhs},
      },
      {.description = "Flag to (de)activate Reynolds-stress term -> residual-based VMM.",
          .default_value = reynolds_stress_stab_none}));

  {
    std::map<std::string, TauType> map{
        {"Taylor_Hughes_Zarins", tau_taylor_hughes_zarins},
        {"Taylor_Hughes_Zarins_wo_dt", tau_taylor_hughes_zarins_wo_dt},
        {"Taylor_Hughes_Zarins_Whiting_Jansen", tau_taylor_hughes_zarins_whiting_jansen},
        {"Taylor_Hughes_Zarins_Whiting_Jansen_wo_dt",
            tau_taylor_hughes_zarins_whiting_jansen_wo_dt},
        {"Taylor_Hughes_Zarins_scaled", tau_taylor_hughes_zarins_scaled},
        {"Taylor_Hughes_Zarins_scaled_wo_dt", tau_taylor_hughes_zarins_scaled_wo_dt},
        {"Franca_Barrenechea_Valentin_Frey_Wall", tau_franca_barrenechea_valentin_frey_wall},
        {"Franca_Barrenechea_Valentin_Frey_Wall_wo_dt",
            tau_franca_barrenechea_valentin_frey_wall_wo_dt},
        {"Shakib_Hughes_Codina", tau_shakib_hughes_codina},
        {"Shakib_Hughes_Codina_wo_dt", tau_shakib_hughes_codina_wo_dt},
        {"Codina", tau_codina},
        {"Codina_wo_dt", tau_codina_wo_dt},
        {"Codina_convscaled", tau_codina_convscaled},
        {"Franca_Madureira_Valentin_Badia_Codina", tau_franca_madureira_valentin_badia_codina},
        {"Franca_Madureira_Valentin_Badia_Codina_wo_dt",
            tau_franca_madureira_valentin_badia_codina_wo_dt},
        {"Hughes_Franca_Balestra_wo_dt", tau_hughes_franca_balestra_wo_dt},
    };

    fdyn_stab.specs.emplace_back(deprecated_selection<TauType>("DEFINITION_TAU", map,
        {.description = "Definition of tau",
            .default_value = tau_franca_barrenechea_valentin_frey_wall}));
  }

  // this parameter selects the characteristic element length for tau_Mu for all
  // stabilization parameter definitions requiring such a length
  fdyn_stab.specs.emplace_back(deprecated_selection<CharEleLengthU>("CHARELELENGTH_U",
      {
          {"streamlength", streamlength_u},
          {"volume_equivalent_diameter", volume_equivalent_diameter_u},
          {"root_of_volume", root_of_volume_u},
      },
      {.description = "Characteristic element length for tau_Mu",
          .default_value = streamlength_u}));

  // this parameter selects the characteristic element length for tau_Mp and tau_C for
  // all stabilization parameter definitions requiring such a length
  fdyn_stab.specs.emplace_back(deprecated_selection<CharEleLengthPC>("CHARELELENGTH_PC",
      {
          {"streamlength", streamlength_pc},
          {"volume_equivalent_diameter", volume_equivalent_diameter_pc},
          {"root_of_volume", root_of_volume_pc},
      },
      {.description = "Characteristic element length for tau_Mp/tau_C",
          .default_value = volume_equivalent_diameter_pc}));

  // this parameter selects the location where tau is evaluated

  std::vector<std::string> evaluation_tau_valid_input = {"element_center", "integration_point"};
  fdyn_stab.specs.emplace_back(
      deprecated_selection<std::string>("EVALUATION_TAU", evaluation_tau_valid_input,
          {.description = "Location where tau is evaluated", .default_value = "element_center"}));

  // this parameter selects the location where the material law is evaluated
  // (does not fit here very well, but parameter transfer is easier)

  std::vector<std::string> evaluation_mat_valid_input = {"element_center", "integration_point"};
  fdyn_stab.specs.emplace_back(
      deprecated_selection<std::string>("EVALUATION_MAT", evaluation_mat_valid_input,
          {.description = "Location where material law is evaluated",
              .default_value = "element_center"}));

  // these parameters active additional terms in loma continuity equation
  // which might be identified as SUPG-/cross- and Reynolds-stress term
  fdyn_stab.specs.emplace_back(parameter<bool>("LOMA_CONTI_SUPG",
      {.description = "Flag to (de)activate SUPG stabilization in loma continuity equation.",
          .default_value = false}));

  fdyn_stab.specs.emplace_back(deprecated_selection<CrossStress>("LOMA_CONTI_CROSS_STRESS",
      {
          {"no_cross", cross_stress_stab_none},
          {"yes_cross", cross_stress_stab},
          {"cross_rhs", cross_stress_stab_only_rhs},
      },
      {.description = "Flag to (de)activate cross-stress term loma continuity equation-> "
                      "residual-based VMM.",
          .default_value = cross_stress_stab_none}));

  fdyn_stab.specs.emplace_back(deprecated_selection<ReynoldsStress>("LOMA_CONTI_REYNOLDS_STRESS",
      {
          {"no_reynolds", reynolds_stress_stab_none},
          {"yes_reynolds", reynolds_stress_stab},
          {"reynolds_rhs", reynolds_stress_stab_only_rhs},
      },
      {.description = "Flag to (de)activate Reynolds-stress term loma continuity equation-> "
                      "residual-based VMM.",
          .default_value = reynolds_stress_stab_none}));

  fdyn_stab.move_into_collection(list);

  /*----------------------------------------------------------------------*/
  Core::Utils::SectionSpecs fdyn_edge_based_stab{fdyn, "EDGE-BASED STABILIZATION"};

  //! Flag to (de)activate edge-based (EOS) pressure stabilization
  fdyn_edge_based_stab.specs.emplace_back(deprecated_selection<EosPres>("EOS_PRES",
      {
          {"none", EOS_PRES_none},
          {"std_eos", EOS_PRES_std_eos},
          {"xfem_gp", EOS_PRES_xfem_gp},
      },
      {.description =
              "Flag to (de)activate pressure edge-based stabilization. Options: do not use "
              "pressure edge-based stabilization, or, use pressure edge-based stabilization as "
              "standard edge-based stabilization on the entire domain, or, use pressure edge-based "
              "stabilization as xfem ghost-penalty stabilization just around cut elements",
          .default_value = EOS_PRES_none}));

  //! Flag to (de)activate edge-based (EOS) convective streamline stabilization
  fdyn_edge_based_stab.specs.emplace_back(deprecated_selection<EosConvStream>("EOS_CONV_STREAM",
      {
          {"none", EOS_CONV_STREAM_none},
          {"std_eos", EOS_CONV_STREAM_std_eos},
          {"xfem_gp", EOS_CONV_STREAM_xfem_gp},
      },
      {.description =
              "Flag to (de)activate convective streamline edge-based stabilization. Options: do "
              "not use convective streamline edge-based stabilization, or, use convective "
              "streamline edge-based stabilization as standard edge-based stabilization on the "
              "entire domain, or, use convective streamline edge-based stabilization as xfem "
              "ghost-penalty stabilization just around cut elements",
          .default_value = EOS_CONV_STREAM_none}));

  //! Flag to (de)activate edge-based (EOS) convective crosswind stabilization
  fdyn_edge_based_stab.specs.emplace_back(deprecated_selection<EosConvCross>("EOS_CONV_CROSS",
      {
          {"none", EOS_CONV_CROSS_none},
          {"std_eos", EOS_CONV_CROSS_std_eos},
          {"xfem_gp", EOS_CONV_CROSS_xfem_gp},
      },
      {.description = "Flag to (de)activate convective crosswind edge-based stabilization. "
                      "Options:do not use convective crosswind edge-based stabilization, or, use "
                      "convective crosswind edge-based stabilization as standard edge-based "
                      "stabilization on the entire domain, or, use convective crosswind edge-based "
                      "stabilization as xfem ghost-penalty stabilization just around cut elements",
          .default_value = EOS_CONV_CROSS_none}));

  //! Flag to (de)activate edge-based (EOS) divergence stabilization
  fdyn_edge_based_stab.specs.emplace_back(deprecated_selection<EosDiv>("EOS_DIV",
      {
          {"none", EOS_DIV_none},
          {"vel_jump_std_eos", EOS_DIV_vel_jump_std_eos},
          {"vel_jump_xfem_gp", EOS_DIV_vel_jump_xfem_gp},
          {"div_jump_std_eos", EOS_DIV_div_jump_std_eos},
          {"div_jump_xfem_gp", EOS_DIV_div_jump_xfem_gp},
      },
      {.description =
              "Flag to (de)activate divergence edge-based stabilization. Options: do not use "
              "divergence edge-based stabilization, or, divergence edge-based stabilization based "
              "on velocity jump on the entire domain, or, divergence edge-based stabilization "
              "based on divergence jump just around cut elements, or, divergence edge-based "
              "stabilization based on velocity jump on the entire domain, or, divergence "
              "edge-based stabilization based on divergence jump just around cut elements",
          .default_value = EOS_DIV_none}));

  //! special least-squares condition for pseudo 2D examples where pressure level is determined via
  //! Krylov-projection
  fdyn_edge_based_stab.specs.emplace_back(parameter<bool>(
      "PRES_KRYLOV_2Dz", {.description = "residual based without second derivatives (i.e. only "
                                         "consistent for tau->0, but faster)",
                             .default_value = false}));

  //! this parameter selects the definition of Edge-based stabilization parameter
  fdyn_edge_based_stab.specs.emplace_back(deprecated_selection<EosTauType>("EOS_DEFINITION_TAU",
      {
          {"Burman_Fernandez_Hansbo", Inpar::FLUID::EOS_tau_burman_fernandez_hansbo},
          {"Burman_Fernandez_Hansbo_wo_dt", Inpar::FLUID::EOS_tau_burman_fernandez_hansbo_wo_dt},
          {"Braack_Burman_John_Lube", Inpar::FLUID::EOS_tau_braack_burman_john_lube},
          {"Braack_Burman_John_Lube_wo_divjump",
              Inpar::FLUID::EOS_tau_braack_burman_john_lube_wo_divjump},
          {"Franca_Barrenechea_Valentin_Wall",
              Inpar::FLUID::EOS_tau_franca_barrenechea_valentin_wall},
          {"Burman_Fernandez", Inpar::FLUID::EOS_tau_burman_fernandez},
          {"Burman_Hansbo_DAngelo_Zunino", Inpar::FLUID::EOS_tau_burman_hansbo_dangelo_zunino},
          {"Burman_Hansbo_DAngelo_Zunino_wo_dt",
              Inpar::FLUID::EOS_tau_burman_hansbo_dangelo_zunino_wo_dt},
          {"Schott_Massing_Burman_DAngelo_Zunino",
              Inpar::FLUID::EOS_tau_schott_massing_burman_dangelo_zunino},
          {"Schott_Massing_Burman_DAngelo_Zunino_wo_dt",
              Inpar::FLUID::EOS_tau_schott_massing_burman_dangelo_zunino_wo_dt},
          {"Burman", Inpar::FLUID::EOS_tau_burman},
          {"Taylor_Hughes_Zarins_Whiting_Jansen_Codina_scaling",
              Inpar::FLUID::EOS_tau_Taylor_Hughes_Zarins_Whiting_Jansen_Codina_scaling},
          {"tau_not_defined", Inpar::FLUID::EOS_tau_not_defined},
      },
      {.description = "Definition of stabilization parameter for edge-based stabilization",
          .default_value = Inpar::FLUID::EOS_tau_burman_hansbo_dangelo_zunino}));

  //! this parameter selects how the element length of Edge-based stabilization is defined
  fdyn_edge_based_stab.specs.emplace_back(parameter<EosElementLength>("EOS_H_DEFINITION",
      {.description =
              "Definition of element length for edge-based stabilization. Options:take the maximal "
              "(nsd-1)D diameter of faces that connect the internal face to its opposite faces, "
              "or, take the maximal 1D distance along 1D edge to opposite surface for both parent "
              "elements, or, take the maximal (nsd-1)D face diameter of all faces for both parent "
              "elements, or, maximal nD diameter of the neighboring elements, or, maximal (n-1)D "
              "diameter of the internal face/edge, or, take the maximal volume equivalent diameter "
              "of adjacent elements",
          .default_value = EOS_he_max_diameter_to_opp_surf}));

  fdyn_edge_based_stab.move_into_collection(list);

  /*----------------------------------------------------------------------*/
  Core::Utils::SectionSpecs fdyn_porostab{fdyn, "POROUS-FLOW STABILIZATION"};

  fdyn_porostab.specs.emplace_back(parameter<bool>("STAB_BIOT",
      {.description = "Flag to (de)activate BIOT stabilization.", .default_value = false}));
  fdyn_porostab.specs.emplace_back(parameter<double>("STAB_BIOT_SCALING",
      {.description =
              "Scaling factor for stabilization parameter for biot stabilization of porous flow.",
          .default_value = 1.0}));

  // this parameter defines various stabilized methods
  fdyn_porostab.specs.emplace_back(deprecated_selection<Inpar::FLUID::StabType>("STABTYPE",
      {
          {"no_stabilization", stabtype_nostab},
          {"residual_based", stabtype_residualbased},
          {"edge_based", stabtype_edgebased},
      },
      {.description = "Apply (un)stabilized fluid formulation. No stabilization is only possible "
                      "for inf-sup stable elements. Use a residual-based stabilization or, more "
                      "generally, a stabilization \nbased on the concept of the residual-based "
                      "variational multiscale method...\nExpecting additional inputUse an "
                      "edge-based stabilization, especially for XFEM",
          .default_value = stabtype_residualbased}));

  fdyn_porostab.specs.emplace_back(parameter<bool>(
      "INCONSISTENT", {.description = "residual based without second derivatives (i.e. only "
                                      "consistent for tau->0, but faster)",
                          .default_value = false}));

  fdyn_porostab.specs.emplace_back(parameter<bool>("Reconstruct_Sec_Der",
      {.description = "residual computed with a reconstruction of the second derivatives via "
                      "projection or superconvergent patch recovery",
          .default_value = false}));

  // the following parameters are necessary only if a residual based stabilized method is applied
  fdyn_porostab.specs.emplace_back(deprecated_selection<SubscalesTD>("TDS",
      {
          {"quasistatic", subscales_quasistatic},
          {"time_dependent", subscales_time_dependent},
      },
      {.description =
              "Flag to allow time dependency of subscales for residual-based stabilization.",
          .default_value = subscales_quasistatic}));

  fdyn_porostab.specs.emplace_back(deprecated_selection<Transient>("TRANSIENT",
      {
          {"no_transient", inertia_stab_drop},
          {"yes_transient", inertia_stab_keep},
          {"transient_complete", inertia_stab_keep_complete},
      },
      {.description = "Specify how to treat the transient term.",
          .default_value = inertia_stab_drop}));

  fdyn_porostab.specs.emplace_back(parameter<bool>(
      "PSPG", {.description = "Flag to (de)activate PSPG stabilization.", .default_value = true}));
  fdyn_porostab.specs.emplace_back(parameter<bool>(
      "SUPG", {.description = "Flag to (de)activate SUPG stabilization.", .default_value = true}));
  fdyn_porostab.specs.emplace_back(parameter<bool>(
      "GRAD_DIV", {.description = "Flag to (de)activate grad-div term.", .default_value = true}));

  fdyn_porostab.specs.emplace_back(deprecated_selection<VStab>("VSTAB",
      {
          {"no_vstab", viscous_stab_none},
          {"vstab_gls", viscous_stab_gls},
          {"vstab_gls_rhs", viscous_stab_gls_only_rhs},
          {"vstab_usfem", viscous_stab_usfem},
          {"vstab_usfem_rhs", viscous_stab_usfem_only_rhs},
      },
      {.description = "Flag to (de)activate viscous term in residual-based stabilization.",
          .default_value = viscous_stab_none}));

  fdyn_porostab.specs.emplace_back(deprecated_selection<RStab>("RSTAB",
      {
          {"no_rstab", reactive_stab_none},
          {"rstab_gls", reactive_stab_gls},
          {"rstab_usfem", reactive_stab_usfem},
      },
      {.description = "Flag to (de)activate reactive term in residual-based stabilization.",
          .default_value = reactive_stab_none}));

  fdyn_porostab.specs.emplace_back(deprecated_selection<CrossStress>("CROSS-STRESS",
      {
          {"no_cross", cross_stress_stab_none},
          {"yes_cross", cross_stress_stab},
          {"cross_rhs", cross_stress_stab_only_rhs},
      },
      {.description = "Flag to (de)activate cross-stress term -> residual-based VMM.",
          .default_value = cross_stress_stab_none}));

  fdyn_porostab.specs.emplace_back(deprecated_selection<ReynoldsStress>("REYNOLDS-STRESS",
      {
          {"no_reynolds", reynolds_stress_stab_none},
          {"yes_reynolds", reynolds_stress_stab},
          {"reynolds_rhs", reynolds_stress_stab_only_rhs},
      },
      {.description = "Flag to (de)activate Reynolds-stress term -> residual-based VMM.",
          .default_value = reynolds_stress_stab_none}));

  // this parameter selects the tau definition applied
  fdyn_porostab.specs.emplace_back(deprecated_selection<TauType>("DEFINITION_TAU",
      {
          {"Taylor_Hughes_Zarins", tau_taylor_hughes_zarins},
          {"Taylor_Hughes_Zarins_wo_dt", tau_taylor_hughes_zarins_wo_dt},
          {"Taylor_Hughes_Zarins_Whiting_Jansen", tau_taylor_hughes_zarins_whiting_jansen},
          {"Taylor_Hughes_Zarins_Whiting_Jansen_wo_dt",
              tau_taylor_hughes_zarins_whiting_jansen_wo_dt},
          {"Taylor_Hughes_Zarins_scaled", tau_taylor_hughes_zarins_scaled},
          {"Taylor_Hughes_Zarins_scaled_wo_dt", tau_taylor_hughes_zarins_scaled_wo_dt},
          {"Franca_Barrenechea_Valentin_Frey_Wall", tau_franca_barrenechea_valentin_frey_wall},
          {"Franca_Barrenechea_Valentin_Frey_Wall_wo_dt",
              tau_franca_barrenechea_valentin_frey_wall_wo_dt},
          {"Shakib_Hughes_Codina", tau_shakib_hughes_codina},
          {"Shakib_Hughes_Codina_wo_dt", tau_shakib_hughes_codina_wo_dt},
          {"Codina", tau_codina},
          {"Codina_wo_dt", tau_codina_wo_dt},
          {"Franca_Madureira_Valentin_Badia_Codina", tau_franca_madureira_valentin_badia_codina},
          {"Franca_Madureira_Valentin_Badia_Codina_wo_dt",
              tau_franca_madureira_valentin_badia_codina_wo_dt},
      },
      {.description = "Definition of tau_M and Tau_C",
          .default_value = tau_franca_barrenechea_valentin_frey_wall}));

  // this parameter selects the characteristic element length for tau_Mu for all
  // stabilization parameter definitions requiring such a length
  fdyn_porostab.specs.emplace_back(deprecated_selection<CharEleLengthU>("CHARELELENGTH_U",
      {
          {"streamlength", streamlength_u},
          {"volume_equivalent_diameter", volume_equivalent_diameter_u},
          {"root_of_volume", root_of_volume_u},
      },
      {.description = "Characteristic element length for tau_Mu",
          .default_value = streamlength_u}));

  // this parameter selects the characteristic element length for tau_Mp and tau_C for
  // all stabilization parameter definitions requiring such a length
  fdyn_porostab.specs.emplace_back(deprecated_selection<CharEleLengthPC>("CHARELELENGTH_PC",
      {
          {"streamlength", streamlength_pc},
          {"volume_equivalent_diameter", volume_equivalent_diameter_pc},
          {"root_of_volume", root_of_volume_pc},
      },
      {.description = "Characteristic element length for tau_Mp/tau_C",
          .default_value = volume_equivalent_diameter_pc}));

  // this parameter selects the location where tau is evaluated
  evaluation_tau_valid_input = {"element_center", "integration_point"};
  fdyn_porostab.specs.emplace_back(
      deprecated_selection<std::string>("EVALUATION_TAU", evaluation_tau_valid_input,
          {.description = "Location where tau is evaluated", .default_value = "element_center"}));

  // this parameter selects the location where the material law is evaluated
  // (does not fit here very well, but parameter transfer is easier)
  evaluation_mat_valid_input = {"element_center", "integration_point"};
  fdyn_porostab.specs.emplace_back(
      deprecated_selection<std::string>("EVALUATION_MAT", evaluation_mat_valid_input,
          {.description = "Location where material law is evaluated",
              .default_value = "element_center"}));


  // these parameters active additional terms in loma continuity equation
  // which might be identified as SUPG-/cross- and Reynolds-stress term
  fdyn_porostab.specs.emplace_back(parameter<bool>("LOMA_CONTI_SUPG",
      {.description = "Flag to (de)activate SUPG stabilization in loma continuity equation.",
          .default_value = false}));

  fdyn_porostab.specs.emplace_back(deprecated_selection<CrossStress>("LOMA_CONTI_CROSS_STRESS",
      {
          {"no_cross", cross_stress_stab_none},
          {"yes_cross", cross_stress_stab},
          {"cross_rhs", cross_stress_stab_only_rhs},
      },
      {.description = "Flag to (de)activate cross-stress term loma continuity equation-> "
                      "residual-based VMM.",
          .default_value = cross_stress_stab_none}));

  fdyn_porostab.specs.emplace_back(
      deprecated_selection<ReynoldsStress>("LOMA_CONTI_REYNOLDS_STRESS",
          {
              {"no_reynolds", reynolds_stress_stab_none},
              {"yes_reynolds", reynolds_stress_stab},
              {"reynolds_rhs", reynolds_stress_stab_only_rhs},
          },
          {.description = "Flag to (de)activate Reynolds-stress term loma continuity equation-> "
                          "residual-based VMM.",
              .default_value = reynolds_stress_stab_none}));

  fdyn_porostab.move_into_collection(list);

  /*----------------------------------------------------------------------*/
  Core::Utils::SectionSpecs fdyn_turbu{fdyn, "TURBULENCE MODEL"};

  //----------------------------------------------------------------------
  // modeling strategies
  //----------------------------------------------------------------------

  std::vector<std::string> turbulence_approach_valid_input = {"DNS_OR_RESVMM_LES", "CLASSICAL_LES"};
  std::string turbulence_approach_doc_string =
      "Try to solve flow as an underresolved DNS. Mind that your stabilisation already acts as a "
      "kind of turbulence model! Perform a classical Large Eddy Simulation adding addititional "
      "turbulent viscosity. This may be based on various physical models.)";

  fdyn_turbu.specs.emplace_back(
      deprecated_selection<std::string>("TURBULENCE_APPROACH", turbulence_approach_valid_input,
          {.description = turbulence_approach_doc_string, .default_value = "DNS_OR_RESVMM_LES"}));

  std::vector<std::string> physical_model_valid_input = {"no_model", "Smagorinsky",
      "Smagorinsky_with_van_Driest_damping", "Dynamic_Smagorinsky", "Multifractal_Subgrid_Scales",
      "Vreman", "Dynamic_Vreman"};
  std::string physical_model_doc_string =
      "If classical LES is our turbulence approach, this is a contradiction and should cause a "
      "FOUR_C_THROW. Classical constant coefficient Smagorinsky model. Be careful if you have a "
      "wall bounded flow domain! Use an exponential damping function for the turbulent viscosity "
      "close to the wall. This is only implemented for a channel geometry of height 2 in y "
      "direction. The viscous lengthscale l_tau is required as additional input. The solution is "
      "filtered and by comparison of the filtered velocity field with the real solution, the "
      "Smagorinsky constant is estimated in each step --- mind that this procedure includes an "
      "averaging in the xz plane, hence this implementation will only work for a channel flow. "
      "Multifractal Subgrid-Scale Modeling based on the work of burton. Vremans constant model. "
      "Dynamic Vreman model according to You and Moin (2007)";
  fdyn_turbu.specs.emplace_back(
      deprecated_selection<std::string>("PHYSICAL_MODEL", physical_model_valid_input,
          {.description = physical_model_doc_string, .default_value = "no_model"}));

  fdyn_turbu.specs.emplace_back(deprecated_selection<Inpar::FLUID::FineSubgridVisc>("FSSUGRVISC",
      {
          {"No", no_fssgv},
          {"Smagorinsky_all", smagorinsky_all},
          {"Smagorinsky_small", smagorinsky_small},
      },
      {.description = "fine-scale subgrid viscosity", .default_value = no_fssgv}));

  //----------------------------------------------------------------------
  // turbulence specific output and statistics
  //----------------------------------------------------------------------

  fdyn_turbu.specs.emplace_back(parameter<int>(
      "SAMPLING_START", {.description = "Time step after when sampling shall be started",
                            .default_value = 10000000}));
  fdyn_turbu.specs.emplace_back(parameter<int>("SAMPLING_STOP",
      {.description = "Time step when sampling shall be stopped", .default_value = 1}));
  fdyn_turbu.specs.emplace_back(parameter<int>("DUMPING_PERIOD",
      {.description = "Period of time steps after which statistical data shall be dumped",
          .default_value = 1}));
  fdyn_turbu.specs.emplace_back(parameter<bool>(
      "SUBGRID_DISSIPATION", {.description = "Flag to (de)activate estimation of subgrid-scale "
                                             "dissipation (only for seclected flows).",
                                 .default_value = false}));
  fdyn_turbu.specs.emplace_back(parameter<bool>("OUTMEAN",
      {.description = "Flag to (de)activate averaged paraview output", .default_value = false}));
  fdyn_turbu.specs.emplace_back(parameter<bool>(
      "TURBMODEL_LS", {.description = "Flag to (de)activate turbulence model in level-set equation",
                          .default_value = true}));

  //----------------------------------------------------------------------
  // turbulent flow problem and general characteristics
  //----------------------------------------------------------------------

  {
    std::vector<std::string> canonical_flow_valid_input = {"no", "time_averaging",
        "channel_flow_of_height_2", "lid_driven_cavity", "backward_facing_step", "square_cylinder",
        "square_cylinder_nurbs", "rotating_circular_cylinder_nurbs",
        "rotating_circular_cylinder_nurbs_scatra", "loma_channel_flow_of_height_2",
        "loma_lid_driven_cavity", "loma_backward_facing_step", "combust_oracles",
        "bubbly_channel_flow", "scatra_channel_flow_of_height_2",
        "decaying_homogeneous_isotropic_turbulence", "forced_homogeneous_isotropic_turbulence",
        "scatra_forced_homogeneous_isotropic_turbulence", "taylor_green_vortex", "periodic_hill",
        "blood_fda_flow", "backward_facing_step2"};

    std::string canonical_flow_doc =
        ""
        "Sampling is different for different canonical flows - so specify what kind of flow you've "
        "got \n\n"
        "no: The flow is not further specified, so spatial averaging and hence the standard "
        "sampling procedure is not possible\n"
        "time_averaging: The flow is not further specified, but time averaging of velocity and "
        "pressure field is performed\n"
        "channel_flow_of_height_2: For this flow, all statistical data could be averaged in the "
        "homogeneous planes - it is essentially a statistically one dimensional flow.\n"
        "lid_driven_cavity: For this flow, all statistical data are evaluated on the center lines "
        "of the xy-midplane, averaged only over time.\n"
        "backward_facing_step: For this flow, statistical data are evaluated on various lines, "
        "averaged over time and z.\n"
        "square_cylinder: For this flow, statistical data are evaluated on various lines of the "
        "xy-midplane, averaged only over time.\n"
        "square_cylinder_nurbs: For this flow, statistical data are evaluated on various lines of "
        "the xy-midplane, averaged over time and eventually in one home.direction.\n"
        "rotating_circular_cylinder_nurbs: For this flow, statistical data is computed in "
        "concentric surfaces and averaged. in time and in one home. direction\n"
        "rotating_circular_cylinder_nurbs_scatra: For this flow with mass transport, statistical "
        "data is computed in concentric surfaces and averaged. in time and in one home. direction\n"
        "loma_channel_flow_of_height_2: For this low-Mach-number flow, all statistical data could "
        "be averaged in the homogeneous planes - it is essentially a statistically one dimensional "
        "flow.\n"
        "loma_lid_driven_cavity: For this low-Mach-number flow, all statistical data are evaluated "
        "on the center lines of the xy-midplane, averaged only over time.\n"
        "loma_backward_facing_step: For this low-Mach-number flow, statistical data are evaluated "
        "on various lines, averaged over time and z.\n"
        "combust_oracles: ORACLES test rig for turbulent premixed combustion.\n"
        "bubbly_channel_flow: Turbulent two-phase flow: bubbly channel flow, statistical data are "
        "averaged in homogeneous planes and over time.\n"
        "scatra_channel_flow_of_height_2: For this flow, all statistical data could be averaged in "
        "the homogeneous planes - it is essentially a statistically one dimensional flow.\n"
        "decaying_homogeneous_isotropic_turbulence: For this flow, all statistical data could be "
        "averaged in the in all homogeneous directions  - it is essentially a statistically zero "
        "dimensional flow.\n"
        "forced_homogeneous_isotropic_turbulence: For this flow, all statistical data could be "
        "averaged in the in all homogeneous directions  - it is essentially a statistically zero "
        "dimensional flow.\n"
        "scatra_forced_homogeneous_isotropic_turbulence: For this flow, all statistical data could "
        "be averaged in the in all homogeneous directions  - it is essentially a statistically "
        "zero dimensional flow.\n"
        "taylor_green_vortex: For this flow, dissipation rate could be averaged in the in all "
        "homogeneous directions  - it is essentially a statistically zero dimensional flow.\n"
        "periodic_hill: For this flow, statistical data is evaluated on various lines, averaged "
        "over time and z.\n"
        "blood_fda_flow: For this flow, statistical data is evaluated on various planes.\n"
        "backward_facing_step2: For this flow, statistical data is evaluated on various planes.\n";

    fdyn_turbu.specs.emplace_back(deprecated_selection<std::string>("CANONICAL_FLOW",
        canonical_flow_valid_input, {.description = canonical_flow_doc, .default_value = "no"}));
  }

  std::vector<std::string> homdir_valid_input = {
      "not_specified", "x", "y", "z", "xy", "xz", "yz", "xyz"};
  std::string homdir_doc =
      "Specify the homogeneous direction(s) of a flow.\n"
      "not_specified: no homogeneous directions available, averaging is restricted to time "
      "averaging\n"
      "x: average along x-direction\n"
      "y: average along y-direction\n"
      "z: average along z-direction\n"
      "xy: Wall normal direction is z, average in x and y direction\n"
      "xz: Wall normal direction is y, average in x and z direction\n"
      "yz: Wall normal direction is x, average in y and z direction\n"
      "xyz: Averaging in all directions\n";
  fdyn_turbu.specs.emplace_back(deprecated_selection<std::string>(
      "HOMDIR", homdir_valid_input, {.description = homdir_doc, .default_value = "not_specified"}));

  //---------------------------------------
  // further problem-specific parameters

  // CHANNEL FLOW
  //--------------

  fdyn_turbu.specs.emplace_back(parameter<double>(
      "CHAN_AMPL_INIT_DIST", {.description = "Max. amplitude of the random disturbance in percent "
                                             "of the initial value in mean flow direction.",
                                 .default_value = 0.1}));

  fdyn_turbu.specs.emplace_back(parameter<ForcingType>(
      "FORCING_TYPE", {.description = "forcing strategy",
                          .default_value = linear_compensation_from_intermediate_spectrum}));

  fdyn_turbu.specs.emplace_back(parameter<int>("CHA_NUMSUBDIVISIONS",
      {.description = "Number of homogeneous sampling planes in element", .default_value = 5}));

  // HIT
  //--------------

  fdyn_turbu.specs.emplace_back(parameter<int>("FORCING_TIME_STEPS",
      {.description = "Number of time steps during which forcing is applied", .default_value = 0}));

  fdyn_turbu.specs.emplace_back(parameter<double>(
      "THRESHOLD_WAVENUMBER", {.description = "Forcing is only applied to wave numbers lower or "
                                              "equal than the given threshold wave number.",
                                  .default_value = 0.0}));

  fdyn_turbu.specs.emplace_back(
      parameter<double>("POWER_INPUT", {.description = "power of forcing", .default_value = 0.0}));

  std::vector<std::string> scalar_forcing_valid_input = {"no", "isotropic", "mean_scalar_gradient"};
  std::string scalar_forcing_doc =
      "no: Do not force the scalar field\n"
      "isotropic: Force scalar field isotropically such as the fluid field.\n"
      "mean_scalar_gradient: Force scalar field by imposed mean-scalar gradient.\n";
  fdyn_turbu.specs.emplace_back(deprecated_selection<std::string>("SCALAR_FORCING",
      scalar_forcing_valid_input, {.description = scalar_forcing_doc, .default_value = "no"}));

  fdyn_turbu.specs.emplace_back(parameter<double>("MEAN_SCALAR_GRADIENT",
      {.description = "Value of imposed mean-scalar gradient to force scalar field.",
          .default_value = 0.0}));

  // filtering with xfem
  //--------------

  fdyn_turbu.specs.emplace_back(parameter<bool>("EXCLUDE_XFEM",
      {.description = "Flag to (de)activate XFEM dofs in calculation of fine-scale velocity.",
          .default_value = false}));

  fdyn_turbu.move_into_collection(list);

  /*----------------------------------------------------------------------*/
  // sublist with additional input parameters for Smagorinsky model
  Core::Utils::SectionSpecs fdyn_turbsgv{fdyn, "SUBGRID VISCOSITY"};

  fdyn_turbsgv.specs.emplace_back(parameter<double>("C_SMAGORINSKY",
      {.description =
              "Constant for the Smagorinsky model. Something between 0.1 to 0.24. Vreman constant "
              "if the constant vreman model is applied (something between 0.07 and 0.01).",
          .default_value = 0.0}));
  fdyn_turbsgv.specs.emplace_back(parameter<double>("C_YOSHIZAWA",
      {.description = "Constant for the compressible Smagorinsky model: isotropic part of "
                      "subgrid-stress tensor. About 0.09 or 0.0066. Ci will not be squared!",
          .default_value = -1.0}));
  fdyn_turbsgv.specs.emplace_back(parameter<bool>("C_SMAGORINSKY_AVERAGED",
      {.description = "Flag to (de)activate averaged Smagorinksy constant",
          .default_value = false}));
  fdyn_turbsgv.specs.emplace_back(parameter<bool>("C_INCLUDE_CI",
      {.description = "Flag to (de)inclusion of Yoshizawa model", .default_value = false}));
  // remark: following Moin et al. 1991, the extension of the dynamic Smagorinsky model to
  // compressibel flow
  //        also contains a model for the isotropic part of the subgrid-stress tensor according to
  //        Yoshizawa 1989 although used in literature for turbulent variable-density flow at low
  //        Mach number, this additional term seems to destabilize the simulation when the flow is
  //        only weakly compressible therefore C_INCLUDE_CI allows to exclude this term if
  //        C_SMAGORINSKY_AVERAGED == true
  //           if C_INCLUDE_CI==true and C_YOSHIZAWA>=0.0 then the given value C_YOSHIZAWA is used
  //           if C_INCLUDE_CI==true and C_YOSHIZAWA<0.0 then C_YOSHIZAWA is determined dynamically
  //        else all values are taken from input

  fdyn_turbsgv.specs.emplace_back(parameter<double>("CHANNEL_L_TAU",
      {.description = "Used for normalisation of the wall normal distance in the Van \nDriest "
                      "Damping function. May be taken from the output of \nthe "
                      "apply_mesh_stretching.pl preprocessing script.",
          .default_value = 0.0}));

  fdyn_turbsgv.specs.emplace_back(parameter<double>("C_TURBPRANDTL",
      {.description =
              "(Constant) turbulent Prandtl number for the Smagorinsky model in scalar transport.",
          .default_value = 1.0}));

  fdyn_turbsgv.specs.emplace_back(deprecated_selection<VremanFiMethod>("FILTER_WIDTH",
      {
          {"CubeRootVol", cuberootvol},
          {"Direction_dependent", dir_dep},
          {"Minimum_length", min_len},
      },
      {.description = "The Vreman model requires a filter width.", .default_value = cuberootvol}));

  fdyn_turbsgv.move_into_collection(list);

  /*----------------------------------------------------------------------*/
  // sublist with additional input parameters for Smagorinsky model
  Core::Utils::SectionSpecs fdyn_wallmodel{fdyn, "WALL MODEL"};

  fdyn_wallmodel.specs.emplace_back(parameter<bool>(
      "X_WALL", {.description = "Flag to switch on the xwall model", .default_value = false}));


  std::vector<std::string> tauw_type_valid_input = {"constant", "between_steps"};
  std::string tauw_type_doc =
      "constant: Use the constant wall shear stress given in the input file for the whole "
      "simulation.\n"
      "between_steps: Calculate wall shear stress in between time steps.\n";
  fdyn_wallmodel.specs.emplace_back(deprecated_selection<std::string>("Tauw_Type",
      tauw_type_valid_input, {.description = tauw_type_doc, .default_value = "constant"}));

  std::vector<std::string> tauw_calc_type_valid_input = {
      "residual", "gradient", "gradient_to_residual"};
  std::string tauw_calc_type_doc =
      "residual: Residual (force) divided by area.\n"
      "gradient: Gradient via shape functions and nodal values.\n"
      "gradient_to_residual: First gradient, then residual.\n";
  fdyn_wallmodel.specs.emplace_back(
      deprecated_selection<std::string>("Tauw_Calc_Type", tauw_calc_type_valid_input,
          {.description = tauw_calc_type_doc, .default_value = "residual"}));


  fdyn_wallmodel.specs.emplace_back(parameter<int>("Switch_Step",
      {.description = "Switch from gradient to residual based tauw.", .default_value = -1}));

  std::vector<std::string> projection_valid_input = {
      "No", "onlyl2projection", "l2projectionwithcontinuityconstraint"};
  std::string projection_doc =
      "Flag to switch projection of the enriched dofs after updating tauw, alternatively with or "
      "without continuity constraint.";
  fdyn_wallmodel.specs.emplace_back(deprecated_selection<std::string>("Projection",
      projection_valid_input, {.description = projection_doc, .default_value = "No"}));

  fdyn_wallmodel.specs.emplace_back(parameter<double>(
      "C_Tauw", {.description = "Constant wall shear stress for Spalding's law, if applicable",
                    .default_value = 1.0}));

  fdyn_wallmodel.specs.emplace_back(parameter<double>(
      "Min_Tauw", {.description = "Minimum wall shear stress preventing system to become singular",
                      .default_value = 2.0e-9}));

  fdyn_wallmodel.specs.emplace_back(parameter<double>(
      "Inc_Tauw", {.description = "Increment of Tauw of full step, between 0.0 and 1.0",
                      .default_value = 1.0}));

  std::vector<std::string> blending_type_valid_input = {"none", "ramp_function"};
  std::string blending_type_doc =
      "Methods for blending the enrichment space.\n"
      "none: No ramp function, does not converge!\n"
      "ramp_function: Enrichment is multiplied with linear ramp function resulting in zero "
      "enrichment at the interface.\n";
  fdyn_wallmodel.specs.emplace_back(deprecated_selection<std::string>("Blending_Type",
      blending_type_valid_input, {.description = blending_type_doc, .default_value = "none"}));


  fdyn_wallmodel.specs.emplace_back(parameter<int>("GP_Wall_Normal",
      {.description = "Gauss points in wall normal direction", .default_value = 3}));
  fdyn_wallmodel.specs.emplace_back(parameter<int>("GP_Wall_Normal_Off_Wall",
      {.description = "Gauss points in wall normal direction, off-wall elements",
          .default_value = 3}));
  fdyn_wallmodel.specs.emplace_back(parameter<int>("GP_Wall_Parallel",
      {.description = "Gauss points in wall parallel direction", .default_value = 3}));

  fdyn_wallmodel.specs.emplace_back(parameter<bool>("Treat_Tauw_on_Dirichlet_Inflow",
      {.description = "Flag to treat residual on Dirichlet inflow nodes for calculation of wall "
                      "shear stress",
          .default_value = false}));

  fdyn_wallmodel.specs.emplace_back(parameter<int>("PROJECTION_SOLVER",
      {.description = "Set solver number for l2-projection.", .default_value = -1}));

  fdyn_wallmodel.move_into_collection(list);

  /*----------------------------------------------------------------------*/
  // sublist with additional input parameters for multifractal subgrid-scales
  Core::Utils::SectionSpecs fdyn_turbmfs{fdyn, "MULTIFRACTAL SUBGRID SCALES"};

  fdyn_turbmfs.specs.emplace_back(parameter<double>("CSGS",
      {.description = "Modelparameter of multifractal subgrid-scales.", .default_value = 0.0}));

  std::vector<std::string> scale_separation_valid_input = {
      "no_scale_sep", "box_filter", "algebraic_multigrid_operator"};
  std::string scale_separation_doc =
      "Specify the filter type for scale separation in LES.\n"
      "no_scale_sep: no scale separation.\n"
      "box_filter: classical box filter.\n"
      "algebraic_multigrid_operator: scale separation by algebraic multigrid operator.\n";
  fdyn_turbmfs.specs.emplace_back(
      deprecated_selection<std::string>("SCALE_SEPARATION", scale_separation_valid_input,
          {.description = scale_separation_doc, .default_value = "no_scale_sep"}));


  fdyn_turbmfs.specs.emplace_back(parameter<int>("ML_SOLVER",
      {.description =
              "Set solver number for scale separation via level set transfer operators from "
              "plain aggregation.",
          .default_value = -1}));

  fdyn_turbmfs.specs.emplace_back(parameter<bool>(
      "CALC_N", {.description = "Flag to (de)activate calculation of N from the Reynolds number.",
                    .default_value = false}));

  fdyn_turbmfs.specs.emplace_back(parameter<double>(
      "N", {.description = "Set grid to viscous scale ratio.", .default_value = 1.0}));

  std::vector<std::string> ref_length_valid_input = {
      "cube_edge", "sphere_diameter", "streamlength", "gradient_based", "metric_tensor"};
  std::string ref_length_doc =
      "Specify the reference length for Re-dependent N.\n"
      "cube_edge: edge length of volume equivalent cube.\n"
      "sphere_diameter: diameter of volume equivalent sphere.\n"
      "streamlength: streamlength taken from stabilization.\n"
      "gradient_based: gradient based length taken from stabilization.\n"
      "metric_tensor: metric tensor taken from stabilization.\n";
  fdyn_turbmfs.specs.emplace_back(deprecated_selection<std::string>("REF_LENGTH",
      ref_length_valid_input, {.description = ref_length_doc, .default_value = "cube_edge"}));

  std::vector<std::string> ref_velocity_valid_input = {"strainrate", "resolved", "fine_scale"};
  std::string ref_velocity_doc =
      "Specify the reference velocity for Re-dependent N.\n"
      "strainrate: norm of strain rate.\n"
      "resolved: resolved velocity.\n"
      "fine_scale: fine-scale velocity.\n";
  fdyn_turbmfs.specs.emplace_back(deprecated_selection<std::string>("REF_VELOCITY",
      ref_velocity_valid_input, {.description = ref_velocity_doc, .default_value = "strainrate"}));


  fdyn_turbmfs.specs.emplace_back(parameter<double>("C_NU",
      {.description =
              "Proportionality constant between Re and ratio viscous scale to element length.",
          .default_value = 1.0}));

  fdyn_turbmfs.specs.emplace_back(parameter<bool>("NEAR_WALL_LIMIT",
      {.description = "Flag to (de)activate near-wall limit.", .default_value = false}));

  std::vector<std::string> evaluation_b_valid_input = {"element_center", "integration_point"};
  std::string evaluation_b_doc =
      "Location where B is evaluated\n"
      "element_center: evaluate B at element center.\n"
      "integration_point: evaluate B at integration point.\n";
  fdyn_turbmfs.specs.emplace_back(
      deprecated_selection<std::string>("EVALUATION_B", evaluation_b_valid_input,
          {.description = evaluation_b_doc, .default_value = "element_center"}));


  fdyn_turbmfs.specs.emplace_back(parameter<double>(
      "BETA", {.description = "Cross- and Reynolds-stress terms only on right-hand-side.",
                  .default_value = 0.0}));

  convform_valid_input = {"convective", "conservative"};
  std::string convform_doc =
      "form of convective term\n"
      "convective: Use the convective form.\n"
      "conservative: Use the conservative form.\n";
  fdyn_turbmfs.specs.emplace_back(deprecated_selection<std::string>("CONVFORM",
      convform_valid_input, {.description = convform_doc, .default_value = "convective"}));

  fdyn_turbmfs.specs.emplace_back(parameter<double>("CSGS_PHI",
      {.description = "Modelparameter of multifractal subgrid-scales for scalar transport.",
          .default_value = 0.0}));

  fdyn_turbmfs.specs.emplace_back(parameter<bool>("ADAPT_CSGS_PHI",
      {.description = "Flag to (de)activate adaption of CsgsD to CsgsB.", .default_value = false}));

  fdyn_turbmfs.specs.emplace_back(parameter<bool>("NEAR_WALL_LIMIT_CSGS_PHI",
      {.description = "Flag to (de)activate near-wall limit for scalar field.",
          .default_value = false}));

  fdyn_turbmfs.specs.emplace_back(parameter<bool>("CONSISTENT_FLUID_RESIDUAL",
      {.description = "Flag to (de)activate the consistency term for residual-based stabilization.",
          .default_value = false}));

  fdyn_turbmfs.specs.emplace_back(parameter<double>(
      "C_DIFF", {.description = "Proportionality constant between Re*Pr and ratio dissipative "
                                "scale to element length. Usually equal cnu.",
                    .default_value = 1.0}));

  fdyn_turbmfs.specs.emplace_back(parameter<bool>("SET_FINE_SCALE_VEL",
      {.description = "Flag to set fine-scale velocity for parallel nightly tests.",
          .default_value = false}));

  // activate cross- and Reynolds-stress terms in loma continuity equation
  fdyn_turbmfs.specs.emplace_back(parameter<bool>("LOMA_CONTI",
      {.description =
              "Flag to (de)activate cross- and Reynolds-stress terms in loma continuity equation.",
          .default_value = false}));

  fdyn_turbmfs.move_into_collection(list);

  /*----------------------------------------------------------------------*/
  Core::Utils::SectionSpecs fdyn_turbinf{fdyn, "TURBULENT INFLOW"};

  fdyn_turbinf.specs.emplace_back(parameter<bool>("TURBULENTINFLOW",
      {.description = "Flag to (de)activate potential separate turbulent inflow section",
          .default_value = false}));

  fdyn_turbinf.specs.emplace_back(deprecated_selection<InitialField>("INITIALINFLOWFIELD",
      {
          {"zero_field", initfield_zero_field},
          {"field_by_function", initfield_field_by_function},
          {"disturbed_field_from_function", initfield_disturbed_field_from_function},
      },
      {.description = "Initial field for inflow section", .default_value = initfield_zero_field}));

  fdyn_turbinf.specs.emplace_back(parameter<int>(
      "INFLOWFUNC", {.description = "Function number for initial flow field in inflow section",
                        .default_value = -1}));

  fdyn_turbinf.specs.emplace_back(parameter<double>(
      "INFLOW_INIT_DIST", {.description = "Max. amplitude of the random disturbance in percent of "
                                          "the initial value in mean flow direction.",
                              .default_value = 0.1}));

  fdyn_turbinf.specs.emplace_back(parameter<int>("NUMINFLOWSTEP",
      {.description = "Total number of time steps for development of turbulent flow",
          .default_value = 1}));

  std::vector<std::string> canonical_inflow_valid_input = {"no", "time_averaging",
      "channel_flow_of_height_2", "loma_channel_flow_of_height_2",
      "scatra_channel_flow_of_height_2"};
  std::string canonical_inflow_doc =
      "Sampling is different for different canonical flows \n--- so specify what kind of flow "
      "you've got\n"
      "no: The flow is not further specified, so spatial averaging \nand hence the standard "
      "sampling procedure is not possible.\n"
      "time_averaging: The flow is not further specified, but time averaging of velocity and "
      "pressure field is performed.\n"
      "channel_flow_of_height_2: For this flow, all statistical data could be averaged in \nthe "
      "homogeneous planes --- it is essentially a statistically one dimensional flow.\n"
      "loma_channel_flow_of_height_2: For this low-Mach-number flow, all statistical data could be "
      "averaged in \nthe homogeneous planes --- it is essentially a statistically one dimensional "
      "flow.\n"
      "scatra_channel_flow_of_height_2: For this flow, all statistical data could be averaged in "
      "\nthe homogeneous planes --- it is essentially a statistically one dimensional flow.\n";
  fdyn_turbinf.specs.emplace_back(deprecated_selection<std::string>("CANONICAL_INFLOW",
      canonical_inflow_valid_input, {.description = canonical_inflow_doc, .default_value = "no"}));


  fdyn_turbinf.specs.emplace_back(parameter<double>("INFLOW_CHA_SIDE",
      {.description = "Most right side of inflow channel. Necessary to define sampling domain.",
          .default_value = 0.0}));

  std::vector<std::string> inflow_homdir_valid_input = {
      "not_specified", "x", "y", "z", "xy", "xz", "yz"};
  std::string inflow_homdir_doc =
      "Specify the homogeneous direction(s) of a flow\n"
      "not_specified: no homogeneous directions available, averaging is restricted to time "
      "averaging.\n"
      "x: average along x-direction.\n"
      "y: average along y-direction.\n"
      "z: average along z-direction.\n"
      "xy: Wall normal direction is z, average in x and y direction.\n"
      "xz: Wall normal direction is y, average in x and z direction (standard case).\n"
      "yz: Wall normal direction is x, average in y and z direction.\n";
  fdyn_turbinf.specs.emplace_back(
      deprecated_selection<std::string>("INFLOW_HOMDIR", inflow_homdir_valid_input,
          {.description = inflow_homdir_doc, .default_value = "not_specified"}));

  fdyn_turbinf.specs.emplace_back(parameter<int>(
      "INFLOW_SAMPLING_START", {.description = "Time step after when sampling shall be started",
                                   .default_value = 10000000}));
  fdyn_turbinf.specs.emplace_back(parameter<int>("INFLOW_SAMPLING_STOP",
      {.description = "Time step when sampling shall be stopped", .default_value = 1}));
  fdyn_turbinf.specs.emplace_back(parameter<int>("INFLOW_DUMPING_PERIOD",
      {.description = "Period of time steps after which statistical data shall be dumped",
          .default_value = 1}));

  fdyn_turbinf.move_into_collection(list);

  /*----------------------------------------------------------------------*/
  // sublist with additional input parameters for time adaptivity in fluid/ coupled problems
  Core::Utils::SectionSpecs fdyn_timintada{fdyn, "TIMEADAPTIVITY"};
  fdyn_timintada.specs.emplace_back(
      deprecated_selection<AdaptiveTimeStepEstimator>("ADAPTIVE_TIME_STEP_ESTIMATOR",
          {
              {"none", const_dt},
              {"cfl_number", cfl_number},
              {"only_print_cfl_number", only_print_cfl_number},
          },
          {.description = "Method used to determine adaptive time step size.",
              .default_value = const_dt}));

  fdyn_timintada.specs.emplace_back(parameter<double>(
      "CFL_NUMBER", {.description = "CFL number for adaptive time step", .default_value = -1.0}));
  fdyn_timintada.specs.emplace_back(parameter<int>("FREEZE_ADAPTIVE_DT_AT",
      {.description = "keep time step constant after this step, otherwise turbulence statistics "
                      "sampling is not consistent",
          .default_value = 1000000}));
  fdyn_timintada.specs.emplace_back(parameter<double>("ADAPTIVE_DT_INC",
      {.description = "Increment of whole step for adaptive dt via CFL", .default_value = 0.8}));

  fdyn_timintada.move_into_collection(list);
}



void Inpar::LowMach::set_valid_parameters(std::map<std::string, Core::IO::InputSpec>& list)
{
  using Teuchos::tuple;
  using namespace Core::IO::InputSpecBuilders;

  Core::Utils::SectionSpecs lomacontrol{"LOMA CONTROL"};

  lomacontrol.specs.emplace_back(
      parameter<bool>("MONOLITHIC", {.description = "monolithic solver", .default_value = false}));
  lomacontrol.specs.emplace_back(parameter<int>(
      "NUMSTEP", {.description = "Total number of time steps", .default_value = 24}));
  lomacontrol.specs.emplace_back(
      parameter<double>("TIMESTEP", {.description = "Time increment dt", .default_value = 0.1}));
  lomacontrol.specs.emplace_back(parameter<double>(
      "MAXTIME", {.description = "Total simulation time", .default_value = 1000.0}));
  lomacontrol.specs.emplace_back(parameter<int>(
      "ITEMAX", {.description = "Maximum number of outer iterations", .default_value = 10}));
  lomacontrol.specs.emplace_back(parameter<int>("ITEMAX_BEFORE_SAMPLING",
      {.description =
              "Maximum number of outer iterations before sampling (for turbulent flows only)",
          .default_value = 1}));
  lomacontrol.specs.emplace_back(parameter<double>(
      "CONVTOL", {.description = "Tolerance for convergence check", .default_value = 1e-6}));
  lomacontrol.specs.emplace_back(parameter<int>(
      "RESULTSEVERY", {.description = "Increment for writing solution", .default_value = 1}));
  lomacontrol.specs.emplace_back(parameter<int>(
      "RESTARTEVERY", {.description = "Increment for writing restart", .default_value = 1}));

  std::vector<std::string> constthermpress_valid_input = {"No_energy", "No_mass", "Yes"};
  lomacontrol.specs.emplace_back(
      deprecated_selection<std::string>("CONSTHERMPRESS", constthermpress_valid_input,
          {.description = "treatment of thermodynamic pressure in time", .default_value = "Yes"}));

  lomacontrol.specs.emplace_back(parameter<bool>(
      "SGS_MATERIAL_UPDATE", {.description = "update material by adding subgrid-scale scalar field",
                                 .default_value = false}));

  // number of linear solver used for LOMA solver
  lomacontrol.specs.emplace_back(parameter<int>("LINEAR_SOLVER",
      {.description = "number of linear solver used for LOMA problem", .default_value = -1}));

  lomacontrol.move_into_collection(list);
}


void Inpar::FLUID::set_valid_conditions(
    std::vector<Core::Conditions::ConditionDefinition>& condlist)
{
  using namespace Core::IO::InputSpecBuilders;

  /*--------------------------------------------------------------------*/
  // transfer boundary condition for turbulent inflow

  Core::Conditions::ConditionDefinition tbc_turb_inflow("DESIGN SURF TURBULENT INFLOW TRANSFER",
      "TransferTurbulentInflow", "TransferTurbulentInflow",
      Core::Conditions::TransferTurbulentInflow, true, Core::Conditions::geometry_type_surface);

  tbc_turb_inflow.add_component(parameter<int>("ID", {.description = ""}));
  tbc_turb_inflow.add_component(
      deprecated_selection<std::string>("toggle", {"master", "slave"}, {.description = "toggle"}));
  tbc_turb_inflow.add_component(parameter<TurbInflowDirection>(
      "DIRECTION", {.description = "transfer direction", .default_value = TurbInflowDirection::x}));
  tbc_turb_inflow.add_component(
      parameter<std::optional<int>>("curve", {.description = "curve id"}));

  condlist.push_back(tbc_turb_inflow);

  /*--------------------------------------------------------------------*/
  // separate domain for turbulent inflow generation

  Core::Conditions::ConditionDefinition turbulentinflowgeneration("FLUID TURBULENT INFLOW VOLUME",
      "TurbulentInflowSection", "TurbulentInflowSection", Core::Conditions::TurbulentInflowSection,
      true, Core::Conditions::geometry_type_volume);

  condlist.push_back(turbulentinflowgeneration);


  /*--------------------------------------------------------------------*/
  // flow-dependent pressure conditions

  Core::Conditions::ConditionDefinition surfflowdeppressure(
      "DESIGN SURFACE FLOW-DEPENDENT PRESSURE CONDITIONS", "SurfaceFlowDepPressure",
      "SurfaceFlowDepPressure", Core::Conditions::SurfaceFlowDepPressure, true,
      Core::Conditions::geometry_type_surface);

  // flow-dependent pressure conditions can be imposed either based on
  // (out)flow rate or (out)flow volume (e.g., for air-cushion condition)
  surfflowdeppressure.add_component(deprecated_selection<std::string>("TYPE_OF_FLOW_DEPENDENCE",
      {"flow_rate", "flow_volume", "fixed_pressure"}, {.description = "type of flow dependence"}));
  surfflowdeppressure.add_component(parameter<double>("ConstCoeff",
      {.description = "constant coefficient for (linear) flow-rate-based condition"}));
  surfflowdeppressure.add_component(parameter<double>(
      "LinCoeff", {.description = "linear coefficient for (linear) flow-rate-based condition"}));
  surfflowdeppressure.add_component(parameter<double>(
      "InitialVolume", {.description = "initial (air-cushion) volume outside of boundary"}));
  surfflowdeppressure.add_component(parameter<double>(
      "ReferencePressure", {.description = " reference pressure outside of boundary"}));
  surfflowdeppressure.add_component(
      parameter<double>("AdiabaticExponent", {.description = "adiabatic exponent"}));
  surfflowdeppressure.add_component(
      parameter<std::optional<int>>("curve", {.description = "curve id"}));

  condlist.emplace_back(surfflowdeppressure);


  /*--------------------------------------------------------------------*/
  // Slip Supplemental Curved Boundary conditions

  Core::Conditions::ConditionDefinition lineslipsupp(
      "DESIGN LINE SLIP SUPPLEMENTAL CURVED BOUNDARY CONDITIONS", "LineSlipSupp", "LineSlipSupp",
      Core::Conditions::LineSlipSupp, true, Core::Conditions::geometry_type_line);

  Core::Conditions::ConditionDefinition surfslipsupp(
      "DESIGN SURFACE SLIP SUPPLEMENTAL CURVED BOUNDARY CONDITIONS", "SurfaceSlipSupp",
      "SurfaceSlipSupp", Core::Conditions::SurfaceSlipSupp, true,
      Core::Conditions::geometry_type_surface);

  const auto make_slip_supp = [&](Core::Conditions::ConditionDefinition& cond)
  {
    cond.add_component(parameter<double>("USEUPDATEDNODEPOS"));
    condlist.emplace_back(cond);
  };

  make_slip_supp(lineslipsupp);
  make_slip_supp(surfslipsupp);

  /*--------------------------------------------------------------------*/
  // Navier-slip boundary conditions

  Core::Conditions::ConditionDefinition linenavierslip(
      "DESIGN LINE NAVIER-SLIP BOUNDARY CONDITIONS", "LineNavierSlip", "LineNavierSlip",
      Core::Conditions::LineNavierSlip, true, Core::Conditions::geometry_type_line);

  Core::Conditions::ConditionDefinition surfnavierslip(
      "DESIGN SURF NAVIER-SLIP BOUNDARY CONDITIONS", "SurfNavierSlip", "SurfNavierSlip",
      Core::Conditions::SurfNavierSlip, true, Core::Conditions::geometry_type_surface);

  const auto make_navierslip = [&](Core::Conditions::ConditionDefinition& cond)
  {
    cond.add_component(parameter<double>("SLIPCOEFFICIENT"));
    condlist.emplace_back(cond);
  };

  make_navierslip(linenavierslip);
  make_navierslip(surfnavierslip);

  /*--------------------------------------------------------------------*/
  // consistent outflow bcs for conservative element formulations

  Core::Conditions::ConditionDefinition surfconsistentoutflowconsistency(
      "DESIGN SURFACE CONSERVATIVE OUTFLOW CONSISTENCY", "SurfaceConservativeOutflowConsistency",
      "SurfaceConservativeOutflowConsistency",
      Core::Conditions::SurfaceConservativeOutflowConsistency, true,
      Core::Conditions::geometry_type_surface);

  condlist.push_back(surfconsistentoutflowconsistency);

  /*--------------------------------------------------------------------*/
  // Neumann condition for fluid that can handle inflow/backflow

  Core::Conditions::ConditionDefinition linefluidneumanninflow(
      "FLUID NEUMANN INFLOW LINE CONDITIONS", "FluidNeumannInflow", "Line Fluid Neumann Inflow",
      Core::Conditions::FluidNeumannInflow, true, Core::Conditions::geometry_type_line);
  Core::Conditions::ConditionDefinition surffluidneumanninflow(
      "FLUID NEUMANN INFLOW SURF CONDITIONS", "FluidNeumannInflow", "Surface Fluid Neumann Inflow",
      Core::Conditions::FluidNeumannInflow, true, Core::Conditions::geometry_type_surface);

  condlist.push_back(linefluidneumanninflow);
  condlist.push_back(surffluidneumanninflow);

  /*--------------------------------------------------------------------*/
  // mixed/hybrid Dirichlet conditions

  Core::Conditions::ConditionDefinition linemixhybDirichlet(
      "DESIGN LINE MIXED/HYBRID DIRICHLET CONDITIONS", "LineMixHybDirichlet", "LineMixHybDirichlet",
      Core::Conditions::LineMixHybDirichlet, true, Core::Conditions::geometry_type_line);
  Core::Conditions::ConditionDefinition surfmixhybDirichlet(
      "DESIGN SURFACE MIXED/HYBRID DIRICHLET CONDITIONS", "SurfaceMixHybDirichlet",
      "SurfaceMixHybDirichlet", Core::Conditions::SurfaceMixHybDirichlet, true,
      Core::Conditions::geometry_type_surface);

  // we attach all the components of this condition to this condition
  const auto make_mixhybDirichlet = [&](Core::Conditions::ConditionDefinition& cond)
  {
    // we provide a vector of 3 values for velocities
    cond.add_component(
        parameter<std::vector<double>>("val", {.description = "velocity", .size = 3}));

    // and optional spatial functions
    cond.add_component(parameter<std::vector<std::optional<int>>>(
        "funct", {.description = "spatial function",
                     .default_value = std::vector(3, std::optional<int>{}),
                     .size = 3}));

    // characteristic velocity
    cond.add_component(parameter<double>("u_C"));

    // the penalty parameter could be computed dynamically (using Spaldings
    // law of the wall) or using a fixed value (1)
    cond.add_component(deprecated_selection<std::string>(
        "PENTYPE", {"constant", "Spalding"}, {.description = "how penalty parameter is computed"}));

    // scaling factor for penalty parameter tauB
    cond.add_component(parameter<double>("hB_divided_by"));

    // if Spaldings law is used, this defines the way how the traction at y is computed from utau
    cond.add_component(
        deprecated_selection<std::string>("utau_computation", {"at_wall", "viscous_tangent"},
            {.description = "how traction at y is computed from utau"}));

    // we append it to the list of all conditions
    condlist.push_back(cond);
  };

  make_mixhybDirichlet(linemixhybDirichlet);
  make_mixhybDirichlet(surfmixhybDirichlet);

  /*--------------------------------------------------------------------*/
  // fluid stress

  Core::Conditions::ConditionDefinition linefluidstress("DESIGN FLUID STRESS CALC LINE CONDITIONS",
      "FluidStressCalc", "Line Fluid Stress Calculation", Core::Conditions::FluidStressCalc, true,
      Core::Conditions::geometry_type_line);
  Core::Conditions::ConditionDefinition surffluidstress("DESIGN FLUID STRESS CALC SURF CONDITIONS",
      "FluidStressCalc", "Surf Fluid Stress Calculation", Core::Conditions::FluidStressCalc, true,
      Core::Conditions::geometry_type_surface);

  condlist.push_back(linefluidstress);
  condlist.push_back(surffluidstress);

  /*--------------------------------------------------------------------*/
  // lift & drag
  Core::Conditions::ConditionDefinition surfliftdrag("DESIGN FLUID SURF LIFT&DRAG", "LIFTDRAG",
      "Surface LIFTDRAG", Core::Conditions::SurfLIFTDRAG, true,
      Core::Conditions::geometry_type_surface);

  surfliftdrag.add_component(parameter<int>("label"));
  surfliftdrag.add_component(
      parameter<std::vector<double>>("CENTER", {.description = "", .size = 3}));
  surfliftdrag.add_component(parameter<std::vector<double>>(
      "AXIS", {.description = "", .default_value = std::vector<double>{0.0, 0.0, 0.0}, .size = 3}));

  condlist.push_back(surfliftdrag);

  /*--------------------------------------------------------------------*/
  // flow rate through line

  Core::Conditions::ConditionDefinition lineflowrate("DESIGN FLOW RATE LINE CONDITIONS",
      "LineFlowRate", "Line Flow Rate", Core::Conditions::FlowRateThroughLine_2D, true,
      Core::Conditions::geometry_type_line);

  lineflowrate.add_component(parameter<int>("ConditionID"));

  condlist.push_back(lineflowrate);

  /*--------------------------------------------------------------------*/
  // flow rate through surface

  Core::Conditions::ConditionDefinition surfflowrate("DESIGN FLOW RATE SURF CONDITIONS",
      "SurfFlowRate", "Surface Flow Rate", Core::Conditions::FlowRateThroughSurface_3D, true,
      Core::Conditions::geometry_type_surface);

  surfflowrate.add_component(parameter<int>("ConditionID"));

  condlist.push_back(surfflowrate);

  /*--------------------------------------------------------------------*/
  // Volumetric surface flow profile condition
  Core::Conditions::ConditionDefinition volumetric_surface_flow_cond(
      "DESIGN SURF VOLUMETRIC FLOW CONDITIONS", "VolumetricSurfaceFlowCond",
      "volumetric surface flow condition", Core::Conditions::VolumetricSurfaceFlowCond, true,
      Core::Conditions::geometry_type_surface);

  volumetric_surface_flow_cond.add_component(parameter<int>("ConditionID"));

  volumetric_surface_flow_cond.add_component(
      deprecated_selection<std::string>("ConditionType", {"POLYNOMIAL", "WOMERSLEY"},
          {.description = "condition type", .default_value = "POLYNOMIAL"}));

  volumetric_surface_flow_cond.add_component(
      deprecated_selection<std::string>("prebiased", {"NOTPREBIASED", "PREBIASED", "FORCED"},
          {.description = "prebiased", .default_value = "NOTPREBIASED"}));

  volumetric_surface_flow_cond.add_component(deprecated_selection<std::string>(
      "FlowType", {"InFlow", "OutFlow"}, {.description = "flow type", .default_value = "InFlow"}));
  volumetric_surface_flow_cond.add_component(
      deprecated_selection<std::string>("CorrectionFlag", {"WithOutCorrection", "WithCorrection"},
          {.description = "correction flag", .default_value = "WithOutCorrection"}));

  volumetric_surface_flow_cond.add_component(parameter<double>("Period"));
  volumetric_surface_flow_cond.add_component(parameter<int>("Order"));
  volumetric_surface_flow_cond.add_component(parameter<int>("Harmonics"));
  volumetric_surface_flow_cond.add_component(parameter<double>("Val"));
  volumetric_surface_flow_cond.add_component(parameter<int>("Funct"));

  volumetric_surface_flow_cond.add_component(
      deprecated_selection<std::string>("NORMAL", {"SelfEvaluateNormal", "UsePrescribedNormal"},
          {.description = "normal", .default_value = "SelfEvaluateNormal"}));
  volumetric_surface_flow_cond.add_component(parameter<double>("n1"));
  volumetric_surface_flow_cond.add_component(parameter<double>("n2"));
  volumetric_surface_flow_cond.add_component(parameter<double>("n3"));

  volumetric_surface_flow_cond.add_component(deprecated_selection<std::string>("CenterOfMass",
      {"SelfEvaluateCenterOfMass", "UsePrescribedCenterOfMass"},
      {.description = "center of mass", .default_value = "SelfEvaluateCenterOfMass"}));
  volumetric_surface_flow_cond.add_component(parameter<double>("c1"));
  volumetric_surface_flow_cond.add_component(parameter<double>("c2"));
  volumetric_surface_flow_cond.add_component(parameter<double>("c3"));

  condlist.push_back(volumetric_surface_flow_cond);



  /*--------------------------------------------------------------------*/
  // Volumetric flow border nodes condition

  Core::Conditions::ConditionDefinition volumetric_border_nodes_cond(
      "DESIGN LINE VOLUMETRIC FLOW BORDER NODES", "VolumetricFlowBorderNodesCond",
      "volumetric flow border nodes condition", Core::Conditions::VolumetricFlowBorderNodes, true,
      Core::Conditions::geometry_type_line);

  volumetric_border_nodes_cond.add_component(parameter<int>("ConditionID"));

  condlist.push_back(volumetric_border_nodes_cond);

  /*--------------------------------------------------------------------*/
  // Volumetric surface total traction corrector
  Core::Conditions::ConditionDefinition total_traction_correction_cond(
      "DESIGN SURF TOTAL TRACTION CORRECTION CONDITIONS", "TotalTractionCorrectionCond",
      "total traction correction condition", Core::Conditions::TotalTractionCorrectionCond, true,
      Core::Conditions::geometry_type_surface);

  total_traction_correction_cond.add_component(parameter<int>("ConditionID"));
  total_traction_correction_cond.add_component(
      deprecated_selection<std::string>("ConditionType", {"POLYNOMIAL", "WOMERSLEY"},
          {.description = "condition type", .default_value = "POLYNOMIAL"}));

  total_traction_correction_cond.add_component(
      deprecated_selection<std::string>("prebiased", {"NOTPREBIASED", "PREBIASED", "FORCED"},
          {.description = "prebiased", .default_value = "NOTPREBIASED"}));

  total_traction_correction_cond.add_component(deprecated_selection<std::string>(
      "FlowType", {"InFlow", "OutFlow"}, {.description = "flow type", .default_value = "InFlow"}));
  total_traction_correction_cond.add_component(
      deprecated_selection<std::string>("CorrectionFlag", {"WithOutCorrection", "WithCorrection"},
          {.description = "correction flag", .default_value = "WithOutCorrection"}));

  total_traction_correction_cond.add_component(parameter<double>("Period"));
  total_traction_correction_cond.add_component(parameter<int>("Order"));
  total_traction_correction_cond.add_component(parameter<int>("Harmonics"));
  total_traction_correction_cond.add_component(parameter<double>("Val"));
  total_traction_correction_cond.add_component(parameter<int>("Funct"));

  total_traction_correction_cond.add_component(
      deprecated_selection<std::string>("NORMAL", {"SelfEvaluateNormal", "UsePrescribedNormal"},
          {.description = "normal", .default_value = "SelfEvaluateNormal"}));
  total_traction_correction_cond.add_component(parameter<double>("n1"));
  total_traction_correction_cond.add_component(parameter<double>("n2"));
  total_traction_correction_cond.add_component(parameter<double>("n3"));

  total_traction_correction_cond.add_component(deprecated_selection<std::string>("CenterOfMass",
      {"SelfEvaluateCenterOfMass", "UsePrescribedCenterOfMass"},
      {.description = "center of mass", .default_value = "SelfEvaluateCenterOfMass"}));
  total_traction_correction_cond.add_component(parameter<double>("c1"));
  total_traction_correction_cond.add_component(parameter<double>("c2"));
  total_traction_correction_cond.add_component(parameter<double>("c3"));

  condlist.push_back(total_traction_correction_cond);

  /*--------------------------------------------------------------------*/
  // Volumetric flow traction correction border nodes condition

  Core::Conditions::ConditionDefinition traction_corrector_border_nodes_cond(
      "DESIGN LINE TOTAL TRACTION CORRECTION BORDER NODES",
      "TotalTractionCorrectionBorderNodesCond", "total traction correction border nodes condition",
      Core::Conditions::TotalTractionCorrectionBorderNodes, true,
      Core::Conditions::geometry_type_line);

  traction_corrector_border_nodes_cond.add_component(parameter<int>("ConditionID"));

  condlist.push_back(traction_corrector_border_nodes_cond);


  /*--------------------------------------------------------------------*/
  // no penetration for darcy flow in porous media

  Core::Conditions::ConditionDefinition nopenetration_surf(
      "DESIGN SURFACE NORMAL NO PENETRATION CONDITION", "no_penetration", "No Penetration",
      Core::Conditions::no_penetration, true, Core::Conditions::geometry_type_surface);

  condlist.push_back(nopenetration_surf);

  /*--------------------------------------------------------------------*/
  // no penetration for darcy flow in porous media

  Core::Conditions::ConditionDefinition nopenetration_line(
      "DESIGN LINE NORMAL NO PENETRATION CONDITION", "no_penetration", "No Penetration",
      Core::Conditions::no_penetration, true, Core::Conditions::geometry_type_line);

  condlist.push_back(nopenetration_line);

  /*--------------------------------------------------------------------*/
  // condition for evaluation of coupling terms in porous media

  Core::Conditions::ConditionDefinition porocoupling_vol("DESIGN VOLUME POROCOUPLING CONDITION",
      "PoroCoupling", "Poro Coupling", Core::Conditions::PoroCoupling, true,
      Core::Conditions::geometry_type_volume);

  condlist.push_back(porocoupling_vol);

  /*--------------------------------------------------------------------*/
  // condition for evaluation of coupling terms in porous media

  Core::Conditions::ConditionDefinition porocoupling_surf("DESIGN SURFACE POROCOUPLING CONDITION",
      "PoroCoupling", "Poro Coupling", Core::Conditions::PoroCoupling, true,
      Core::Conditions::geometry_type_surface);

  condlist.push_back(porocoupling_surf);

  /*--------------------------------------------------------------------*/
  // condition for evaluation of boundary terms in porous media problems

  Core::Conditions::ConditionDefinition poropartint_surf("DESIGN SURFACE PORO PARTIAL INTEGRATION",
      "PoroPartInt", "Poro Partial Integration", Core::Conditions::PoroPartInt, true,
      Core::Conditions::geometry_type_surface);

  condlist.push_back(poropartint_surf);

  /*--------------------------------------------------------------------*/
  // condition for evaluation of boundary terms in porous media problems

  Core::Conditions::ConditionDefinition poropartint_line("DESIGN LINE PORO PARTIAL INTEGRATION",
      "PoroPartInt", "Poro Partial Integration", Core::Conditions::PoroPartInt, true,
      Core::Conditions::geometry_type_line);

  condlist.push_back(poropartint_line);

  /*--------------------------------------------------------------------*/
  // condition for evaluation of boundary terms in porous media problems

  Core::Conditions::ConditionDefinition poropresint_surf("DESIGN SURFACE PORO PRESSURE INTEGRATION",
      "PoroPresInt", "Poro Pressure Integration", Core::Conditions::PoroPresInt, true,
      Core::Conditions::geometry_type_surface);

  condlist.push_back(poropresint_surf);

  /*--------------------------------------------------------------------*/
  // condition for evaluation of boundary terms in porous media problems

  Core::Conditions::ConditionDefinition poropresint_line("DESIGN LINE PORO PRESSURE INTEGRATION",
      "PoroPresInt", "Poro Pressure Integration", Core::Conditions::PoroPresInt, true,
      Core::Conditions::geometry_type_line);

  condlist.push_back(poropresint_line);
}

std::string Inpar::FLUID::to_string(Inpar::FLUID::StabType stabtype)
{
  switch (stabtype)
  {
    case stabtype_nostab:
      return "no_stabilization";
    case stabtype_residualbased:
      return "residual_based";
    case stabtype_edgebased:
      return "edge_based";
    case stabtype_pressureprojection:
      return "pressure_projection";
    default:
      return "unknown";
  }
}

std::string Inpar::FLUID::to_string(Inpar::FLUID::TauType tau)
{
  switch (tau)
  {
    case tau_taylor_hughes_zarins:
      return "Taylor_Hughes_Zarins";
    case tau_taylor_hughes_zarins_wo_dt:
      return "Taylor_Hughes_Zarins_wo_dt";
    case tau_taylor_hughes_zarins_whiting_jansen:
      return "Taylor_Hughes_Zarins_Whiting_Jansen";
    case tau_taylor_hughes_zarins_whiting_jansen_wo_dt:
      return "Taylor_Hughes_Zarins_Whiting_Jansen_wo_dt";
    case tau_taylor_hughes_zarins_scaled:
      return "Taylor_Hughes_Zarins_scaled";
    case tau_taylor_hughes_zarins_scaled_wo_dt:
      return "Taylor_Hughes_Zarins_scaled_wo_dt";
    case tau_franca_barrenechea_valentin_frey_wall:
      return "Franca_Barrenechea_Valentin_Frey_Wall";
    case tau_franca_barrenechea_valentin_frey_wall_wo_dt:
      return "Franca_Barrenechea_Valentin_Frey_Wall_wo_dt";
    case tau_shakib_hughes_codina:
      return "Shakib_Hughes_Codina";
    case tau_shakib_hughes_codina_wo_dt:
      return "Shakib_Hughes_Codina_wo_dt";
    case tau_codina:
      return "Codina";
    case tau_codina_wo_dt:
      return "Codina_wo_dt";
    case tau_codina_convscaled:
      return "Codina_convscaled";
    case tau_franca_madureira_valentin_badia_codina:
      return "Franca_Madureira_Valentin_Badia_Codina";
    case tau_franca_madureira_valentin_badia_codina_wo_dt:
      return "Franca_Madureira_Valentin_Badia_Codina_wo_dt";
    case tau_hughes_franca_balestra_wo_dt:
      return "Hughes_Franca_Balestra_wo_dt";
    case tau_not_defined:
      return "not_defined";
  }
  return "not_defined";
}

std::string Inpar::FLUID::to_string(Inpar::FLUID::CrossStress crossstress)
{
  switch (crossstress)
  {
    case cross_stress_stab_none:
      return "no_cross";
    case cross_stress_stab:
      return "yes_cross";
    case cross_stress_stab_only_rhs:
      return "cross_rhs";
    default:
      return "unknown";
  }
}

std::string Inpar::FLUID::to_string(Inpar::FLUID::VStab vstab)
{
  switch (vstab)
  {
    case viscous_stab_none:
      return "no_viscous";
    case viscous_stab_gls:
      return "yes_viscous";
    case viscous_stab_gls_only_rhs:
      return "viscous_rhs";
    case viscous_stab_usfem:
      return "usfem_viscous";
    case viscous_stab_usfem_only_rhs:
      return "usfem_viscous_rhs";
    default:
      return "unknown";
  }
}

std::string Inpar::FLUID::to_string(Inpar::FLUID::RStab rstab)
{
  switch (rstab)
  {
    case reactive_stab_none:
      return "no_reactive";
    case reactive_stab_gls:
      return "yes_reactive";
    case reactive_stab_usfem:
      return "usfem_reactive";
    default:
      return "unknown";
  }
}

std::string Inpar::FLUID::to_string(Inpar::FLUID::ReynoldsStress reynoldsstress)
{
  switch (reynoldsstress)
  {
    case reynolds_stress_stab_none:
      return "no_reynolds";
    case reynolds_stress_stab:
      return "yes_reynolds";
    case reynolds_stress_stab_only_rhs:
      return "reynolds_rhs";
    default:
      return "unknown";
  }
}

std::string Inpar::FLUID::to_string(Inpar::FLUID::Transient transient)
{
  switch (transient)
  {
    case inertia_stab_drop:
      return "no_transient";
    case inertia_stab_keep:
      return "yes_transient";
    case inertia_stab_keep_complete:
      return "transient_complete";
    default:
      return "unknown";
  }
}

std::string Inpar::FLUID::to_string(Inpar::FLUID::EosPres eospres)
{
  switch (eospres)
  {
    case EOS_PRES_none:
      return "none";
    case EOS_PRES_std_eos:
      return "std_eos";
    case EOS_PRES_xfem_gp:
      return "xfem_gp";
    default:
      return "unknown";
  }
}

std::string Inpar::FLUID::to_string(Inpar::FLUID::EosConvStream eosconvstream)
{
  switch (eosconvstream)
  {
    case EOS_CONV_STREAM_none:
      return "none";
    case EOS_CONV_STREAM_std_eos:
      return "std_eos";
    case EOS_CONV_STREAM_xfem_gp:
      return "xfem_gp";
    default:
      return "unknown";
  }
}

std::string Inpar::FLUID::to_string(Inpar::FLUID::EosConvCross eosconvcross)
{
  switch (eosconvcross)
  {
    case EOS_CONV_CROSS_none:
      return "none";
    case EOS_CONV_CROSS_std_eos:
      return "std_eos";
    case EOS_CONV_CROSS_xfem_gp:
      return "xfem_gp";
    default:
      return "unknown";
  }
}

std::string Inpar::FLUID::to_string(Inpar::FLUID::EosDiv eosdiv)
{
  switch (eosdiv)
  {
    case EOS_DIV_none:
      return "none";
    case EOS_DIV_vel_jump_std_eos:
      return "vel_jump_std_eos";
    case EOS_DIV_vel_jump_xfem_gp:
      return "vel_jump_xfem_gp";
    case EOS_DIV_div_jump_std_eos:
      return "div_jump_std_eos";
    case EOS_DIV_div_jump_xfem_gp:
      return "div_jump_xfem_gp";
    default:
      return "unknown";
  }
}

std::string Inpar::FLUID::to_string(Inpar::FLUID::EosTauType eostautype)
{
  switch (eostautype)
  {
    case EOS_tau_burman:
      return "Burman";
    case EOS_tau_burman_fernandez_hansbo:
      return "Burman_Fernandez_Hansbo";
    case EOS_tau_burman_fernandez_hansbo_wo_dt:
      return "Burman_Fernandez_Hansbo_wo_dt";
    case EOS_tau_braack_burman_john_lube:
      return "Braack_Burman_John_Lube";
    case EOS_tau_braack_burman_john_lube_wo_divjump:
      return "Braack_Burman_John_Lube_wo_divjump";
    case EOS_tau_franca_barrenechea_valentin_wall:
      return "Franca_Barrenechea_Valentin_Wall";
    case EOS_tau_burman_fernandez:
      return "Burman_Fernandez";
    case EOS_tau_burman_hansbo_dangelo_zunino:
      return "Burman_Hansbo_DAngelo_Zunino";
    case EOS_tau_burman_hansbo_dangelo_zunino_wo_dt:
      return "Burman_Hansbo_DAngelo_Zunino_wo_dt";
    case EOS_tau_schott_massing_burman_dangelo_zunino:
      return "Schott_Massing_Burman_DAngelo_Zunino";
    case EOS_tau_schott_massing_burman_dangelo_zunino_wo_dt:
      return "Schott_Massing_Burman_DAngelo_Zunino_wo_dt";
    case EOS_tau_Taylor_Hughes_Zarins_Whiting_Jansen_Codina_scaling:
      return "Taylor_Hughes_Zarins_Whiting_Jansen_Codina_scaling";
    case EOS_tau_poroelast_fluid:
      return "Poroelast_Fluid";
    case EOS_tau_not_defined:
      return "not_defined";
    default:
      return "unknown";
  }
}

std::string Inpar::FLUID::to_string(Inpar::FLUID::EosElementLength eoselementlength)
{
  switch (eoselementlength)
  {
    case EOS_he_max_diameter_to_opp_surf:
      return "EOS_he_max_diameter_to_opp_surf";
    case EOS_he_max_dist_to_opp_surf:
      return "EOS_he_max_dist_to_opp_surf";
    case EOS_he_surf_with_max_diameter:
      return "EOS_he_surf_with_max_diameter";
    case EOS_hk_max_diameter:
      return "EOS_hk_max_diameter";
    case EOS_he_surf_diameter:
      return "EOS_he_surf_diameter";
    case EOS_he_vol_eq_diameter:
      return "EOS_he_vol_eq_diameter";
    default:
      return "unknown";
  }
}

std::string Inpar::FLUID::to_string(Inpar::FLUID::VremanFiMethod vremanfimethod)
{
  switch (vremanfimethod)
  {
    case cuberootvol:
      return "CubeRootVol";
    case dir_dep:
      return "Direction_dependent";
    case min_len:
      return "Minimum_length";
    default:
      return "unknown";
  }
}

FOUR_C_NAMESPACE_CLOSE
