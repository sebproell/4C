// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_inpar_particle.hpp"

#include "4C_fem_condition_definition.hpp"
#include "4C_io_input_spec_builders.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | set the particle parameters                                               |
 *---------------------------------------------------------------------------*/
void Inpar::PARTICLE::set_valid_parameters(std::map<std::string, Core::IO::InputSpec>& list)
{
  using Teuchos::tuple;
  using namespace Core::IO::InputSpecBuilders;

  /*-------------------------------------------------------------------------*
   | general control parameters for particle simulations                     |
   *-------------------------------------------------------------------------*/
  Core::Utils::SectionSpecs particledyn{"PARTICLE DYNAMIC"};

  // type of particle time integration
  particledyn.specs.emplace_back(deprecated_selection<DynamicType>("DYNAMICTYPE",
      {
          {"SemiImplicitEuler", Inpar::PARTICLE::dyna_semiimpliciteuler},
          {"VelocityVerlet", Inpar::PARTICLE::dyna_velocityverlet},
      },
      {.description = "type of particle time integration",
          .default_value = Inpar::PARTICLE::dyna_velocityverlet}));

  // type of particle interaction
  particledyn.specs.emplace_back(deprecated_selection<InteractionType>("INTERACTION",
      {
          {"None", Inpar::PARTICLE::interaction_none},
          {"SPH", Inpar::PARTICLE::interaction_sph},
          {"DEM", Inpar::PARTICLE::interaction_dem},
      },
      {.description = "type of particle interaction",
          .default_value = Inpar::PARTICLE::interaction_none}));

  // output type
  particledyn.specs.emplace_back(parameter<int>(
      "RESULTSEVERY", {.description = "write particle runtime output every RESULTSEVERY steps",
                          .default_value = 1}));
  particledyn.specs.emplace_back(parameter<int>("RESTARTEVERY",
      {.description = "write restart possibility every RESTARTEVERY steps", .default_value = 1}));

  // write ghosted particles
  particledyn.specs.emplace_back(parameter<bool>("WRITE_GHOSTED_PARTICLES",
      {.description = "write ghosted particles (debug feature)", .default_value = false}));

  // time loop control
  particledyn.specs.emplace_back(
      parameter<double>("TIMESTEP", {.description = "time step size", .default_value = 0.01}));
  particledyn.specs.emplace_back(
      parameter<int>("NUMSTEP", {.description = "maximum number of steps", .default_value = 100}));
  particledyn.specs.emplace_back(
      parameter<double>("MAXTIME", {.description = "maximum time", .default_value = 1.0}));

  // gravity acceleration control
  particledyn.specs.emplace_back(parameter<std::string>("GRAVITY_ACCELERATION",
      {.description = "acceleration due to gravity", .default_value = "0.0 0.0 0.0"}));
  particledyn.specs.emplace_back(parameter<int>("GRAVITY_RAMP_FUNCT",
      {.description = "number of function governing gravity ramp", .default_value = -1}));

  // viscous damping factor
  particledyn.specs.emplace_back(parameter<double>("VISCOUS_DAMPING",
      {.description = "apply viscous damping force to determine static equilibrium solutions",
          .default_value = -1.0}));

  // transfer particles to new bins every time step
  particledyn.specs.emplace_back(parameter<bool>("TRANSFER_EVERY",
      {.description = "transfer particles to new bins every time step", .default_value = false}));

  // considered particle phases with dynamic load balance weighting factor
  particledyn.specs.emplace_back(parameter<std::string>("PHASE_TO_DYNLOADBALFAC",
      {.description = "considered particle phases with dynamic load balance weighting factor",
          .default_value = "none"}));

  // relate particle phase to material id
  particledyn.specs.emplace_back(parameter<std::string>("PHASE_TO_MATERIAL_ID",
      {.description = "relate particle phase to material id", .default_value = "none"}));

  // amplitude of noise added to initial position for each spatial direction
  particledyn.specs.emplace_back(parameter<std::string>("INITIAL_POSITION_AMPLITUDE",
      {.description = "amplitude of noise added to initial position for each spatial direction",
          .default_value = "0.0 0.0 0.0"}));

  // type of particle wall source
  particledyn.specs.emplace_back(deprecated_selection<ParticleWallSource>("PARTICLE_WALL_SOURCE",
      {
          {"NoParticleWall", Inpar::PARTICLE::NoParticleWall},
          {"DiscretCondition", Inpar::PARTICLE::DiscretCondition},
          {"BoundingBox", Inpar::PARTICLE::BoundingBox},
      },
      {.description = "type of particle wall source",
          .default_value = Inpar::PARTICLE::NoParticleWall}));

  // material id for particle wall from bounding box source
  particledyn.specs.emplace_back(parameter<int>(
      "PARTICLE_WALL_MAT", {.description = "material id for particle wall from bounding box source",
                               .default_value = -1}));

  // flags defining considered states of particle wall
  particledyn.specs.emplace_back(parameter<bool>("PARTICLE_WALL_MOVING",
      {.description = "consider a moving particle wall", .default_value = false}));
  particledyn.specs.emplace_back(parameter<bool>("PARTICLE_WALL_LOADED",
      {.description = "consider loading on particle wall", .default_value = false}));

  // consider rigid body motion
  particledyn.specs.emplace_back(parameter<bool>(
      "RIGID_BODY_MOTION", {.description = "consider rigid body motion", .default_value = false}));

  particledyn.specs.emplace_back(parameter<double>("RIGID_BODY_PHASECHANGE_RADIUS",
      {.description = "search radius for neighboring rigid bodies in case of phase change",
          .default_value = -1.0}));

  particledyn.move_into_collection(list);

  /*-------------------------------------------------------------------------*
   | control parameters for initial/boundary conditions                      |
   *-------------------------------------------------------------------------*/
  Core::Utils::SectionSpecs particledynconditions{particledyn, "INITIAL AND BOUNDARY CONDITIONS"};

  // initial temperature field of particle phase given by function
  particledynconditions.specs.emplace_back(parameter<std::string>("INITIAL_TEMP_FIELD",
      {.description =
              "Refer to the function ID describing the initial temperature field of particle phase",
          .default_value = "none"}));

  // initial velocity field of particle phase given by function
  particledynconditions.specs.emplace_back(parameter<std::string>("INITIAL_VELOCITY_FIELD",
      {.description =
              "Refer to the function ID describing the initial velocity field of particle phase",
          .default_value = "none"}));

  // initial angular velocity field of particle phase given by function
  particledynconditions.specs.emplace_back(parameter<std::string>("INITIAL_ANGULAR_VELOCITY_FIELD",
      {.description = "Refer to the function ID describing the initial angular velocity field of "
                      "rigid body  phase/DEM particle",
          .default_value = "none"}));


  // initial acceleration field of particle phase given by function
  particledynconditions.specs.emplace_back(parameter<std::string>(
      "INITIAL_ACCELERATION_FIELD", {.description = "Refer to the function ID describing the "
                                                    "initial acceleration field of particle phase",
                                        .default_value = "none"}));

  // initial angular acceleration field of particle phase given by function
  particledynconditions.specs.emplace_back(
      parameter<std::string>("INITIAL_ANGULAR_ACCELERATION_FIELD",
          {.description = "Refer to the function ID describing the initial angular acceleration "
                          "field of rigid body  phase/DEM particle",
              .default_value = "none"}));


  // dirichlet boundary condition of particle phase given by function
  particledynconditions.specs.emplace_back(parameter<std::string>("DIRICHLET_BOUNDARY_CONDITION",
      {.description = "Refer to the function ID describing the dirichlet boundary condition of "
                      "particle phase",
          .default_value = "none"}));

  // temperature boundary condition of particle phase given by function
  particledynconditions.specs.emplace_back(parameter<std::string>("TEMPERATURE_BOUNDARY_CONDITION",
      {.description = "Refer to the function ID describing the temperature boundary condition of "
                      "particle phase",
          .default_value = "none"}));

  particledynconditions.move_into_collection(list);

  /*-------------------------------------------------------------------------*
   | smoothed particle hydrodynamics (SPH) specific control parameters       |
   *-------------------------------------------------------------------------*/
  Core::Utils::SectionSpecs particledynsph{particledyn, "SPH"};

  // write particle-wall interaction output
  particledynsph.specs.emplace_back(parameter<bool>("WRITE_PARTICLE_WALL_INTERACTION",
      {.description = "write particle-wall interaction output", .default_value = false}));

  // type of smoothed particle hydrodynamics kernel
  particledynsph.specs.emplace_back(deprecated_selection<KernelType>("KERNEL",
      {
          {"CubicSpline", Inpar::PARTICLE::CubicSpline},
          {"QuinticSpline", Inpar::PARTICLE::QuinticSpline},
      },
      {.description = "type of smoothed particle hydrodynamics kernel",
          .default_value = Inpar::PARTICLE::CubicSpline}));

  // kernel space dimension number
  particledynsph.specs.emplace_back(deprecated_selection<KernelSpaceDimension>("KERNEL_SPACE_DIM",
      {
          {"Kernel1D", Inpar::PARTICLE::Kernel1D},
          {"Kernel2D", Inpar::PARTICLE::Kernel2D},
          {"Kernel3D", Inpar::PARTICLE::Kernel3D},
      },
      {.description = "kernel space dimension number",
          .default_value = Inpar::PARTICLE::Kernel3D}));

  particledynsph.specs.emplace_back(parameter<double>("INITIALPARTICLESPACING",
      {.description = "initial spacing of particles", .default_value = 0.0}));

  // type of smoothed particle hydrodynamics equation of state
  particledynsph.specs.emplace_back(deprecated_selection<EquationOfStateType>("EQUATIONOFSTATE",
      {
          {"GenTait", Inpar::PARTICLE::GenTait},
          {"IdealGas", Inpar::PARTICLE::IdealGas},
      },
      {.description = "type of smoothed particle hydrodynamics equation of state",
          .default_value = Inpar::PARTICLE::GenTait}));

  // type of smoothed particle hydrodynamics momentum formulation
  particledynsph.specs.emplace_back(
      deprecated_selection<MomentumFormulationType>("MOMENTUMFORMULATION",
          {
              {"AdamiMomentumFormulation", Inpar::PARTICLE::AdamiMomentumFormulation},
              {"MonaghanMomentumFormulation", Inpar::PARTICLE::MonaghanMomentumFormulation},
          },
          {.description = "type of smoothed particle hydrodynamics momentum formulation",
              .default_value = Inpar::PARTICLE::AdamiMomentumFormulation}));

  // type of density evaluation scheme
  particledynsph.specs.emplace_back(
      deprecated_selection<DensityEvaluationScheme>("DENSITYEVALUATION",
          {
              {"DensitySummation", Inpar::PARTICLE::DensitySummation},
              {"DensityIntegration", Inpar::PARTICLE::DensityIntegration},
              {"DensityPredictCorrect", Inpar::PARTICLE::DensityPredictCorrect},
          },
          {.description = "type of density evaluation scheme",
              .default_value = Inpar::PARTICLE::DensitySummation}));

  // type of density correction scheme
  particledynsph.specs.emplace_back(
      deprecated_selection<DensityCorrectionScheme>("DENSITYCORRECTION",
          {
              {"NoCorrection", Inpar::PARTICLE::NoCorrection},
              {"InteriorCorrection", Inpar::PARTICLE::InteriorCorrection},
              {"NormalizedCorrection", Inpar::PARTICLE::NormalizedCorrection},
              {"RandlesCorrection", Inpar::PARTICLE::RandlesCorrection},
          },
          {.description = "type of density correction scheme",
              .default_value = Inpar::PARTICLE::NoCorrection}));

  // type of boundary particle formulation
  particledynsph.specs.emplace_back(
      deprecated_selection<BoundaryParticleFormulationType>("BOUNDARYPARTICLEFORMULATION",
          {
              {"NoBoundaryFormulation", Inpar::PARTICLE::NoBoundaryFormulation},
              {"AdamiBoundaryFormulation", Inpar::PARTICLE::AdamiBoundaryFormulation},
          },
          {.description = "type of boundary particle formulation",
              .default_value = Inpar::PARTICLE::NoBoundaryFormulation}));

  // type of boundary particle interaction
  particledynsph.specs.emplace_back(
      deprecated_selection<BoundaryParticleInteraction>("BOUNDARYPARTICLEINTERACTION",
          {
              {"NoSlipBoundaryParticle", Inpar::PARTICLE::NoSlipBoundaryParticle},
              {"FreeSlipBoundaryParticle", Inpar::PARTICLE::FreeSlipBoundaryParticle},
          },
          {.description = "type of boundary particle interaction",
              .default_value = Inpar::PARTICLE::NoSlipBoundaryParticle}));

  // type of wall formulation
  particledynsph.specs.emplace_back(deprecated_selection<WallFormulationType>("WALLFORMULATION",
      {
          {"NoWallFormulation", Inpar::PARTICLE::NoWallFormulation},
          {"VirtualParticleWallFormulation", Inpar::PARTICLE::VirtualParticleWallFormulation},
      },
      {.description = "type of wall formulation",
          .default_value = Inpar::PARTICLE::NoWallFormulation}));

  // type of transport velocity formulation
  particledynsph.specs.emplace_back(
      deprecated_selection<TransportVelocityFormulation>("TRANSPORTVELOCITYFORMULATION",
          {
              {"NoTransportVelocity", Inpar::PARTICLE::NoTransportVelocity},
              {"StandardTransportVelocity", Inpar::PARTICLE::StandardTransportVelocity},
              {"GeneralizedTransportVelocity", Inpar::PARTICLE::GeneralizedTransportVelocity},
          },
          {.description = "type of transport velocity formulation",
              .default_value = Inpar::PARTICLE::NoTransportVelocity}));

  // type of temperature evaluation scheme
  particledynsph.specs.emplace_back(
      deprecated_selection<TemperatureEvaluationScheme>("TEMPERATUREEVALUATION",
          {
              {"NoTemperatureEvaluation", Inpar::PARTICLE::NoTemperatureEvaluation},
              {"TemperatureIntegration", Inpar::PARTICLE::TemperatureIntegration},
          },
          {.description = "type of temperature evaluation scheme",
              .default_value = Inpar::PARTICLE::NoTemperatureEvaluation}));

  particledynsph.specs.emplace_back(parameter<bool>("TEMPERATUREGRADIENT",
      {.description = "evaluate temperature gradient", .default_value = false}));

  // type of heat source
  particledynsph.specs.emplace_back(deprecated_selection<HeatSourceType>("HEATSOURCETYPE",
      {
          {"NoHeatSource", Inpar::PARTICLE::NoHeatSource},
          {"VolumeHeatSource", Inpar::PARTICLE::VolumeHeatSource},
          {"SurfaceHeatSource", Inpar::PARTICLE::SurfaceHeatSource},
      },
      {.description = "type of heat source", .default_value = Inpar::PARTICLE::NoHeatSource}));

  particledynsph.specs.emplace_back(parameter<int>("HEATSOURCE_FUNCT",
      {.description = "number of function governing heat source", .default_value = -1}));

  particledynsph.specs.emplace_back(parameter<std::string>("HEATSOURCE_DIRECTION",
      {.description = "direction of surface heat source", .default_value = "0.0 0.0 0.0"}));

  // evaporation induced heat loss
  particledynsph.specs.emplace_back(parameter<bool>("VAPOR_HEATLOSS",
      {.description = "evaluate evaporation induced heat loss", .default_value = false}));
  particledynsph.specs.emplace_back(parameter<double>("VAPOR_HEATLOSS_LATENTHEAT",
      {.description = "latent heat in heat loss formula", .default_value = 0.0}));
  particledynsph.specs.emplace_back(parameter<double>("VAPOR_HEATLOSS_ENTHALPY_REFTEMP",
      {.description = "enthalpy reference temperature in heat loss formula",
          .default_value = 0.0}));
  particledynsph.specs.emplace_back(parameter<double>("VAPOR_HEATLOSS_PFAC",
      {.description = "pressure factor in heat loss formula", .default_value = 0.0}));
  particledynsph.specs.emplace_back(parameter<double>("VAPOR_HEATLOSS_TFAC",
      {.description = "temperature factor in heat loss formula", .default_value = 0.0}));

  // evaporation induced recoil pressure
  particledynsph.specs.emplace_back(parameter<bool>("VAPOR_RECOIL",
      {.description = "evaluate evaporation induced recoil pressure", .default_value = false}));
  particledynsph.specs.emplace_back(parameter<double>("VAPOR_RECOIL_BOILINGTEMPERATURE",
      {.description = "boiling temperature in recoil pressure formula", .default_value = 0.0}));
  particledynsph.specs.emplace_back(parameter<double>("VAPOR_RECOIL_PFAC",
      {.description = "pressure factor in recoil pressure formula", .default_value = 0.0}));
  particledynsph.specs.emplace_back(parameter<double>("VAPOR_RECOIL_TFAC",
      {.description = "temperature factor in recoil pressure formula", .default_value = 0.0}));

  // type of surface tension formulation
  particledynsph.specs.emplace_back(
      deprecated_selection<SurfaceTensionFormulation>("SURFACETENSIONFORMULATION",
          {
              {"NoSurfaceTension", Inpar::PARTICLE::NoSurfaceTension},
              {"ContinuumSurfaceForce", Inpar::PARTICLE::ContinuumSurfaceForce},
          },
          {.description = "type of surface tension formulation",
              .default_value = Inpar::PARTICLE::NoSurfaceTension}));

  particledynsph.specs.emplace_back(parameter<int>("SURFACETENSION_RAMP_FUNCT",
      {.description = "number of function governing surface tension ramp", .default_value = -1}));

  particledynsph.specs.emplace_back(parameter<double>("SURFACETENSIONCOEFFICIENT",
      {.description = "constant part of surface tension coefficient", .default_value = -1.0}));
  particledynsph.specs.emplace_back(parameter<double>("SURFACETENSIONMINIMUM",
      {.description = "minimum surface tension coefficient in case of temperature dependence",
          .default_value = 0.0}));
  particledynsph.specs.emplace_back(parameter<double>("SURFACETENSIONTEMPFAC",
      {.description = "factor of dependence of surface tension coefficient on temperature",
          .default_value = 0.0}));
  particledynsph.specs.emplace_back(parameter<double>("SURFACETENSIONREFTEMP",
      {.description = "reference temperature for surface tension coefficient",
          .default_value = 0.0}));

  // wetting
  particledynsph.specs.emplace_back(parameter<double>(
      "STATICCONTACTANGLE", {.description = "static contact angle in degree with wetting effects",
                                .default_value = 0.0}));
  particledynsph.specs.emplace_back(parameter<double>("TRIPLEPOINTNORMAL_CORR_CF_LOW",
      {.description = "triple point normal correction wall color field low",
          .default_value = 0.0}));
  particledynsph.specs.emplace_back(parameter<double>("TRIPLEPOINTNORMAL_CORR_CF_UP",
      {.description = "triple point normal correction wall color field up", .default_value = 0.0}));

  // interface viscosity
  particledynsph.specs.emplace_back(parameter<bool>("INTERFACE_VISCOSITY",
      {.description = "evaluate interface viscosity", .default_value = false}));
  particledynsph.specs.emplace_back(parameter<double>("INTERFACE_VISCOSITY_LIQUIDGAS",
      {.description = "artificial viscosity on liquid-gas interface", .default_value = 0.0}));
  particledynsph.specs.emplace_back(parameter<double>("INTERFACE_VISCOSITY_SOLIDLIQUID",
      {.description = "artificial viscosity on solid-liquid interface", .default_value = 0.0}));

  // barrier force
  particledynsph.specs.emplace_back(parameter<bool>(
      "BARRIER_FORCE", {.description = "evaluate barrier force", .default_value = false}));
  particledynsph.specs.emplace_back(parameter<double>(
      "BARRIER_FORCE_DISTANCE", {.description = "barrier force distance", .default_value = 0.0}));
  particledynsph.specs.emplace_back(parameter<double>("BARRIER_FORCE_TEMPSCALE",
      {.description = "barrier force temperature scaling", .default_value = 0.0}));
  particledynsph.specs.emplace_back(parameter<double>("BARRIER_FORCE_STIFF_HEAVY",
      {.description = "barrier force stiffness of heavy phase", .default_value = -1.0}));
  particledynsph.specs.emplace_back(parameter<double>("BARRIER_FORCE_DAMP_HEAVY",
      {.description = "barrier force damping parameter of heavy phase", .default_value = 0.0}));
  particledynsph.specs.emplace_back(parameter<double>("BARRIER_FORCE_STIFF_GAS",
      {.description = "barrier force stiffness of gas phase", .default_value = -1.0}));
  particledynsph.specs.emplace_back(parameter<double>("BARRIER_FORCE_DAMP_GAS",
      {.description = "barrier force damping parameter of gas phase", .default_value = 0.0}));

  // linear transition in surface tension evaluation
  particledynsph.specs.emplace_back(parameter<double>("TRANS_REF_TEMPERATURE",
      {.description = "transition reference temperature", .default_value = 0.0}));
  particledynsph.specs.emplace_back(parameter<double>("TRANS_DT_SURFACETENSION",
      {.description = "transition temperature difference for surface tension evaluation",
          .default_value = 0.0}));
  particledynsph.specs.emplace_back(parameter<double>("TRANS_DT_MARANGONI",
      {.description = "transition temperature difference for marangoni evaluation",
          .default_value = 0.0}));
  particledynsph.specs.emplace_back(parameter<double>("TRANS_DT_CURVATURE",
      {.description = "transition temperature difference for curvature evaluation",
          .default_value = 0.0}));
  particledynsph.specs.emplace_back(parameter<double>("TRANS_DT_WETTING",
      {.description = "transition temperature difference for wetting evaluation",
          .default_value = 0.0}));
  particledynsph.specs.emplace_back(parameter<double>("TRANS_DT_INTVISC",
      {.description = "transition temperature difference for interface viscosity evaluation",
          .default_value = 0.0}));
  particledynsph.specs.emplace_back(parameter<double>("TRANS_DT_BARRIER",
      {.description = "transition temperature difference for barrier force evaluation",
          .default_value = 0.0}));

  // type of dirichlet open boundary
  particledynsph.specs.emplace_back(
      deprecated_selection<DirichletOpenBoundaryType>("DIRICHLETBOUNDARYTYPE",
          {
              {"NoDirichletOpenBoundary", Inpar::PARTICLE::NoDirichletOpenBoundary},
              {"DirichletNormalToPlane", Inpar::PARTICLE::DirichletNormalToPlane},
          },
          {.description = "type of dirichlet open boundary",
              .default_value = Inpar::PARTICLE::NoDirichletOpenBoundary}));

  particledynsph.specs.emplace_back(parameter<int>("DIRICHLET_FUNCT",
      {.description = "number of function governing velocity condition on dirichlet open boundary",
          .default_value = -1}));

  particledynsph.specs.emplace_back(parameter<std::string>("DIRICHLET_OUTWARD_NORMAL",
      {.description = "direction of outward normal on dirichlet open boundary",
          .default_value = "0.0 0.0 0.0"}));
  particledynsph.specs.emplace_back(parameter<std::string>("DIRICHLET_PLANE_POINT",
      {.description = "point on dirichlet open boundary plane", .default_value = "0.0 0.0 0.0"}));

  // type of neumann open boundary
  particledynsph.specs.emplace_back(
      deprecated_selection<NeumannOpenBoundaryType>("NEUMANNBOUNDARYTYPE",
          {
              {"NoNeumannOpenBoundary", Inpar::PARTICLE::NoNeumannOpenBoundary},
              {"NeumannNormalToPlane", Inpar::PARTICLE::NeumannNormalToPlane},
          },
          {.description = "type of neumann open boundary",
              .default_value = Inpar::PARTICLE::NoNeumannOpenBoundary}));

  particledynsph.specs.emplace_back(parameter<int>("NEUMANN_FUNCT",
      {.description = "number of function governing pressure condition on neumann open boundary",
          .default_value = -1}));

  particledynsph.specs.emplace_back(parameter<std::string>("NEUMANN_OUTWARD_NORMAL",
      {.description = "direction of outward normal on neumann open boundary",
          .default_value = "0.0 0.0 0.0"}));
  particledynsph.specs.emplace_back(parameter<std::string>("NEUMANN_PLANE_POINT",
      {.description = "point on neumann open boundary plane", .default_value = "0.0 0.0 0.0"}));

  // type of phase change
  particledynsph.specs.emplace_back(deprecated_selection<PhaseChangeType>("PHASECHANGETYPE",
      {
          {"NoPhaseChange", Inpar::PARTICLE::NoPhaseChange},
          {"OneWayScalarBelowToAbovePhaseChange",
              Inpar::PARTICLE::OneWayScalarBelowToAbovePhaseChange},
          {"OneWayScalarAboveToBelowPhaseChange",
              Inpar::PARTICLE::OneWayScalarAboveToBelowPhaseChange},
          {"TwoWayScalarPhaseChange", Inpar::PARTICLE::TwoWayScalarPhaseChange},
      },
      {.description = "type of phase change", .default_value = Inpar::PARTICLE::NoPhaseChange}));

  // definition of phase change
  particledynsph.specs.emplace_back(parameter<std::string>("PHASECHANGEDEFINITION",
      {.description = "phase change definition", .default_value = "none"}));

  // type of rigid particle contact
  particledynsph.specs.emplace_back(
      deprecated_selection<RigidParticleContactType>("RIGIDPARTICLECONTACTTYPE",
          {
              {"NoRigidParticleContact", Inpar::PARTICLE::NoRigidParticleContact},
              {"ElasticRigidParticleContact", Inpar::PARTICLE::ElasticRigidParticleContact},
          },
          {.description = "type of rigid particle contact",
              .default_value = Inpar::PARTICLE::NoRigidParticleContact}));

  particledynsph.specs.emplace_back(parameter<double>("RIGIDPARTICLECONTACTSTIFF",
      {.description = "rigid particle contact stiffness", .default_value = -1.0}));
  particledynsph.specs.emplace_back(parameter<double>("RIGIDPARTICLECONTACTDAMP",
      {.description = "rigid particle contact damping parameter", .default_value = 0.0}));

  particledynsph.move_into_collection(list);

  /*-------------------------------------------------------------------------*
   | discrete element method (DEM) specific control parameters               |
   *-------------------------------------------------------------------------*/
  Core::Utils::SectionSpecs particledyndem{particledyn, "DEM"};

  // write particle energy output
  particledyndem.specs.emplace_back(parameter<bool>("WRITE_PARTICLE_ENERGY",
      {.description = "write particle energy output", .default_value = false}));

  // write particle-wall interaction output
  particledyndem.specs.emplace_back(parameter<bool>("WRITE_PARTICLE_WALL_INTERACTION",
      {.description = "write particle-wall interaction output", .default_value = false}));

  // type of normal contact law
  particledyndem.specs.emplace_back(deprecated_selection<NormalContact>("NORMALCONTACTLAW",
      {
          {"NormalLinearSpring", Inpar::PARTICLE::NormalLinSpring},
          {"NormalLinearSpringDamp", Inpar::PARTICLE::NormalLinSpringDamp},
          {"NormalHertz", Inpar::PARTICLE::NormalHertz},
          {"NormalLeeHerrmann", Inpar::PARTICLE::NormalLeeHerrmann},
          {"NormalKuwabaraKono", Inpar::PARTICLE::NormalKuwabaraKono},
          {"NormalTsuji", Inpar::PARTICLE::NormalTsuji},
      },
      {.description = "normal contact law for particles",
          .default_value = Inpar::PARTICLE::NormalLinSpring}));

  // type of tangential contact law
  particledyndem.specs.emplace_back(deprecated_selection<TangentialContact>("TANGENTIALCONTACTLAW",
      {
          {"NoTangentialContact", Inpar::PARTICLE::NoTangentialContact},
          {"TangentialLinSpringDamp", Inpar::PARTICLE::TangentialLinSpringDamp},
      },
      {.description = "tangential contact law for particles",
          .default_value = Inpar::PARTICLE::NoTangentialContact}));

  // type of rolling contact law
  particledyndem.specs.emplace_back(deprecated_selection<RollingContact>("ROLLINGCONTACTLAW",
      {
          {"NoRollingContact", Inpar::PARTICLE::NoRollingContact},
          {"RollingViscous", Inpar::PARTICLE::RollingViscous},
          {"RollingCoulomb", Inpar::PARTICLE::RollingCoulomb},
      },
      {.description = "rolling contact law for particles",
          .default_value = Inpar::PARTICLE::NoRollingContact}));

  // type of normal adhesion law
  particledyndem.specs.emplace_back(deprecated_selection<AdhesionLaw>("ADHESIONLAW",
      {
          {"NoAdhesion", Inpar::PARTICLE::NoAdhesion},
          {"AdhesionVdWDMT", Inpar::PARTICLE::AdhesionVdWDMT},
          {"AdhesionRegDMT", Inpar::PARTICLE::AdhesionRegDMT},
      },
      {.description = "type of adhesion law for particles",
          .default_value = Inpar::PARTICLE::NoAdhesion}));

  // type of (random) surface energy distribution
  particledyndem.specs.emplace_back(
      deprecated_selection<SurfaceEnergyDistribution>("ADHESION_SURFACE_ENERGY_DISTRIBUTION",
          {
              {"ConstantSurfaceEnergy", Inpar::PARTICLE::ConstantSurfaceEnergy},
              {"NormalSurfaceEnergyDistribution", Inpar::PARTICLE::NormalSurfaceEnergyDistribution},
              {"LogNormalSurfaceEnergyDistribution",
                  Inpar::PARTICLE::LogNormalSurfaceEnergyDistribution},
          },
          {.description = "type of (random) surface energy distribution",
              .default_value = Inpar::PARTICLE::ConstantSurfaceEnergy}));

  particledyndem.specs.emplace_back(parameter<double>(
      "MIN_RADIUS", {.description = "minimum allowed particle radius", .default_value = 0.0}));
  particledyndem.specs.emplace_back(parameter<double>(
      "MAX_RADIUS", {.description = "maximum allowed particle radius", .default_value = 0.0}));
  particledyndem.specs.emplace_back(parameter<double>("MAX_VELOCITY",
      {.description = "maximum expected particle velocity", .default_value = -1.0}));

  // type of initial particle radius assignment
  particledyndem.specs.emplace_back(deprecated_selection<InitialRadiusAssignment>("INITIAL_RADIUS",
      {
          {"RadiusFromParticleMaterial", Inpar::PARTICLE::RadiusFromParticleMaterial},
          {"RadiusFromParticleInput", Inpar::PARTICLE::RadiusFromParticleInput},
          {"NormalRadiusDistribution", Inpar::PARTICLE::NormalRadiusDistribution},
          {"LogNormalRadiusDistribution", Inpar::PARTICLE::LogNormalRadiusDistribution},
      },
      {.description = "type of initial particle radius assignment",
          .default_value = Inpar::PARTICLE::RadiusFromParticleMaterial}));

  particledyndem.specs.emplace_back(parameter<double>("RADIUSDISTRIBUTION_SIGMA",
      {.description = "sigma of random particle radius distribution", .default_value = -1.0}));

  particledyndem.specs.emplace_back(parameter<double>("REL_PENETRATION",
      {.description = "maximum allowed relative penetration", .default_value = -1.0}));
  particledyndem.specs.emplace_back(parameter<double>(
      "NORMAL_STIFF", {.description = "normal contact stiffness", .default_value = -1.0}));
  particledyndem.specs.emplace_back(parameter<double>(
      "NORMAL_DAMP", {.description = "normal contact damping parameter", .default_value = -1.0}));
  particledyndem.specs.emplace_back(parameter<double>(
      "COEFF_RESTITUTION", {.description = "coefficient of restitution", .default_value = -1.0}));
  particledyndem.specs.emplace_back(parameter<double>(
      "DAMP_REG_FAC", {.description = "linearly regularized damping normal force in the interval "
                                      "$|g| < (\\text{DAMP_REG_FAC} \\cdot r_{\\min})$",
                          .default_value = -1.0}));
  particledyndem.specs.emplace_back(parameter<bool>("TENSION_CUTOFF",
      {.description = "evaluate tension cutoff of normal contact force", .default_value = true}));

  particledyndem.specs.emplace_back(
      parameter<double>("POISSON_RATIO", {.description = "poisson ratio", .default_value = -1.0}));
  particledyndem.specs.emplace_back(parameter<double>(
      "YOUNG_MODULUS", {.description = "young's modulus", .default_value = -1.0}));

  particledyndem.specs.emplace_back(parameter<double>("FRICT_COEFF_TANG",
      {.description = "friction coefficient for tangential contact", .default_value = -1.0}));
  particledyndem.specs.emplace_back(parameter<double>("FRICT_COEFF_ROLL",
      {.description = "friction coefficient for rolling contact", .default_value = -1.0}));

  particledyndem.specs.emplace_back(parameter<double>("ADHESION_DISTANCE",
      {.description = "adhesion distance between interacting surfaces", .default_value = -1.0}));

  particledyndem.specs.emplace_back(parameter<double>("ADHESION_MAX_CONTACT_PRESSURE",
      {.description = "adhesion maximum contact pressure", .default_value = 0.0}));
  particledyndem.specs.emplace_back(parameter<double>("ADHESION_MAX_CONTACT_FORCE",
      {.description = "adhesion maximum contact force", .default_value = 0.0}));
  particledyndem.specs.emplace_back(parameter<bool>("ADHESION_USE_MAX_CONTACT_FORCE",
      {.description = "use maximum contact force instead of maximum contact pressure",
          .default_value = false}));

  particledyndem.specs.emplace_back(parameter<bool>("ADHESION_VDW_CURVE_SHIFT",
      {.description = "shifts van-der-Waals-curve to g = 0", .default_value = false}));

  particledyndem.specs.emplace_back(parameter<double>("ADHESION_HAMAKER",
      {.description = "hamaker constant of van-der-Waals interaction", .default_value = -1.0}));

  particledyndem.specs.emplace_back(parameter<double>("ADHESION_SURFACE_ENERGY",
      {.description = "adhesion surface energy for the calculation of the pull-out force",
          .default_value = -1.0}));
  particledyndem.specs.emplace_back(parameter<double>("ADHESION_SURFACE_ENERGY_DISTRIBUTION_VAR",
      {.description = "variance of adhesion surface energy distribution", .default_value = -1.0}));
  particledyndem.specs.emplace_back(
      parameter<double>("ADHESION_SURFACE_ENERGY_DISTRIBUTION_CUTOFF_FACTOR",
          {.description = "adhesion surface energy distribution limited by multiple of variance",
              .default_value = -1.0}));
  particledyndem.specs.emplace_back(parameter<double>("ADHESION_SURFACE_ENERGY_FACTOR",
      {.description = "factor to calculate minimum adhesion surface energy",
          .default_value = 1.0}));

  particledyndem.move_into_collection(list);
}

/*---------------------------------------------------------------------------*
 | set the particle conditions                                               |
 *---------------------------------------------------------------------------*/
void Inpar::PARTICLE::set_valid_conditions(
    std::vector<Core::Conditions::ConditionDefinition>& condlist)
{
  using namespace Core::IO::InputSpecBuilders;

  /*-------------------------------------------------------------------------*
   | particle wall condition                                                 |
   *-------------------------------------------------------------------------*/
  Core::Conditions::ConditionDefinition surfpartwall("DESIGN SURFACE PARTICLE WALL", "ParticleWall",
      "Wall for particle interaction with (optional) material definition",
      Core::Conditions::ParticleWall, true, Core::Conditions::geometry_type_surface);

  surfpartwall.add_component(parameter<int>("MAT"));

  condlist.push_back(surfpartwall);
}

FOUR_C_NAMESPACE_CLOSE
