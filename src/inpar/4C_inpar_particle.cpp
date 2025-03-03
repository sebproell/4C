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
  Core::Utils::string_to_integral_parameter<DynamicType>("DYNAMICTYPE", "VelocityVerlet",
      "type of particle time integration",
      tuple<std::string>("SemiImplicitEuler", "VelocityVerlet"),
      tuple<DynamicType>(
          Inpar::PARTICLE::dyna_semiimpliciteuler, Inpar::PARTICLE::dyna_velocityverlet),
      particledyn);

  // type of particle interaction
  Core::Utils::string_to_integral_parameter<InteractionType>("INTERACTION", "None",
      "type of particle interaction", tuple<std::string>("None", "SPH", "DEM"),
      tuple<InteractionType>(Inpar::PARTICLE::interaction_none, Inpar::PARTICLE::interaction_sph,
          Inpar::PARTICLE::interaction_dem),
      particledyn);

  // output type
  Core::Utils::int_parameter(
      "RESULTSEVERY", 1, "write particle runtime output every RESULTSEVERY steps", particledyn);
  Core::Utils::int_parameter(
      "RESTARTEVERY", 1, "write restart possibility every RESTARTEVERY steps", particledyn);

  // write ghosted particles
  particledyn.specs.emplace_back(parameter<bool>("WRITE_GHOSTED_PARTICLES",
      {.description = "write ghosted particles (debug feature)", .default_value = false}));

  // time loop control
  Core::Utils::double_parameter("TIMESTEP", 0.01, "time step size", particledyn);
  Core::Utils::int_parameter("NUMSTEP", 100, "maximum number of steps", particledyn);
  Core::Utils::double_parameter("MAXTIME", 1.0, "maximum time", particledyn);

  // gravity acceleration control
  Core::Utils::string_parameter(
      "GRAVITY_ACCELERATION", "0.0 0.0 0.0", "acceleration due to gravity", particledyn);
  Core::Utils::int_parameter(
      "GRAVITY_RAMP_FUNCT", -1, "number of function governing gravity ramp", particledyn);

  // viscous damping factor
  Core::Utils::double_parameter("VISCOUS_DAMPING", -1.0,
      "apply viscous damping force to determine static equilibrium solutions", particledyn);

  // transfer particles to new bins every time step
  particledyn.specs.emplace_back(parameter<bool>("TRANSFER_EVERY",
      {.description = "transfer particles to new bins every time step", .default_value = false}));

  // considered particle phases with dynamic load balance weighting factor
  Core::Utils::string_parameter("PHASE_TO_DYNLOADBALFAC", "none",
      "considered particle phases with dynamic load balance weighting factor", particledyn);

  // relate particle phase to material id
  Core::Utils::string_parameter(
      "PHASE_TO_MATERIAL_ID", "none", "relate particle phase to material id", particledyn);

  // amplitude of noise added to initial position for each spatial direction
  Core::Utils::string_parameter("INITIAL_POSITION_AMPLITUDE", "0.0 0.0 0.0",
      "amplitude of noise added to initial position for each spatial direction", particledyn);

  // type of particle wall source
  Core::Utils::string_to_integral_parameter<ParticleWallSource>("PARTICLE_WALL_SOURCE",
      "NoParticleWall", "type of particle wall source",
      tuple<std::string>("NoParticleWall", "DiscretCondition", "BoundingBox"),
      tuple<ParticleWallSource>(Inpar::PARTICLE::NoParticleWall, Inpar::PARTICLE::DiscretCondition,
          Inpar::PARTICLE::BoundingBox),
      particledyn);

  // material id for particle wall from bounding box source
  Core::Utils::int_parameter("PARTICLE_WALL_MAT", -1,
      "material id for particle wall from bounding box source", particledyn);

  // flags defining considered states of particle wall
  particledyn.specs.emplace_back(parameter<bool>("PARTICLE_WALL_MOVING",
      {.description = "consider a moving particle wall", .default_value = false}));
  particledyn.specs.emplace_back(parameter<bool>("PARTICLE_WALL_LOADED",
      {.description = "consider loading on particle wall", .default_value = false}));

  // consider rigid body motion
  particledyn.specs.emplace_back(parameter<bool>(
      "RIGID_BODY_MOTION", {.description = "consider rigid body motion", .default_value = false}));

  Core::Utils::double_parameter("RIGID_BODY_PHASECHANGE_RADIUS", -1.0,
      "search radius for neighboring rigid bodies in case of phase change", particledyn);

  particledyn.move_into_collection(list);

  /*-------------------------------------------------------------------------*
   | control parameters for initial/boundary conditions                      |
   *-------------------------------------------------------------------------*/
  Core::Utils::SectionSpecs particledynconditions{particledyn, "INITIAL AND BOUNDARY CONDITIONS"};

  // initial temperature field of particle phase given by function
  Core::Utils::string_parameter("INITIAL_TEMP_FIELD", "none",
      "Refer to the function ID describing the initial temperature field of particle phase",
      particledynconditions);

  // initial velocity field of particle phase given by function
  Core::Utils::string_parameter("INITIAL_VELOCITY_FIELD", "none",
      "Refer to the function ID describing the initial velocity field of particle phase",
      particledynconditions);

  // initial angular velocity field of particle phase given by function
  Core::Utils::string_parameter("INITIAL_ANGULAR_VELOCITY_FIELD", "none",
      "Refer to the function ID describing the initial angular velocity field of rigid body "
      "phase/DEM particle",
      particledynconditions);

  // initial acceleration field of particle phase given by function
  Core::Utils::string_parameter("INITIAL_ACCELERATION_FIELD", "none",
      "Refer to the function ID describing the initial acceleration field of particle phase",
      particledynconditions);

  // initial angular acceleration field of particle phase given by function
  Core::Utils::string_parameter("INITIAL_ANGULAR_ACCELERATION_FIELD", "none",
      "Refer to the function ID describing the initial angular acceleration field of rigid body "
      "phase/DEM particle",
      particledynconditions);

  // dirichlet boundary condition of particle phase given by function
  Core::Utils::string_parameter("DIRICHLET_BOUNDARY_CONDITION", "none",
      "Refer to the function ID describing the dirichlet boundary condition of particle phase",
      particledynconditions);

  // temperature boundary condition of particle phase given by function
  Core::Utils::string_parameter("TEMPERATURE_BOUNDARY_CONDITION", "none",
      "Refer to the function ID describing the temperature boundary condition of particle phase",
      particledynconditions);

  particledynconditions.move_into_collection(list);

  /*-------------------------------------------------------------------------*
   | smoothed particle hydrodynamics (SPH) specific control parameters       |
   *-------------------------------------------------------------------------*/
  Core::Utils::SectionSpecs particledynsph{particledyn, "SPH"};

  // write particle-wall interaction output
  particledynsph.specs.emplace_back(parameter<bool>("WRITE_PARTICLE_WALL_INTERACTION",
      {.description = "write particle-wall interaction output", .default_value = false}));

  // type of smoothed particle hydrodynamics kernel
  Core::Utils::string_to_integral_parameter<KernelType>("KERNEL", "CubicSpline",
      "type of smoothed particle hydrodynamics kernel",
      tuple<std::string>("CubicSpline", "QuinticSpline"),
      tuple<KernelType>(Inpar::PARTICLE::CubicSpline, Inpar::PARTICLE::QuinticSpline),
      particledynsph);

  // kernel space dimension number
  Core::Utils::string_to_integral_parameter<KernelSpaceDimension>("KERNEL_SPACE_DIM", "Kernel3D",
      "kernel space dimension number", tuple<std::string>("Kernel1D", "Kernel2D", "Kernel3D"),
      tuple<KernelSpaceDimension>(
          Inpar::PARTICLE::Kernel1D, Inpar::PARTICLE::Kernel2D, Inpar::PARTICLE::Kernel3D),
      particledynsph);

  Core::Utils::double_parameter(
      "INITIALPARTICLESPACING", 0.0, "initial spacing of particles", particledynsph);

  // type of smoothed particle hydrodynamics equation of state
  Core::Utils::string_to_integral_parameter<EquationOfStateType>("EQUATIONOFSTATE", "GenTait",
      "type of smoothed particle hydrodynamics equation of state",
      tuple<std::string>("GenTait", "IdealGas"),
      tuple<EquationOfStateType>(Inpar::PARTICLE::GenTait, Inpar::PARTICLE::IdealGas),
      particledynsph);

  // type of smoothed particle hydrodynamics momentum formulation
  Core::Utils::string_to_integral_parameter<MomentumFormulationType>("MOMENTUMFORMULATION",
      "AdamiMomentumFormulation", "type of smoothed particle hydrodynamics momentum formulation",
      tuple<std::string>("AdamiMomentumFormulation", "MonaghanMomentumFormulation"),
      tuple<MomentumFormulationType>(
          Inpar::PARTICLE::AdamiMomentumFormulation, Inpar::PARTICLE::MonaghanMomentumFormulation),
      particledynsph);

  // type of density evaluation scheme
  Core::Utils::string_to_integral_parameter<DensityEvaluationScheme>("DENSITYEVALUATION",
      "DensitySummation", "type of density evaluation scheme",
      tuple<std::string>("DensitySummation", "DensityIntegration", "DensityPredictCorrect"),
      tuple<DensityEvaluationScheme>(Inpar::PARTICLE::DensitySummation,
          Inpar::PARTICLE::DensityIntegration, Inpar::PARTICLE::DensityPredictCorrect),
      particledynsph);

  // type of density correction scheme
  Core::Utils::string_to_integral_parameter<DensityCorrectionScheme>("DENSITYCORRECTION",
      "NoCorrection", "type of density correction scheme",
      tuple<std::string>(
          "NoCorrection", "InteriorCorrection", "NormalizedCorrection", "RandlesCorrection"),
      tuple<DensityCorrectionScheme>(Inpar::PARTICLE::NoCorrection,
          Inpar::PARTICLE::InteriorCorrection, Inpar::PARTICLE::NormalizedCorrection,
          Inpar::PARTICLE::RandlesCorrection),
      particledynsph);

  // type of boundary particle formulation
  Core::Utils::string_to_integral_parameter<BoundaryParticleFormulationType>(
      "BOUNDARYPARTICLEFORMULATION", "NoBoundaryFormulation",
      "type of boundary particle formulation",
      tuple<std::string>("NoBoundaryFormulation", "AdamiBoundaryFormulation"),
      tuple<BoundaryParticleFormulationType>(
          Inpar::PARTICLE::NoBoundaryFormulation, Inpar::PARTICLE::AdamiBoundaryFormulation),
      particledynsph);

  // type of boundary particle interaction
  Core::Utils::string_to_integral_parameter<BoundaryParticleInteraction>(
      "BOUNDARYPARTICLEINTERACTION", "NoSlipBoundaryParticle",
      "type of boundary particle interaction",
      tuple<std::string>("NoSlipBoundaryParticle", "FreeSlipBoundaryParticle"),
      tuple<BoundaryParticleInteraction>(
          Inpar::PARTICLE::NoSlipBoundaryParticle, Inpar::PARTICLE::FreeSlipBoundaryParticle),
      particledynsph);

  // type of wall formulation
  Core::Utils::string_to_integral_parameter<WallFormulationType>("WALLFORMULATION",
      "NoWallFormulation", "type of wall formulation",
      tuple<std::string>("NoWallFormulation", "VirtualParticleWallFormulation"),
      tuple<WallFormulationType>(
          Inpar::PARTICLE::NoWallFormulation, Inpar::PARTICLE::VirtualParticleWallFormulation),
      particledynsph);

  // type of transport velocity formulation
  Core::Utils::string_to_integral_parameter<TransportVelocityFormulation>(
      "TRANSPORTVELOCITYFORMULATION", "NoTransportVelocity",
      "type of transport velocity formulation",
      tuple<std::string>(
          "NoTransportVelocity", "StandardTransportVelocity", "GeneralizedTransportVelocity"),
      tuple<TransportVelocityFormulation>(Inpar::PARTICLE::NoTransportVelocity,
          Inpar::PARTICLE::StandardTransportVelocity,
          Inpar::PARTICLE::GeneralizedTransportVelocity),
      particledynsph);

  // type of temperature evaluation scheme
  Core::Utils::string_to_integral_parameter<TemperatureEvaluationScheme>("TEMPERATUREEVALUATION",
      "NoTemperatureEvaluation", "type of temperature evaluation scheme",
      tuple<std::string>("NoTemperatureEvaluation", "TemperatureIntegration"),
      tuple<TemperatureEvaluationScheme>(
          Inpar::PARTICLE::NoTemperatureEvaluation, Inpar::PARTICLE::TemperatureIntegration),
      particledynsph);

  particledynsph.specs.emplace_back(parameter<bool>("TEMPERATUREGRADIENT",
      {.description = "evaluate temperature gradient", .default_value = false}));

  // type of heat source
  Core::Utils::string_to_integral_parameter<HeatSourceType>("HEATSOURCETYPE", "NoHeatSource",
      "type of heat source",
      tuple<std::string>("NoHeatSource", "VolumeHeatSource", "SurfaceHeatSource"),
      tuple<HeatSourceType>(Inpar::PARTICLE::NoHeatSource, Inpar::PARTICLE::VolumeHeatSource,
          Inpar::PARTICLE::SurfaceHeatSource),
      particledynsph);

  Core::Utils::int_parameter(
      "HEATSOURCE_FUNCT", -1, "number of function governing heat source", particledynsph);

  Core::Utils::string_parameter(
      "HEATSOURCE_DIRECTION", "0.0 0.0 0.0", "direction of surface heat source", particledynsph);

  // evaporation induced heat loss
  particledynsph.specs.emplace_back(parameter<bool>("VAPOR_HEATLOSS",
      {.description = "evaluate evaporation induced heat loss", .default_value = false}));
  Core::Utils::double_parameter(
      "VAPOR_HEATLOSS_LATENTHEAT", 0.0, "latent heat in heat loss formula", particledynsph);
  Core::Utils::double_parameter("VAPOR_HEATLOSS_ENTHALPY_REFTEMP", 0.0,
      "enthalpy reference temperature in heat loss formula", particledynsph);
  Core::Utils::double_parameter(
      "VAPOR_HEATLOSS_PFAC", 0.0, "pressure factor in heat loss formula", particledynsph);
  Core::Utils::double_parameter(
      "VAPOR_HEATLOSS_TFAC", 0.0, "temperature factor in heat loss formula", particledynsph);

  // evaporation induced recoil pressure
  particledynsph.specs.emplace_back(parameter<bool>("VAPOR_RECOIL",
      {.description = "evaluate evaporation induced recoil pressure", .default_value = false}));
  Core::Utils::double_parameter("VAPOR_RECOIL_BOILINGTEMPERATURE", 0.0,
      "boiling temperature in recoil pressure formula", particledynsph);
  Core::Utils::double_parameter(
      "VAPOR_RECOIL_PFAC", 0.0, "pressure factor in recoil pressure formula", particledynsph);
  Core::Utils::double_parameter(
      "VAPOR_RECOIL_TFAC", 0.0, "temperature factor in recoil pressure formula", particledynsph);

  // type of surface tension formulation
  Core::Utils::string_to_integral_parameter<SurfaceTensionFormulation>("SURFACETENSIONFORMULATION",
      "NoSurfaceTension", "type of surface tension formulation",
      tuple<std::string>("NoSurfaceTension", "ContinuumSurfaceForce"),
      tuple<SurfaceTensionFormulation>(
          Inpar::PARTICLE::NoSurfaceTension, Inpar::PARTICLE::ContinuumSurfaceForce),
      particledynsph);

  Core::Utils::int_parameter("SURFACETENSION_RAMP_FUNCT", -1,
      "number of function governing surface tension ramp", particledynsph);

  Core::Utils::double_parameter("SURFACETENSIONCOEFFICIENT", -1.0,
      "constant part of surface tension coefficient", particledynsph);
  Core::Utils::double_parameter("SURFACETENSIONMINIMUM", 0.0,
      "minimum surface tension coefficient in case of temperature dependence", particledynsph);
  Core::Utils::double_parameter("SURFACETENSIONTEMPFAC", 0.0,
      "factor of dependence of surface tension coefficient on temperature", particledynsph);
  Core::Utils::double_parameter("SURFACETENSIONREFTEMP", 0.0,
      "reference temperature for surface tension coefficient", particledynsph);

  // wetting
  Core::Utils::double_parameter("STATICCONTACTANGLE", 0.0,
      "static contact angle in degree with wetting effects", particledynsph);
  Core::Utils::double_parameter("TRIPLEPOINTNORMAL_CORR_CF_LOW", 0.0,
      "triple point normal correction wall color field low", particledynsph);
  Core::Utils::double_parameter("TRIPLEPOINTNORMAL_CORR_CF_UP", 0.0,
      "triple point normal correction wall color field up", particledynsph);

  // interface viscosity
  particledynsph.specs.emplace_back(parameter<bool>("INTERFACE_VISCOSITY",
      {.description = "evaluate interface viscosity", .default_value = false}));
  Core::Utils::double_parameter("INTERFACE_VISCOSITY_LIQUIDGAS", 0.0,
      "artificial viscosity on liquid-gas interface", particledynsph);
  Core::Utils::double_parameter("INTERFACE_VISCOSITY_SOLIDLIQUID", 0.0,
      "artificial viscosity on solid-liquid interface", particledynsph);

  // barrier force
  particledynsph.specs.emplace_back(parameter<bool>(
      "BARRIER_FORCE", {.description = "evaluate barrier force", .default_value = false}));
  Core::Utils::double_parameter(
      "BARRIER_FORCE_DISTANCE", 0.0, "barrier force distance", particledynsph);
  Core::Utils::double_parameter(
      "BARRIER_FORCE_TEMPSCALE", 0.0, "barrier force temperature scaling", particledynsph);
  Core::Utils::double_parameter(
      "BARRIER_FORCE_STIFF_HEAVY", -1.0, "barrier force stiffness of heavy phase", particledynsph);
  Core::Utils::double_parameter("BARRIER_FORCE_DAMP_HEAVY", 0.0,
      "barrier force damping parameter of heavy phase", particledynsph);
  Core::Utils::double_parameter(
      "BARRIER_FORCE_STIFF_GAS", -1.0, "barrier force stiffness of gas phase", particledynsph);
  Core::Utils::double_parameter("BARRIER_FORCE_DAMP_GAS", 0.0,
      "barrier force damping parameter of gas phase", particledynsph);

  // linear transition in surface tension evaluation
  Core::Utils::double_parameter(
      "TRANS_REF_TEMPERATURE", 0.0, "transition reference temperature", particledynsph);
  Core::Utils::double_parameter("TRANS_DT_SURFACETENSION", 0.0,
      "transition temperature difference for surface tension evaluation", particledynsph);
  Core::Utils::double_parameter("TRANS_DT_MARANGONI", 0.0,
      "transition temperature difference for marangoni evaluation", particledynsph);
  Core::Utils::double_parameter("TRANS_DT_CURVATURE", 0.0,
      "transition temperature difference for curvature evaluation", particledynsph);
  Core::Utils::double_parameter("TRANS_DT_WETTING", 0.0,
      "transition temperature difference for wetting evaluation", particledynsph);
  Core::Utils::double_parameter("TRANS_DT_INTVISC", 0.0,
      "transition temperature difference for interface viscosity evaluation", particledynsph);
  Core::Utils::double_parameter("TRANS_DT_BARRIER", 0.0,
      "transition temperature difference for barrier force evaluation", particledynsph);

  // type of dirichlet open boundary
  Core::Utils::string_to_integral_parameter<DirichletOpenBoundaryType>("DIRICHLETBOUNDARYTYPE",
      "NoDirichletOpenBoundary", "type of dirichlet open boundary",
      tuple<std::string>("NoDirichletOpenBoundary", "DirichletNormalToPlane"),
      tuple<DirichletOpenBoundaryType>(
          Inpar::PARTICLE::NoDirichletOpenBoundary, Inpar::PARTICLE::DirichletNormalToPlane),
      particledynsph);

  Core::Utils::int_parameter("DIRICHLET_FUNCT", -1,
      "number of function governing velocity condition on dirichlet open boundary", particledynsph);

  Core::Utils::string_parameter("DIRICHLET_OUTWARD_NORMAL", "0.0 0.0 0.0",
      "direction of outward normal on dirichlet open boundary", particledynsph);
  Core::Utils::string_parameter("DIRICHLET_PLANE_POINT", "0.0 0.0 0.0",
      "point on dirichlet open boundary plane", particledynsph);

  // type of neumann open boundary
  Core::Utils::string_to_integral_parameter<NeumannOpenBoundaryType>("NEUMANNBOUNDARYTYPE",
      "NoNeumannOpenBoundary", "type of neumann open boundary",
      tuple<std::string>("NoNeumannOpenBoundary", "NeumannNormalToPlane"),
      tuple<NeumannOpenBoundaryType>(
          Inpar::PARTICLE::NoNeumannOpenBoundary, Inpar::PARTICLE::NeumannNormalToPlane),
      particledynsph);

  Core::Utils::int_parameter("NEUMANN_FUNCT", -1,
      "number of function governing pressure condition on neumann open boundary", particledynsph);

  Core::Utils::string_parameter("NEUMANN_OUTWARD_NORMAL", "0.0 0.0 0.0",
      "direction of outward normal on neumann open boundary", particledynsph);
  Core::Utils::string_parameter(
      "NEUMANN_PLANE_POINT", "0.0 0.0 0.0", "point on neumann open boundary plane", particledynsph);

  // type of phase change
  Core::Utils::string_to_integral_parameter<PhaseChangeType>("PHASECHANGETYPE", "NoPhaseChange",
      "type of phase change",
      tuple<std::string>("NoPhaseChange", "OneWayScalarBelowToAbovePhaseChange",
          "OneWayScalarAboveToBelowPhaseChange", "TwoWayScalarPhaseChange"),
      tuple<PhaseChangeType>(Inpar::PARTICLE::NoPhaseChange,
          Inpar::PARTICLE::OneWayScalarBelowToAbovePhaseChange,
          Inpar::PARTICLE::OneWayScalarAboveToBelowPhaseChange,
          Inpar::PARTICLE::TwoWayScalarPhaseChange),
      particledynsph);

  // definition of phase change
  Core::Utils::string_parameter(
      "PHASECHANGEDEFINITION", "none", "phase change definition", particledynsph);

  // type of rigid particle contact
  Core::Utils::string_to_integral_parameter<RigidParticleContactType>("RIGIDPARTICLECONTACTTYPE",
      "NoRigidParticleContact", "type of rigid particle contact",
      tuple<std::string>("NoRigidParticleContact", "ElasticRigidParticleContact"),
      tuple<RigidParticleContactType>(
          Inpar::PARTICLE::NoRigidParticleContact, Inpar::PARTICLE::ElasticRigidParticleContact),
      particledynsph);

  Core::Utils::double_parameter(
      "RIGIDPARTICLECONTACTSTIFF", -1.0, "rigid particle contact stiffness", particledynsph);
  Core::Utils::double_parameter(
      "RIGIDPARTICLECONTACTDAMP", 0.0, "rigid particle contact damping parameter", particledynsph);

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
  Core::Utils::string_to_integral_parameter<NormalContact>("NORMALCONTACTLAW", "NormalLinearSpring",
      "normal contact law for particles",
      tuple<std::string>("NormalLinearSpring", "NormalLinearSpringDamp", "NormalHertz",
          "NormalLeeHerrmann", "NormalKuwabaraKono", "NormalTsuji"),
      tuple<NormalContact>(Inpar::PARTICLE::NormalLinSpring, Inpar::PARTICLE::NormalLinSpringDamp,
          Inpar::PARTICLE::NormalHertz, Inpar::PARTICLE::NormalLeeHerrmann,
          Inpar::PARTICLE::NormalKuwabaraKono, Inpar::PARTICLE::NormalTsuji),
      particledyndem);

  // type of tangential contact law
  Core::Utils::string_to_integral_parameter<TangentialContact>("TANGENTIALCONTACTLAW",
      "NoTangentialContact", "tangential contact law for particles",
      tuple<std::string>("NoTangentialContact", "TangentialLinSpringDamp"),
      tuple<TangentialContact>(
          Inpar::PARTICLE::NoTangentialContact, Inpar::PARTICLE::TangentialLinSpringDamp),
      particledyndem);

  // type of rolling contact law
  Core::Utils::string_to_integral_parameter<RollingContact>("ROLLINGCONTACTLAW", "NoRollingContact",
      "rolling contact law for particles",
      tuple<std::string>("NoRollingContact", "RollingViscous", "RollingCoulomb"),
      tuple<RollingContact>(Inpar::PARTICLE::NoRollingContact, Inpar::PARTICLE::RollingViscous,
          Inpar::PARTICLE::RollingCoulomb),
      particledyndem);

  // type of normal adhesion law
  Core::Utils::string_to_integral_parameter<AdhesionLaw>("ADHESIONLAW", "NoAdhesion",
      "type of adhesion law for particles",
      tuple<std::string>("NoAdhesion", "AdhesionVdWDMT", "AdhesionRegDMT"),
      tuple<AdhesionLaw>(Inpar::PARTICLE::NoAdhesion, Inpar::PARTICLE::AdhesionVdWDMT,
          Inpar::PARTICLE::AdhesionRegDMT),
      particledyndem);

  // type of (random) surface energy distribution
  Core::Utils::string_to_integral_parameter<SurfaceEnergyDistribution>(
      "ADHESION_SURFACE_ENERGY_DISTRIBUTION", "ConstantSurfaceEnergy",
      "type of (random) surface energy distribution",
      tuple<std::string>("ConstantSurfaceEnergy", "NormalSurfaceEnergyDistribution",
          "LogNormalSurfaceEnergyDistribution"),
      tuple<SurfaceEnergyDistribution>(Inpar::PARTICLE::ConstantSurfaceEnergy,
          Inpar::PARTICLE::NormalSurfaceEnergyDistribution,
          Inpar::PARTICLE::LogNormalSurfaceEnergyDistribution),
      particledyndem);

  Core::Utils::double_parameter(
      "MIN_RADIUS", 0.0, "minimum allowed particle radius", particledyndem);
  Core::Utils::double_parameter(
      "MAX_RADIUS", 0.0, "maximum allowed particle radius", particledyndem);
  Core::Utils::double_parameter(
      "MAX_VELOCITY", -1.0, "maximum expected particle velocity", particledyndem);

  // type of initial particle radius assignment
  Core::Utils::string_to_integral_parameter<InitialRadiusAssignment>("INITIAL_RADIUS",
      "RadiusFromParticleMaterial", "type of initial particle radius assignment",
      tuple<std::string>("RadiusFromParticleMaterial", "RadiusFromParticleInput",
          "NormalRadiusDistribution", "LogNormalRadiusDistribution"),
      tuple<InitialRadiusAssignment>(Inpar::PARTICLE::RadiusFromParticleMaterial,
          Inpar::PARTICLE::RadiusFromParticleInput, Inpar::PARTICLE::NormalRadiusDistribution,
          Inpar::PARTICLE::LogNormalRadiusDistribution),
      particledyndem);

  Core::Utils::double_parameter("RADIUSDISTRIBUTION_SIGMA", -1.0,
      "sigma of random particle radius distribution", particledyndem);

  Core::Utils::double_parameter(
      "REL_PENETRATION", -1.0, "maximum allowed relative penetration", particledyndem);
  Core::Utils::double_parameter("NORMAL_STIFF", -1.0, "normal contact stiffness", particledyndem);
  Core::Utils::double_parameter(
      "NORMAL_DAMP", -1.0, "normal contact damping parameter", particledyndem);
  Core::Utils::double_parameter(
      "COEFF_RESTITUTION", -1.0, "coefficient of restitution", particledyndem);
  Core::Utils::double_parameter("DAMP_REG_FAC", -1.0,
      "linearly regularized damping normal force in the interval "
      "$|g| < (\\text{DAMP_REG_FAC} \\cdot r_{\\min})$",
      particledyndem);
  particledyndem.specs.emplace_back(parameter<bool>("TENSION_CUTOFF",
      {.description = "evaluate tension cutoff of normal contact force", .default_value = true}));

  Core::Utils::double_parameter("POISSON_RATIO", -1.0, "poisson ratio", particledyndem);
  Core::Utils::double_parameter("YOUNG_MODULUS", -1.0, "young's modulus", particledyndem);

  Core::Utils::double_parameter(
      "FRICT_COEFF_TANG", -1.0, "friction coefficient for tangential contact", particledyndem);
  Core::Utils::double_parameter(
      "FRICT_COEFF_ROLL", -1.0, "friction coefficient for rolling contact", particledyndem);

  Core::Utils::double_parameter(
      "ADHESION_DISTANCE", -1.0, "adhesion distance between interacting surfaces", particledyndem);

  Core::Utils::double_parameter(
      "ADHESION_MAX_CONTACT_PRESSURE", 0.0, "adhesion maximum contact pressure", particledyndem);
  Core::Utils::double_parameter(
      "ADHESION_MAX_CONTACT_FORCE", 0.0, "adhesion maximum contact force", particledyndem);
  particledyndem.specs.emplace_back(parameter<bool>("ADHESION_USE_MAX_CONTACT_FORCE",
      {.description = "use maximum contact force instead of maximum contact pressure",
          .default_value = false}));

  particledyndem.specs.emplace_back(parameter<bool>("ADHESION_VDW_CURVE_SHIFT",
      {.description = "shifts van-der-Waals-curve to g = 0", .default_value = false}));

  Core::Utils::double_parameter(
      "ADHESION_HAMAKER", -1.0, "hamaker constant of van-der-Waals interaction", particledyndem);

  Core::Utils::double_parameter("ADHESION_SURFACE_ENERGY", -1.0,
      "adhesion surface energy for the calculation of the pull-out force", particledyndem);
  Core::Utils::double_parameter("ADHESION_SURFACE_ENERGY_DISTRIBUTION_VAR", -1.0,
      "variance of adhesion surface energy distribution", particledyndem);
  Core::Utils::double_parameter("ADHESION_SURFACE_ENERGY_DISTRIBUTION_CUTOFF_FACTOR", -1.0,
      "adhesion surface energy distribution limited by multiple of variance", particledyndem);
  Core::Utils::double_parameter("ADHESION_SURFACE_ENERGY_FACTOR", 1.0,
      "factor to calculate minimum adhesion surface energy", particledyndem);

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
