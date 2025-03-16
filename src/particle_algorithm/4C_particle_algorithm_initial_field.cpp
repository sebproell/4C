// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_particle_algorithm_initial_field.hpp"

#include "4C_global_data.hpp"
#include "4C_particle_algorithm_utils.hpp"
#include "4C_particle_engine_container.hpp"
#include "4C_particle_engine_container_bundle.hpp"
#include "4C_particle_engine_enums.hpp"
#include "4C_particle_engine_interface.hpp"
#include "4C_utils_function.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEALGORITHM::InitialFieldHandler::InitialFieldHandler(const Teuchos::ParameterList& params)
    : params_(params)
{
  // empty constructor
}

void PARTICLEALGORITHM::InitialFieldHandler::init()
{
  // get control parameters for initial/boundary conditions
  const Teuchos::ParameterList& params_conditions =
      params_.sublist("INITIAL AND BOUNDARY CONDITIONS");

  // relate particle state to input name
  std::map<std::string, PARTICLEENGINE::StateEnum> initialfieldtostateenum = {
      std::make_pair("INITIAL_TEMP_FIELD", PARTICLEENGINE::Temperature),
      std::make_pair("INITIAL_VELOCITY_FIELD", PARTICLEENGINE::Velocity),
      std::make_pair("INITIAL_ANGULAR_VELOCITY_FIELD", PARTICLEENGINE::AngularVelocity),
      std::make_pair("INITIAL_ACCELERATION_FIELD", PARTICLEENGINE::Acceleration),
      std::make_pair("INITIAL_ANGULAR_ACCELERATION_FIELD", PARTICLEENGINE::AngularAcceleration)};

  // iterate over particle states
  for (auto& stateIt : initialfieldtostateenum)
  {
    // get reference to sub-map
    std::map<PARTICLEENGINE::TypeEnum, int>& currentstatetypetofunctidmap =
        statetotypetofunctidmap_[stateIt.second];

    // read parameters relating particle types to values
    PARTICLEALGORITHM::Utils::read_params_types_related_to_values(
        params_conditions, stateIt.first, currentstatetypetofunctidmap);
  }
}

void PARTICLEALGORITHM::InitialFieldHandler::setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;
}

void PARTICLEALGORITHM::InitialFieldHandler::set_initial_fields()
{
  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->get_particle_container_bundle();

  for (auto& stateIt : statetotypetofunctidmap_)
  {
    // get state of particles
    PARTICLEENGINE::StateEnum particleState = stateIt.first;

    // iterate over particle types
    for (auto& initialFieldIt : stateIt.second)
    {
      // get type of particles
      PARTICLEENGINE::TypeEnum particleType = initialFieldIt.first;

      // get container of owned particles of current particle type
      PARTICLEENGINE::ParticleContainer* container =
          particlecontainerbundle->get_specific_container(particleType, PARTICLEENGINE::Owned);

      // get number of particles stored in container
      const int particlestored = container->particles_stored();

      // no owned particles of current particle type
      if (particlestored <= 0) continue;

      if (not container->have_stored_state(particleState)) continue;

      // get id of function
      const int functid = initialFieldIt.second;

      // get reference to function
      const auto& function =
          Global::Problem::instance()->function_by_id<Core::Utils::FunctionOfSpaceTime>(functid);

      // get pointer to particle states
      const double* pos = container->get_ptr_to_state(PARTICLEENGINE::Position, 0);
      double* state = container->get_ptr_to_state(particleState, 0);

      // get particle state dimensions
      int posstatedim = container->get_state_dim(PARTICLEENGINE::Position);
      int statedim = container->get_state_dim(particleState);

      // safety check
      if (static_cast<std::size_t>(statedim) != function.number_components())
        FOUR_C_THROW(
            "dimensions of function defining initial field and of state '{}' not matching!",
            PARTICLEENGINE::enum_to_state_name(particleState).c_str());

      // iterate over owned particles of current type
      for (int i = 0; i < particlestored; ++i)
      {
        // evaluate function to set initial field
        for (int dim = 0; dim < statedim; ++dim)
          state[statedim * i + dim] = function.evaluate(&(pos[posstatedim * i]), 0.0, dim);
      }
    }
  }
}

FOUR_C_NAMESPACE_CLOSE
