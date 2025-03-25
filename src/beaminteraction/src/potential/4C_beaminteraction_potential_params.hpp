// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_BEAMINTERACTION_POTENTIAL_PARAMS_HPP
#define FOUR_C_BEAMINTERACTION_POTENTIAL_PARAMS_HPP

#include "4C_config.hpp"

#include "4C_beaminteraction_potential_input.hpp"
#include "4C_beaminteraction_potential_runtime_visualization_output_params.hpp"

#include <unordered_map>

FOUR_C_NAMESPACE_OPEN

namespace BeamInteraction
{

  struct BeamPotentialParams
  {
    //! potential law parameters
    std::vector<double> potential_law_exponents;
    std::vector<double> potential_law_prefactors;

    //! potential strategy
    BeamPotential::Strategy strategy = BeamPotential::Strategy::vague;

    //! potential type
    BeamPotential::Type potential_type = BeamPotential::Type::vague;

    //! cutoff radius
    std::optional<double> cutoff_radius;

    //! regularization
    BeamPotential::RegularizationType regularization_type = BeamPotential::RegularizationType::none;
    double regularization_separation;

    //! potential reduction length
    std::optional<double> potential_reduction_length = std::nullopt;

    //! integration parameters
    int n_integration_segments;
    int n_gauss_points;

    //! automatic differentiation
    bool use_fad;

    //! master/slave choice
    BeamPotential::MasterSlaveChoice choice_master_slave;

    //! visualization output
    bool write_visualization_output;

    BeamInteraction::BeamToBeamPotentialRuntimeOutputParams
        params_runtime_visualization_output_btb_potential =
            BeamToBeamPotentialRuntimeOutputParams(0.0);

    //! data container for prior element lengths for potential reduction strategy
    //! first entry is left prior length and second entry is right prior length
    // this is stored in the beam potential params for easy access during evaluation
    std::unordered_map<int, std::pair<double, double>> ele_gid_prior_length_map;
  };

  void initialize_validate_beam_potential_params(
      BeamPotentialParams& beam_potential_params, const double restart_time);

}  // namespace BeamInteraction

FOUR_C_NAMESPACE_CLOSE

#endif
