/*----------------------------------------------------------------------*/
/*! \file

\brief Data container holding all beam to solid surface contact input parameters.

\level 3
*/


#include "4C_beaminteraction_beam_to_solid_surface_contact_params.hpp"

#include "4C_beaminteraction_beam_to_solid_surface_visualization_output_params.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_geometry_pair.hpp"

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
BEAMINTERACTION::BeamToSolidSurfaceContactParams::BeamToSolidSurfaceContactParams()
    : BeamToSolidParamsBase(),
      contact_type_(INPAR::BEAMTOSOLID::BeamToSolidSurfaceContact::gap_variation),
      penalty_law_(INPAR::BEAMTOSOLID::BeamToSolidSurfaceContactPenaltyLaw::none),
      penalty_parameter_g0_(0.0),
      output_params_ptr_(Teuchos::null)
{
  // Empty Constructor.
}


/**
 *
 */
void BEAMINTERACTION::BeamToSolidSurfaceContactParams::Init()
{
  // Teuchos parameter list for beam contact
  const Teuchos::ParameterList& beam_to_solid_contact_params_list =
      GLOBAL::Problem::Instance()->BeamInteractionParams().sublist("BEAM TO SOLID SURFACE CONTACT");

  // Set the common beam-to-solid parameters.
  SetBaseParams(beam_to_solid_contact_params_list);

  // Get parameters form input file.
  {
    contact_type_ = Teuchos::getIntegralValue<INPAR::BEAMTOSOLID::BeamToSolidSurfaceContact>(
        beam_to_solid_contact_params_list, "CONTACT_TYPE");

    penalty_law_ =
        Teuchos::getIntegralValue<INPAR::BEAMTOSOLID::BeamToSolidSurfaceContactPenaltyLaw>(
            beam_to_solid_contact_params_list, "PENALTY_LAW");

    penalty_parameter_g0_ = beam_to_solid_contact_params_list.get<double>("PENALTY_PARAMETER_G0");
  }

  // Setup the output parameter object.
  {
    output_params_ptr_ = Teuchos::rcp<BeamToSolidSurfaceVisualizationOutputParams>(
        new BeamToSolidSurfaceVisualizationOutputParams());
    output_params_ptr_->Init();
    output_params_ptr_->Setup();
  }

  isinit_ = true;
}


/**
 *
 */
int BEAMINTERACTION::BeamToSolidSurfaceContactParams::GetFADOrder() const

{
  switch (GetContactType())
  {
    case INPAR::BEAMTOSOLID::BeamToSolidSurfaceContact::gap_variation:
      return 1;
      break;
    case INPAR::BEAMTOSOLID::BeamToSolidSurfaceContact::potential:
      return 2;
      break;
    default:
      FOUR_C_THROW("Got unexpected contact type.");
      return 0;
  }
}

FOUR_C_NAMESPACE_CLOSE