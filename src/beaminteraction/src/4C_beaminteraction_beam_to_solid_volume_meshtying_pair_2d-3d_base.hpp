// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_BEAMINTERACTION_BEAM_TO_SOLID_VOLUME_MESHTYING_PAIR_2D_3D_BASE_HPP
#define FOUR_C_BEAMINTERACTION_BEAM_TO_SOLID_VOLUME_MESHTYING_PAIR_2D_3D_BASE_HPP


#include "4C_config.hpp"

#include "4C_beam3_triad_interpolation_local_rotation_vectors.hpp"
#include "4C_beaminteraction_beam_to_solid_volume_meshtying_pair_base.hpp"

FOUR_C_NAMESPACE_OPEN

// Forward declaration.
namespace GeometryPair
{
  template <typename ScalarType, typename Line, typename Volume>
  class GeometryPairLineToVolumeGaussPointProjectionCrossSection;
}  // namespace GeometryPair


namespace BeamInteraction
{
  /**
   * \brief Base class for 2D-3D beam to solid volume mesh tying
   * @tparam ScalarType Scalar FAD type to be used in this pair.
   * @tparam beam Type from GeometryPair::ElementDiscretization... representing the beam.
   * @tparam solid Type from GeometryPair::ElementDiscretization... representing the solid.
   */
  template <typename ScalarType, typename Beam, typename Solid>
  class BeamToSolidVolumeMeshtyingPair2D3DBase
      : public BeamToSolidVolumeMeshtyingPairBase<ScalarType, Beam, Solid>
  {
   protected:
    //! Shortcut to the base class.
    using base_class = BeamToSolidVolumeMeshtyingPairBase<ScalarType, Beam, Solid>;

    //! Type to be used for scalar AD variables. This can not be inherited from the base class.
    using scalar_type = typename base_class::scalar_type;

   public:
    /**
     * \brief Standard Constructor
     */
    BeamToSolidVolumeMeshtyingPair2D3DBase() = default;


    /**
     * \brief Create the geometry pair for this contact pair. We overload that function because this
     * pair requires explicitly that a cross section projection pair is created.
     * @param element1 Pointer to the first element
     * @param element2 Pointer to the second element
     * @param geometry_evaluation_data_ptr Evaluation data that will be linked to the pair.
     */
    void create_geometry_pair(const Core::Elements::Element* element1,
        const Core::Elements::Element* element2,
        const std::shared_ptr<GeometryPair::GeometryEvaluationDataBase>&
            geometry_evaluation_data_ptr) override;

   protected:
    /**
     * \brief Calculate the position on the beam, also taking into account parameter coordinates on
     * the cross section.
     * @param integration_point (in) Integration where the position should be evaluated.
     * @param r_beam (out) Position on the beam.
     * @param reference (in) True -> the reference position is calculated, False -> the current
     * position is calculated.
     */
    void evaluate_beam_position_double(
        const GeometryPair::ProjectionPoint1DTo3D<double>& integration_point,
        Core::LinAlg::Matrix<3, 1, double>& r_beam, bool reference) const override;

    /**
     * \brief Return a cast of the geometry pair to the type for this contact pair.
     * @return RPC with the type of geometry pair for this beam contact pair.
     */
    inline std::shared_ptr<
        GeometryPair::GeometryPairLineToVolumeGaussPointProjectionCrossSection<double, Beam, Solid>>
    cast_geometry_pair() const
    {
      return std::dynamic_pointer_cast<GeometryPair::
              GeometryPairLineToVolumeGaussPointProjectionCrossSection<double, Beam, Solid>>(
          this->geometry_pair_);
    };

    /**
     * \brief Get the triad of the beam at the parameter coordinate xi
     * @param xi (in) Parameter coordinate on the beam
     * @param triad (out) Beam cross section triad
     * @param reference (in) If the triad in the reference or current configuration should be
     * returned
     */
    virtual void get_triad_at_xi_double(
        const double xi, Core::LinAlg::Matrix<3, 3, double>& triad, const bool reference) const = 0;
  };
}  // namespace BeamInteraction

FOUR_C_NAMESPACE_CLOSE

#endif
