// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_GEOMETRY_PAIR_LINE_TO_VOLUME_GAUSS_POINT_PROJECTION_HPP
#define FOUR_C_GEOMETRY_PAIR_LINE_TO_VOLUME_GAUSS_POINT_PROJECTION_HPP


#include "4C_config.hpp"

#include "4C_geometry_pair_line_to_volume.hpp"

FOUR_C_NAMESPACE_OPEN

namespace GeometryPair
{
  /**
   * \brief Class that handles the geometrical interactions of a line and a volume by projecting
   * Gauss points from the line to the volume. In case a line pokes out of the volumes it is
   * segmented.
   * @param scalar_type Type that will be used for scalar values.
   * @param line Type of line element.
   * @param volume Type of volume element.
   */
  template <typename ScalarType, typename Line, typename Volume>
  class GeometryPairLineToVolumeGaussPointProjection
      : public GeometryPairLineToVolume<ScalarType, Line, Volume>
  {
   public:
    //! Public alias for scalar type so that other classes can use this type.
    using t_scalar_type = ScalarType;

    //! Public alias for line type so that other classes can use this type.
    using t_line = Line;

    //! Public alias for volume type so that other classes can use this type.
    using t_other = Volume;

   public:
    /**
     * \brief Constructor.
     */
    GeometryPairLineToVolumeGaussPointProjection(const Core::Elements::Element* element1,
        const Core::Elements::Element* element2,
        const std::shared_ptr<GeometryPair::LineTo3DEvaluationData>& evaluation_data);


    /**
     * \brief Try to project all Gauss points to the volume.
     *
     * Only points are checked that do not
     * already have a valid projection in the projection tracker of the evaluation data container.
     * Eventually needed segmentation at lines poking out of the volume is done in the Evaluate
     * method.
     *
     * @param element_data_line (in) Degrees of freedom for the line.
     * @param element_data_volume (in) Degrees of freedom for the volume.
     * @param segments (out) Vector with the segments of this line to volume pair.
     */
    void pre_evaluate(const ElementData<Line, ScalarType>& element_data_line,
        const ElementData<Volume, ScalarType>& element_data_volume,
        std::vector<LineSegment<ScalarType>>& segments) const override;

    /**
     * \brief Check if a Gauss point projected valid for this pair in pre_evaluate.
     *
     * If so, all Gauss points have to project valid (in the tracker, since some can be valid on
     * other pairs). If not all project, the beam pokes out of the volumes and in this method
     * segmentation is performed.
     *
     * @param element_data_line (in) Degrees of freedom for the line.
     * @param element_data_volume (in) Degrees of freedom for the volume.
     * @param segments (out) Vector with the segments of this line to volume pair.
     */
    void evaluate(const ElementData<Line, ScalarType>& element_data_line,
        const ElementData<Volume, ScalarType>& element_data_volume,
        std::vector<LineSegment<ScalarType>>& segments) const override;
  };
}  // namespace GeometryPair

FOUR_C_NAMESPACE_CLOSE

#endif
