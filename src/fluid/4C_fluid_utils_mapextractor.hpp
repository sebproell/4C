// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FLUID_UTILS_MAPEXTRACTOR_HPP
#define FOUR_C_FLUID_UTILS_MAPEXTRACTOR_HPP


#include "4C_config.hpp"

#include "4C_linalg_mapextractor.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace FLD
{
  namespace Utils
  {
    /// specific MultiMapExtractor to handle the fluid field
    class MapExtractor : public Core::LinAlg::MultiMapExtractor
    {
     public:
      enum
      {
        cond_other = 0,
        cond_fsi = 1,
        cond_lung_asi = 3,
        cond_mortar = 4,
        cond_au = 5
      };

      /// setup the whole thing
      void setup(const Core::FE::Discretization& dis, bool withpressure = false,
          bool overlapping = false, const int nds_master = 0);

      /*!
       * \brief setup from an existing extractor
       * By calling this setup version we create a map extractor from
       * (1) an existing map extractor and
       * (2) a DOF-map from another discretization, which is appended to othermap.
       * We need this in the context of XFFSI.
       * \param (in) additionalothermap : map of additional unconditioned DOF
       * \param (in) extractor : extractor, from which the conditions are cloned

       */
      void setup(std::shared_ptr<const Core::LinAlg::Map>& additionalothermap,
          const FLD::Utils::MapExtractor& extractor);

      /// get all element gids those nodes are touched by any condition
      std::shared_ptr<std::set<int>> conditioned_element_map(
          const Core::FE::Discretization& dis) const;

      std::shared_ptr<Core::LinAlg::Vector<double>> extract_other_vector(
          const Core::LinAlg::Vector<double>& full) const
      {
        return MultiMapExtractor::extract_vector(full, cond_other);
      }
      void extract_other_vector(
          const Core::LinAlg::Vector<double>& full, Core::LinAlg::Vector<double>& cond) const
      {
        extract_vector(full, cond_other, cond);
      }
      std::shared_ptr<Core::LinAlg::Vector<double>> insert_other_vector(
          const Core::LinAlg::Vector<double>& cond) const
      {
        return insert_vector(cond, cond_other);
      }
      void insert_other_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        insert_vector(cond, cond_other, full);
      }
      void add_other_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_other, full);
      }
      void add_other_vector(double scale, const Core::LinAlg::Vector<double>& cond,
          Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_other, full, scale);
      }
      const std::shared_ptr<const Core::LinAlg::Map>& other_map() const { return Map(cond_other); }
      bool other_relevant() const { return other_map()->NumGlobalElements() != 0; }
      void other_put_scalar(Core::LinAlg::Vector<double>& full, double scalar) const
      {
        put_scalar(full, cond_other, scalar);
      }
      std::shared_ptr<Core::LinAlg::Vector<double>> extract_fsi_cond_vector(
          const Core::LinAlg::Vector<double>& full) const
      {
        return MultiMapExtractor::extract_vector(full, cond_fsi);
      }
      void extract_fsi_cond_vector(
          const Core::LinAlg::Vector<double>& full, Core::LinAlg::Vector<double>& cond) const
      {
        extract_vector(full, cond_fsi, cond);
      }
      std::shared_ptr<Core::LinAlg::Vector<double>> insert_fsi_cond_vector(
          const Core::LinAlg::Vector<double>& cond) const
      {
        return insert_vector(cond, cond_fsi);
      }
      void insert_fsi_cond_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        insert_vector(cond, cond_fsi, full);
      }
      void add_fsi_cond_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_fsi, full);
      }
      void add_fsi_cond_vector(double scale, const Core::LinAlg::Vector<double>& cond,
          Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_fsi, full, scale);
      }
      const std::shared_ptr<const Core::LinAlg::Map>& fsi_cond_map() const { return Map(cond_fsi); }
      bool fsi_cond_relevant() const { return fsi_cond_map()->NumGlobalElements() != 0; }
      void fsi_cond_put_scalar(Core::LinAlg::Vector<double>& full, double scalar) const
      {
        put_scalar(full, cond_fsi, scalar);
      }
      std::shared_ptr<Core::LinAlg::Vector<double>> extract_lung_asi_cond_vector(
          const Core::LinAlg::Vector<double>& full) const
      {
        return MultiMapExtractor::extract_vector(full, cond_lung_asi);
      }
      void extract_lung_asi_cond_vector(
          const Core::LinAlg::Vector<double>& full, Core::LinAlg::Vector<double>& cond) const
      {
        extract_vector(full, cond_lung_asi, cond);
      }
      std::shared_ptr<Core::LinAlg::Vector<double>> insert_lung_asi_cond_vector(
          const Core::LinAlg::Vector<double>& cond) const
      {
        return insert_vector(cond, cond_lung_asi);
      }
      void insert_lung_asi_cond_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        insert_vector(cond, cond_lung_asi, full);
      }
      void add_lung_asi_cond_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_lung_asi, full);
      }
      void add_lung_asi_cond_vector(double scale, const Core::LinAlg::Vector<double>& cond,
          Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_lung_asi, full, scale);
      }
      const std::shared_ptr<const Core::LinAlg::Map>& lung_asi_cond_map() const
      {
        return Map(cond_lung_asi);
      }
      bool lung_asi_cond_relevant() const { return lung_asi_cond_map()->NumGlobalElements() != 0; }
      void lung_asi_cond_put_scalar(Core::LinAlg::Vector<double>& full, double scalar) const
      {
        put_scalar(full, cond_lung_asi, scalar);
      }
      std::shared_ptr<Core::LinAlg::Vector<double>> extract_mortar_cond_vector(
          const Core::LinAlg::Vector<double>& full) const
      {
        return MultiMapExtractor::extract_vector(full, cond_mortar);
      }
      void extract_mortar_cond_vector(
          const Core::LinAlg::Vector<double>& full, Core::LinAlg::Vector<double>& cond) const
      {
        extract_vector(full, cond_mortar, cond);
      }
      std::shared_ptr<Core::LinAlg::Vector<double>> insert_mortar_cond_vector(
          const Core::LinAlg::Vector<double>& cond) const
      {
        return insert_vector(cond, cond_mortar);
      }
      void insert_mortar_cond_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        insert_vector(cond, cond_mortar, full);
      }
      void add_mortar_cond_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_mortar, full);
      }
      void add_mortar_cond_vector(double scale, const Core::LinAlg::Vector<double>& cond,
          Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_mortar, full, scale);
      }
      const std::shared_ptr<const Core::LinAlg::Map>& mortar_cond_map() const
      {
        return Map(cond_mortar);
      }
      bool mortar_cond_relevant() const { return mortar_cond_map()->NumGlobalElements() != 0; }
      void mortar_cond_put_scalar(Core::LinAlg::Vector<double>& full, double scalar) const
      {
        put_scalar(full, cond_mortar, scalar);
      }
      std::shared_ptr<Core::LinAlg::Vector<double>> extract_au_cond_vector(
          const Core::LinAlg::Vector<double>& full) const
      {
        return MultiMapExtractor::extract_vector(full, cond_au);
      }
      void extract_au_cond_vector(
          const Core::LinAlg::Vector<double>& full, Core::LinAlg::Vector<double>& cond) const
      {
        extract_vector(full, cond_au, cond);
      }
      std::shared_ptr<Core::LinAlg::Vector<double>> insert_au_cond_vector(
          const Core::LinAlg::Vector<double>& cond) const
      {
        return insert_vector(cond, cond_au);
      }
      void insert_au_cond_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        insert_vector(cond, cond_au, full);
      }
      void add_au_cond_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_au, full);
      }
      void add_au_cond_vector(double scale, const Core::LinAlg::Vector<double>& cond,
          Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_au, full, scale);
      }
      const std::shared_ptr<const Core::LinAlg::Map>& au_cond_map() const { return Map(cond_au); }
      bool au_cond_relevant() const { return au_cond_map()->NumGlobalElements() != 0; }
      void au_cond_put_scalar(Core::LinAlg::Vector<double>& full, double scalar) const
      {
        put_scalar(full, cond_au, scalar);
      }
    };

    /// specific MultiMapExtractor to handle the part of fluid with volumetric surface flow
    /// condition
    class VolumetricFlowMapExtractor : public Core::LinAlg::MultiMapExtractor
    {
     public:
      enum
      {
        cond_other = 0,
        cond_vol_surf_flow = 1
      };

      /// setup the whole thing
      void setup(const Core::FE::Discretization& dis);

      std::shared_ptr<Core::LinAlg::Vector<double>> extract_other_vector(
          const Core::LinAlg::Vector<double>& full) const
      {
        return MultiMapExtractor::extract_vector(full, cond_other);
      }
      void extract_other_vector(
          const Core::LinAlg::Vector<double>& full, Core::LinAlg::Vector<double>& cond) const
      {
        extract_vector(full, cond_other, cond);
      }
      std::shared_ptr<Core::LinAlg::Vector<double>> insert_other_vector(
          const Core::LinAlg::Vector<double>& cond) const
      {
        return insert_vector(cond, cond_other);
      }
      void insert_other_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        insert_vector(cond, cond_other, full);
      }
      void add_other_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_other, full);
      }
      void add_other_vector(double scale, const Core::LinAlg::Vector<double>& cond,
          Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_other, full, scale);
      }
      const std::shared_ptr<const Core::LinAlg::Map>& other_map() const { return Map(cond_other); }
      bool other_relevant() const { return other_map()->NumGlobalElements() != 0; }
      void other_put_scalar(Core::LinAlg::Vector<double>& full, double scalar) const
      {
        put_scalar(full, cond_other, scalar);
      }
      std::shared_ptr<Core::LinAlg::Vector<double>> extract_volumetric_surface_flow_cond_vector(
          const Core::LinAlg::Vector<double>& full) const
      {
        return MultiMapExtractor::extract_vector(full, cond_vol_surf_flow);
      }
      void extract_volumetric_surface_flow_cond_vector(
          const Core::LinAlg::Vector<double>& full, Core::LinAlg::Vector<double>& cond) const
      {
        extract_vector(full, cond_vol_surf_flow, cond);
      }
      std::shared_ptr<Core::LinAlg::Vector<double>> insert_volumetric_surface_flow_cond_vector(
          const Core::LinAlg::Vector<double>& cond) const
      {
        return insert_vector(cond, cond_vol_surf_flow);
      }
      void insert_volumetric_surface_flow_cond_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        insert_vector(cond, cond_vol_surf_flow, full);
      }
      void add_volumetric_surface_flow_cond_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_vol_surf_flow, full);
      }
      void add_volumetric_surface_flow_cond_vector(double scale,
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_vol_surf_flow, full, scale);
      }
      const std::shared_ptr<const Core::LinAlg::Map>& volumetric_surface_flow_cond_map() const
      {
        return Map(cond_vol_surf_flow);
      }
      bool volumetric_surface_flow_cond_relevant() const
      {
        return volumetric_surface_flow_cond_map()->NumGlobalElements() != 0;
      }
      void volumetric_surface_flow_cond_put_scalar(
          Core::LinAlg::Vector<double>& full, double scalar) const
      {
        put_scalar(full, cond_vol_surf_flow, scalar);
      }
    };

    /// specific MultiMapExtractor to handle the part of fluid with Krylov space projection
    class KSPMapExtractor : public Core::LinAlg::MultiMapExtractor
    {
     public:
      enum
      {
        cond_other = 0,
        cond_ksp = 1
      };

      /// setup the whole thing
      void setup(const Core::FE::Discretization& dis);

      /// get all element gids those nodes are touched by any condition
      std::shared_ptr<std::set<int>> conditioned_element_map(
          const Core::FE::Discretization& dis) const;

      std::shared_ptr<Core::LinAlg::Vector<double>> extract_other_vector(
          const Core::LinAlg::Vector<double>& full) const
      {
        return MultiMapExtractor::extract_vector(full, cond_other);
      }
      void extract_other_vector(
          const Core::LinAlg::Vector<double>& full, Core::LinAlg::Vector<double>& cond) const
      {
        extract_vector(full, cond_other, cond);
      }
      std::shared_ptr<Core::LinAlg::Vector<double>> insert_other_vector(
          const Core::LinAlg::Vector<double>& cond) const
      {
        return insert_vector(cond, cond_other);
      }
      void insert_other_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        insert_vector(cond, cond_other, full);
      }
      void add_other_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_other, full);
      }
      void add_other_vector(double scale, const Core::LinAlg::Vector<double>& cond,
          Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_other, full, scale);
      }
      const std::shared_ptr<const Core::LinAlg::Map>& other_map() const { return Map(cond_other); }
      bool other_relevant() const { return other_map()->NumGlobalElements() != 0; }
      void other_put_scalar(Core::LinAlg::Vector<double>& full, double scalar) const
      {
        put_scalar(full, cond_other, scalar);
      }
      std::shared_ptr<Core::LinAlg::Vector<double>> extract_ksp_cond_vector(
          const Core::LinAlg::Vector<double>& full) const
      {
        return MultiMapExtractor::extract_vector(full, cond_ksp);
      }
      void extract_ksp_cond_vector(
          const Core::LinAlg::Vector<double>& full, Core::LinAlg::Vector<double>& cond) const
      {
        extract_vector(full, cond_ksp, cond);
      }
      std::shared_ptr<Core::LinAlg::Vector<double>> insert_ksp_cond_vector(
          const Core::LinAlg::Vector<double>& cond) const
      {
        return insert_vector(cond, cond_ksp);
      }
      void insert_ksp_cond_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        insert_vector(cond, cond_ksp, full);
      }
      void add_ksp_cond_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_ksp, full);
      }
      void add_ksp_cond_vector(double scale, const Core::LinAlg::Vector<double>& cond,
          Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_ksp, full, scale);
      }
      const std::shared_ptr<const Core::LinAlg::Map>& ksp_cond_map() const { return Map(cond_ksp); }
      bool ksp_cond_relevant() const { return ksp_cond_map()->NumGlobalElements() != 0; }
      void ksp_cond_put_scalar(Core::LinAlg::Vector<double>& full, double scalar) const
      {
        put_scalar(full, cond_ksp, scalar);
      }
    };

    /// specific MultiMapExtractor to handle the velocity-pressure split
    class VelPressExtractor : public Core::LinAlg::MultiMapExtractor
    {
     public:
      /// setup the whole thing
      void setup(const Core::FE::Discretization& dis);

      std::shared_ptr<Core::LinAlg::Vector<double>> extract_velocity_vector(
          const Core::LinAlg::Vector<double>& full) const
      {
        return MultiMapExtractor::extract_vector(full, 0);
      }
      void extract_velocity_vector(
          const Core::LinAlg::Vector<double>& full, Core::LinAlg::Vector<double>& cond) const
      {
        extract_vector(full, 0, cond);
      }
      std::shared_ptr<Core::LinAlg::Vector<double>> insert_velocity_vector(
          const Core::LinAlg::Vector<double>& cond) const
      {
        return insert_vector(cond, 0);
      }
      void insert_velocity_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        insert_vector(cond, 0, full);
      }
      void add_velocity_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, 0, full);
      }
      void add_velocity_vector(double scale, const Core::LinAlg::Vector<double>& cond,
          Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, 0, full, scale);
      }
      const std::shared_ptr<const Core::LinAlg::Map>& velocity_map() const { return Map(0); }
      bool velocity_relevant() const { return velocity_map()->NumGlobalElements() != 0; }
      void velocity_put_scalar(Core::LinAlg::Vector<double>& full, double scalar) const
      {
        put_scalar(full, 0, scalar);
      }
      std::shared_ptr<Core::LinAlg::Vector<double>> extract_pressure_vector(
          const Core::LinAlg::Vector<double>& full) const
      {
        return MultiMapExtractor::extract_vector(full, 1);
      }
      void extract_pressure_vector(
          const Core::LinAlg::Vector<double>& full, Core::LinAlg::Vector<double>& cond) const
      {
        extract_vector(full, 1, cond);
      }
      std::shared_ptr<Core::LinAlg::Vector<double>> insert_pressure_vector(
          const Core::LinAlg::Vector<double>& cond) const
      {
        return insert_vector(cond, 1);
      }
      void insert_pressure_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        insert_vector(cond, 1, full);
      }
      void add_pressure_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, 1, full);
      }
      void add_pressure_vector(double scale, const Core::LinAlg::Vector<double>& cond,
          Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, 1, full, scale);
      }
      const std::shared_ptr<const Core::LinAlg::Map>& pressure_map() const { return Map(1); }
      bool pressure_relevant() const { return pressure_map()->NumGlobalElements() != 0; }
      void pressure_put_scalar(Core::LinAlg::Vector<double>& full, double scalar) const
      {
        put_scalar(full, 1, scalar);
      }
    };

    /// specific MultiMapExtractor to handle the fsi and ale meshtying at the same time
    class FsiMapExtractor : public Core::LinAlg::MultiMapExtractor
    {
     public:
      enum
      {
        cond_other = 0,
        cond_fsi = 1
      };

      /// setup the whole thing
      void setup(const Core::FE::Discretization& dis);

      void setup(std::shared_ptr<const Core::LinAlg::Map>& additionalothermap,
          const FLD::Utils::FsiMapExtractor& extractor);

      std::shared_ptr<Core::LinAlg::Vector<double>> extract_other_vector(
          const Core::LinAlg::Vector<double>& full) const
      {
        return MultiMapExtractor::extract_vector(full, cond_other);
      }
      void extract_other_vector(
          const Core::LinAlg::Vector<double>& full, Core::LinAlg::Vector<double>& cond) const
      {
        extract_vector(full, cond_other, cond);
      }
      std::shared_ptr<Core::LinAlg::Vector<double>> insert_other_vector(
          const Core::LinAlg::Vector<double>& cond) const
      {
        return insert_vector(cond, cond_other);
      }
      void insert_other_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        insert_vector(cond, cond_other, full);
      }
      void add_other_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_other, full);
      }
      void add_other_vector(double scale, const Core::LinAlg::Vector<double>& cond,
          Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_other, full, scale);
      }
      const std::shared_ptr<const Core::LinAlg::Map>& other_map() const { return Map(cond_other); }
      bool other_relevant() const { return other_map()->NumGlobalElements() != 0; }
      void other_put_scalar(Core::LinAlg::Vector<double>& full, double scalar) const
      {
        put_scalar(full, cond_other, scalar);
      }
      std::shared_ptr<Core::LinAlg::Vector<double>> extract_fsi_vector(
          const Core::LinAlg::Vector<double>& full) const
      {
        return MultiMapExtractor::extract_vector(full, cond_fsi);
      }
      void extract_fsi_vector(
          const Core::LinAlg::Vector<double>& full, Core::LinAlg::Vector<double>& cond) const
      {
        extract_vector(full, cond_fsi, cond);
      }
      std::shared_ptr<Core::LinAlg::Vector<double>> insert_fsi_vector(
          const Core::LinAlg::Vector<double>& cond) const
      {
        return insert_vector(cond, cond_fsi);
      }
      void insert_fsi_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        insert_vector(cond, cond_fsi, full);
      }
      void add_fsi_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_fsi, full);
      }
      void add_fsi_vector(double scale, const Core::LinAlg::Vector<double>& cond,
          Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_fsi, full, scale);
      }
      const std::shared_ptr<const Core::LinAlg::Map>& fsi_map() const { return Map(cond_fsi); }
      bool fsi_relevant() const { return fsi_map()->NumGlobalElements() != 0; }
      void fsi_put_scalar(Core::LinAlg::Vector<double>& full, double scalar) const
      {
        put_scalar(full, cond_fsi, scalar);
      }
    };

    /// specific MultiMapExtractor to handle the fluid field
    class XFluidFluidMapExtractor : public Core::LinAlg::MultiMapExtractor
    {
     public:
      enum
      {
        cond_fluid = 0,
        cond_xfluid = 1,
      };

      /// setup the whole thing
      void setup(const Core::LinAlg::Map& fullmap,
          std::shared_ptr<const Core::LinAlg::Map> fluidmap,
          std::shared_ptr<const Core::LinAlg::Map> xfluidmap);

      std::shared_ptr<Core::LinAlg::Vector<double>> extract_fluid_vector(
          const Core::LinAlg::Vector<double>& full) const
      {
        return MultiMapExtractor::extract_vector(full, cond_fluid);
      }
      void extract_fluid_vector(
          const Core::LinAlg::Vector<double>& full, Core::LinAlg::Vector<double>& cond) const
      {
        extract_vector(full, cond_fluid, cond);
      }
      std::shared_ptr<Core::LinAlg::Vector<double>> insert_fluid_vector(
          const Core::LinAlg::Vector<double>& cond) const
      {
        return insert_vector(cond, cond_fluid);
      }
      void insert_fluid_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        insert_vector(cond, cond_fluid, full);
      }
      void add_fluid_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_fluid, full);
      }
      void add_fluid_vector(double scale, const Core::LinAlg::Vector<double>& cond,
          Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_fluid, full, scale);
      }
      const std::shared_ptr<const Core::LinAlg::Map>& fluid_map() const { return Map(cond_fluid); }
      bool fluid_relevant() const { return fluid_map()->NumGlobalElements() != 0; }
      void fluid_put_scalar(Core::LinAlg::Vector<double>& full, double scalar) const
      {
        put_scalar(full, cond_fluid, scalar);
      }
      std::shared_ptr<Core::LinAlg::Vector<double>> extract_x_fluid_vector(
          const Core::LinAlg::Vector<double>& full) const
      {
        return MultiMapExtractor::extract_vector(full, cond_xfluid);
      }
      void extract_x_fluid_vector(
          const Core::LinAlg::Vector<double>& full, Core::LinAlg::Vector<double>& cond) const
      {
        extract_vector(full, cond_xfluid, cond);
      }
      std::shared_ptr<Core::LinAlg::Vector<double>> insert_x_fluid_vector(
          const Core::LinAlg::Vector<double>& cond) const
      {
        return insert_vector(cond, cond_xfluid);
      }
      void insert_x_fluid_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        insert_vector(cond, cond_xfluid, full);
      }
      void add_x_fluid_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_xfluid, full);
      }
      void add_x_fluid_vector(double scale, const Core::LinAlg::Vector<double>& cond,
          Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_xfluid, full, scale);
      }
      const std::shared_ptr<const Core::LinAlg::Map>& x_fluid_map() const
      {
        return Map(cond_xfluid);
      }
      bool x_fluid_relevant() const { return x_fluid_map()->NumGlobalElements() != 0; }
      void x_fluid_put_scalar(Core::LinAlg::Vector<double>& full, double scalar) const
      {
        put_scalar(full, cond_xfluid, scalar);
      }
    };

  }  // namespace Utils
}  // namespace FLD

FOUR_C_NAMESPACE_CLOSE

#endif
