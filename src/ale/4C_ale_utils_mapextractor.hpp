// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ALE_UTILS_MAPEXTRACTOR_HPP
#define FOUR_C_ALE_UTILS_MAPEXTRACTOR_HPP

#include "4C_config.hpp"

#include "4C_linalg_mapextractor.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace ALE
{
  namespace Utils
  {
    /// specific MultiMapExtractor to handle the ale field
    class MapExtractor : public Core::LinAlg::MultiMapExtractor
    {
     public:
      enum
      {
        cond_other = 0,
        cond_fsi = 1,
        cond_fs = 2,
        cond_lung_asi = 3,
        cond_ale_wear = 4,
        cond_bio_gr = 5,
        cond_au = 6,
        cond_fpsi = 7,
        cond_mortar = 8
      };

      /// application-specific types of Dirichlet sets
      enum AleDBCSetType
      {
        dbc_set_std = 0,       ///< type of Dirichlet set in standard ALE-FSI
        dbc_set_x_ff = 1,      ///< Dirichlet sets include fluid-fluid-interface DOF (during Newton
                               ///< iteration of XFFSI)
        dbc_set_x_fsi = 2,     ///< Dirichlet sets include ALE-sided FSI interface DOF (in
                               ///< ALE-relaxation step of XFFSI)
        dbc_set_biofilm = 3,   ///< type of Dirichlet set for biofilm applications
        dbc_set_part_fsi = 4,  ///< type of Dirichlet set for partitioned FSI
        dbc_set_wear = 5,      ///< type of Dirichlet set for wear application
        dbc_set_count = 6      ///< total number of types
      };

      /// setup the whole thing
      void setup(const Core::FE::Discretization& dis, bool overlapping = false);

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
      const std::shared_ptr<const Core::LinAlg::Map>& other_map() const { return map(cond_other); }
      bool other_relevant() const { return other_map()->num_global_elements() != 0; }
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
      const std::shared_ptr<const Core::LinAlg::Map>& fsi_cond_map() const { return map(cond_fsi); }
      bool fsi_cond_relevant() const { return fsi_cond_map()->num_global_elements() != 0; }
      void fsi_cond_put_scalar(Core::LinAlg::Vector<double>& full, double scalar) const
      {
        put_scalar(full, cond_fsi, scalar);
      }
      std::shared_ptr<Core::LinAlg::Vector<double>> extract_fs_cond_vector(
          const Core::LinAlg::Vector<double>& full) const
      {
        return MultiMapExtractor::extract_vector(full, cond_fs);
      }
      void extract_fs_cond_vector(
          const Core::LinAlg::Vector<double>& full, Core::LinAlg::Vector<double>& cond) const
      {
        extract_vector(full, cond_fs, cond);
      }
      std::shared_ptr<Core::LinAlg::Vector<double>> insert_fs_cond_vector(
          const Core::LinAlg::Vector<double>& cond) const
      {
        return insert_vector(cond, cond_fs);
      }
      void insert_fs_cond_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        insert_vector(cond, cond_fs, full);
      }
      void add_fs_cond_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_fs, full);
      }
      void add_fs_cond_vector(double scale, const Core::LinAlg::Vector<double>& cond,
          Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_fs, full, scale);
      }
      const std::shared_ptr<const Core::LinAlg::Map>& fs_cond_map() const { return map(cond_fs); }
      bool fs_cond_relevant() const { return fs_cond_map()->num_global_elements() != 0; }
      void fs_cond_put_scalar(Core::LinAlg::Vector<double>& full, double scalar) const
      {
        put_scalar(full, cond_fs, scalar);
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
        return map(cond_lung_asi);
      }
      bool lung_asi_cond_relevant() const
      {
        return lung_asi_cond_map()->num_global_elements() != 0;
      }
      void lung_asi_cond_put_scalar(Core::LinAlg::Vector<double>& full, double scalar) const
      {
        put_scalar(full, cond_lung_asi, scalar);
      }
      std::shared_ptr<Core::LinAlg::Vector<double>> extract_ale_wear_cond_vector(
          const Core::LinAlg::Vector<double>& full) const
      {
        return MultiMapExtractor::extract_vector(full, cond_ale_wear);
      }
      void extract_ale_wear_cond_vector(
          const Core::LinAlg::Vector<double>& full, Core::LinAlg::Vector<double>& cond) const
      {
        extract_vector(full, cond_ale_wear, cond);
      }
      std::shared_ptr<Core::LinAlg::Vector<double>> insert_ale_wear_cond_vector(
          const Core::LinAlg::Vector<double>& cond) const
      {
        return insert_vector(cond, cond_ale_wear);
      }
      void insert_ale_wear_cond_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        insert_vector(cond, cond_ale_wear, full);
      }
      void add_ale_wear_cond_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_ale_wear, full);
      }
      void add_ale_wear_cond_vector(double scale, const Core::LinAlg::Vector<double>& cond,
          Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_ale_wear, full, scale);
      }
      const std::shared_ptr<const Core::LinAlg::Map>& ale_wear_cond_map() const
      {
        return map(cond_ale_wear);
      }
      bool ale_wear_cond_relevant() const
      {
        return ale_wear_cond_map()->num_global_elements() != 0;
      }
      void ale_wear_cond_put_scalar(Core::LinAlg::Vector<double>& full, double scalar) const
      {
        put_scalar(full, cond_ale_wear, scalar);
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
      const std::shared_ptr<const Core::LinAlg::Map>& au_cond_map() const { return map(cond_au); }
      bool au_cond_relevant() const { return au_cond_map()->num_global_elements() != 0; }
      void au_cond_put_scalar(Core::LinAlg::Vector<double>& full, double scalar) const
      {
        put_scalar(full, cond_au, scalar);
      }
      std::shared_ptr<Core::LinAlg::Vector<double>> extract_fpsi_cond_vector(
          const Core::LinAlg::Vector<double>& full) const
      {
        return MultiMapExtractor::extract_vector(full, cond_fpsi);
      }
      void extract_fpsi_cond_vector(
          const Core::LinAlg::Vector<double>& full, Core::LinAlg::Vector<double>& cond) const
      {
        extract_vector(full, cond_fpsi, cond);
      }
      std::shared_ptr<Core::LinAlg::Vector<double>> insert_fpsi_cond_vector(
          const Core::LinAlg::Vector<double>& cond) const
      {
        return insert_vector(cond, cond_fpsi);
      }
      void insert_fpsi_cond_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        insert_vector(cond, cond_fpsi, full);
      }
      void add_fpsi_cond_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_fpsi, full);
      }
      void add_fpsi_cond_vector(double scale, const Core::LinAlg::Vector<double>& cond,
          Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_fpsi, full, scale);
      }
      const std::shared_ptr<const Core::LinAlg::Map>& fpsi_cond_map() const
      {
        return map(cond_fpsi);
      }
      bool fpsi_cond_relevant() const { return fpsi_cond_map()->num_global_elements() != 0; }
      void fpsi_cond_put_scalar(Core::LinAlg::Vector<double>& full, double scalar) const
      {
        put_scalar(full, cond_fpsi, scalar);
      }
      std::shared_ptr<Core::LinAlg::Vector<double>> extract_mortar_vector(
          const Core::LinAlg::Vector<double>& full) const
      {
        return MultiMapExtractor::extract_vector(full, cond_mortar);
      }
      void extract_mortar_vector(
          const Core::LinAlg::Vector<double>& full, Core::LinAlg::Vector<double>& cond) const
      {
        extract_vector(full, cond_mortar, cond);
      }
      std::shared_ptr<Core::LinAlg::Vector<double>> insert_mortar_vector(
          const Core::LinAlg::Vector<double>& cond) const
      {
        return insert_vector(cond, cond_mortar);
      }
      void insert_mortar_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        insert_vector(cond, cond_mortar, full);
      }
      void add_mortar_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_mortar, full);
      }
      void add_mortar_vector(double scale, const Core::LinAlg::Vector<double>& cond,
          Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_mortar, full, scale);
      }
      const std::shared_ptr<const Core::LinAlg::Map>& mortar_map() const
      {
        return map(cond_mortar);
      }
      bool mortar_relevant() const { return mortar_map()->num_global_elements() != 0; }
      void mortar_put_scalar(Core::LinAlg::Vector<double>& full, double scalar) const
      {
        put_scalar(full, cond_mortar, scalar);
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
      const std::shared_ptr<const Core::LinAlg::Map>& other_map() const { return map(cond_other); }
      bool other_relevant() const { return other_map()->num_global_elements() != 0; }
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
    };

    /// specific MultiMapExtractor to handle the fluid_fluid_Coupling
    class XFluidFluidMapExtractor : public Core::LinAlg::MultiMapExtractor
    {
     public:
      enum
      {
        cond_other = 0,
        cond_xfluidfluid = 1
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
      const std::shared_ptr<const Core::LinAlg::Map>& other_map() const { return map(cond_other); }
      bool other_relevant() const { return other_map()->num_global_elements() != 0; }
      void other_put_scalar(Core::LinAlg::Vector<double>& full, double scalar) const
      {
        put_scalar(full, cond_other, scalar);
      }
      std::shared_ptr<Core::LinAlg::Vector<double>> extract_xfluid_fluid_cond_vector(
          const Core::LinAlg::Vector<double>& full) const
      {
        return MultiMapExtractor::extract_vector(full, cond_xfluidfluid);
      }
      void extract_xfluid_fluid_cond_vector(
          const Core::LinAlg::Vector<double>& full, Core::LinAlg::Vector<double>& cond) const
      {
        extract_vector(full, cond_xfluidfluid, cond);
      }
      std::shared_ptr<Core::LinAlg::Vector<double>> insert_xfluid_fluid_cond_vector(
          const Core::LinAlg::Vector<double>& cond) const
      {
        return insert_vector(cond, cond_xfluidfluid);
      }
      void insert_xfluid_fluid_cond_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        insert_vector(cond, cond_xfluidfluid, full);
      }
      void add_xfluid_fluid_cond_vector(
          const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_xfluidfluid, full);
      }
      void add_xfluid_fluid_cond_vector(double scale, const Core::LinAlg::Vector<double>& cond,
          Core::LinAlg::Vector<double>& full) const
      {
        add_vector(cond, cond_xfluidfluid, full, scale);
      }
      const std::shared_ptr<const Core::LinAlg::Map>& xfluid_fluid_cond_map() const
      {
        return map(cond_xfluidfluid);
      }
      bool xfluid_fluid_cond_relevant() const
      {
        return xfluid_fluid_cond_map()->num_global_elements() != 0;
      }
      void xfluid_fluid_cond_put_scalar(Core::LinAlg::Vector<double>& full, double scalar) const
      {
        put_scalar(full, cond_xfluidfluid, scalar);
      }
    };
  }  // namespace Utils
}  // namespace ALE

FOUR_C_NAMESPACE_CLOSE

#endif
