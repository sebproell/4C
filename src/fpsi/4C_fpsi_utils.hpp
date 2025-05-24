// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FPSI_UTILS_HPP
#define FOUR_C_FPSI_UTILS_HPP

/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_linalg_mapextractor.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
namespace FPSI
{
  class FpsiBase;

  class InterfaceUtils
  {
   public:
    //! Singleton access method
    static InterfaceUtils* instance();

    //! singleton object
    static std::shared_ptr<FPSI::InterfaceUtils> instance_;

    //! Setup Discretizations for FPSI problem (clone ALE and porofluid and setup interfaces)
    std::shared_ptr<FPSI::FpsiBase> setup_discretizations(MPI_Comm comm,
        const Teuchos::ParameterList& fpsidynparams,
        const Teuchos::ParameterList& poroelastdynparams);

    //! redistribute interface for parallel computations
    void redistribute_interface(Core::FE::Discretization& masterdis, const std::string& condname,
        std::map<int, int>& interfacefacingelementmap);

    //! build map for fpsi interface
    void setup_interface_map(MPI_Comm comm, Core::FE::Discretization& structdis,
        std::shared_ptr<Core::FE::Discretization> porofluiddis,
        std::shared_ptr<Core::FE::Discretization> fluiddis, Core::FE::Discretization& aledis);

    //! Fills a map that matches the global id of an interface element on the slave side to the
    //! global id of the opposing bulk element. This is done processor locally. Works only for
    //! matching grids.
    /*!
      \param In
             masterdis - Reference of discretization of master field.
      \param In
             slavedis  - Reference of discretization of slave field.
      \param In
             condname  - String with name of condition on interface to be considered.
      \param Out
             interfacefacingelementmap - processor local map to be filled

       See Detailed Description section for further discussion.
    */
    void setup_local_interface_facing_element_map(Core::FE::Discretization& masterdis,
        const Core::FE::Discretization& slavedis, const std::string& condname,
        std::map<int, int>& interfacefacingelementmap);

    //! access methods
    std::shared_ptr<std::map<int, int>> get_fluid_poro_fluid_interface_map()
    {
      return fluid_poro_fluid_interface_map_;
    };
    std::shared_ptr<std::map<int, int>> get_poro_fluid_fluid_interface_map()
    {
      return poro_fluid_fluid_interface_map_;
    };

   private:
    InterfaceUtils() = default;

    //! interface maps
    std::shared_ptr<std::map<int, int>> fluid_poro_fluid_interface_map_;
    std::shared_ptr<std::map<int, int>> poro_fluid_fluid_interface_map_;

  };  // class Utils


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
        cond_fpsi = 2
      };

      /// setup the whole thing
      void setup(
          const Core::FE::Discretization& dis, bool withpressure = false, bool overlapping = false);

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
          const FPSI::Utils::MapExtractor& extractor);

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
    };

  }  // namespace Utils
}  // namespace FPSI

FOUR_C_NAMESPACE_CLOSE

#endif
