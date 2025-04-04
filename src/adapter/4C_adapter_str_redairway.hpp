// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ADAPTER_STR_REDAIRWAY_HPP
#define FOUR_C_ADAPTER_STR_REDAIRWAY_HPP
/*----------------------------------------------------------------------*/
/* headers */
#include "4C_config.hpp"

#include "4C_adapter_str_wrapper.hpp"
#include "4C_fem_condition.hpp"
#include "4C_linalg_vector.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  class StructureRedAirway : public StructureWrapper
  {
   public:
    /// Constructor
    StructureRedAirway(std::shared_ptr<Structure> stru);

    /// set pressure calculated from reduced-d airway tree
    void set_pressure(Core::LinAlg::Vector<double>& couppres);

    /// calculate outlet fluxes for reduced-d airway tree
    void calc_flux(
        Core::LinAlg::Vector<double>& coupflux, Core::LinAlg::Vector<double>& coupvol, double dt);

    /// calculate volume
    void calc_vol(std::map<int, double>& V);

    /// calculate initial volume
    void init_vol();

    //! (derived)
    void update() override;

   private:
    /// map between coupling ID and conditions on structure
    std::map<int, Core::Conditions::Condition*> coupcond_;

    /// map of coupling IDs
    std::shared_ptr<Core::LinAlg::Map> coupmap_;

    std::map<int, double> vn_;
    std::map<int, double> vnp_;
  };

}  // namespace Adapter
FOUR_C_NAMESPACE_CLOSE

#endif
