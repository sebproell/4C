// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SCATRA_ELE_CALC_NO_PHYSICS_HPP
#define FOUR_C_SCATRA_ELE_CALC_NO_PHYSICS_HPP

#include "4C_config.hpp"

#include "4C_scatra_ele_calc.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace Elements
  {
    /**
     * \brief ScaTra ImplType containing no physics.
     *
     * Evaluation of a scatra element that does not contain any physics. Currently, this class only
     * implements the minimal set of actions needed for reading the scatra results from a restart
     * file and simulating a one-way coupling to the structure. This ImplType can currently not be
     * used in solving the scatra equations, as the needed actions are not implemented yet.
     *
     * @tparam distype Element shape type
     * @tparam probdim Number of space dimensions of the problem
     */
    template <Core::FE::CellType distype, int probdim>
    class ScaTraEleCalcNoPhysics : public ScaTraEleCalc<distype, probdim>
    {
     public:
      //! abbreviation
      using my = ScaTraEleCalc<distype, probdim>;

      //! singleton access method
      static ScaTraEleCalcNoPhysics<distype, probdim>* instance(
          int numdofpernode, int numscal, const std::string& disname);

      //! evaluate the element
      int evaluate_action(Core::Elements::Element* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, const ScaTra::Action& action,
          Core::Elements::LocationArray& la, Core::LinAlg::SerialDenseMatrix& elemat1,
          Core::LinAlg::SerialDenseMatrix& elemat2, Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseVector& elevec2,
          Core::LinAlg::SerialDenseVector& elevec3) override;

     protected:
      //! protected constructor for singletons
      ScaTraEleCalcNoPhysics(int numdofpernode, int numscal, const std::string& disname);
    };
  }  // namespace Elements
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
