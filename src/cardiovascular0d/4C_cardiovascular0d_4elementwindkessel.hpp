// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_CARDIOVASCULAR0D_4ELEMENTWINDKESSEL_HPP
#define FOUR_C_CARDIOVASCULAR0D_4ELEMENTWINDKESSEL_HPP

#include "4C_config.hpp"

#include "4C_cardiovascular0d.hpp"
#include "4C_fem_condition.hpp"
#include "4C_fem_general_utils_integration.hpp"
#include "4C_inpar_cardiovascular0d.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::LinAlg
{
  class SparseMatrix;
  class SparseOperator;
}  // namespace Core::LinAlg

namespace Utils
{
  /*! \brief Four-element Windkessel model without valves
   *
   * Boundary condition: DESIGN SURF CARDIOVASCULAR 0D WINDKESSEL ONLY CONDITIONS
   *
   * Original equation is:
   * C * dp/dt + (p - p_ref)/R_p - (1 + Z_c/R_p) * Q - (C Z_c  + * L/R_p) * dQ/dt - L C * d2Q/dt2 =
   * 0 where Q := -dV/dt
   *
   * -> we reformulate the ODE to three first-order ODEs with variables p, q, s:
   *
   *       [C * dp/dt + (p - p_ref)/R_p + (1 + Z_c/R_p) * q + (C Z_c  + L/R_p) * s + L C * ds/dt] [
   * 0 ] Res = [dV/dt - q ] = [ 0 ] [dq/dt - s ]   [ 0 ]
   *
   * The classical 3- or 2-element windkessel models are reproduced by setting L, or L and Z_c to
   * zero, respectively
   *
   */
  class Cardiovascular0D4ElementWindkessel : public Cardiovascular0D

  {
   public:
    /*!
    \brief Constructor of a Cardiovascular0D based on conditions with a given name. It also
    takes care of the Cardiovascular0D IDs.
    */

    Cardiovascular0D4ElementWindkessel(
        std::shared_ptr<Core::FE::Discretization>
            discr,                         ///< discretization where Cardiovascular0D lives on
        const std::string& conditionname,  ///< Name of condition to create Cardiovascular0D from
        std::vector<int>& curID            ///< current ID
    );



    /// initialization routine called by the manager ctor to get correct reference base values and
    /// activating the right conditions at the beginning
    void initialize(
        Teuchos::ParameterList&
            params,  ///< parameter list to communicate between elements and discretization
        std::shared_ptr<Core::LinAlg::Vector<double>>
            sysvec1,  ///< distributed vector that may be filled by
                      ///< assembly of element contributions
        std::shared_ptr<Core::LinAlg::Vector<double>>
            sysvec2  ///< distributed vector that may be filled by assembly of element contributions
        ) override;

    //! Evaluate routine to call from outside. In here the right action is determined and the
    //! #EvaluateCardiovascular0D routine is called
    void evaluate(
        Teuchos::ParameterList&
            params,  ///< parameter list to communicate between elements and discretization
        std::shared_ptr<Core::LinAlg::SparseMatrix> sysmat1,  ///< Cardiovascular0D stiffness matrix
        std::shared_ptr<Core::LinAlg::SparseOperator>
            sysmat2,  ///< Cardiovascular0D offdiagonal matrix dV/dd
        std::shared_ptr<Core::LinAlg::SparseOperator>
            sysmat3,  ///< Cardiovascular0D offdiagonal matrix dfext/dp
        std::shared_ptr<Core::LinAlg::Vector<double>>
            sysvec1,  ///< distributed vectors that may be filled by
                      ///< assembly of element contributions
        std::shared_ptr<Core::LinAlg::Vector<double>> sysvec2,
        std::shared_ptr<Core::LinAlg::Vector<double>> sysvec3,
        const std::shared_ptr<Core::LinAlg::Vector<double>> sysvec4,
        std::shared_ptr<Core::LinAlg::Vector<double>> sysvec5) override;

   private:
    // don't want = operator, cctor and destructor

    Cardiovascular0D4ElementWindkessel operator=(const Cardiovascular0D4ElementWindkessel& old);
    Cardiovascular0D4ElementWindkessel(const Cardiovascular0D4ElementWindkessel& old);



  };  // class
}  // namespace Utils

FOUR_C_NAMESPACE_CLOSE

#endif
