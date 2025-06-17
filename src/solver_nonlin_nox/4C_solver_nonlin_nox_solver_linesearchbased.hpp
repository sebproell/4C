// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLVER_NONLIN_NOX_SOLVER_LINESEARCHBASED_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_SOLVER_LINESEARCHBASED_HPP

#include "4C_config.hpp"

#include "4C_solver_nonlin_nox_forward_decl.hpp"

#include <NOX_Solver_LineSearchBased.H>  // base class

FOUR_C_NAMESPACE_OPEN

namespace NOX
{
  namespace Nln
  {
    namespace StatusTest
    {
      enum QuantityType : int;
    }  // namespace StatusTest
    namespace Inner
    {
      namespace StatusTest
      {
        class Generic;
      }  // namespace StatusTest
    }  // namespace Inner
    namespace Solver
    {
      class LineSearchBased : public ::NOX::Solver::LineSearchBased
      {
       public:
        //! Constructor
        /*!
          See reset(::NOX::Abstract::Group&, ::NOX::StatusTest::Generic&, Teuchos::ParameterList&)
          for description
         */
        LineSearchBased(const Teuchos::RCP<::NOX::Abstract::Group>& grp,
            const Teuchos::RCP<::NOX::StatusTest::Generic>& outerTests,
            const Teuchos::RCP<NOX::Nln::Inner::StatusTest::Generic>& innerTests,
            const Teuchos::RCP<Teuchos::ParameterList>& params);

        //! Returns the stopping test.
        virtual const ::NOX::StatusTest::Generic& get_outer_status_test() const;

        //! \brief Returns a pointer to one specific stopping test %T.
        /** If there is no outer test of type T, a nullptr pointer will be returned.
         *
         *  */
        template <class T>
        ::NOX::StatusTest::Generic* get_outer_status_test() const;

        //! Returns the ::NOX::Utils object
        [[nodiscard]] const ::NOX::Utils& get_utils() const;

        //! Returns the global status of the given test name
        template <class T>
        [[nodiscard]] ::NOX::StatusTest::StatusType get_status() const;

       protected:
        //! initialize additional variables or overwrite base class initialization
        virtual void init(const Teuchos::RCP<NOX::Nln::Inner::StatusTest::Generic>& innerTests);

        void printUpdate() override;
      };  // class LineSearchBased
    }  // namespace Solver
  }  // namespace Nln
}  // namespace NOX

FOUR_C_NAMESPACE_CLOSE

#endif
