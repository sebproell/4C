// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_UTILS_SYMBOLIC_EXPRESSION_HPP
#define FOUR_C_UTILS_SYMBOLIC_EXPRESSION_HPP

#include "4C_config.hpp"

#include "4C_utils_symbolic_expression_details.hpp"

#include <memory>
#include <numeric>
#include <string>

FOUR_C_NAMESPACE_OPEN

namespace Core::Utils
{
  /*!
   *  @brief The SymbolicExpression class evaluates and forms the first and second derivatives of
   * arbitrary symbolic expressions.
   *
   * The class constructs a SymbolicExpression from an expression string. The expression must only
   * contain supported functions ("acos", "asin", "atan", "cos", "sin", "tan", "cosh", "sinh",
   * "tanh", "exp", "log", "log10", "sqrt", "heaviside", "fabs", "atan2") literals
   * ('1.0', 'pi', 'e', 'E', etc) and supported operators ("+", "-", "*", "/", "^", ".", ","). In
   * addition, an arbitrary number of variables can be contained. Any substring that is not a number
   * or supported function is parsed as a variable. When calling Value(), FirstDerivative() or
   * SecondDerivative(), the variables that have been parsed need to be supplied with a value.
   *
   * \note If you want to evaluate the same expression more than once, it is better to reuse that
   * object of the SymbolicExpression instead of creating a new object of that class with the same
   * expression so that the expression only needs to be parsed once.
   *
   * @tparam Number: Only an arithmetic type is allowed for template parameter. So far only double
   * is supported.
   */

  template <typename Number>
  class SymbolicExpression
  {
    static_assert(std::is_arithmetic_v<Number>, "Need an arithmetic type.");

   public:
    //! Type returned by the Value() function
    using ValueType = Number;
    //! Type returned by the FirstDerivative() function
    using FirstDerivativeType = Sacado::Fad::DFad<Number>;
    //! Type returned by the SecondDerivative() function
    using SecondDerivativeType = Sacado::Fad::DFad<Sacado::Fad::DFad<Number>>;

    //! Construct a SymbolicExpression from the given @p expression string. The expression must only
    //! contain supported functions, literals and operators, as well as arbitrary number of
    //! variables. See the class documentation for more details.
    SymbolicExpression(const std::string& expression);

    //! destructor
    ~SymbolicExpression();

    //! copy constructor
    SymbolicExpression(const SymbolicExpression& other);

    //! copy assignment operator
    SymbolicExpression& operator=(const SymbolicExpression& other);

    //! move constructor
    SymbolicExpression(SymbolicExpression&& other) noexcept;

    //! move assignment operator
    SymbolicExpression& operator=(SymbolicExpression&& other) noexcept;


    /*!
     * @brief Evaluates the parsed expression for a given set of variables
     *
     * The variables can be either passed as key-value pairs (preferred) or as a map of keys and
     * values (slower).
     *
     * @code
     * symbolic_expression.value("x", 1.0, "y", 2.0);
     * // or
     * symbolic_expression.value({{"x", 1.0}, {"y", 2.0}});
     * @endcode
     *
     * Failing to provide a value for a variable that is parsed in the expression will result in an
     * exception being thrown. It is *not* an error to provide a value for a variable that is
     * not parsed in the expression; the variable is simply ignored.
     *
     * @return A ValueType with the value of the expression evaluated for the given variables.
     */
    template <typename... Args>
    ValueType value(Args&&... args) const;


    /*!
     * @brief Evaluates the first derivative of the parsed expression.
     *
     * This functions makes use of the FirstDerivativeType to compute the first derivative. The
     * variables can be either passed as key-value pairs (preferred) or as a map of keys and
     * values (slower). The values need to be of FirstDerivativeType. See value() for the
     * general syntax.
     *
     * @return A FirstDerivativeType object that contains the first derivatives of the expression
     * with respect to the variables that were passed.
     */
    template <typename... Args>
    FirstDerivativeType first_derivative(Args&&... args) const;


    /*!
     * @brief Evaluates the second derivative of the parsed expression.
     *
     * This functions makes use of the SecondDerivativeType to compute the second derivative. The
     * variables can be either passed as key-value pairs (preferred) or as a map of keys and
     * values (slower). The values need to be of SecondDerivativeType. See value() for the
     * general syntax.
     *
     * @returns A SecondDerivativeType object that contains the second derivatives of the expression
     * with respect to the variables that were passed.
     */
    template <typename... Args>
    SecondDerivativeType second_derivative(Args&&... args) const;

   private:
    //! Parser for the symbolic expression evaluation
    std::unique_ptr<Core::Utils::SymbolicExpressionDetails::Parser<ValueType>> parser_for_value_;

    //! Parser for the symbolic expression first derivative evaluation
    std::unique_ptr<Core::Utils::SymbolicExpressionDetails::Parser<FirstDerivativeType>>
        parser_for_firstderivative_;

    //! Parser for the symbolic expression second derivative evaluation
    std::unique_ptr<Core::Utils::SymbolicExpressionDetails::Parser<SecondDerivativeType>>
        parser_for_secondderivative_;
  };

  template <typename T>
  template <typename... Args>
  auto SymbolicExpression<T>::value(Args&&... args) const -> ValueType
  {
    return parser_for_value_->evaluate(std::forward<Args>(args)...);
  }


  template <typename T>
  template <typename... Args>
  auto SymbolicExpression<T>::first_derivative(Args&&... args) const -> FirstDerivativeType
  {
    return parser_for_firstderivative_->evaluate(std::forward<Args>(args)...);
  }


  template <typename T>
  template <typename... Args>
  auto SymbolicExpression<T>::second_derivative(Args&&... args) const -> SecondDerivativeType
  {
    return parser_for_secondderivative_->evaluate(std::forward<Args>(args)...);
  }

}  // namespace Core::Utils

FOUR_C_NAMESPACE_CLOSE

#endif
