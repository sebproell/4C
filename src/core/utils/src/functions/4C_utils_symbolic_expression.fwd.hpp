// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_UTILS_SYMBOLIC_EXPRESSION_FWD_HPP
#define FOUR_C_UTILS_SYMBOLIC_EXPRESSION_FWD_HPP

#include "4C_config.hpp"

#include <algorithm>
#include <cstddef>

FOUR_C_NAMESPACE_OPEN

namespace Core::Utils
{
  /**
   * A compile-time string class that can be used to represent strings as template parameters.
   */
  template <std::size_t n>
  struct CompileTimeString
  {
    char value[n];

    //! Implicit conversion from a string literal to a CompileTimeString.
    constexpr CompileTimeString(const char (&str)[n]) { std::copy_n(str, n, value); }

    template <std::size_t m>
    constexpr bool operator==(const CompileTimeString<m>& other) const
    {
      if (n != m) return false;
      for (std::size_t i = 0; i < n; ++i)
      {
        if (value[i] != other.value[i]) return false;
      }
      return true;
    }
  };

  // Deduction guide.
  template <size_t n>
  CompileTimeString(const char (&)[n]) -> CompileTimeString<n>;

  /**
   * Get index of CompileTimeString @p s in the list of CompileTimeStrings @p strings. Returns -1 if
   * @p s is not found in @p strings.
   */
  template <CompileTimeString s, CompileTimeString... strings>
  consteval int index_of()
  {
    int index = 0;

    bool found = ((s == strings ? true : (++index, false)) || ...);
    if (!found) return -1;
    return index;
  }


  /**
   * @brief A variable in a symbolic expression.
   *
   * Represents a compile-time variable @p name and a runtime @p value.
   */
  template <CompileTimeString name, typename Number>
  struct VarWrapper
  {
    Number value;
  };

  /**
   * Create a variable for passing it to SymbolicExpression.
   */
  template <CompileTimeString name, typename Number>
  constexpr VarWrapper<name, Number> var(Number value)
  {
    return VarWrapper<name, Number>{value};
  }


  template <typename Number, CompileTimeString... variables>
  class SymbolicExpression;
}  // namespace Core::Utils

FOUR_C_NAMESPACE_CLOSE

#endif
