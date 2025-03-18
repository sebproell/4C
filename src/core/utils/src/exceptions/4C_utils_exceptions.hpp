// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_UTILS_EXCEPTIONS_HPP
#define FOUR_C_UTILS_EXCEPTIONS_HPP

#include "4C_config.hpp"

#include <magic_enum/magic_enum.hpp>

#include <format>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>

// Specialize formatting for any enum type.
template <typename T>
  requires std::is_enum_v<T>
struct std::formatter<T> : std::formatter<std::string_view>
{
  auto format(T e, std::format_context& ctx) const
  {
    return std::formatter<std::string_view>::format(magic_enum::enum_name(e), ctx);
  }
};

FOUR_C_NAMESPACE_OPEN

namespace Core
{
  namespace Internal
  {
    class ExceptionImplementation;

    /**
     * A helper struct taking the file name and line number from the error macro.
     */
    struct ErrorHelper
    {
      const char* file_name;
      int line_number;
      const char* failed_assertion_string = nullptr;
    };

    /**
     * Check if a string literal contains C-style format specifiers.
     */
    [[nodiscard]] consteval bool contains_c_format_specifiers(std::string_view str)
    {
      for (size_t i = 0; i < str.size(); ++i)
      {
        // A rough check for C-style format specifiers
        if (str[i] == '%' && i + 1 < str.size() && str[i + 1] != '%') return true;
      }
      return false;
    }

    template <typename... Args>
    [[noreturn]] void constexpr format_and_throw_error(
        ErrorHelper error_helper, std::format_string<Args...> fmt, Args&&... args)
    {
      throw_error(error_helper, std::format(std::move(fmt), std::forward<Args>(args)...));
    }

    [[noreturn]] void throw_error(
        const ErrorHelper& error_helper, const std::string& formatted_message);
  }  // namespace Internal

  /**
   * @brief Base class for all 4C exceptions.
   *
   * Any exceptions generated directly by 4C will have this or a derived type. This allows to
   * catch 4C exceptions specifically by using this type in the `catch` clause.
   */
  class Exception : public std::exception
  {
   public:
    /**
     * Generate an Exception with the given message. A stacktrace is automatically attached to this
     * Exception.
     */
    explicit Exception(std::string message);

    /**
     * Destructor.
     */
    ~Exception() override;

    /**
     * Return a message that describes what happened.
     */
    [[nodiscard]] const char* what() const noexcept override;

    /**
     * Return a message that describes what happened and includes a stack trace.
     *
     * @note Calling this function can be a lot more expensive than the what() function because the
     * stacktrace needs to be symbolized.
     */
    [[nodiscard]] std::string what_with_stacktrace() const noexcept;

   private:
    /**
     * Pointer to implementation. This technique is used to minimize the footprint of the exception
     * class that is put on the stack.
     */
    std::unique_ptr<Internal::ExceptionImplementation> pimpl_;
  };
}  // namespace Core


/**
 * Throw an error in the form of a Core::Exception.
 *
 * @note Consider using the more expressive FOUR_C_ASSERT and FOUR_C_ASSERT_ALWAYS macros which
 * take a violated assertion as an argument and print it in the error message.
 *
 * This macro takes an error message, which may contain replacement fields for formatting.
 * The format arguments are passed as additional arguments. For example:
 *
 * @code
 *    FOUR_C_THROW("An error occurred in iteration {}.", iter);
 * @endcode
 */
#define FOUR_C_THROW(fmt, ...)                                                              \
  do                                                                                        \
  {                                                                                         \
    static_assert(!FourC::Core::Internal::contains_c_format_specifiers(fmt),                \
        "Use C++20 formatting with braces {} instead of C-style formatting with %.");       \
    FourC::Core::Internal::format_and_throw_error(                                          \
        FourC::Core::Internal::ErrorHelper{.file_name = __FILE__, .line_number = __LINE__}, \
        fmt __VA_OPT__(, ) __VA_ARGS__);                                                    \
  } while (0)

/**
 * Assert that @p test is `true`. If not issue an error in the form of a Core::Exception.
 *
 * This macro takes an error message, which may contain replacement fields for formatting.
 * The format arguments are passed as additional arguments. For example:
 *
 * @code
 *    FOUR_C_ASSERT_ALWAYS(vector.size() == dim, "Vector size {} does not equal dimension {}.",
 *      vector.size(), dim);
 * @endcode
 */
#define FOUR_C_ASSERT_ALWAYS(test, fmt, ...)                                                     \
  do                                                                                             \
  {                                                                                              \
    static_assert(!FourC::Core::Internal::contains_c_format_specifiers(fmt),                     \
        "Use C++20 formatting with braces {} instead of C-style formatting with %.");            \
    if (!(test))                                                                                 \
    {                                                                                            \
      FourC::Core::Internal::format_and_throw_error(                                             \
          FourC::Core::Internal::ErrorHelper{                                                    \
              .file_name = __FILE__, .line_number = __LINE__, .failed_assertion_string = #test}, \
          fmt __VA_OPT__(, ) __VA_ARGS__);                                                       \
    }                                                                                            \
  } while (0)


#ifdef FOUR_C_ENABLE_ASSERTIONS

/**
 * Assert that @p test is `true`. If not issue an error in the form of a Core::Exception.
 * This macro is only active if FOUR_C_ENABLE_ASSERTIONS is set and may therefore be used for
 * expensive checks. Use FOUR_C_ASSERT_ALWAYS if you want to evaluate the test in any case.
 *
 * This macro takes an error message, which may contain replacement fields for formatting.
 * The format arguments are passed as additional arguments. For example:
 *
 * @code
 *   FOUR_C_ASSERT(vector.size() == dim, "Vector size {} does not equal dimension {}.",
 *     vector.size(), dim);
 * @endcode
 */
#define FOUR_C_ASSERT FOUR_C_ASSERT_ALWAYS

#else

/**
 * This macro would assert that @p test is true, but only if FOUR_C_ENABLE_ASSERTIONS is set.
 */
#define FOUR_C_ASSERT(test, args...) static_assert(true, "Terminate with a comma.")

#endif


FOUR_C_NAMESPACE_CLOSE

#endif
