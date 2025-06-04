// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_IO_INPUT_SPEC_VALIDATORS_HPP
#define FOUR_C_IO_INPUT_SPEC_VALIDATORS_HPP

#include "4C_config.hpp"

#include "4C_utils_enum.hpp"
#include "4C_utils_exceptions.hpp"

#include <concepts>
#include <functional>
#include <limits>
#include <memory>
#include <set>
#include <sstream>

FOUR_C_NAMESPACE_OPEN

/**
 * This namespace provides a set of basic validators that can be used to validate input values. The
 * set of validators is deliberately kept small to avoid confusion and to keep you from doing
 * too specific checks while input is being parsed.
 */
namespace Core::IO::InputSpecBuilders::Validators
{
  /**
   * A generic validator providing the interface to check an input value against a predicate.
   * This class is heavily tied to the InputSpec system and is used to validate input values in a
   * type-erased manner. Usually, there is nothing you should do with this class directly. Use the
   * various factory functions to create instances of this class that can then be passed to one of
   * the InputSpecBuilders functions.
   */
  template <typename T>
  class Validator
  {
    using Predicate = std::function<bool(const T&)>;
    using Describe = std::function<void(std::ostream&)>;

    struct Concept
    {
      Predicate pred;
      Describe desc;

      Concept(Predicate p, Describe d) : pred(std::move(p)), desc(std::move(d)) {}

      bool operator()(const T& v) const { return pred(v); }
      void describe(std::ostream& os) const { desc(os); }
    };

    std::shared_ptr<const Concept> impl_;

   public:
    //! Type-erasing constructor.
    template <typename FPredicate, typename FDescribe>
    Validator(FPredicate p, FDescribe d)
        : impl_(std::make_shared<Concept>(Predicate(std::move(p)), Describe(std::move(d))))
    {
    }

    /**
     * The main validation function. Returns true if the value is valid according to the
     * predicate stored inside this validator. Returns false otherwise.
     *
     * @note This operator can only be called with a value that is exactly of type T without any
     * implicit conversions.
     */
    template <typename U>
      requires(std::same_as<T, std::decay_t<U>>)
    [[nodiscard]] bool operator()(const U& v) const
    {
      return (*impl_)(v);
    }

    /**
     * Describes what this validator expects from the value in human-readable form .
     */
    void describe(std::ostream& os) const { impl_->describe(os); }

    friend std::ostream& operator<<(std::ostream& os, const Validator& v)
    {
      v.describe(os);
      return os;
    }
  };

  namespace Internal
  {
    //! The types of numbers that we support in the input.
    template <typename T>
    concept Numeric = std::is_same_v<T, int> || std::is_same_v<T, double>;

    enum class InclExclType
    {
      incl,
      excl,
    };

    //! Mark values for inclusion or exclusion for the in_range() function.
    template <Numeric T>
    struct InclExclTag
    {
      T value;
      InclExclType incl_excl;

      //! Deliberately an implicit constructor defaulting to inclusion.
      InclExclTag(T v, InclExclType ie = InclExclType::incl) : value(v), incl_excl(ie) {}
    };

    template <typename T>
    struct InclExclTypeExtractor
    {
      using type = T;
    };

    template <Numeric T>
    struct InclExclTypeExtractor<InclExclTag<T>>
    {
      using type = T;
    };

    template <typename T1, typename T2>
    struct LowHighCommonType
    {
      using type = std::enable_if_t<std::is_same_v<typename InclExclTypeExtractor<T1>::type,
                                        typename InclExclTypeExtractor<T2>::type>,
          typename InclExclTypeExtractor<T1>::type>;
    };


  }  // namespace Internal

  /**
   * @brief Create a Validator that checks if a value is within a specified range.
   *
   * The range is defined by two values, @p low and @p high which can be either inclusive or
   * exclusive. By default, both values are inclusive. Whether a bound is inclusive or exclusive
   * can be specified using the incl() and excl() functions.
   *
   * @code
   * auto validator = in_range(incl(0), incl(10)); // value must be in [0, 10]
   * auto validator = in_range(0, 10); // same as above
   *
   * auto validator = in_range(excl(0.0), incl(1.23)); // value must be in (0.0, 1.23]
   *
   * // This will not compile because the types (int, double) are not identical:
   * // auto validator = in_range(0, 10.0);
   * @endcode
   *
   * @note You can only use this function for types @p int and @p double.
   *
   */
  template <typename Low, typename High,
      typename T = typename Internal::LowHighCommonType<Low, High>::type>
    requires((std::integral<T> && !std::same_as<T, bool>) || std::floating_point<T>)
  [[nodiscard]] auto in_range(const Low& low, const High& high);

  /**
   * Create a Validator that checks if an enum value is in a @p set of values.
   *
   * @note This only works for enum types to keep you from constructing input that is hard to
   * understand and use. Use the in_range() function for numeric types. If you find this comment
   * while you wanted to check a string value, this means you should really be using an enum type
   * with the enum values being the strings you wanted to check.
   */
  template <typename T>
    requires(std::is_enum_v<T>)
  [[nodiscard]] auto in_set(const std::set<T>& set);

  /**
   * Create a Validator that checks if an enum value is in a @p set of values.
   *
   * @note This only works for enum types to keep you from constructing input that is hard to
   * understand and use. Use the in_range() function for numeric types. If you find this comment
   * while you wanted to check a string value, this means you should really be using an enum type
   * with the enum values being the strings you wanted to check.
   */
  template <typename T>
    requires(std::is_enum_v<T>)
  [[nodiscard]] auto in_set(std::initializer_list<T> set);

  /**
   * Helper to mark a range value as inclusive in the in_range() function.
   *
   * @note This is not a Validator!
   *
   * @see in_range()
   */
  template <Internal::Numeric T>
  [[nodiscard]] auto incl(const T& value);

  /**
   * Helper to mark a range value as exclusive in the in_range() function.
   *
   * @note This is not a Validator!
   *
   * @see in_range()
   */
  template <Internal::Numeric T>
  [[nodiscard]] auto excl(const T& value);

  /**
   * A shorthand for creating a Validator that checks if a value is in the range (0, max].
   *
   * @note You need to specify the type of the value to be checked explicitly, since there are no
   * arguments to deduce the type from.
   */
  template <Internal::Numeric T>
  [[nodiscard]] auto positive();
}  // namespace Core::IO::InputSpecBuilders::Validators

// --- template definitions --- //


template <typename Low, typename High, typename T>
  requires((std::integral<T> && !std::same_as<T, bool>) || std::floating_point<T>)
[[nodiscard]] auto Core::IO::InputSpecBuilders::Validators::in_range(
    const Low& low, const High& high)
{
  const auto [low_val, low_type] = Internal::InclExclTag<T>(low);
  const auto [high_val, high_type] = Internal::InclExclTag<T>(high);
  FOUR_C_ASSERT_ALWAYS(low_val < high_val,
      "Invalid range: low value {} is not less than high value {}.", low_val, high_val);

  const auto description =
      std::format("in_range{}{},{}{}", (low_type == Internal::InclExclType::incl) ? "[" : "(",
          low_val, high_val, (high_type == Internal::InclExclType::incl) ? "]" : ")");


  return Validator<T>(
      // Rename the variables so Clang OpenMP does not complain about capturing structured bindings.
      [lv = low_val, hv = high_val, lt = low_type, ht = high_type](const T& v)
      {
        return ((lt == Internal::InclExclType::incl) ? (v >= lv) : (v > lv)) &&
               ((ht == Internal::InclExclType::incl) ? (v <= hv) : (v < hv));
      },
      [description](std::ostream& os) { os << description; });
}

template <Core::IO::InputSpecBuilders::Validators::Internal::Numeric T>
[[nodiscard]] auto Core::IO::InputSpecBuilders::Validators::incl(const T& value)
{
  return Internal::InclExclTag<T>{value, Internal::InclExclType::incl};
}

template <Core::IO::InputSpecBuilders::Validators::Internal::Numeric T>
[[nodiscard]] auto Core::IO::InputSpecBuilders::Validators::excl(const T& value)
{
  return Internal::InclExclTag<T>{value, Internal::InclExclType::excl};
}

template <Core::IO::InputSpecBuilders::Validators::Internal::Numeric T>
[[nodiscard]] auto Core::IO::InputSpecBuilders::Validators::positive()
{
  return in_range(excl(static_cast<T>(0)), std::numeric_limits<T>::max());
}

template <typename T>
  requires(std::is_enum_v<T>)
[[nodiscard]] auto Core::IO::InputSpecBuilders::Validators::in_set(const std::set<T>& set)
{
  FOUR_C_ASSERT_ALWAYS(!set.empty(), "You passed an empty set to in_set().");

  using namespace EnumTools;
  std::ostringstream set_description;
  set_description << "in_set{";
  for (const auto& val : set) set_description << val << ",";
  auto set_description_str = set_description.str();
  set_description_str.pop_back();
  set_description_str += "}";

  return Validator<T>([set](const T& v) { return set.contains(v); },
      [set_description_str](std::ostream& os) { os << set_description_str; });
}

template <typename T>
  requires(std::is_enum_v<T>)
[[nodiscard]] auto Core::IO::InputSpecBuilders::Validators::in_set(std::initializer_list<T> set)
{
  return in_set(std::set<T>(set.begin(), set.end()));
}


FOUR_C_NAMESPACE_CLOSE

#endif
