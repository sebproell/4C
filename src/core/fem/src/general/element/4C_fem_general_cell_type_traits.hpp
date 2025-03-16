// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FEM_GENERAL_CELL_TYPE_TRAITS_HPP
#define FOUR_C_FEM_GENERAL_CELL_TYPE_TRAITS_HPP

#include "4C_config.hpp"

#include "4C_fem_general_cell_type.hpp"
#include "4C_fem_general_cell_type_traits_def.hpp"
#include "4C_fem_general_cell_type_traits_details.hpp"
#include "4C_utils_exceptions.hpp"

#include <algorithm>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <utility>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  //! @name Type traits for cell shape based on cell types
  /// @{
  template <Core::FE::CellType celltype>
  inline static constexpr bool is_tet = Internal::is_tet<celltype>::value;

  template <Core::FE::CellType celltype>
  inline static constexpr bool is_hex = Internal::is_hex<celltype>::value;

  template <Core::FE::CellType celltype>
  inline static constexpr bool is_nurbs = Internal::is_nurbs<celltype>::value;

  template <Core::FE::CellType celltype>
  inline static constexpr bool is_wedge = Internal::is_wedge<celltype>::value;

  template <Core::FE::CellType celltype>
  inline static constexpr bool is_pyramid = Internal::is_pyramid<celltype>::value;

  template <Core::FE::CellType celltype>
  inline static constexpr bool is_quad = Internal::is_quad<celltype>::value;

  template <Core::FE::CellType celltype>
  inline static constexpr bool is_tri = Internal::is_tri<celltype>::value;

  template <Core::FE::CellType celltype>
  inline static constexpr bool is_line = Internal::is_line<celltype>::value;

  template <Core::FE::CellType celltype>
  inline static constexpr bool use_lagrange_shapefnct =
      Internal::use_lagrange_shapefnct<celltype>::value;
  /// @}


  //! @name Compile time lists of types with cell type as template argument
  ///@{
  template <typename... Ts>
  using BaseTypeList = std::tuple<Ts...>;

  /*!
   * @brief Join two compile time base type lists
   *
   * Example:
   * @code{.cpp}
   * using BaseHex = BaseTypeList<Base<CellType::hex8>, Base<CellType::hex27>>
   * using BaseTet = BaseTypeList<Base<CellType::tet4>, Base<CellType::tet10>>
   *
   * // result in BaseTypeList<Base<CellType::hex8>, Base<CellType::hex27>, Base<CellType::tet4>,
   * //   Base<CellType::tet10>>
   * using BaseHexAndTet = Join<BaseHex, BaseTet>;
   * @endcode
   *
   * @tparam t
   */
  template <typename... T>
  using Join = typename Internal::Join<T...>::type;
  ///@}

  /*!
   * @brief Apply cell type list to a class that takes the cell type as template parameter and
   * return a compile time list with the base type for each cell type in the cell type list.
   *
   * Suppose you have a template class @c Base<CellType> and you want to generate a compile time
   * list for different cell types. You can do it like so:
   *
   * @code{.cpp}
   * using CellTypes = Core::FE::celltype_sequence<CellType, CellType::hex8, CellType::hex27>
   *
   * // results in BaseTypeList<Base<CellType::hex8>, Base<CellType::hex27>>
   * using BaseTypeList = apply_celltype_sequence<Base, CellTypes>;
   *
   * @endcode
   *
   * You could use this BaseTypeList to create a variant of a the template class @c Base for each
   * cell type
   *
   * @code{.cpp}
   * namespace Internal
   * {
   *   template <typename Tuple>
   *   struct CreateVariant;
   *
   *   template <typename... Ts>
   *   struct CreateVariant<BaseTypeList<Ts...>>
   *   {
   *     using type = std::variant<Ts...>;
   *   };
   * }
   *
   * template <typename... Ts>
   * using CreateVariantType = typename Internal::CreateVariant<Ts...>::type;
   *
   * // results in std::variant<Base<CellType::hex8>, Base<CellType::hex27>>
   * using Variant = CreateVariantType<BaseTypeList>;
   * @endcode
   *
   * @tparam Base
   * @tparam CellTypeList
   */
  template <template <Core::FE::CellType> typename Base, typename CellTypeSequence>
  using apply_celltype_sequence =
      typename Internal::apply_celltype_sequence<BaseTypeList, Base, CellTypeSequence>::type;

  /*!
   * @brief Returns a std::array of celltypes defined in the given integer sequence
   *
   * @tparam celltypes Integer sequence of cell types
   */
  template <typename Celltypes>
  static constexpr std::array celltype_array =
      Internal::celltype_sequence_to_array<Celltypes>::value;

  /*!
   * @brief a Core::FE::celltype_sequence holding all @p CellTypes including @p none and @p max
   */
  using all_celltypes = Internal::make_celltype_sequence<CellType::dis_none, CellType::max_distype>;
  /*!
   * @brief a Core::FE::celltype_sequence holding all physical @p CellTypes , i.e. all CellTypes
   * except @p none and @p max
   */
  using all_physical_celltypes = Internal::make_celltype_sequence<static_cast<CellType>(1),
      static_cast<CellType>(static_cast<std::size_t>(CellType::max_distype) - 1)>;
  ///@}

  /*!
   * @brief Call function with compile time celltype template argument with runtime celltype.
   *
   * Suppose you want to have a compile time celltype from a runtime celltype. Typically, you
   * would create a switch statement with all supported cell types. This function generates the
   * switch-case statements at compile time and inserts the respective template parameter.
   *
   * @code{.cpp}
   * auto result = CellTypeSwitch(celltype, [](auto celltype_t) {
   *   // celltype_t is an std::integral_constant<CellType, celltype>
   *   // you could use it as a template parameter of a function
   *   return do_sth<celltype_t()>();
   * });
   * @endcode
   *
   * @note by default, only the cases for all physical celltypes are defined, i.e. all celltypes
   * except @p none and @p max . If you want to define the result for all, you have to specify @p
   * all_celltypes as template parameter:
   *
   * @code{.cpp}
   * auto result = CellTypeSwitch<all_celltypes>(celltype, [](auto celltype_t) {
   *   // do_sth<...>(); must be defined for all celltypes including none and max in order to
   *   // compile
   *   return do_sth<celltype_t()>();
   * });
   * @endcode
   *
   * You can also pass a sequence of supported celltypes
   *
   * @code{.cpp}
   * using hex_celltypes = Core::FE::celltype_sequence<CellType::hex8, CellType::hex16,
   * CellType::hex18, CellType::hex20, CellType::hex27>; auto result =
   * CellTypeSwitch<hex_celltypes>(celltype, [](auto celltype_t) {
   *   // do_sth<...>(); must be defined for all celltypes listed in hex_celltypes
   *   return do_sth<celltype_t()>();
   * });
   * @endcode
   *
   * Any of the methods above also support passing a callable with the default operation. Per
   * default, an error message is thrown with the calling celltype and the supported celltypes. To
   * override the default operation, you have to provide a callable object (e.g. a lambda function)
   * that either returns a convertible type (for some default operation) or throws an error.
   *
   * In order to throw a custom error message, do it like so:
   * @code{.cpp}
   * using hex_celltypes = Core::FE::celltype_sequence<CellType::hex8, CellType::hex16,
   * CellType::hex18, CellType::hex20, CellType::hex27>; auto result =
   * CellTypeSwitch<hex_celltypes>(celltype, [](auto celltype_t) {
   *   // do_sth<...>(); must be defined for all celltypes listed in hex_celltypes
   *   return do_sth<celltype_t()>();
   * }, [](auto celltype_t) {
   *   FOUR_C_THROW("unsupported celltype"); // new error message.
   * });
   * @endcode
   *
   * In order to do a default operation, do it like so:
   * @code{.cpp}
   * using hex_celltypes = Core::FE::celltype_sequence<CellType::hex8, CellType::hex16,
   * CellType::hex18, CellType::hex20, CellType::hex27>; auto result =
   * CellTypeSwitch<hex_celltypes>(celltype, [](auto celltype_t) {
   *   // do_sth<...>(); must be defined for all celltypes listed in hex_celltypes
   *   return do_sth<celltype_t()>();
   * }, [](auto celltype_t) {
   *   return do_sth_default<celltype_t()>();
   * });
   * @endcode
   *
   * @note the callable @p fct does not have to return anything. If it returns @p void , @p
   * CellTypeSwitch will also return @p void .
   *
   * @param celltype (in) : the runtime celltype
   * @param fct (in) : A callable taking an integral_constant with the celltype as argument. The
   * compile time celltype can be obtained with the @p operator() . Typically, a lambda function
   * with @p auto parameter is passed.
   * @param unsupported_celltype_callable (in) : A callable taking an integral_constant with the
   * celltype as argument for the default operation. The compile time celltype can be obtained with
   * the @p operator() . Typically, a lambda function with @p auto parameter is passed. either throw
   * an error or return a convertible type to @p fct return type.
   * @tparam celltype_sequence : a Core::FE::celltype_sequence holding a list of all implemented
   * celltypes (default: all physical celltypes, i.e. all except none and max)
   * @return returns the invoked result of @p fct
   */
  template <typename CelltypeSequence = all_physical_celltypes, typename Function,
      typename UnsupportedCellTypeCallable =
          Internal::ThrowUnsupportedCellTypeError<CelltypeSequence>>
  auto cell_type_switch(CellType celltype, Function fct,
      UnsupportedCellTypeCallable unsupported_celltype_callable = {})
  {
    switch (celltype)
    {
      case CellType::dis_none:
        return Internal::cell_type_switch_item<CellType::dis_none, CelltypeSequence>(
            fct, unsupported_celltype_callable);
      case CellType::quad4:
        return Internal::cell_type_switch_item<CellType::quad4, CelltypeSequence>(
            fct, unsupported_celltype_callable);
      case CellType::quad6:
        return Internal::cell_type_switch_item<CellType::quad6, CelltypeSequence>(
            fct, unsupported_celltype_callable);
      case CellType::quad8:
        return Internal::cell_type_switch_item<CellType::quad8, CelltypeSequence>(
            fct, unsupported_celltype_callable);
      case CellType::quad9:
        return Internal::cell_type_switch_item<CellType::quad9, CelltypeSequence>(
            fct, unsupported_celltype_callable);
      case CellType::tri3:
        return Internal::cell_type_switch_item<CellType::tri3, CelltypeSequence>(
            fct, unsupported_celltype_callable);
      case CellType::tri6:
        return Internal::cell_type_switch_item<CellType::tri6, CelltypeSequence>(
            fct, unsupported_celltype_callable);
      case CellType::hex8:
        return Internal::cell_type_switch_item<CellType::hex8, CelltypeSequence>(
            fct, unsupported_celltype_callable);
      case CellType::hex16:
        return Internal::cell_type_switch_item<CellType::hex16, CelltypeSequence>(
            fct, unsupported_celltype_callable);
      case CellType::hex18:
        return Internal::cell_type_switch_item<CellType::hex18, CelltypeSequence>(
            fct, unsupported_celltype_callable);
      case CellType::hex20:
        return Internal::cell_type_switch_item<CellType::hex20, CelltypeSequence>(
            fct, unsupported_celltype_callable);
      case CellType::hex27:
        return Internal::cell_type_switch_item<CellType::hex27, CelltypeSequence>(
            fct, unsupported_celltype_callable);
      case CellType::tet4:
        return Internal::cell_type_switch_item<CellType::tet4, CelltypeSequence>(
            fct, unsupported_celltype_callable);
      case CellType::tet10:
        return Internal::cell_type_switch_item<CellType::tet10, CelltypeSequence>(
            fct, unsupported_celltype_callable);
      case CellType::wedge6:
        return Internal::cell_type_switch_item<CellType::wedge6, CelltypeSequence>(
            fct, unsupported_celltype_callable);
      case CellType::wedge15:
        return Internal::cell_type_switch_item<CellType::wedge15, CelltypeSequence>(
            fct, unsupported_celltype_callable);
      case CellType::pyramid5:
        return Internal::cell_type_switch_item<CellType::pyramid5, CelltypeSequence>(
            fct, unsupported_celltype_callable);
      case CellType::line2:
        return Internal::cell_type_switch_item<CellType::line2, CelltypeSequence>(
            fct, unsupported_celltype_callable);
      case CellType::line3:
        return Internal::cell_type_switch_item<CellType::line3, CelltypeSequence>(
            fct, unsupported_celltype_callable);
      case CellType::line4:
        return Internal::cell_type_switch_item<CellType::line4, CelltypeSequence>(
            fct, unsupported_celltype_callable);
      case CellType::line5:
        return Internal::cell_type_switch_item<CellType::line5, CelltypeSequence>(
            fct, unsupported_celltype_callable);
      case CellType::line6:
        return Internal::cell_type_switch_item<CellType::line6, CelltypeSequence>(
            fct, unsupported_celltype_callable);
      case CellType::point1:
        return Internal::cell_type_switch_item<CellType::point1, CelltypeSequence>(
            fct, unsupported_celltype_callable);
      case CellType::nurbs2:
        return Internal::cell_type_switch_item<CellType::nurbs2, CelltypeSequence>(
            fct, unsupported_celltype_callable);
      case CellType::nurbs3:
        return Internal::cell_type_switch_item<CellType::nurbs3, CelltypeSequence>(
            fct, unsupported_celltype_callable);
      case CellType::nurbs4:
        return Internal::cell_type_switch_item<CellType::nurbs4, CelltypeSequence>(
            fct, unsupported_celltype_callable);
      case CellType::nurbs9:
        return Internal::cell_type_switch_item<CellType::nurbs9, CelltypeSequence>(
            fct, unsupported_celltype_callable);
      case CellType::nurbs8:
        return Internal::cell_type_switch_item<CellType::nurbs8, CelltypeSequence>(
            fct, unsupported_celltype_callable);
      case CellType::nurbs27:
        return Internal::cell_type_switch_item<CellType::nurbs27, CelltypeSequence>(
            fct, unsupported_celltype_callable);
      case CellType::max_distype:
        return Internal::cell_type_switch_item<CellType::max_distype, CelltypeSequence>(
            fct, unsupported_celltype_callable);
    }
    FOUR_C_THROW(
        "You are calling this celltype_switch with an unknown cell type. Please add"
        "your cell type to this list above so that the function works as expected.");
  }

  /*!
   * @brief Compile time string (precisely const char * const) representation of a celltype
   *
   * @tparam celltype
   */
  template <CellType celltype>
  static constexpr auto celltype_string = Internal::CellTypeInformation<celltype>::name;

  /*!
   * @brief Returns a string representation of the celltype
   *
   * @param celltype
   * @return std::string
   */
  inline std::string cell_type_to_string(CellType celltype)
  {
    return cell_type_switch<all_celltypes>(
        celltype, [](auto celltype_t) { return celltype_string<celltype_t()>; });
  }

  inline CellType string_to_cell_type(const std::string& celltype_str)
  {
    constexpr std::array all_names =
        Internal::celltype_sequence_to_string_array<all_celltypes>::value;

    const auto* found_position =
        std::find(std::begin(all_names), std::end(all_names), celltype_str);

    if (found_position == std::end(all_names))
    {
      FOUR_C_THROW("Unknown celltype {}", celltype_str);
    }

    constexpr std::array all_celltypes_enum = celltype_array<all_celltypes>;
    return all_celltypes_enum[std::distance(std::begin(all_names), found_position)];
  }


  /*!
   * @brief Compile time mapping from celltype to dimension of the element
   *
   * @tparam celltype
   */
  template <CellType celltype>
  static constexpr int dim = Internal::CellTypeInformation<celltype>::dim;

  /*!
   * @brief Compile time mapping from celltype to number of nodes of the element
   *
   * @tparam celltype
   */
  template <CellType celltype>
  static constexpr int num_nodes = Internal::CellTypeInformation<celltype>::num_nodes;

  /*!
   * @brief Compile time mapping from celltype to number of element faces
   *
   * @tparam celltype
   */
  template <CellType celltype>
  static constexpr int num_faces = Internal::CellTypeInformation<celltype>::num_faces;

  /*!
   * @brief Compile time mapping from celltype to order of the shape functions
   *
   * @tparam celltype
   */
  template <CellType celltype>
  static constexpr CellType order = Internal::CellTypeInformation<celltype>::order;

  /*!
   * @brief Check whether the celltype is a nurbs cell type
   *
   * @param celltype
   * @return true
   * @return false
   */
  inline bool is_nurbs_celltype(CellType celltype)
  {
    return cell_type_switch<all_physical_celltypes>(
        celltype, [&](auto celltype_t) { return is_nurbs<celltype_t()>; });
  }
}  // namespace Core::FE

FOUR_C_NAMESPACE_CLOSE

#endif