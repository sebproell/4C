// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_IO_INPUT_SPEC_HPP
#define FOUR_C_IO_INPUT_SPEC_HPP

#include "4C_config.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace Core::IO
{
  class InputParameterContainer;
  class ValueParser;

  namespace Internal
  {
    class InputSpecTypeErasedBase;
  }  // namespace Internal

  /**
   * Objects of this class encapsulate knowledge about the input. You can create objects using the
   * helper functions in the InputSpecBuilders namespace. See the function
   * InputSpecBuilders::entry() for more information on how to create InputSpecs.
   */
  class InputSpec
  {
   public:
    InputSpec() = default;

    ~InputSpec();

    InputSpec(std::unique_ptr<Internal::InputSpecTypeErasedBase> pimpl);

    InputSpec(const InputSpec& other);
    InputSpec& operator=(const InputSpec& other);

    InputSpec(InputSpec&&) noexcept = default;
    InputSpec& operator=(InputSpec&&) noexcept = default;

    /**
     * Use the @p parser to parse whatever this InputSpec expects. The results are stored in the
     * @p container. If parsing fails, an exception is thrown.
     */
    void fully_parse(ValueParser& parser, InputParameterContainer& container) const;

    /**
     * Print the expected input format of this InputSpec to @p stream in dat format. The @p
     * container is used to print the default values of the input. If a value is not present in the
     * container, the type of the value is printed instead.
     */
    void print_as_dat(std::ostream& stream, const InputParameterContainer& container) const;

    /**
     * Access the opaque implementation class. This is used in the implementation files where the
     * definition is known. There is nothing you can or should do with this function in the user
     * code.
     */
    Internal::InputSpecTypeErasedBase& impl();

    /**
     * Access the opaque implementation class. This is used in the implementation files where the
     * definition is known. There is nothing you can or should do with this function in the user
     * code.
     */
    [[nodiscard]] const Internal::InputSpecTypeErasedBase& impl() const;

   private:
    std::unique_ptr<Internal::InputSpecTypeErasedBase> pimpl_;
  };
}  // namespace Core::IO

FOUR_C_NAMESPACE_CLOSE

#endif
