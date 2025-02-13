// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_io_input_spec.hpp"

#include "4C_io_input_spec_builders.hpp"
#include "4C_io_yaml.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

Core::IO::InputSpec::InputSpec(std::unique_ptr<Internal::InputSpecTypeErasedBase> pimpl)
    : pimpl_(std::move(pimpl))
{
}

Core::IO::InputSpec::~InputSpec() = default;

Core::IO::InputSpec::InputSpec(const InputSpec& other)
    : pimpl_(other.pimpl_ ? other.pimpl_->clone() : nullptr)
{
}

Core::IO::InputSpec& Core::IO::InputSpec::operator=(const InputSpec& other)
{
  if (this == &other) return *this;

  pimpl_ = other.pimpl_ ? other.pimpl_->clone() : nullptr;

  return *this;
}

void Core::IO::InputSpec::fully_parse(
    ValueParser& parser, Core::IO::InputParameterContainer& container) const
{
  FOUR_C_ASSERT(pimpl_, "InputSpec is empty.");

  pimpl_->parse(parser, container);
  if (!parser.at_end())
  {
    std::stringstream ss;
    container.print(ss);
    std::string remainder(parser.get_unparsed_remainder());
    FOUR_C_THROW("After parsing, the line still contains '%s'.\nParsed parameters: %s",
        remainder.c_str(), ss.str().c_str());
  }
}

void Core::IO::InputSpec::match(ConstYamlNodeRef yaml, InputParameterContainer& container) const
{
  FOUR_C_ASSERT(pimpl_, "InputSpec is empty.");

  Internal::MatchTree match_tree{*this, yaml};
  pimpl_->match(yaml, container, match_tree.root());

  match_tree.assert_match();
}

void Core::IO::InputSpec::print_as_dat(std::ostream& stream) const
{
  FOUR_C_ASSERT(pimpl_, "InputSpec is empty.");

  pimpl_->print(stream, 0u);
}

void Core::IO::InputSpec::emit_metadata(YamlNodeRef yaml) const
{
  FOUR_C_ASSERT(pimpl_, "InputSpec is empty.");

  auto root = yaml.node;
  pimpl_->emit_metadata(root);
}

Core::IO::Internal::InputSpecTypeErasedBase& Core::IO::InputSpec::impl()
{
  FOUR_C_ASSERT(pimpl_, "InputSpec is empty.");
  return *pimpl_;
}

const Core::IO::Internal::InputSpecTypeErasedBase& Core::IO::InputSpec::impl() const
{
  FOUR_C_ASSERT(pimpl_, "InputSpec is empty.");
  return *pimpl_;
}

FOUR_C_NAMESPACE_CLOSE
