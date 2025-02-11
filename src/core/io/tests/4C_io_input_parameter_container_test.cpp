// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_io_input_parameter_container.hpp"

namespace
{
  using namespace FourC;

  TEST(InputParameterContainerTest, ToTeuchosParameterList)
  {
    Core::IO::InputParameterContainer container;
    container.add<int>("a", 1);
    container.add<double>("b", 2.0);
    container.add<std::string>("c", "string");
    container.add<bool>("d", true);
    container.add("v", std::vector<int>{1, 2, 3});
    container.add("n", Core::IO::Noneable<int>{});
    container.group("group").add<std::string>("e", "group string");
    container.group("group").group("deeply").group("nested").add<int>("f", 42);

    Teuchos::ParameterList list("dummy");
    container.to_teuchos_parameter_list(list);
    EXPECT_EQ(list.name(), "dummy");
    EXPECT_EQ(list.get<int>("a"), 1);
    EXPECT_EQ(list.get<double>("b"), 2.0);
    EXPECT_EQ(list.get<std::string>("c"), "string");
    EXPECT_EQ(list.get<bool>("d"), true);
    EXPECT_EQ(list.get<std::vector<int>>("v"), (std::vector<int>{1, 2, 3}));
    EXPECT_EQ(list.get<Core::IO::Noneable<int>>("n"), Core::IO::none<int>);
    EXPECT_EQ(list.sublist("group").get<std::string>("e"), "group string");
    EXPECT_EQ(list.sublist("group").sublist("deeply").sublist("nested").get<int>("f"), 42);
  }
}  // namespace