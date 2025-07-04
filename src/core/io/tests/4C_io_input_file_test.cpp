// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_io_input_file.hpp"

#include "4C_io_input_file_utils.hpp"
#include "4C_io_input_spec_builders.hpp"
#include "4C_unittest_utils_assertions_test.hpp"
#include "4C_unittest_utils_support_files_test.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_ParameterList.hpp>

namespace
{
  using namespace FourC;

  void check_section_rank_0_only(
      Core::IO::InputFile& input, const std::string& section, const std::vector<std::string>& lines)
  {
    SCOPED_TRACE("Checking section " + section);
    ASSERT_TRUE(input.has_section(section));
    auto section_lines = input.in_section_rank_0_only(section);
    std::vector<std::string> section_lines_str;
    std::ranges::copy(section_lines | std::views::transform([](const auto& line)
                                          { return std::string{line.get_as_dat_style_string()}; }),
        std::back_inserter(section_lines_str));
    EXPECT_EQ(lines.size(), section_lines_str.size());
    for (std::size_t i = 0; i < lines.size(); ++i)
    {
      EXPECT_EQ(section_lines_str[i], lines[i]);
    }
  }

  TEST(InputFile, BasicYaml)
  {
    const std::string input_file_name = TESTING::get_support_file_path("test_files/yaml/basic.yml");

    MPI_Comm comm(MPI_COMM_WORLD);
    Core::IO::InputFile input{{},
        {
            "EMPTY",
            "SECTION1",
            "SECTION2",
            "SECTION WITH LINES",
        },
        comm};
    input.read(input_file_name);

    EXPECT_FALSE(input.has_section("EMPTY"));
    EXPECT_FALSE(input.has_section("NONEXISTENT SECTION"));

    check_section_rank_0_only(
        input, "SECTION WITH LINES", {"first line", "second line", "third line"});
  }

  TEST(InputFile, YamlIncludes)
  {
    using namespace Core::IO::InputSpecBuilders;
    const std::string input_file_name =
        TESTING::get_support_file_path("test_files/yaml_includes/main.yaml");

    MPI_Comm comm(MPI_COMM_WORLD);
    Core::IO::InputFile input(
        {{"INCLUDED SECTION 2", group("INCLUDED SECTION 2",
                                    {
                                        parameter<int>("a"),
                                        parameter<double>("b"),
                                        parameter<bool>("c"),
                                    })},
            {"SECTION WITH SUBSTRUCTURE", list("SECTION WITH SUBSTRUCTURE",
                                              all_of({
                                                  parameter<int>("MAT"),
                                                  group("THERMO",
                                                      {
                                                          parameter<std::vector<double>>("COND"),
                                                          parameter<double>("CAPA"),
                                                      }),
                                              }))}},
        {
            "MAIN SECTION",
        },
        comm);
    input.read(input_file_name);

    EXPECT_EQ(input.file_for_section("INCLUDED SECTION 2").filename(), "included.yaml");

    Core::IO::InputParameterContainer container;
    input.match_section("INCLUDED SECTION 2", container);
    const auto& included_section_2 = container.group("INCLUDED SECTION 2");
    EXPECT_EQ(included_section_2.get<int>("a"), 1);
    EXPECT_EQ(included_section_2.get<double>("b"), 2.0);
    EXPECT_EQ(included_section_2.get<bool>("c"), true);

    container.clear();
    input.match_section("SECTION WITH SUBSTRUCTURE", container);
    EXPECT_EQ(container.get_list("SECTION WITH SUBSTRUCTURE")[0].get<int>("MAT"), 1);
  }

}  // namespace