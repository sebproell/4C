// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_io_input_file.hpp"

#include "4C_io_input_file_utils.hpp"
#include "4C_unittest_utils_assertions_test.hpp"
#include "4C_unittest_utils_support_files_test.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_ParameterList.hpp>

namespace
{
  using namespace FourC;



  void check_section(
      Core::IO::InputFile& input, const std::string& section, const std::vector<std::string>& lines)
  {
    SCOPED_TRACE("Checking section " + section);
    ASSERT_TRUE(input.has_section(section));
    auto section_lines = input.in_section(section);
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

  TEST(InputFile, Test1)
  {
    const std::string input_file_name = TESTING::get_support_file_path("test_files/test1.dat");

    MPI_Comm comm(MPI_COMM_WORLD);
    Core::IO::InputFile input{input_file_name, comm};

    EXPECT_TRUE(input.has_section("EMPTY"));
    EXPECT_FALSE(input.has_section("NONEXISTENT SECTION"));

    check_section(input, "SECTION WITH SPACES", {"line in section with spaces"});
    check_section(input, "SECTION/WITH/SLASHES", {"line in section with slashes"});
    check_section(input, "SHORT SECTION", std::vector<std::string>(3, "line in short section"));
    check_section(input, "PARTICLES", std::vector<std::string>(30, "line in long section"));
  }

  TEST(InputFile, HasIncludes)
  {
    const std::string input_file_name =
        TESTING::get_support_file_path("test_files/has_includes/main.dat");

    MPI_Comm comm(MPI_COMM_WORLD);
    Core::IO::InputFile input{input_file_name, comm};

    EXPECT_EQ(input.file_for_section("INCLUDED SECTION 1a").filename(), "include1a.dat");
    EXPECT_EQ(input.file_for_section("SECTION 1").filename(), "main.dat");

    check_section(input, "INCLUDED SECTION 1a", std::vector<std::string>(2, "line"));
    check_section(input, "INCLUDED SECTION 1b", std::vector<std::string>(2, "line"));
    check_section(input, "INCLUDED SECTION 2", std::vector<std::string>(2, "line"));
    check_section(input, "INCLUDED SECTION 3", std::vector<std::string>(2, "line"));
    // Check that an on-the-fly section can be read from an include
    check_section(input, "PARTICLES", std::vector<std::string>(5, "line"));
  }

  TEST(InputFile, CyclicIncludes)
  {
    const std::string input_file_name =
        TESTING::get_support_file_path("test_files/cyclic_includes/cycle1.dat");

    MPI_Comm comm(MPI_COMM_WORLD);
    FOUR_C_EXPECT_THROW_WITH_MESSAGE(Core::IO::InputFile(input_file_name, comm), Core::Exception,
        "cycle1.dat' was already included before.");
  }

  TEST(InputFile, BasicYaml)
  {
    const std::string input_file_name = TESTING::get_support_file_path("test_files/yaml/basic.yml");

    MPI_Comm comm(MPI_COMM_WORLD);
    Core::IO::InputFile input{input_file_name, comm};

    EXPECT_FALSE(input.has_section("EMPTY"));
    EXPECT_FALSE(input.has_section("NONEXISTENT SECTION"));

    check_section(input, "SECTION WITH LINES", {"first line", "second line", "third line"});
  }

  TEST(InputFile, YamlIncludes)
  {
    const std::string input_file_name =
        TESTING::get_support_file_path("test_files/yaml_includes/main.yaml");

    MPI_Comm comm(MPI_COMM_WORLD);
    Core::IO::InputFile input{input_file_name, comm};

    check_section(input, "INCLUDED SECTION 1", std::vector<std::string>(2, "line"));
    check_section(input, "SECTION WITH SUBSTRUCTURE", {"MAT 1 THERMO COND 1 2 3 CAPA 2"});

    EXPECT_EQ(input.file_for_section("INCLUDED SECTION 2").filename(), "included.yaml");

    Teuchos::ParameterList pl;
    Core::IO::read_parameters_in_section(input, "INCLUDED SECTION 2", pl);
    EXPECT_EQ(pl.sublist("INCLUDED SECTION 2").get<int>("a"), 1);
    EXPECT_EQ(pl.sublist("INCLUDED SECTION 2").get<double>("b"), 2.0);
  }

}  // namespace