// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_utils_function_manager.hpp"

#include "4C_io_input_file.hpp"
#include "4C_io_input_file_utils.hpp"
#include "4C_io_input_spec_builders.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_function.hpp"
#include "4C_utils_function_of_time.hpp"

FOUR_C_NAMESPACE_OPEN

namespace
{
  using InputParameters = std::vector<Core::IO::InputParameterContainer>;

  using TypeErasedFunctionCreator = std::function<std::any(const InputParameters&)>;

  template <typename T>
  using FunctionCreator = std::shared_ptr<T> (*)(const InputParameters&);

  /**
   * Utility function that takes a function object returning a std::shared_ptr<T> and erases its
   * return type via std::any. In addition, if the returned object would be nullptr, discard
   * it and return an empty std::any instead.
   */
  template <typename T>
  TypeErasedFunctionCreator wrap_function(FunctionCreator<T> fun)
  {
    return [fun](const InputParameters& linedefs) -> std::any
    {
      std::shared_ptr<T> created = fun(linedefs);
      if (created == nullptr)
        return {};
      else
        return created;
    };
  }


  std::any create_builtin_function(const std::vector<Core::IO::InputParameterContainer>& parameters)
  {
    // List all known TryCreate functions in a vector, so they can be called with a unified
    // syntax below. Also, erase their exact return type, since we can only store std::any.
    std::vector<TypeErasedFunctionCreator> try_create_function_vector{
        wrap_function(Core::Utils::try_create_symbolic_function_of_anything),
        wrap_function(Core::Utils::try_create_symbolic_function_of_space_time),
        wrap_function(Core::Utils::try_create_function_of_time)};

    for (const auto& try_create_function : try_create_function_vector)
    {
      auto maybe_function = try_create_function(parameters);
      if (maybe_function.has_value()) return maybe_function;
    }

    FOUR_C_THROW("Internal error: could not create a function that I should be able to create.");
  }

}  // namespace


void Core::Utils::add_valid_builtin_functions(Core::Utils::FunctionManager& function_manager)
{
  using namespace IO::InputSpecBuilders;

  auto time_info = all_of({
      entry<int>("NUMPOINTS"),
      one_of({
          group("BYNUM",
              {
                  entry<std::vector<double>>("TIMERANGE", {.size = 2}),
              },
              {.description = "Linearly distribute NUMPOINTS time points in the TIMERANGE."}),
          entry<std::vector<double>>("TIMES", {.size = from_parameter<int>("NUMPOINTS")}),
      }),
  });

  auto spec = one_of({
      all_of({
          entry<int>("COMPONENT", {.required = false}),
          entry<std::string>("SYMBOLIC_FUNCTION_OF_SPACE_TIME"),
      }),

      entry<std::string>("SYMBOLIC_FUNCTION_OF_TIME"),

      all_of({
          entry<int>("VARIABLE"),
          entry<std::string>("NAME"),
          one_of({
              all_of({
                  selection<std::string>("TYPE", {"expression"}),
                  entry<std::string>("DESCRIPTION"),
              }),
              all_of({
                  selection<std::string>("TYPE", {"linearinterpolation", "fourierinterpolation"}),
                  time_info,
                  entry<std::vector<double>>("VALUES", {.size = from_parameter<int>("NUMPOINTS")}),
              }),
              all_of({
                  selection<std::string>("TYPE", {"multifunction"}),
                  time_info,
                  entry<std::vector<std::string>>(
                      "DESCRIPTION", {.size = [](const IO::InputParameterContainer& container)
                                         { return container.get<int>("NUMPOINTS") - 1; }}),
              }),
          }),
          group("PERIODIC",
              {
                  entry<double>("T1"),
                  entry<double>("T2"),
              },
              {.required = false}),
      }),

      all_of({
          entry<std::string>("VARFUNCTION"),
          entry<int>("NUMCONSTANTS", {.required = false}),
          entry<std::vector<std::pair<std::string, double>>>(
              "CONSTANTS", {.required = false, .size = from_parameter<int>("NUMCONSTANTS")}),
      }),
  });

  function_manager.add_function_definition(spec, create_builtin_function);
}


Core::IO::InputSpec Core::Utils::FunctionManager::valid_function_lines()
{
  std::vector<IO::InputSpec> specs;
  for (const auto& [spec, _] : attached_function_data_)
  {
    specs.emplace_back(spec);
  }
  return Core::IO::InputSpecBuilders::one_of(specs);
}


void Core::Utils::FunctionManager::add_function_definition(
    IO::InputSpec spec, FunctionFactory function_factory)
{
  attached_function_data_.emplace_back(std::move(spec), std::move(function_factory));
}


void Core::Utils::FunctionManager::read_input(Core::IO::InputFile& input)
{
  functions_.clear();

  // Read FUNCT sections starting from FUNCT1 until the first empty one is encountered.
  // This implies that the FUNCT sections must form a contiguous range in the input file.
  // Otherwise, the read fails later.
  for (int funct_suffix = 1;; ++funct_suffix)
  {
    const bool stop_parsing = std::invoke(
        [&]()
        {
          for (auto& [spec, function_factory] : attached_function_data_)
          {
            auto [parsed_parameters, unparsed_lines] = Core::IO::read_matching_lines_in_section(
                input, "FUNCT" + std::to_string(funct_suffix), spec);

            // A convoluted way of saying that there are no lines in the section, thus, stop
            // parsing. This can only be refactored if the reading mechanism is overhauled in
            // general.
            if (parsed_parameters.size() + unparsed_lines.size() == 0)
            {
              return true;
            }

            if (parsed_parameters.size() > 0 && unparsed_lines.size() == 0)
            {
              functions_.emplace_back(function_factory(parsed_parameters));
              return false;
            }
          }

          // If we end up here, the current sections function definition could not be parsed.
          {
            std::stringstream ss;
            for (const auto& line : input.in_section("FUNCT" + std::to_string(funct_suffix)))
            {
              ss << '\n' << line.get_as_dat_style_string();
            }

            FOUR_C_THROW("Could not parse the following lines into a Function known to 4C:\n%s",
                ss.str().c_str());
          }
        });

    // Stop reading as soon as the first FUNCT section in the input file is empty
    if (stop_parsing) break;
  }
}

FOUR_C_NAMESPACE_CLOSE
