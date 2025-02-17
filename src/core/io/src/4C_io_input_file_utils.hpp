// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_IO_INPUT_FILE_UTILS_HPP
#define FOUR_C_IO_INPUT_FILE_UTILS_HPP

#include "4C_config.hpp"

#include "4C_io_input_parameter_container.hpp"
#include "4C_io_input_spec.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <functional>
#include <ostream>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;

  namespace Nurbs
  {
    class Knotvector;
  }
}  // namespace Core::FE

namespace Core::IO
{
  class InputFile;
  class InputSpec;
}  // namespace Core::IO

namespace Core::IO
{

  /**
   * Print a section header padded with dashes to 67 characters.
   */
  void print_section_header(std::ostream& out, const std::string& header);

  /**
   * Print the values of parameter entries for a dat file based on the provided parameter list.
   *
   * This function prints key-value pairs in the dat file format, including parameter names and
   * their values. It processes sublists and parameters, formatting them appropriately.
   *
   * @param stream      The output stream to which the information will be printed.
   * @param list        The parameter list containing the parameters and sublists to be printed.
   * @param comment     A flag indicating whether to print comments (default is true).
   */
  void print_dat(std::ostream& stream, const Teuchos::ParameterList& list, bool comment = true);


  void print_dat(std::ostream& stream, const std::map<std::string, Core::IO::InputSpec>& map);

  /**
   * Return true if the @p list contains any parameter that has whitespace in the key name.
   *
   * @note This is needed for the NOX parameters whose keywords and value have white spaces and
   * thus '=' are inserted to distinguish them.
   */
  bool need_to_print_equal_sign(const Teuchos::ParameterList& list);

  /**
   * Print @p spec into a dat file section with given @p header.
   */
  void print_section(std::ostream& out, const std::string& header, const InputSpec& spec);


  /**
   * Read all lines in a @p section of @p input that match the @p spec. Every line in the @p section
   * must match the @p spec. Otherwise, an exception is thrown.
   */
  std::vector<Core::IO::InputParameterContainer> read_all_lines_in_section(
      Core::IO::InputFile& input, const std::string& section, const InputSpec& spec);


  /**
   * Read only lines in a @p section of @p input that match the @p spec. This implies
   * that, potentially, no lines are read at all, resulting in an empty returned vector. In
   * addition to the vector of parsed lines, the second returned value contains all unparsed input
   * lines.
   *
   * @see read_all_lines_in_section()
   */
  std::pair<std::vector<Core::IO::InputParameterContainer>, std::vector<std::string>>
  read_matching_lines_in_section(
      Core::IO::InputFile& input, const std::string& section, const InputSpec& spec);

  /**
   * Split the given @p line into a key-value pair. Key and value are normally separated by
   * whitespace. In case there are multiple distinct whitespace groups in one line, the first of
   * these is assumed to be the separator and all the other whitespace is assumed to be part of
   * the value. Key and value may also be separated by an equals sign "=" and at least one
   * whitespace character on both sides. In this case, key and value may contain spaces
   * internally. Leading and trailing whitespace is trimmed from both key and value.
   *
   * @throws Core::Exception If the @p line cannot be read.
   *
   * @return A pair of key and value.
   */
  std::pair<std::string, std::string> read_key_value(const std::string& line);

  void read_parameters_in_section(
      InputFile& input, const std::string& section_name, Teuchos::ParameterList& list);

  /**
   * Read a node-design topology section. This is a collective call that propagates data that
   * may only be available on rank 0 to all ranks.
   *
   * @param input The input file.
   * @param name Name of the topology to read
   * @param dobj_fenode Resulting collection of all nodes that belong to a design.
   * @param get_discretization Callback to return a discretization by name.
   */
  void read_design(InputFile& input, const std::string& name,
      std::vector<std::vector<int>>& dobj_fenode,
      const std::function<const Core::FE::Discretization&(const std::string& name)>&
          get_discretization);

  /**
   * \brief read the knotvector section (for isogeometric analysis)
   *
   * \param  reader         (in ): InputFile object
   * \param  name           (in ): Name/type of discretisation
   * \param  disknots       (out): node vector coordinates
   *
   */
  void read_knots(InputFile& input, const std::string& name,
      std::shared_ptr<Core::FE::Nurbs::Knotvector>& disknots);

}  // namespace Core::IO


FOUR_C_NAMESPACE_CLOSE

#endif
