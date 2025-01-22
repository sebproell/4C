// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_IO_INPUT_FILE_HPP
#define FOUR_C_IO_INPUT_FILE_HPP

#include "4C_config.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_comm_pack_buffer.hpp"

#include <filesystem>
#include <list>
#include <map>
#include <memory>
#include <ranges>
#include <set>
#include <string>
#include <variant>

FOUR_C_NAMESPACE_OPEN

namespace Core::IO
{
  /**
   * This class encapsulates input files independent of their format.
   *
   * Objects of this class read the content of a file and grant access to it. The input is not
   * yet interpreted in any way. Input files contain different sections. A section either contains
   * key-value pairs or a list of arbitrary lines of text. Sections may appear in an arbitrary order
   * and each section name must be unique. An exception is the special section named "INCLUDES"
   * which can contain a list of other files that should be read in addition to the current file.
   *
   * Three file formats are supported: the custom .dat file format and the standard .yaml (or .yml)
   * and .json formats. The format of a file is detected based on its ending. If the ending is not
   * one of the above mentioned, the file is assumed to be in the .dat format.
   * Included files do not have to use the same format as the including file.
   *
   * The following example shows the structure of a .dat file:
   *
   * @code
   * // A comment. Sections start with at least two dashes in the .dat file format.
   * --INCLUDES
   * include1.yml
   * include2.json
   * // More dashes are allowed to start a section.
   * -------------SECTION1
   * key1 = value1
   * key2 = value2
   * // The "=" is optional and can be replaced by whitespace. This only works if key and value do
   * // not contain whitespace themselves.
   * key3 value3
   * --SECTION2
   * A line with content that is not a key-value pair. It can contain anything.
   * @endcode
   *
   * A similar example looks like this in .yaml format:
   *
   * @code
   * INCLUDES:
   *   - include1.yml
   *   - include2.json
   * # A comment. Note that sections do NOT start with dashes in the .yaml file format.
   * SECTION1:
   *   key1: value1
   *   key2: value2
   * SECTION2:
   *   - A line with content that is not a key-value pair.
   *   - It can contain anything.
   * @endcode
   *
   *
   * @note The file is only read on rank 0 to save memory. Sections that are huge are only
   * distributed to other ranks if they are accessed through line_in_section(). If you only
   * want to read a section on rank 0, use lines_in_section_rank_0_only().
   */
  class InputFile
  {
   public:
    /// Construct a reader for a given file
    InputFile(std::string filename, MPI_Comm comm);

    /**
     * Get the (absolute) file path of the input file that contained a section. If the section is
     * unknown or was not read from any file, an empty path is returned.
     */
    [[nodiscard]] std::filesystem::path file_for_section(const std::string& section_name) const;

    /**
     * Get a a range of lines inside a section that have actual content, i.e., they contain
     * something other than whitespace or comments. Any line returned will have comment stripped
     * and whitespace trimmed. The usual way to do something with the lines is
     *
     * @code
     *   for (const auto& line : input.lines_in_section("section_name"))
     *   {
     *     // do something with line
     *   }
     * @endcode
     *
     * @return A range of string_views to the lines in this section.
     *
     * @note This is a collective call that needs to be called on all MPI ranks in the communicator
     * associated with this object. Depending on the section size, the content might need to be
     * distributed from rank 0 to all other ranks. This happens automatically.
     */
    std::ranges::view auto lines_in_section(const std::string& section_name);

    /**
     * This function is similar to lines_in_section(), but it only returns the lines on rank 0 and
     * returns an empty range on all other ranks. This is useful for sections that might be huge and
     * are not necessary on all ranks.
     */
    std::ranges::view auto lines_in_section_rank_0_only(const std::string& section_name);

    /**
     * Returns true if the input file contains a section with the given name.
     *
     * @note This is a collective call that needs to be called on all MPI ranks in the communicator.
     */
    [[nodiscard]] bool has_section(const std::string& section_name) const;

    /**
     * Access MPI communicator associated with this object.
     */
    [[nodiscard]] MPI_Comm get_comm() const { return comm_; }

    /**
     * Print a list of all sections that are contained in the input file but never
     * accessed through this object.
     *
     * @return True if there were unknown sections, false otherwise.
     */
    bool print_unknown_sections(std::ostream& out) const;

    /**
     * Internal storage for the content of a section.
     */
    struct SectionContent
    {
      //! The raw chars of the input.
      std::vector<char> raw_content;
      //! String views into #raw_content.
      std::vector<std::string_view> lines;
      //! The file the section was read from.
      std::string file;

      void pack(Core::Communication::PackBuffer& data) const;

      void unpack(Core::Communication::UnpackBuffer& buffer);

      SectionContent() = default;

      // Prevent copies: making a copy would invalidate the string views.
      SectionContent(const SectionContent&) = delete;
      SectionContent& operator=(const SectionContent&) = delete;

      SectionContent(SectionContent&&) = default;
      SectionContent& operator=(SectionContent&&) = default;
    };

   private:
    /**
     * The shared part of reading a file and postprocessing its content.
     */
    void read_generic(const std::filesystem::path& top_level_file);

    //! Remember that a section was used.
    void record_section_used(const std::string& section_name);

    /// The communicator associated with this object.
    MPI_Comm comm_;

    std::unordered_map<std::string, SectionContent> content_by_section_;

    /// Protocol of known and unknown section names
    std::map<std::string, bool> knownsections_;
  };


  /// -- template and inline functions --- //

  inline std::ranges::view auto InputFile::lines_in_section(const std::string& section_name)
  {
    record_section_used(section_name);

    static const std::vector<std::string_view> empty;

    const bool known_somewhere = has_section(section_name);
    if (!known_somewhere)
    {
      return std::views::all(empty);
    }

    const bool locally_known = content_by_section_.contains(section_name);
    const bool known_everywhere = Core::Communication::all_reduce<bool>(
        locally_known, [](const bool& r, const bool& in) { return r && in; }, comm_);
    if (known_everywhere)
    {
      // Take a const reference to the section content.
      const auto& lines = content_by_section_.at(section_name).lines;
      return std::views::all(lines);
    }

    // Distribute the content of the section to all ranks.
    {
      FOUR_C_ASSERT((!locally_known && (Core::Communication::my_mpi_rank(comm_) > 0)) ||
                        (locally_known && (Core::Communication::my_mpi_rank(comm_) == 0)),
          "Implementation error: section should be known on rank 0 and unknown on others.");

      auto& content = content_by_section_[section_name];
      Core::Communication::broadcast(content, 0, comm_);
    }

    const auto& lines = content_by_section_.at(section_name).lines;
    return std::views::all(lines);
  }

  inline std::ranges::view auto InputFile::lines_in_section_rank_0_only(
      const std::string& section_name)
  {
    if (Core::Communication::my_mpi_rank(comm_) == 0 && content_by_section_.contains(section_name))
    {
      record_section_used(section_name);
      const auto& lines = content_by_section_.at(section_name).lines;
      return std::views::all(lines);
    }
    else
    {
      static const std::vector<std::string_view> empty;
      return std::views::all(empty);
    }
  }

}  // namespace Core::IO

FOUR_C_NAMESPACE_CLOSE

#endif
