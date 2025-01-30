// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_io_input_file.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_io_input_spec.hpp"
#include "4C_io_value_parser.hpp"
#include "4C_utils_string.hpp"

#include <ryml.hpp>
#include <ryml_std.hpp>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <utility>

FOUR_C_NAMESPACE_OPEN

namespace Core::IO
{
  namespace
  {
    std::string to_string(const ryml::csubstr str) { return std::string(str.data(), str.size()); };
  }  // namespace

  namespace Internal
  {
    struct SectionContent;

    //! Helper for PIMPL idiom.
    class InputFileFragmentImpl
    {
     public:
      // This depends on the source of the fragment. Dat files produce string_views, yaml files
      // can directly store the node with all internal structure.
      std::variant<std::string_view, ryml::ConstNodeRef> fragment;

      /**
       * Store a pointer to the whole section this fragment belongs to.
       */
      SectionContent* section;
    };

    /**
     * Internal storage for the content of a section.
     */
    struct SectionContent
    {
      struct DatContent
      {
        //! The raw chars of the input.
        std::vector<char> raw_content;

        //! String views into #raw_content.
        std::vector<std::string_view> lines;
      };

      struct YamlContent
      {
        //! Node in the tree of the associated InputFileImpl.
        ryml::NodeRef node;

        //! String representation of the section in the .dat format. This is intended as a buffer
        //! to store the content if it should be requested via `get_as_dat_style_string()`. Do not
        //! assume this to contain a particular content.
        std::string dat_style_string{};
      };

      //! Content of the section is either in the .dat format or in the yaml format.
      std::variant<DatContent, YamlContent> content;

      /**
       * The actual content of the section sorted into fragments.
       */
      std::vector<InputFileFragmentImpl> fragments_impl;

      /**
       * The InputFile::Fragment is exposed to the user and contains a pointer to the actual
       * data stored in #fragments_impl.
       */
      std::vector<InputFile::Fragment> fragments;

      //! The file the section was read from.
      std::string file;

      DatContent& as_dat() { return std::get<DatContent>(content); }
      [[nodiscard]] const DatContent& as_dat() const { return std::get<DatContent>(content); }

      YamlContent& as_yaml() { return std::get<YamlContent>(content); }
      [[nodiscard]] const YamlContent& as_yaml() const { return std::get<YamlContent>(content); }

      /**
       * Depending on what is stored in this section, set up the Fragment objects to point to the
       * correct backend data.
       */
      void set_up_fragments()
      {
        if (std::holds_alternative<DatContent>(content))
        {
          auto& lines = as_dat().lines;
          fragments_impl.clear();
          fragments.clear();
          fragments_impl.reserve(lines.size());
          fragments.reserve(lines.size());

          for (const auto& line : lines)
          {
            auto& ref = fragments_impl.emplace_back(InputFileFragmentImpl{
                .fragment = line,
                .section = this,
            });
            fragments.emplace_back(&ref);
          }
        }
        else
        {
          // extract the fragments from the yaml node
          auto& node = as_yaml().node;
          fragments_impl.clear();
          fragments.clear();
          fragments_impl.reserve(node.num_children());
          fragments.reserve(node.num_children());

          FOUR_C_ASSERT(node.is_seq() || node.is_map(),
              "Section '%s' is neither a sequence nor a map.", to_string(node.key()).c_str());

          for (auto child : node.children())
          {
            auto& ref = fragments_impl.emplace_back(InputFileFragmentImpl{
                .fragment = child,
                .section = this,
            });
            fragments.emplace_back(&ref);
          }
        }
      }

      void pack(Communication::PackBuffer& data) const
      {
        Core::Communication::add_to_pack(data, file);
        Core::Communication::add_to_pack(data, content.index());
        if (std::holds_alternative<DatContent>(content))
        {
          const auto& raw_content = as_dat().raw_content;
          auto& lines = as_dat().lines;
          Core::Communication::add_to_pack(data, raw_content);

          // String_views are not packable, so we store offsets.
          std::vector<std::size_t> offsets;
          offsets.reserve(lines.size());
          for (const auto& line : lines)
          {
            FOUR_C_ASSERT(line.data() >= raw_content.data() &&
                              line.data() <= raw_content.data() + raw_content.size(),
                "Line data out of bounds.");
            const std::size_t offset = (line.data() - raw_content.data()) / sizeof(char);
            FOUR_C_ASSERT(offset <= raw_content.size(), "Offset out of bounds.");
            offsets.push_back(offset);
          }
          FOUR_C_ASSERT(
              offsets.empty() || offsets.back() + lines.back().size() == raw_content.size(),
              "Offset out of bounds.");
          Core::Communication::add_to_pack(data, offsets);
        }
        else if (std::holds_alternative<YamlContent>(content))
        {
          // Nothing to communicate for yaml content.
        }
        else
        {
          FOUR_C_THROW("Unknown content type.");
        }
      }


      void unpack(Communication::UnpackBuffer& buffer)
      {
        Core::Communication::extract_from_pack(buffer, file);
        std::size_t content_index;
        Core::Communication::extract_from_pack(buffer, content_index);

        if (content_index == 0)
        {
          auto& raw_content = std::get<DatContent>(content).raw_content;
          auto& lines = std::get<DatContent>(content).lines;

          Core::Communication::extract_from_pack(buffer, raw_content);

          std::vector<std::size_t> offsets;
          Core::Communication::extract_from_pack(buffer, offsets);
          lines.clear();
          for (std::size_t i = 0; i < offsets.size(); ++i)
          {
            const char* start = raw_content.data() + offsets[i];
            const std::size_t length = (i + 1 < offsets.size() ? (offsets[i + 1] - offsets[i])
                                                               : (raw_content.size()) - offsets[i]);
            FOUR_C_ASSERT(start >= raw_content.data() &&
                              start + length <= raw_content.data() + raw_content.size(),
                "Line data out of bounds.");
            lines.emplace_back(start, length);
          }
        }
        else if (content_index == 1)
        {
          // Nothing to communicate for yaml content.
        }
        else
        {
          FOUR_C_THROW("Unknown content index %d", content_index);
        }
      }

      SectionContent() = default;

      // Prevent copies: making a copy would invalidate the string views.
      SectionContent(const SectionContent&) = delete;
      SectionContent& operator=(const SectionContent&) = delete;

      SectionContent(SectionContent&&) = default;
      SectionContent& operator=(SectionContent&&) = default;
    };


    //! Helper for PIMPL idiom.
    class InputFileImpl
    {
     public:
      InputFileImpl(MPI_Comm comm) : comm_(comm)
      {
        // Root node is a map that stores the file paths as keys and the file content as values.
        yaml_tree_.rootref() |= ryml::MAP;
      }

      MPI_Comm comm_;

      std::unordered_map<std::string, SectionContent> content_by_section_;

      /**
       * This is the merged tree of all yaml input files combined, excluding any "INCLUDES"
       * sections. The data from different files will be stored under top-level keys corresponding
       * to the file paths.
       */
      ryml::Tree yaml_tree_;

      std::map<std::string, bool> knownsections_;
    };

  }  // namespace Internal

  namespace
  {
    /**
     * Sections that contain at least this number of entries are considered huge and are only
     * available on rank 0.
     */
    constexpr std::size_t huge_section_threshold = 10'000;

    //! The different ways we want to handle sections in the input file.
    enum class SectionType
    {
      //! A section that is read directly.
      normal,
      //! A section that mentions other files that are included and need to be read.
      include,
    };


    std::filesystem::path get_include_path(
        const std::string& include_line, const std::filesystem::path& current_file)
    {
      // Interpret the path as relative to the currently read file, if is not absolute
      std::filesystem::path included_file(include_line);
      if (!included_file.is_absolute())
      {
        included_file = current_file.parent_path() / included_file;
      }
      FOUR_C_ASSERT_ALWAYS(
          std::filesystem::status(included_file).type() == std::filesystem::file_type::regular,
          "Included file '%s' is not a regular file. Does the file exist?", included_file.c_str());
      return included_file;
    }

    void join_lines(std::list<std::string>& list_of_lines, std::vector<char>& raw_content,
        std::vector<std::string_view>& lines)
    {
      FOUR_C_ASSERT(raw_content.empty() && lines.empty(),
          "Implementation error: raw_content and lines must be empty.");

      // Sum up the length of all lines to reserve the memory for the raw content.
      const std::size_t raw_content_size =
          std::accumulate(list_of_lines.begin(), list_of_lines.end(), std::size_t{0},
              [](std::size_t sum, const auto& line) { return sum + line.size(); });

      raw_content.reserve(raw_content_size);
      lines.reserve(list_of_lines.size());

      for (const auto& line : list_of_lines)
      {
        FOUR_C_ASSERT(raw_content.data(), "Implementation error: raw_content must be allocated.");
        const auto* start_of_line = raw_content.data() + raw_content.size();
        raw_content.insert(raw_content.end(), line.begin(), line.end());
        lines.emplace_back(start_of_line, line.size());
      }
    }

    std::vector<std::filesystem::path> read_dat_content(const std::filesystem::path& file_path,
        std::unordered_map<std::string, Internal::SectionContent>& content_by_section)
    {
      const auto name_of_section = [](const std::string& section_header)
      {
        auto pos = section_header.rfind("--");
        if (pos == std::string::npos) return std::string{};
        return Core::Utils::trim(section_header.substr(pos + 2));
      };

      std::ifstream file(file_path);
      if (not file) FOUR_C_THROW("Unable to open file: %s", file_path.c_str());

      // Tracking variables while walking through the file
      std::vector<std::filesystem::path> included_files;
      SectionType current_section_type = SectionType::normal;
      std::list<std::string> list_of_lines;
      Internal::SectionContent* current_section_content = nullptr;
      std::string line;

      // Loop over all input lines. This reads the actual file contents and determines whether a
      // line is to be read immediately or should be excluded because it is in one of the excluded
      // sections.
      while (getline(file, line))
      {
        // In case we are reading an include section, a comment needs to be preceded by
        // whitespace. Otherwise, we would treat double slashes as comments, although they are
        // part of the file path.
        if (current_section_type == SectionType::include)
        {
          // Take care to remove comments only if they are preceded by whitespace.
          line = Core::Utils::strip_comment(line, " //");
          if (line.empty()) continue;

          // Additionally check if the first token is a comment to handle the case where the
          // comment starts at the beginning of the line.
          if (line.starts_with("//")) continue;
        }
        // Remove comments, trailing and leading whitespaces, compact internal whitespaces
        else
        {
          line = Core::Utils::strip_comment(line);
        }

        // line is now empty
        if (line.size() == 0) continue;

        // This line starts a new section
        if (line.starts_with("--"))
        {
          // Finish the current section.
          if (current_section_content)
          {
            join_lines(list_of_lines, current_section_content->as_dat().raw_content,
                current_section_content->as_dat().lines);
          }
          list_of_lines.clear();

          const auto name = name_of_section(line);
          FOUR_C_ASSERT(name.size() > 0, "Section name must not be empty.");

          // Determine what kind of new section we started.
          if (line.rfind("--INCLUDES") != std::string::npos)
          {
            current_section_type = SectionType::include;
            current_section_content = nullptr;
          }
          else
          {
            current_section_type = SectionType::normal;
            FOUR_C_ASSERT_ALWAYS(content_by_section.find(name) == content_by_section.end(),
                "Section '%s' is defined again in file '%s'.", name.c_str(), file_path.c_str());

            content_by_section[name] = {};
            current_section_content = &content_by_section[name];
            current_section_content->file = file_path;
          }
        }
        // The line is part of a section.
        else
        {
          switch (current_section_type)
          {
            case SectionType::normal:
            {
              list_of_lines.emplace_back(line);
              break;
            }
            case SectionType::include:
            {
              if (!line.starts_with("--"))
              {
                included_files.emplace_back(get_include_path(line, file_path));
              }
              break;
            }
          }
        }
      }

      // Finish the current section.
      if (current_section_content)
      {
        join_lines(list_of_lines, current_section_content->as_dat().raw_content,
            current_section_content->as_dat().lines);
      }

      return included_files;
    }


    std::vector<std::filesystem::path> read_yaml_content(
        const std::filesystem::path& file_path, Internal::InputFileImpl& input_file_impl)
    {
      std::vector<std::filesystem::path> included_files;

      // Read the whole file into a string ...
      std::ifstream file(file_path);
      std::ostringstream ss;
      ss << file.rdbuf();
      const std::string file_content = ss.str();

      // ... and parse it into a new node in the yaml tree.
      ryml::NodeRef file_node = input_file_impl.yaml_tree_.rootref().append_child();
      {
        std::string file_path_str = file_path.string();
        file_node << ryml::key(file_path_str);
        ryml::parse_in_arena(ryml::to_csubstr(file_content), file_node);
      }

      // Treat the special section "INCLUDES" to find more included files. Remove the section once
      // processed.
      if (file_node.has_child("INCLUDES"))
      {
        ryml::ConstNodeRef node = file_node["INCLUDES"];
        if (node.has_val())
        {
          included_files.emplace_back(get_include_path(to_string(node.val()), file_path));
        }
        else if (node.is_seq())
        {
          for (const auto& include_node : node)
          {
            included_files.emplace_back(get_include_path(to_string(include_node.val()), file_path));
          }
        }
        else
        {
          FOUR_C_THROW("INCLUDES section must contain a single file or a sequence.");
        }

        // Now we can drop the node containing the INCLUDES from the tree.
        input_file_impl.yaml_tree_.remove(node.id());
      }

      return included_files;
    }
  }  // namespace


  std::string_view InputFile::Fragment::get_as_dat_style_string() const
  {
    FOUR_C_ASSERT(pimpl_, "Implementation error: fragment is not initialized.");
    if (std::holds_alternative<std::string_view>(pimpl_->fragment))
    {
      return std::get<std::string_view>(pimpl_->fragment);
    }
    else
    {
      const auto& node = std::get<ryml::ConstNodeRef>(pimpl_->fragment);
      if (node.is_val())
      {
        return std::string_view(node.val().data(), node.val().size());
      }
      else
      {
        // flatten the yaml node into a string
        std::ostringstream ss;
        ss << node;
        // replace all the yaml markup characters with spaces to receive a dat-style string
        auto& str = pimpl_->section->as_yaml().dat_style_string;
        str = ss.str();
        std::replace(str.begin(), str.end(), '\n', ' ');
        std::replace(str.begin(), str.end(), ':', ' ');
        std::replace(str.begin(), str.end(), ',', ' ');
        std::replace(str.begin(), str.end(), '[', ' ');
        std::replace(str.begin(), str.end(), ']', ' ');
        str = Core::Utils::trim(str);

        return std::string_view(str);
      }
    }
  }


  std::optional<InputParameterContainer> InputFile::Fragment::match(const InputSpec& spec) const
  {
    FOUR_C_ASSERT(pimpl_, "Implementation error: fragment is not initialized.");

    Core::IO::ValueParser parser{get_as_dat_style_string(),
        {.base_path = std::filesystem::path(pimpl_->section->file).parent_path()}};
    InputParameterContainer container;
    spec.fully_parse(parser, container);
    return container;
  }


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  InputFile::InputFile(std::string filename, MPI_Comm comm)
      : pimpl_(std::make_unique<Internal::InputFileImpl>(comm))
  {
    read_generic(filename);
  }


  // Note: defaulted in implementation file to allow for use of incomplete type in PIMPL unique_ptr.
  InputFile::~InputFile() = default;


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  std::filesystem::path InputFile::file_for_section(const std::string& section_name) const
  {
    auto it = pimpl_->content_by_section_.find(section_name);
    if (it == pimpl_->content_by_section_.end()) return std::filesystem::path{};
    return std::filesystem::path{it->second.file};
  }


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  bool InputFile::has_section(const std::string& section_name) const
  {
    const bool known_somewhere = Core::Communication::all_reduce<bool>(
        pimpl_->content_by_section_.contains(section_name),
        [](const bool& r, const bool& in) { return r || in; }, pimpl_->comm_);
    return known_somewhere;
  }



  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void InputFile::read_generic(const std::filesystem::path& top_level_file)
  {
    if (Core::Communication::my_mpi_rank(pimpl_->comm_) == 0)
    {
      // Start by "including" the top-level file.
      std::list<std::filesystem::path> included_files{top_level_file};

      // We use a hand-rolled loop here because the list keeps growing; thus we need to
      // continuously re-evaluate where the end of the list is.
      for (auto file_it = included_files.begin(); file_it != included_files.end(); ++file_it)
      {
        std::vector<std::filesystem::path> new_include_files = std::invoke(
            [&]
            {
              const auto file_extension = file_it->extension().string();
              // Note that json is valid yaml and we can read it with the yaml parser.
              if (file_extension == ".yaml" || file_extension == ".yml" ||
                  file_extension == ".json")
              {
                return read_yaml_content(*file_it, *pimpl_);
              }
              else
              {
                return read_dat_content(*file_it, pimpl_->content_by_section_);
              }
            });

        // Check that the file is not included twice
        for (const auto& file : new_include_files)
        {
          if (std::ranges::find(included_files, file) != included_files.end())
          {
            FOUR_C_THROW(
                "File '%s' was already included before.\n Cycles are not allowed.", file.c_str());
          }
          else
          {
            included_files.emplace_back(file);
          }
        }
      }
    }

    // Communicate the dat content
    if (Core::Communication::my_mpi_rank(pimpl_->comm_) == 0)
    {
      // Temporarily move the sections that are not huge into a separate map.
      std::unordered_map<std::string, Internal::SectionContent> non_huge_sections;

      for (auto&& [section_name, content] : pimpl_->content_by_section_)
      {
        if (std::holds_alternative<Internal::SectionContent::DatContent>(content.content))
        {
          if (content.as_dat().lines.size() < huge_section_threshold)
          {
            non_huge_sections[section_name] = std::move(content);
          }
        }
      }

      Core::Communication::broadcast(non_huge_sections, 0, pimpl_->comm_);

      // Move the non-huge sections back into the main map.
      for (auto&& [section_name, content] : non_huge_sections)
      {
        pimpl_->content_by_section_[section_name] = std::move(content);
      }
    }
    else
    {
      // Other ranks receive the non-huge sections.
      Core::Communication::broadcast(pimpl_->content_by_section_, 0, pimpl_->comm_);
    }

    // Communicate the yaml content
    {
      if (Core::Communication::my_mpi_rank(pimpl_->comm_) == 0)
      {
        std::stringstream ss;
        ss << pimpl_->yaml_tree_;
        std::string yaml_str = ss.str();
        Core::Communication::broadcast(yaml_str, /*root*/ 0, pimpl_->comm_);
      }
      else
      {
        std::string yaml_str;
        Core::Communication::broadcast(yaml_str, /*root*/ 0, pimpl_->comm_);
        pimpl_->yaml_tree_.reserve_arena(yaml_str.size());
        ryml::parse_in_arena(ryml::to_csubstr(yaml_str), pimpl_->yaml_tree_.rootref());
      }

      // Now all ranks agree on the yaml tree. Parse out the sections and fill the SectionContent.
      for (auto file_node : pimpl_->yaml_tree_.rootref())
      {
        for (auto node : file_node.children())
        {
          const std::string section_name = to_string(node.key());
          FOUR_C_ASSERT_ALWAYS(!pimpl_->content_by_section_.contains(section_name),
              "Section '%s' is defined more than once.", section_name.c_str());

          Internal::SectionContent& content = pimpl_->content_by_section_[section_name];
          content.content = Internal::SectionContent::YamlContent{.node = node};
          content.file = to_string(file_node.key());
        }
      }
    }

    for (auto& [name, content] : pimpl_->content_by_section_)
    {
      content.set_up_fragments();
    }

    // the following section names are always regarded as valid
    record_section_used("TITLE");
    record_section_used("FUNCT1");
    record_section_used("FUNCT2");
    record_section_used("FUNCT3");
    record_section_used("FUNCT4");
    record_section_used("FUNCT5");
    record_section_used("FUNCT6");
    record_section_used("FUNCT7");
    record_section_used("FUNCT8");
    record_section_used("FUNCT9");
    record_section_used("FUNCT10");
    record_section_used("FUNCT11");
    record_section_used("FUNCT12");
    record_section_used("FUNCT13");
    record_section_used("FUNCT14");
    record_section_used("FUNCT15");
    record_section_used("FUNCT16");
    record_section_used("FUNCT17");
    record_section_used("FUNCT18");
    record_section_used("FUNCT19");
    record_section_used("FUNCT20");
  }

  InputFile::FragmentIteratorRange InputFile::in_section(const std::string& section_name)
  {
    static const std::vector<Fragment> empty;

    // Early return in case the section does not exist at all.
    const bool known_somewhere = has_section(section_name);
    if (!known_somewhere)
    {
      return std::views::all(empty);
    }

    record_section_used(section_name);
    const bool locally_known = pimpl_->content_by_section_.contains(section_name);
    const bool known_everywhere = Core::Communication::all_reduce<bool>(
        locally_known, [](const bool& r, const bool& in) { return r && in; }, pimpl_->comm_);

    if (known_everywhere)
    {
      // Take a const reference to the section content to match the return type.
      const auto& lines = pimpl_->content_by_section_.at(section_name).fragments;
      return std::views::all(lines);
    }
    else
    // Distribute the content of the section to all ranks.
    {
      FOUR_C_ASSERT((!locally_known && (Core::Communication::my_mpi_rank(pimpl_->comm_) > 0)) ||
                        (locally_known && (Core::Communication::my_mpi_rank(pimpl_->comm_) == 0)),
          "Implementation error: section should be known on rank 0 and unknown on others.");

      auto& content = pimpl_->content_by_section_[section_name];
      Core::Communication::broadcast(content, 0, pimpl_->comm_);
      content.set_up_fragments();

      const auto& lines = pimpl_->content_by_section_.at(section_name).fragments;
      return std::views::all(lines);
    }
  }



  InputFile::FragmentIteratorRange InputFile::in_section_rank_0_only(
      const std::string& section_name)
  {
    if (Core::Communication::my_mpi_rank(pimpl_->comm_) == 0 &&
        pimpl_->content_by_section_.contains(section_name))
    {
      record_section_used(section_name);
      const auto& lines = pimpl_->content_by_section_.at(section_name).fragments;
      return std::views::all(lines);
    }
    else
    {
      static const std::vector<Fragment> empty;
      return std::views::all(empty);
    }
  }


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  bool InputFile::print_unknown_sections(std::ostream& out) const
  {
    using MapType = decltype(pimpl_->knownsections_);
    const auto merged_map = Core::Communication::all_reduce<MapType>(
        pimpl_->knownsections_,
        [](const MapType& r, const MapType& in)
        {
          MapType result = r;
          for (const auto& [key, value] : in)
          {
            result[key] |= value;
          }
          return result;
        },
        pimpl_->comm_);
    const bool printout = std::any_of(
        merged_map.begin(), merged_map.end(), [](const auto& kv) { return !kv.second; });

    if (printout and (Core::Communication::my_mpi_rank(get_comm()) == 0))
    {
      out << "\nERROR!"
          << "\n--------"
          << "\nThe following input file sections remained unused (obsolete or typo?):\n";
      for (const auto& [section_name, known] : pimpl_->knownsections_)
      {
        if (!known) out << section_name << '\n';
      }
      out << '\n';
    }

    return printout;
  }


  void InputFile::record_section_used(const std::string& section_name)
  {
    pimpl_->knownsections_[section_name] = true;
  }


  MPI_Comm InputFile::get_comm() const { return pimpl_->comm_; }

}  // namespace Core::IO

FOUR_C_NAMESPACE_CLOSE
