// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_io_input_file.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_fem_nurbs_discretization_knotvector.hpp"
#include "4C_io_pstream.hpp"
#include "4C_linalg_utils_densematrix_communication.hpp"
#include "4C_utils_string.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Time.hpp>
#include <yaml-cpp/yaml.h>

#include <filesystem>
#include <sstream>
#include <utility>

#ifdef FOUR_C_ENABLE_FE_TRAPPING
#include <cfenv>
#include <format>
#endif


FOUR_C_NAMESPACE_OPEN

namespace Core::IO
{
  namespace
  {
    constexpr double tolerance_n = 1.0e-14;

    /**
     * Sections that contain at least this number of entries are considered huge and are only
     * available on rank 0.
     */
    constexpr std::size_t huge_section_threshold = 10'000;

    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
    Teuchos::ParameterList& find_sublist(std::string name, Teuchos::ParameterList& list)
    {
      Teuchos::ParameterList* sublist = &list;

      for (std::string::size_type pos = name.find('/'); pos != std::string::npos;
          pos = name.find('/'))
      {
        sublist = &sublist->sublist(name.substr(0, pos));
        name = name.substr(pos + 1);
      }

      return sublist->sublist(name);
    }

    void add_entry(const std::string& key, const std::string& value, Teuchos::ParameterList& list)
    {
      // safety check: Is there a duplicate of the same parameter?
      if (list.isParameter(key))
        FOUR_C_THROW("Duplicate parameter '%s' in sublist '%s'", key.c_str(), list.name().c_str());

      if (key.empty()) FOUR_C_THROW("Internal error: missing key.", key.c_str());
      // safety check: Is the parameter without any specified value?
      if (value.empty())
        FOUR_C_THROW("Missing value for parameter %s. Fix your input file!", key.c_str());

      {  // try to find an int
        std::stringstream ssi;
        int iv;

        ssi << value;
        ssi >> iv;

        if (ssi.eof())
        {
          list.set(key, iv);
          return;
        }
      }

#ifdef FOUR_C_ENABLE_FE_TRAPPING
      // somehow the following test whether we have a double or not
      // creates always an internal floating point exception (FE_INVALID). An alternative
      // implementation using boost::lexical_cast<double> does not solve this problem!
      // Better temporarily disable this floating point exception in the following,
      // so that we can go on.
      feclearexcept(FE_INVALID);
      /*feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW);*/
      fedisableexcept(FE_INVALID);
#endif

      {  // try to find a double
        std::stringstream ssd;
        double dv;

        ssd << value;
        ssd >> dv;

#ifdef FOUR_C_ENABLE_FE_TRAPPING
        feclearexcept(FE_INVALID);
        /*feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW);*/
        feenableexcept(FE_INVALID | FE_DIVBYZERO);
#endif

        if (ssd.eof())
        {
          list.set(key, dv);
          return;
        }
      }

      // if it is not an int or a double it must be a string
      list.set(key, value);
    }

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
        std::unordered_map<std::string, InputFile::SectionContent>& content_by_section)
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
      InputFile::SectionContent* current_section_content = nullptr;
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
            join_lines(list_of_lines, current_section_content->raw_content,
                current_section_content->lines);
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
        join_lines(
            list_of_lines, current_section_content->raw_content, current_section_content->lines);
      }

      return included_files;
    }

    std::vector<std::filesystem::path> read_yaml_content(const std::filesystem::path& file_path,
        std::unordered_map<std::string, InputFile::SectionContent>& content_by_section)
    {
      std::vector<std::filesystem::path> included_files;

      // In this first iteration of the YAML support, we map the constructs from a YAML file back
      // to constructs in a dat file. This means that top-level sections are pre-fixed with "--" and
      // the key-value pairs are mapped to "key = value" lines.
      //
      // caveats:
      // - YAML files can only be read in full.
      YAML::Node config = YAML::LoadFile(file_path);
      for (const auto& section : config)
      {
        const std::string& section_name = section.first.as<std::string>();
        const auto& entries = section.second;

        // If this is the special section "INCLUDES", we need to handle it differently.
        if (section_name == "INCLUDES")
        {
          if (entries.IsScalar())
            included_files.emplace_back(get_include_path(entries.as<std::string>(), file_path));
          else if (entries.IsSequence())
          {
            for (const auto& entry : entries)
            {
              FOUR_C_ASSERT_ALWAYS(
                  entry.IsScalar(), "Only scalar entries are supported in the INCLUDES section.");
              const auto line = entry.as<std::string>();
              included_files.emplace_back(get_include_path(line, file_path));
            }
          }
          else
            FOUR_C_THROW("INCLUDES section must contain a single file or a sequence.");

          continue;
        }

        FOUR_C_ASSERT_ALWAYS(content_by_section.find(section_name) == content_by_section.end(),
            "Section '%s' is defined again in file '%s'.", section_name.c_str(), file_path.c_str());

        auto& current_content = content_by_section[section_name];
        current_content.file = file_path;
        std::list<std::string> list_of_lines;

        const auto read_flat_sequence = [&](const YAML::Node& node)
        {
          for (const auto& entry : node)
          {
            FOUR_C_ASSERT_ALWAYS(entry.IsScalar(),
                "While reading section '%s': "
                "only scalar entries are supported in sequences.",
                section_name.c_str());
            const auto line = entry.as<std::string>();
            list_of_lines.emplace_back(line);
          }
        };

        const auto read_map = [&](const YAML::Node& node)
        {
          for (const auto& entry : node)
          {
            FOUR_C_ASSERT_ALWAYS(entry.first.IsScalar() && entry.second.IsScalar(),
                "While reading section '%s': "
                "only scalar key-value pairs are supported in maps.",
                section_name.c_str());
            const auto line =
                entry.first.as<std::string>() + " = " + entry.second.as<std::string>();
            list_of_lines.emplace_back(line);
          }
        };

        switch (entries.Type())
        {
          case YAML::NodeType::Undefined:
          case YAML::NodeType::Null:
          case YAML::NodeType::Scalar:
            FOUR_C_THROW(
                "Entries in section %s must either form a map or a sequence", section_name.c_str());
          case YAML::NodeType::Sequence:
            read_flat_sequence(entries);
            break;
          case YAML::NodeType::Map:
            read_map(entries);
            break;
        }

        // Finish the current section by condensing the lines into the content.
        join_lines(list_of_lines, current_content.raw_content, current_content.lines);
      }

      return included_files;
    }
  }  // namespace

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  InputFile::InputFile(std::string filename, MPI_Comm comm) : comm_(comm)
  {
    read_generic(filename);
  }


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  std::filesystem::path InputFile::file_for_section(const std::string& section_name) const
  {
    auto it = content_by_section_.find(section_name);
    if (it == content_by_section_.end()) return std::filesystem::path{};
    return std::filesystem::path{it->second.file};
  }


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  bool InputFile::has_section(const std::string& section_name) const
  {
    const bool known_somewhere = Core::Communication::all_reduce<bool>(
        content_by_section_.contains(section_name),
        [](const bool& r, const bool& in) { return r || in; }, comm_);
    return known_somewhere;
  }


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  bool read_parameters_in_section(
      InputFile& input, const std::string& section_name, Teuchos::ParameterList& list)
  {
    if (section_name.empty()) FOUR_C_THROW("Empty section name given.");

    Teuchos::ParameterList& sublist = find_sublist(section_name, list);

    for (const auto& line : input.lines_in_section(section_name))
    {
      const auto& [key, value] = read_key_value(std::string(line));

      add_entry(key, value, sublist);
    }

    return true;
  }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void read_design(InputFile& input, const std::string& name,
      std::vector<std::vector<int>>& dobj_fenode,
      const std::function<const Core::FE::Discretization&(const std::string& name)>&
          get_discretization)
  {
    std::map<int, std::set<int>> topology;

    std::string sectionname = name + "-NODE TOPOLOGY";
    std::string marker = sectionname;

    // Store lines that need special treatment
    std::vector<std::string> box_face_conditions;

    for (const auto& l : input.lines_in_section_rank_0_only(marker))
    {
      int dobj;
      int nodeid;
      std::string nname;
      std::string dname;
      std::istringstream stream{std::string(l)};
      stream >> nname;
      if (not stream) FOUR_C_THROW("Illegal line in section '%s': '%s'", marker.c_str(), l.data());

      if (nname == "NODE")  // plain old reading of the design nodes from the .dat-file
      {
        stream >> nodeid >> dname >> dobj;
        topology[dobj - 1].insert(nodeid - 1);
      }
      else  // fancy specification of the design nodes by specifying min or max of the domain
      {     // works best on rectangular domains ;)
        // Store the specification and broadcast it to all processors
        box_face_conditions.emplace_back(l);
      }
    }

    // Broadcast what we have read on rank 0 to all other ranks
    Core::Communication::broadcast(topology, 0, input.get_comm());
    Core::Communication::broadcast(box_face_conditions, 0, input.get_comm());

    // Special treatment if we have any box face conditions
    {
      for (const auto& l : box_face_conditions)
      {
        int dobj;
        std::string nname;
        std::string dname;
        std::string disname;
        std::array<int, 3> dir = {0, 0, 0};

        std::istringstream stream{l};
        stream >> nname;
        if (not stream)
          FOUR_C_THROW("Illegal line in section '%s': '%s'", marker.c_str(), l.data());

        if (nname == "CORNER" && name == "DNODE")
        {
          std::string tmp;
          stream >> disname;
          for (int i = 0; i < 3; ++i)
          {
            stream >> tmp;
            if (tmp.size() != 2 || tmp[0] < 'x' || tmp[0] > 'z' || (tmp[1] != '+' && tmp[1] != '-'))
              FOUR_C_THROW("Illegal design node definition.");
            dir[tmp[0] - 'x'] = (tmp[1] == '+') ? 1 : -1;
          }
          stream >> dname >> dobj;
        }
        else if (nname == "EDGE" && name == "DLINE")
        {
          std::string tmp;
          stream >> disname;
          for (int i = 0; i < 2; ++i)
          {
            stream >> tmp;
            if (tmp.size() != 2 || tmp[0] < 'x' || tmp[0] > 'z' || (tmp[1] != '+' && tmp[1] != '-'))
              FOUR_C_THROW("Illegal design node definition.");
            dir[tmp[0] - 'x'] = (tmp[1] == '+') ? 1 : -1;
          }
          stream >> dname >> dobj;
        }
        else if (nname == "SIDE" && name == "DSURF")
        {
          std::string tmp;
          stream >> disname;
          stream >> tmp;
          if (tmp.size() != 2 || tmp[0] < 'x' || tmp[0] > 'z' || (tmp[1] != '+' && tmp[1] != '-'))
            FOUR_C_THROW("Illegal design node definition.");
          dir[tmp[0] - 'x'] = (tmp[1] == '+') ? 1 : -1;
          stream >> dname >> dobj;
        }
        else if (nname == "VOLUME" && name == "DVOL")
        {
          stream >> disname;
          stream >> dname >> dobj;
        }
        else
        {
          FOUR_C_THROW("Illegal line in section '%s': '%s'", marker.c_str(), l.data());
        }

        const Core::FE::Discretization& actdis = get_discretization(disname);

        std::vector<double> box_specifications;
        {
          for (int init = 0; init < 9; ++init) box_specifications.push_back(0.0);
          if (Core::Communication::my_mpi_rank(input.get_comm()) == 0)  // Reading is done by proc 0
          {
            // get original domain section from the *.dat-file
            std::string dommarker = disname + " DOMAIN";
            std::transform(dommarker.begin(), dommarker.end(), dommarker.begin(), ::toupper);

            for (const auto& line : input.lines_in_section_rank_0_only(dommarker))
            {
              std::istringstream t{std::string{line}};
              std::string key;
              t >> key;

              if (key == "LOWER_BOUND")
              {
                t >> box_specifications[0] >> box_specifications[1] >> box_specifications[2];
              }
              else if (key == "UPPER_BOUND")
              {
                t >> box_specifications[3] >> box_specifications[4] >> box_specifications[5];
              }
              else if (key == "ROTATION")
              {
                t >> box_specifications[6] >> box_specifications[7] >> box_specifications[8];
              }
            }
          }
          // All other processors get this info broadcasted
          Core::Communication::broadcast(box_specifications.data(),
              static_cast<int>(box_specifications.size()), 0, input.get_comm());
        }

        // determine the active discretizations bounding box
        std::array<double, 6> bbox;
        for (size_t i = 0; i < sizeof(bbox) / sizeof(bbox[0]); ++i) bbox[i] = box_specifications[i];

        // manipulate the bounding box according to the specified condition
        for (size_t i = 0; i < 3; ++i)
        {
          switch (dir[i])
          {
            case 0:
              bbox[i + 0] = std::numeric_limits<double>::max();
              bbox[i + 3] = -std::numeric_limits<double>::max();
              break;
            case -1:
              bbox[i] += tolerance_n;
              bbox[i + 3] = std::numeric_limits<double>::max();
              break;
            case 1:
              bbox[i] = -std::numeric_limits<double>::max();
              bbox[i + 3] -= tolerance_n;
              break;
            default:
              FOUR_C_THROW("Invalid BC specification");
          }
        }

        // collect all nodes which are outside the adapted bounding box
        std::set<int> dnodes;
        for (const auto* node : actdis.my_row_node_range())
        {
          const auto& coord = node->x();
          std::array<double, 3> coords;
          coords[0] = coord[0];
          coords[1] = coord[1];
          coords[2] = coord[2];
          // rotate back to identify condition, if a rotation is defined
          static const int rotoffset = 6;
          for (int rotaxis = 2; rotaxis > -1; --rotaxis)
          {
            if (box_specifications[rotaxis + rotoffset] != 0.0)
            {
              std::array<double, 3> coordm;
              coordm[0] = (box_specifications[0] + box_specifications[3]) / 2.;
              coordm[1] = (box_specifications[1] + box_specifications[4]) / 2.;
              coordm[2] = (box_specifications[2] + box_specifications[5]) / 2.;
              // add rotation around mitpoint here.
              std::array<double, 3> dx;
              dx[0] = coords[0] - coordm[0];
              dx[1] = coords[1] - coordm[1];
              dx[2] = coords[2] - coordm[2];

              double calpha = cos(-box_specifications[rotaxis + rotoffset] * M_PI / 180);
              double salpha = sin(-box_specifications[rotaxis + rotoffset] * M_PI / 180);

              coords[0] = coordm[0];  //+ calpha*dx[0] + salpha*dx[1];
              coords[1] = coordm[1];  //+ -salpha*dx[0] + calpha*dx[1];
              coords[2] = coordm[2];

              coords[(rotaxis + 1) % 3] +=
                  calpha * dx[(rotaxis + 1) % 3] + salpha * dx[(rotaxis + 2) % 3];
              coords[(rotaxis + 2) % 3] +=
                  calpha * dx[(rotaxis + 2) % 3] - salpha * dx[(rotaxis + 1) % 3];
              coords[rotaxis] += dx[rotaxis];
            }
          }

          if ((coords[0] <= bbox[0] || coords[0] >= bbox[3]) &&
              (coords[1] <= bbox[1] || coords[1] >= bbox[4]) &&
              (coords[2] <= bbox[2] || coords[2] >= bbox[5]))
            dnodes.insert(node->id());
        }
        Core::LinAlg::gather_all(dnodes, input.get_comm());
        topology[dobj - 1].insert(dnodes.begin(), dnodes.end());


        if (dname.substr(0, name.length()) != name)
          FOUR_C_THROW("Illegal line in section '%s': '%s'\n%s found, where %s was expected",
              marker.c_str(), l.data(), dname.substr(0, name.length()).c_str(), name.c_str());
      }
    }

    if (topology.size() > 0)
    {
      int max_num_dobj = topology.rbegin()->first;
      if (max_num_dobj >= static_cast<int>(dobj_fenode.size()))
        dobj_fenode.resize(max_num_dobj + 1);
      // copy all design object entries
      for (auto& topo : topology)
      {
        // we copy from a std::set, thus the gids are sorted
        dobj_fenode[topo.first].reserve(topo.second.size());
        dobj_fenode[topo.first].assign(topo.second.begin(), topo.second.end());
      }
    }
  }


  //----------------------------------------------------------------------
  /// read a knotvector section (for isogeometric analysis)
  //----------------------------------------------------------------------
  void read_knots(InputFile& input, const std::string& name,
      std::shared_ptr<Core::FE::Nurbs::Knotvector>& disknots)
  {
    // io to shell
    const int myrank = Core::Communication::my_mpi_rank(input.get_comm());

    Teuchos::Time time("", true);

    // only the knotvector section of this discretisation
    // type is of interest
    std::string field;
    if (name == "fluid" or name == "xfluid" or name == "porofluid")
    {
      field = "FLUID";
    }
    else if (name == "structure")
    {
      field = "STRUCTURE";
    }
    else if (name == "ale")
    {
      field = "ALE";
    }
    else if (name == "scatra")
    {
      field = "TRANSPORT";
    }
    else if (name == "thermo")
    {
      field = "THERMO";
    }
    else if (name == "scatra_micro")
    {
      field = "TRANSPORT2";
    }
    else
    {
      FOUR_C_THROW("Unknown discretization name for knotvector input\n");
    }

    // another valid section name was found
    const std::string sectionname = field + " KNOTVECTORS";

    if (myrank == 0)
    {
      Core::IO::cout << "Reading knot vectors for " << name << " discretization :\n";
      fflush(stdout);
    }

    // number of patches to be determined
    int npatches = 0;

    // dimension of nurbs patches
    int nurbs_dim = 0;

    //--------------------------------------------------------------------
    //--------------------------------------------------------------------
    //      first, determine number of patches and dimension of nurbs
    //--------------------------------------------------------------------
    //--------------------------------------------------------------------
    {
      // temporary string
      std::string tmp;
      // loop lines in file
      for (const auto& line : input.lines_in_section(sectionname))
      {
        // count number of patches in knotvector section of
        // this discretisation
        {
          std::string::size_type loc;
          std::istringstream file{std::string{line}};
          file >> tmp;

          // check for the number of dimensions
          loc = tmp.rfind("NURBS_DIMENSION");
          if (loc != std::string::npos)
          {
            // set number of nurbs dimension
            std::string str_nurbs_dim;
            file >> str_nurbs_dim;
            char* endptr = nullptr;
            nurbs_dim = static_cast<int>(strtol(str_nurbs_dim.c_str(), &endptr, 10));

            continue;
          }

          // check for a new patch
          loc = tmp.rfind("ID");
          if (loc != std::string::npos)
          {
            // increase number of patches
            npatches++;

            continue;
          }
        }
      }  // end loop through file
    }

    if (myrank == 0)
    {
      printf("                        %8d patches", npatches);
      fflush(stdout);
    }


    //--------------------------------------------------------------------
    //--------------------------------------------------------------------
    //                alloc knotvector object to fill
    //--------------------------------------------------------------------
    //--------------------------------------------------------------------

    // allocate knotvector for this dis
    disknots = std::make_shared<Core::FE::Nurbs::Knotvector>(nurbs_dim, npatches);

    // make sure that we have some Knotvector object to fill
    if (disknots == nullptr)
    {
      FOUR_C_THROW("disknots should have been allocated before");
    }

    //--------------------------------------------------------------------
    //--------------------------------------------------------------------
    //                finally read knotvector section
    //--------------------------------------------------------------------
    //--------------------------------------------------------------------
    {
      // this is a pointer to the knots of one patch in one direction
      // we will read them and put them
      std::vector<std::shared_ptr<std::vector<double>>> patch_knots(nurbs_dim);

      // temporary string
      std::string tmp;

      // start to read something when read is true
      bool read = false;

      // index for number of patch
      int npatch = 0;
      // index for u/v/w
      int actdim = -1;
      // ints for the number of knots
      std::vector<int> n_x_m_x_l(nurbs_dim);
      // ints for patches degrees
      std::vector<int> degree(nurbs_dim);
      // a vector of strings holding the knotvectortypes read
      std::vector<std::string> knotvectortype(nurbs_dim);

      // count for sanity check
      int count_read = 0;
      std::vector<int> count_vals(nurbs_dim);

      // loop lines in file
      for (const auto& line : input.lines_in_section(sectionname))
      {
        std::istringstream file{std::string{line}};
        file >> tmp;

        // check for a new patch
        std::string::size_type loc = tmp.rfind("BEGIN");
        if (loc != std::string::npos)
        {
          file >> tmp;

          // activate reading
          read = true;

          actdim = -1;

          // create vectors for knots in this patch
          for (int rr = 0; rr < nurbs_dim; ++rr)
          {
            patch_knots[rr] = std::make_shared<std::vector<double>>();
            (*(patch_knots[rr])).clear();
          }

          // reset counter for knot values
          for (int rr = 0; rr < nurbs_dim; rr++)
          {
            count_vals[rr] = 0;
          }

          continue;
        }

        // get ID of patch we are currently reading
        loc = tmp.rfind("ID");
        if (loc != std::string::npos)
        {
          std::string str_npatch;
          file >> str_npatch;

          char* endptr = nullptr;
          npatch = static_cast<int>(strtol(str_npatch.c_str(), &endptr, 10));
          npatch--;

          continue;
        }

        // get number of knots in the knotvector direction
        // we are currently reading
        loc = tmp.rfind("NUMKNOTS");
        if (loc != std::string::npos)
        {
          std::string str_numknots;
          file >> str_numknots;

          // increase dimesion for knotvector (i.e. next time
          // we'll fill the following knot vector)
          actdim++;
          if (actdim > nurbs_dim)
          {
            FOUR_C_THROW(
                "too many knotvectors, we only need one for each dimension (nurbs_dim = %d)\n",
                nurbs_dim);
          }

          char* endptr = nullptr;
          n_x_m_x_l[actdim] = static_cast<int>(strtol(str_numknots.c_str(), &endptr, 10));

          continue;
        }

        // get number of bspline polinomial associated with
        // knots in this direction
        loc = tmp.rfind("DEGREE");
        if (loc != std::string::npos)
        {
          std::string str_degree;
          file >> str_degree;

          char* endptr = nullptr;
          degree[actdim] = static_cast<int>(strtol(str_degree.c_str(), &endptr, 10));

          continue;
        }

        // get type of knotvector (interpolated or periodic)
        loc = tmp.rfind("TYPE");
        if (loc != std::string::npos)
        {
          std::string type;

          file >> type;
          knotvectortype[actdim] = type;

          continue;
        }

        // locate end of patch
        loc = tmp.rfind("END");
        if (loc != std::string::npos)
        {
          for (int rr = 0; rr < nurbs_dim; ++rr)
          {
            disknots->set_knots(
                rr, npatch, degree[rr], n_x_m_x_l[rr], knotvectortype[rr], patch_knots[rr]);
          }
          file >> tmp;
          // stop reading of knot values if we are here
          read = false;

          for (int rr = 0; rr < nurbs_dim; rr++)
          {
            if (n_x_m_x_l[rr] != count_vals[rr])
            {
              FOUR_C_THROW("not enough knots read in dim %d (%d!=NUMKNOTS=%d), nurbs_dim=%d\n", rr,
                  count_vals[rr], n_x_m_x_l[rr], nurbs_dim);
            }
          }

          // count for sanity check
          count_read++;

          continue;
        }

        //  reading of knot values if read is true and no
        // other keyword was found
        if (read)
        {
          char* endptr = nullptr;

          double dv = strtod(tmp.c_str(), &endptr);

          // count for sanity check
          count_vals[actdim]++;

          (*(patch_knots[actdim])).push_back(dv);
        }
      }  // end loop through file

      if (count_read != npatches)
      {
        FOUR_C_THROW("wasn't able to read enough patches\n");
      }
    }

    if (myrank == 0)
    {
      Core::IO::cout << " in...." << time.totalElapsedTime(true) << " secs\n";

      time.reset();
      fflush(stdout);
    }
  }



  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void InputFile::read_generic(const std::filesystem::path& top_level_file)
  {
    if (Core::Communication::my_mpi_rank(comm_) == 0)
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
                return read_yaml_content(*file_it, content_by_section_);
              }
              else
              {
                return read_dat_content(*file_it, content_by_section_);
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

    if (Core::Communication::my_mpi_rank(comm_) == 0)
    {
      // Temporarily move the sections that are not huge into a separate map.
      std::unordered_map<std::string, SectionContent> non_huge_sections;

      for (auto&& [section_name, content] : content_by_section_)
      {
        if (content.lines.size() < huge_section_threshold)
        {
          non_huge_sections[section_name] = std::move(content);
        }
      }

      Core::Communication::broadcast(non_huge_sections, 0, comm_);

      // Move the non-huge sections back into the main map.
      for (auto&& [section_name, content] : non_huge_sections)
      {
        content_by_section_[section_name] = std::move(content);
      }
    }
    else
    {
      // Other ranks receive the non-huge sections.
      Core::Communication::broadcast(content_by_section_, 0, comm_);
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



  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  bool InputFile::print_unknown_sections(std::ostream& out) const
  {
    using MapType = decltype(knownsections_);
    const MapType merged_map = Core::Communication::all_reduce<MapType>(
        knownsections_,
        [](const MapType& r, const MapType& in)
        {
          MapType result = r;
          for (const auto& [key, value] : in)
          {
            result[key] |= value;
          }
          return result;
        },
        comm_);
    const bool printout = std::any_of(
        merged_map.begin(), merged_map.end(), [](const auto& kv) { return !kv.second; });

    // now it's time to create noise on the screen
    if (printout and (Core::Communication::my_mpi_rank(get_comm()) == 0))
    {
      out << "\nERROR!"
          << "\n--------"
          << "\nThe following input file sections remained unused (obsolete or typo?):\n";
      for (const auto& [section_name, known] : knownsections_)
      {
        if (!known) out << section_name << '\n';
      }
      out << std::endl;
    }

    return printout;
  }


  void InputFile::record_section_used(const std::string& section_name)
  {
    knownsections_[section_name] = true;
  }


  void InputFile::SectionContent::pack(Communication::PackBuffer& data) const
  {
    Core::Communication::add_to_pack(data, file);
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
    FOUR_C_ASSERT(offsets.empty() || offsets.back() + lines.back().size() == raw_content.size(),
        "Offset out of bounds.");
    Core::Communication::add_to_pack(data, offsets);
  }


  void InputFile::SectionContent::unpack(Communication::UnpackBuffer& buffer)
  {
    Core::Communication::extract_from_pack(buffer, file);
    Core::Communication::extract_from_pack(buffer, raw_content);

    std::vector<std::size_t> offsets;
    Core::Communication::extract_from_pack(buffer, offsets);
    lines.clear();
    for (std::size_t i = 0; i < offsets.size(); ++i)
    {
      const char* start = raw_content.data() + offsets[i];
      const std::size_t length = (i + 1 < offsets.size() ? (offsets[i + 1] - offsets[i])
                                                         : (raw_content.size()) - offsets[i]);
      FOUR_C_ASSERT(
          start >= raw_content.data() && start + length <= raw_content.data() + raw_content.size(),
          "Line data out of bounds.");
      lines.emplace_back(start, length);
    }
  }


  std::pair<std::string, std::string> read_key_value(const std::string& line)
  {
    std::string::size_type separator_index = line.find('=');
    // The equals sign is only treated as a separator when surrounded by whitespace.
    if (separator_index != std::string::npos &&
        !(std::isspace(line[separator_index - 1]) && std::isspace(line[separator_index + 1])))
      separator_index = std::string::npos;

    // In case we didn't find an "=" separator, look for a space instead
    if (separator_index == std::string::npos)
    {
      separator_index = line.find(' ');

      if (separator_index == std::string::npos)
        FOUR_C_THROW("Line '%s' with just one word in parameter section", line.c_str());
    }

    std::string key = Core::Utils::trim(line.substr(0, separator_index));
    std::string value = Core::Utils::trim(line.substr(separator_index + 1));

    if (key.empty()) FOUR_C_THROW("Cannot get key from line '%s'", line.c_str());
    if (value.empty()) FOUR_C_THROW("Cannot get value from line '%s'", line.c_str());

    return {std::move(key), std::move(value)};
  }

}  // namespace Core::IO

FOUR_C_NAMESPACE_CLOSE
