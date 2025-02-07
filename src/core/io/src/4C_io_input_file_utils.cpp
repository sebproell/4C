// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_config_revision.hpp"

#include "4C_io_input_file_utils.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_fem_nurbs_discretization_knotvector.hpp"
#include "4C_io_input_file.hpp"
#include "4C_io_input_parameter_container.hpp"
#include "4C_io_input_spec.hpp"
#include "4C_io_pstream.hpp"
#include "4C_io_value_parser.hpp"
#include "4C_io_yaml.hpp"
#include "4C_linalg_utils_densematrix_communication.hpp"
#include "4C_utils_string.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StrUtils.hpp>
#include <Teuchos_Time.hpp>

#include <cfenv>

FOUR_C_NAMESPACE_OPEN

namespace
{
  void print_dat_impl(std::ostream& stream, const Teuchos::ParameterList& list,
      const std::string& parentname, bool comment);

  void print_documentation(std::ostream& stream, const Teuchos::ParameterEntry& entry)
  {
    // Helper function to print documentation
    std::string doc = entry.docString();
    if (!doc.empty())
    {
      Teuchos::StrUtils::printLines(stream, "// ", doc);
    }
  }


  void print_sublist(std::ostream& stream, const std::string& parentname, const std::string& name,
      const Teuchos::ParameterList& list, bool comment)
  {
    // Helper function to print a sublist
    std::string secname = parentname;
    if (!secname.empty()) secname += "/";
    secname += name;
    unsigned l = secname.length();
    stream << "--" << std::string(std::max<int>(65 - l, 0), '-');
    stream << secname << "\n";
    print_dat_impl(stream, list.sublist(name), secname, comment);
  }

  void print_parameter(std::ostream& stream, const Teuchos::ParameterEntry& entry,
      const std::string& name, const Teuchos::ParameterList& list, bool comment)
  {
    // Retrieve the parameter entry's validator (if any)
    Teuchos::RCP<const Teuchos::ParameterEntryValidator> validator = entry.validator();

    // Print comments if requested
    if (comment)
    {
      // Check if the validator has valid string values
      if (validator != Teuchos::null)
      {
        Teuchos::RCP<const Teuchos::Array<std::string>> validValues =
            validator->validStringValues();

        // If valid values exist, print them
        if (validValues != Teuchos::null)
        {
          unsigned totalLength = 0;
          // Calculate the total length of all valid values
          for (const auto& value : *validValues)
          {
            totalLength += value.length() + 1;  // Include space/comma
          }
          // Print valid values in a compact or expanded format based on total length
          if (totalLength < 74)
          {
            // Print all values in a single line, separated by commas
            stream << "//     ";
            for (auto it = validValues->begin(); it != validValues->end(); ++it)
            {
              stream << *it;
              if (std::next(it) != validValues->end())
              {
                stream << ",";  // Add a comma if it's not the last element
              }
            }
            stream << "\n";
          }
          else
          {
            // Print each value on a new line
            for (const auto& value : *validValues)
            {
              stream << "//     " << value << '\n';
            }
          }
        }
      }
    }

    // Print the parameter's name and value
    const Teuchos::any& value = entry.getAny(false);
    stream << name;
    unsigned nameLength = name.length();
    // Ensure proper spacing for alignment
    stream << std::string(std::max<int>(31 - nameLength, 0), ' ');

    // Optionally print an equal sign if needed
    if (Core::IO::need_to_print_equal_sign(list)) stream << " =";

    try
    {
      // print true/false for bool values to distinguish them from type int
      if (value.type() == typeid(bool))
      {
        stream << " " << (Teuchos::any_cast<bool>(value) ? "true" : "false") << "\n";
      }
      else
      {
        // For non-boolean types, print the value directly
        stream << " " << value << "\n";
      }
    }
    catch (const Teuchos::NonprintableTypeException&)
    {
      // Handle non-printable enum class types
      stream << value.typeName() << "\n";
    }
  }

  void print_dat_impl(std::ostream& stream, const Teuchos::ParameterList& list,
      const std::string& parentname, bool comment)
  {
    // Main loop over the parameter list that calls the helper functions to print
    // documentation, sublists or parameters:
    //
    // Iterate through the parameter list in two distinct phases to ensure proper ordering and
    // handling:
    // - **Phase 0**:
    //    Print all parameters that are not sublists. This ensures that top-level parameters
    //    are written to stream first, without any nested content interfering.
    // - **Phase 1**:
    //    Recursively handle and print all sublists. This phase is executed after all non-sublists
    //    have been processed, allowing sublists to be printed in their hierarchical order.
    //
    // By separating the iteration into these phases, we avoid issues with alphabetical
    // ordering that could cause invalid output sequences for nested lists.
    for (int iterationPhase = 0; iterationPhase < 2; ++iterationPhase)
    {
      for (auto paramIter = list.begin(); paramIter != list.end(); ++paramIter)
      {
        const Teuchos::ParameterEntry& entry = list.entry(paramIter);
        const std::string& name = list.name(paramIter);

        if ((entry.isList() && iterationPhase == 0) || (!entry.isList() && iterationPhase == 1))
        {
          continue;
        }
        if (comment)
        {
          stream << "//\n";
          print_documentation(stream, entry);
        }
        if (entry.isList())
        {
          print_sublist(stream, parentname, name, list, comment);
        }
        else
        {
          print_parameter(stream, entry, name, list, comment);
        }
      }
    }
    stream << std::endl;
  }


  void recursively_determine_sublists(const Teuchos::ParameterList& list,
      std::vector<std::pair<std::string, const Teuchos::ParameterList*>>& sublists,
      const std::string& parent_section_name = "")
  {
    for (const auto& key_value : list)
    {
      const Teuchos::ParameterEntry& entry = key_value.second;
      const std::string& name = key_value.first;
      if (entry.isList())
      {
        const std::string current_section_full_name =
            (parent_section_name == "") ? name : parent_section_name + "/" + name;

        sublists.emplace_back(current_section_full_name, &list.sublist(name));
        recursively_determine_sublists(list.sublist(name), sublists, current_section_full_name);
      }
    }
  }


  void print_metadata_yaml_parameter_list(ryml::NodeRef node, const Teuchos::ParameterList& list,
      const std::string& parent_section_name)
  {
    // prevent invalid ordering of parameters caused by alphabetical output:
    // determine all sublists first to pull them out onto the same indentation level
    std::vector<std::pair<std::string, const Teuchos::ParameterList*>> sublists;
    recursively_determine_sublists(list, sublists);



    const auto print_key_value =
        [](ryml::NodeRef parent, const std::string& key, const Teuchos::ParameterEntry& entry)
    {
      const auto to_string = [](const Teuchos::any& any)
      {
        std::stringstream s;
        s << any;
        return s.str();
      };

      auto yaml_entry = parent.append_child();
      // Serialize the key since the ParameterList gives us a temporary copy
      yaml_entry << ryml::key(key);
      yaml_entry |= ryml::MAP;

      const Teuchos::any& v = entry.getAny(false);
      yaml_entry["type"] << v.typeName();
      yaml_entry["default"] << to_string(v);

      const std::string& doc = entry.docString();
      if (doc != "")
      {
        yaml_entry["description"] << doc;
        // Add double quotes to the description to prevent YAML from interpreting special characters
        yaml_entry["description"] |= ryml::VAL_DQUO;
      }

      Teuchos::RCP<const Teuchos::ParameterEntryValidator> validator = entry.validator();
      if (validator != Teuchos::null)
      {
        Teuchos::RCP<const Teuchos::Array<std::string>> values = validator->validStringValues();
        if (values != Teuchos::null)
        {
          auto yaml_values = yaml_entry["valid options"];
          yaml_values |= ryml::SEQ;

          for (int i = 0; i < (int)values->size(); ++i)
          {
            yaml_values[i] << (*values)[i];
          }
        }
      }
    };

    for (const auto& [name, sublist] : sublists)
    {
      auto yaml_sublist = node.append_child();
      // Serialize the name since the ParameterList gives us a temporary copy
      yaml_sublist << ryml::key(name);
      yaml_sublist |= ryml::MAP;

      for (const auto& key_value : *sublist)
      {
        if (!key_value.second.isList())
          print_key_value(yaml_sublist, key_value.first, key_value.second);
      }
    }
  }
}  // namespace

void Core::IO::print_section_header(std::ostream& out, const std::string& header)
{
  constexpr std::size_t max_padding = 65ul;
  const std::size_t padding = (header.length() < max_padding) ? (max_padding - header.length()) : 0;

  out << "--" << std::string(padding, '-') << header << '\n';
}



void Core::IO::print_section(std::ostream& out, const std::string& header, const InputSpec& spec)
{
  print_section_header(out, header);
  spec.print_as_dat(out);
}


void Core::IO::print_dat(std::ostream& stream, const Teuchos::ParameterList& list, bool comment)
{
  print_dat_impl(stream, list, "", comment);
}


void Core::IO::print_metadata_yaml(std::ostream& stream, const Teuchos::ParameterList& list,
    const std::map<std::string, InputSpec>& section_specs)
{
  ryml::Tree tree = init_yaml_tree_with_exceptions();
  ryml::NodeRef root = tree.rootref();
  root |= ryml::MAP;

  {
    auto metadata = root["metadata"];
    metadata |= ryml::MAP;
    metadata["commit_hash"] << VersionControl::git_hash;
  }

  {
    auto parameters = root["parameters"];
    parameters |= ryml::MAP;
    print_metadata_yaml_parameter_list(parameters, list, "");
  }

  {
    auto sections = root["sections"];
    sections |= ryml::MAP;
    for (const auto& [name, spec] : section_specs)
    {
      auto section = sections.append_child();
      section |= ryml::MAP;
      section << ryml::key(name);
      YamlNodeRef spec_emitter{section, ""};
      spec.emit_metadata(spec_emitter);
    }
  }

  stream << tree;
}


bool Core::IO::need_to_print_equal_sign(const Teuchos::ParameterList& list)
{
  // Helper function to check if string contains a space.
  const auto string_has_space = [](const std::string& s)
  { return std::any_of(s.begin(), s.end(), [](unsigned char c) { return std::isspace(c); }); };

  return std::any_of(list.begin(), list.end(),
      [&](const auto& it)
      {
        // skip entries that are lists: they are allowed to have spaces
        if (it.second.isList()) return false;

        const std::string& name = it.key;

        // Only check string values for spaces.
        if (list.isType<std::string>(name))
        {
          const std::string& value = list.get<std::string>(name);
          return string_has_space(value) || string_has_space(name);
        }
        else
        {
          return string_has_space(name);
        }
      });
}



std::pair<std::string, std::string> Core::IO::read_key_value(const std::string& line)
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


std::vector<Core::IO::InputParameterContainer> Core::IO::read_all_lines_in_section(
    Core::IO::InputFile& input, const std::string& section, const InputSpec& spec)
{
  std::vector<Core::IO::InputParameterContainer> parsed_lines;

  for (const auto& line : input.in_section(section))
  {
    ValueParser parser{line.get_as_dat_style_string()};
    Core::IO::InputParameterContainer container;
    spec.fully_parse(parser, container);
    parsed_lines.emplace_back(std::move(container));
  }

  return parsed_lines;
}


std::pair<std::vector<Core::IO::InputParameterContainer>, std::vector<std::string>>
Core::IO::read_matching_lines_in_section(
    Core::IO::InputFile& input, const std::string& section, const IO::InputSpec& spec)
{
  std::vector<std::string> unparsed_lines;
  std::vector<Core::IO::InputParameterContainer> parsed_lines;

  for (const auto& input_line : input.in_section(section))
  {
    try
    {
      ValueParser parser{input_line.get_as_dat_style_string(),
          {.base_path = input.file_for_section(section).parent_path()}};
      InputParameterContainer container;
      spec.fully_parse(parser, container);
      parsed_lines.emplace_back(std::move(container));
    }
    catch (const Core::Exception& e)
    {
      unparsed_lines.emplace_back(input_line.get_as_dat_style_string());
    }
  }

  return {parsed_lines, unparsed_lines};
}

namespace
{
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
}  // namespace



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool Core::IO::read_parameters_in_section(
    InputFile& input, const std::string& section_name, Teuchos::ParameterList& list)
{
  if (section_name.empty()) FOUR_C_THROW("Empty section name given.");

  Teuchos::ParameterList& sublist = find_sublist(section_name, list);

  for (const auto& line : input.in_section(section_name))
  {
    const auto& [key, value] = read_key_value(std::string(line.get_as_dat_style_string()));

    add_entry(key, value, sublist);
  }

  return true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::IO::read_design(InputFile& input, const std::string& name,
    std::vector<std::vector<int>>& dobj_fenode,
    const std::function<const Core::FE::Discretization&(const std::string& name)>&
        get_discretization)
{
  std::map<int, std::set<int>> topology;

  std::string sectionname = name + "-NODE TOPOLOGY";
  std::string marker = sectionname;

  // Store lines that need special treatment
  std::vector<std::string> box_face_conditions;

  for (const auto& l : input.in_section_rank_0_only(marker))
  {
    int dobj;
    int nodeid;
    std::string nname;
    std::string dname;
    std::istringstream stream{std::string(l.get_as_dat_style_string())};
    stream >> nname;
    if (not stream)
    {
      auto s = l.get_as_dat_style_string();
      FOUR_C_THROW("Illegal line in section '%s': '%*s'", marker.c_str(), s.size(), s.data());
    }

    if (nname == "NODE")  // plain old reading of the design nodes from the .dat-file
    {
      stream >> nodeid >> dname >> dobj;
      topology[dobj - 1].insert(nodeid - 1);
    }
    else  // fancy specification of the design nodes by specifying min or max of the domain
    {     // works best on rectangular domains ;)
      // Store the specification and broadcast it to all processors
      box_face_conditions.emplace_back(l.get_as_dat_style_string());
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
      if (not stream) FOUR_C_THROW("Illegal line in section '%s': '%s'", marker.c_str(), l.data());

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

          for (const auto& line : input.in_section_rank_0_only(dommarker))
          {
            std::istringstream t{std::string{line.get_as_dat_style_string()}};
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

      constexpr double tolerance_n = 1.0e-14;
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
    if (max_num_dobj >= static_cast<int>(dobj_fenode.size())) dobj_fenode.resize(max_num_dobj + 1);
    // copy all design object entries
    for (auto& topo : topology)
    {
      // we copy from a std::set, thus the gids are sorted
      dobj_fenode[topo.first].reserve(topo.second.size());
      dobj_fenode[topo.first].assign(topo.second.begin(), topo.second.end());
    }
  }
}


void Core::IO::read_knots(InputFile& input, const std::string& name,
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
    for (const auto& line : input.in_section(sectionname))
    {
      // count number of patches in knotvector section of
      // this discretisation
      {
        std::string::size_type loc;
        std::istringstream file{std::string{line.get_as_dat_style_string()}};
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
    for (const auto& line : input.in_section(sectionname))
    {
      std::istringstream file{std::string{line.get_as_dat_style_string()}};
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

FOUR_C_NAMESPACE_CLOSE
