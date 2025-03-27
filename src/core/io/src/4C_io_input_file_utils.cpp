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

void Core::IO::print_dat(
    std::ostream& stream, const std::map<std::string, Core::IO::InputSpec>& map)
{
  for (const auto& [name, spec] : map)
  {
    print_section_header(stream, name);
    spec.print_as_dat(stream);
  }
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
      FOUR_C_THROW("Line '{}' with just one word in parameter section", line.c_str());
  }

  std::string key = Core::Utils::trim(line.substr(0, separator_index));
  std::string value = Core::Utils::trim(line.substr(separator_index + 1));

  if (key.empty()) FOUR_C_THROW("Cannot get key from line '{}'", line.c_str());
  if (value.empty()) FOUR_C_THROW("Cannot get value from line '{}'", line.c_str());

  return {std::move(key), std::move(value)};
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

}  // namespace



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::IO::read_parameters_in_section(
    InputFile& input, const std::string& section_name, Teuchos::ParameterList& list)
{
  if (section_name.empty()) FOUR_C_THROW("Empty section name given.");

  InputParameterContainer container;
  input.match_section(section_name, container);

  FOUR_C_ASSERT(container.has_group(section_name),
      "Internal error: group '{}' not found in container.", section_name.c_str());

  container.group(section_name).to_teuchos_parameter_list(find_sublist(section_name, list));
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
      FOUR_C_THROW("Illegal line in section '{}': '{}'", marker.c_str(), s);
    }

    if (nname == "NODE")  // plain old reading of the design nodes from the input file
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
      if (not stream) FOUR_C_THROW("Illegal line in section '{}': '{}'", marker.c_str(), l.data());

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
        FOUR_C_THROW("Illegal line in section '{}': '{}'", marker.c_str(), l.data());
      }

      const Core::FE::Discretization& actdis = get_discretization(disname);

      std::vector<double> box_specifications;
      {
        for (int init = 0; init < 9; ++init) box_specifications.push_back(0.0);
        if (Core::Communication::my_mpi_rank(input.get_comm()) == 0)  // Reading is done by proc 0
        {
          // get original domain section from the input file
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
        FOUR_C_THROW("Illegal line in section '{}': '{}'\n{} found, where {} was expected",
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
  const int myrank = Core::Communication::my_mpi_rank(input.get_comm());
  Teuchos::Time time("", true);

  std::string field;
  if (name == "fluid" || name == "xfluid" || name == "porofluid")
    field = "FLUID";
  else if (name == "structure")
    field = "STRUCTURE";
  else if (name == "ale")
    field = "ALE";
  else if (name == "scatra")
    field = "TRANSPORT";
  else if (name == "thermo")
    field = "THERMO";
  else if (name == "scatra_micro")
    field = "TRANSPORT2";
  else
    FOUR_C_THROW("Unknown discretization name for knotvector input.");

  const std::string sectionname = field + " KNOTVECTORS";

  if (myrank == 0)
  {
    Core::IO::cout << "Reading knot vectors for " << name << " discretization :\n";
    fflush(stdout);
  }

  int npatches = 0;
  int nurbs_dim = 0;

  for (const auto& line : input.in_section(sectionname))
  {
    std::string line_str{line.get_as_dat_style_string()};
    std::istringstream file{line_str};
    std::string tmp;
    file >> tmp;

    if (tmp == "NURBS_DIMENSION")
    {
      file >> nurbs_dim;
    }
    else if (tmp == "ID")
    {
      npatches++;
    }
  }

  if (myrank == 0)
  {
    printf("                        %8d patches", npatches);
    fflush(stdout);
  }

  disknots = std::make_shared<Core::FE::Nurbs::Knotvector>(nurbs_dim, npatches);
  if (!disknots) FOUR_C_THROW("disknots should have been allocated before");

  std::vector<std::shared_ptr<std::vector<double>>> patch_knots(nurbs_dim);
  std::vector<int> n_x_m_x_l(nurbs_dim), degree(nurbs_dim), count_vals(nurbs_dim);
  std::vector<std::string> knotvectortype(nurbs_dim);

  bool read = false;
  int npatch = 0, actdim = -1, count_read = 0;

  for (const auto& line : input.in_section(sectionname))
  {
    std::string line_str{line.get_as_dat_style_string()};
    std::istringstream file{line_str};
    std::string tmp;
    file >> tmp;

    if (tmp == "BEGIN")
    {
      read = true;
      actdim = -1;
      for (auto& knots : patch_knots) knots = std::make_shared<std::vector<double>>();
      std::fill(count_vals.begin(), count_vals.end(), 0);
    }
    else if (tmp == "ID")
    {
      file >> npatch;
      npatch--;
    }
    else if (tmp == "NUMKNOTS")
    {
      file >> n_x_m_x_l[++actdim];
    }
    else if (tmp == "DEGREE")
    {
      file >> degree[actdim];
    }
    else if (tmp == "TYPE")
    {
      file >> knotvectortype[actdim];
    }
    else if (tmp == "END")
    {
      for (int rr = 0; rr < nurbs_dim; ++rr)
      {
        disknots->set_knots(
            rr, npatch, degree[rr], n_x_m_x_l[rr], knotvectortype[rr], patch_knots[rr]);
      }
      read = false;
      for (int rr = 0; rr < nurbs_dim; ++rr)
      {
        if (n_x_m_x_l[rr] != count_vals[rr])
        {
          FOUR_C_THROW("not enough knots read in dim {} ({}!=NUMKNOTS={}), nurbs_dim={}\n", rr,
              count_vals[rr], n_x_m_x_l[rr], nurbs_dim);
        }
      }
      count_read++;
    }
    else if (read)
    {
      patch_knots[actdim]->push_back(std::stod(tmp));
      count_vals[actdim]++;
    }
  }

  if (count_read != npatches)
  {
    FOUR_C_THROW("wasn't able to read enough patches\n");
  }

  if (myrank == 0)
  {
    Core::IO::cout << " in...." << time.totalElapsedTime(true) << " secs\n";
    time.reset();
    fflush(stdout);
  }
}

FOUR_C_NAMESPACE_CLOSE
