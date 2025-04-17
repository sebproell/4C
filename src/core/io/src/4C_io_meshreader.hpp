// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_IO_MESHREADER_HPP
#define FOUR_C_IO_MESHREADER_HPP

#include "4C_config.hpp"

#include "4C_io_domainreader.hpp"
#include "4C_io_elementreader.hpp"
#include "4C_io_exodus.hpp"
#include "4C_linalg_graph.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::IO
{
  class InputFile;

  namespace Internal
  {
    /**
     * Support class to read a mesh from an exodus file.
     */
    struct ExodusReader
    {
      /**
       *The discretization that should be filled with the information from the exodus file.
       */
      Core::FE::Discretization& target_discretization;

      /**
       * The section in the input file that has the necessary data for the reader (e.g. the file
       * name).
       */
      std::string section_name;

      /**
       * The actual exodus mesh object. This is only created on rank 0.
       */
      std::unique_ptr<Exodus::Mesh> mesh_on_rank_zero{};
    };
  }  // namespace Internal

  /*!
    \brief helper class to read a mesh

    This is an interface for handling node, element and domain readers
    and to set up a discretization from an input file which is fill_complete().
   */
  class MeshReader
  {
   public:
    /**
     * Additional parameters that governt the reading process.
     */
    struct MeshReaderParameters
    {
      /**
       * How to partition then mesh among processes.
       */
      Teuchos::ParameterList mesh_partitioning_parameters;

      /**
       * Geometric search parameters for certain partitiong methods.
       */
      Teuchos::ParameterList geometric_search_parameters;

      /**
       * General verbosity settings and I/O parameters.
       */
      Teuchos::ParameterList io_parameters;
    };

    /**
     * Construct a mesh reader. Read nodes from the given @p input under section
     * @p node_section_name.
     */
    MeshReader(Core::IO::InputFile& input, std::string node_section_name,
        MeshReaderParameters parameters = {});

    /// add an element reader for each discretization
    /*!
      Each discretization needs its own ElementReader. These readers
      have to be registered at the MeshReader.

      \param er (i) a reader of one discretization that uses (a fraction of) our nodes
     */
    void add_element_reader(const ElementReader& er) { element_readers_.emplace_back(er); }

    /*!
     * \brief Adds the selected reader to this meshreader
     *
     * This is a version without specific elementtypes. It just calls the full
     * version with a dummy set.
     *
     * \param dis            [in] This discretization will be passed on
     * \param input          [in] The input file.
     * \param sectionname    [in] The section name in the input file.
     */
    void add_advanced_reader(
        std::shared_ptr<Core::FE::Discretization> dis, const std::string& sectionname);

    /// do the actual reading
    /*!
      This method contains the whole machinery. The reading consists of
      three steps:

      - Reading and distributing elements using each registered
        ElementReader. This includes creating the connectivity graph,
        building the node row and column maps.

      - Reading and distributing all nodes. Each node gets assigned to
        its discretization.

      - Finalizing the discretizations.

      Actually most of the work gets done by the ElementReader. The
      reading of both elements and nodes happens in blocks on processor
      0. After each block read the discretizations are redistributed.

     */
    void read_and_partition();

    /**
     * Get MPI communicator of this mesh reader.
     */
    MPI_Comm get_comm() const;

    /**
     * Access the exodus mesh on rank 0. This is only available if an exodus file was actually read.
     * On ranks other than 0, this will always return a nullptr.
     */
    const Exodus::Mesh* get_exodus_mesh_on_rank_zero() const;

   private:
    /*!
    \brief Read pre-generated mesh from input file and generate related FillCompleted()
    discretizations

    \param[in/out] max_node_id Maximum node id in a given discretization. To be used as global
                               offset to start node numbering (based on already existing nodes)
    */
    void read_mesh_from_dat_file(int& max_node_id);

    //! Read all mesh information from exodus files.
    void read_mesh_from_exodus();

    /*!
    \brief Rebalance discretizations built in read_mesh_from_dat_file()
    */
    void rebalance() const;

    /*!
    \brief Create inline mesh

    Ask the DomainReader to process input data and create a box shaped mesh at runtime without
    having to read or process nodes and elements from the input file.

    \param[in/out] max_node_id Maximum node id in a given discretization. To be used as global
                               offset to start node numbering (based on already existing nodes)
    */
    void create_inline_mesh(int& max_node_id) const;

    /// my comm
    MPI_Comm comm_;

    //! graphs of each discretization
    std::vector<std::shared_ptr<const Core::LinAlg::Graph>> graph_;

    /// my element readers
    std::vector<ElementReader> element_readers_;

    /// my domain readers
    std::vector<DomainReader> domain_readers_;

    std::vector<Internal::ExodusReader> exodus_readers_;

    /// Input file contents
    Core::IO::InputFile& input_;

    /// The name of the section under which we will read the nodes.
    std::string node_section_name_;

    /// Additional parameters for reading meshes.
    MeshReaderParameters parameters_;

    /// The discretizations to be filled. The key is an identifier for the sections in the input.
    /// Multiple discretizations might be filled from the same section.
    std::vector<std::pair<std::string, std::shared_ptr<Core::FE::Discretization>>>
        target_discretizations_;
  };
}  // namespace Core::IO

FOUR_C_NAMESPACE_CLOSE

#endif
