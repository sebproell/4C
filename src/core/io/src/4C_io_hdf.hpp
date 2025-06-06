// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_IO_HDF_HPP
#define FOUR_C_IO_HDF_HPP

#include "4C_config.hpp"

#include "4C_linalg_map.hpp"
#include "4C_linalg_vector.hpp"

#include <hdf5.h>
#include <hdf5_hl.h>

#include <algorithm>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::IO
{
  /*!
    \brief Helper class that handles the HDF5 files while reading.

    After parallel output we need parallel input. Each processor writes
    its own file. On input therefore each processor needs to read some
    files (zero to many, depending on the number of processors that
    wrote and read.) Additionally each processor might use several files
    (for several time steps, as it happens with restart.) This class
    handles the basic HDF5 file access.

  */
  class HDFReader
  {
   public:
    //! Construct with the directory where the files can be found.
    explicit HDFReader(std::string dir);
    ~HDFReader();

    //! Open a new set of input files.
    /*!
      With num_output_procs==1 this function opens the result data
      file with name basename. If num_output_procs>1 it opens the result
      files of all processors, by appending .p<proc_num> to the
      basename.
    */
    void open(std::string basename, int num_output_procs, int new_proc_num, int my_id);
    //!
    void close();

    /*!
     * \brief read the packed element data from the mesh files
     * \note  this function should only be called while the HDFReader reads
     *        mesh files
     */
    std::shared_ptr<std::vector<char>> read_element_data(
        int step, int new_proc_num, int my_id) const;

    //!
    /*!
     * \brief read the packed nodal data from the mesh files
     * \note: this function should only be called while the HDFReader reads
     *         mesh files
     */
    std::shared_ptr<std::vector<char>> read_node_data(int step, int new_proc_num, int my_id) const;

    /*!
     * \brief reads the packed periodic boundary condition data from the mesh files
     *
     *
     * \note this function should only be called when the HDFReader opened
     *       the mesh files
     */
    std::shared_ptr<std::vector<char>> read_condition(
        const int step, const int new_proc_num, const int my_id, const std::string condname) const;

    /*!
    //  \brief reads the packed knotvector data from the mesh files
    //
    //    \param   step (i)
    //
    //      \return  The whole knotvector data in a char vector
    //      */
    std::shared_ptr<std::vector<char>> read_knotvector(const int step) const;


    //! read an Core::LinAlg::MultiVector<double> from the result files
    /*!
      Right now an Core::LinAlg::MultiVector<double> has to be read along with its map. Thus
      we read an integer and a double array here.

      \note If columns==1, we create an Core::LinAlg::Vector<double>.

      \param id_path      (in): hdf5 path to map array (from control file)
      \param value_path   (in): hdf5 path to value array (from control file)
      \param columns      (in): number of vector columns
      \param Comm         (in): the communicator
     */
    std::shared_ptr<Core::LinAlg::MultiVector<double>> read_result_data(
        std::string id_path, std::string value_path, int columns, MPI_Comm Comm) const;

    //! read a std::vector<char> from the result files
    /*!
      Right now a std::vector<char> has to be read along with its map. Thus
      we read an integer and a double array here.

      \param id_path      (in): hdf5 path to map array (from control file)
      \param value_path   (in): hdf5 path to value array (from control file)
      \param columns      (in): number of vector columns
      \param Comm         (in): the communicator
      \param elemap      (out): element map
     */
    std::shared_ptr<std::vector<char>> read_result_data_vec_char(std::string id_path,
        std::string value_path, int columns, MPI_Comm Comm,
        std::shared_ptr<Core::LinAlg::Map>& elemap) const;

    std::shared_ptr<std::vector<char>> read_char_vector(
        std::string value_path, MPI_Comm Comm) const;

    std::shared_ptr<std::vector<double>> read_double_vector(std::string path) const
    {
      std::vector<int> length;
      std::shared_ptr<std::vector<double>> values = read_double_data(path, 0, 1, length);
      return values;
    };

    std::shared_ptr<std::vector<int>> read_int_vector(std::string path) const
    {
      std::vector<int> length;
      std::shared_ptr<std::vector<int>> values = read_int_data(path, 0, 1);
      return values;
    };

   private:
    /// reads the dataset 'path' in all the files in the range [start,end)
    /*!
      Here we finally loop all the files the local processor has to read.
      returns all the data in one vector. The data is assumed to by of type char.
     */
    std::shared_ptr<std::vector<char>> read_char_data(std::string path, int start, int end) const;

    /// reads the dataset 'path' in all the files in the range [start,end)
    /*!
      Here we finally loop all the files the local processor has to read.
      returns all the data in one vector<int>
    */
    std::shared_ptr<std::vector<int>> read_int_data(std::string path, int start, int end) const;

    /// reads the dataset 'path' in all the files in the range [start,end)
    /*!
      Here we finally loop all the files the local processor has to
      read.
    */
    std::shared_ptr<std::vector<double>> read_double_data(
        std::string path, int start, int end, std::vector<int>& lengths) const;

    //! Figure out which subset of files this process needs to read
    /*!
      In a parallel run each processor writes one file. If we are to
      read this set of files on with a different number of processors,
      we need to know how many files each one has to read. One file will
      always be read by one processor. The reader already knows how many
      files there are (how many processors wrote this set.)

      We assign a consecutive range of files to each processor. Linearly
      distributed. After the read we can redistribute. :)

      \param new_proc_num (in): numbers of processors that read
      \param my_id        (in): my rank. id. processor number.
      \param start       (out): first file I have to read
      \param end         (out): first file I do not read
     */
    void calculate_range(int new_proc_num, int my_id, int& start, int& end) const;

    //! the names of the files opened here
    std::vector<std::string> filenames_;

    //! HDF5 ids of the open files
    std::vector<hid_t> files_;

    //! directory where to find the files
    std::string input_dir_;

    //! number of processors that wrote this set of files
    int num_output_proc_;

    //! file access property list for HDF5 files
    hid_t h5_plist_;
  };
}  // namespace Core::IO


FOUR_C_NAMESPACE_CLOSE

#endif
