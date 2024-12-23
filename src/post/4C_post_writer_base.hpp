// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_POST_WRITER_BASE_HPP
#define FOUR_C_POST_WRITER_BASE_HPP

#include "4C_config.hpp"

#include "4C_post_filter_base.hpp"

#include <fstream>
#include <map>
#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  template <typename T>
  class MultiVector;
}

class PostField;

//! Special writer class that is used to invoke a particular method in the output writer
//! e.g.
struct SpecialFieldInterface
{
  virtual ~SpecialFieldInterface() = default;
  virtual std::vector<int> num_df_map() = 0;

  virtual void operator()(std::vector<std::shared_ptr<std::ofstream>>& files, PostResult& result,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string& groupname, const std::vector<std::string>& names) = 0;
};


//! Base class for various output writers that use generic interfaces (Ensight, VTU)
class PostWriterBase
{
 public:
  //! constructor. initializes the writer to a certain field
  PostWriterBase(PostField* field, const std::string& filename);

  //! destructor
  virtual ~PostWriterBase() = default;
  //! return the field specified at construction
  PostField* get_field()
  {
    FOUR_C_ASSERT(field_ != nullptr, "No field has been set");
    return field_;
  }

  const PostField* get_field() const
  {
    FOUR_C_ASSERT(field_ != nullptr, "No field has been set");
    return field_;
  }

  virtual void write_files(PostFilterBase& filter) = 0;

  /*!
   \brief write all time steps of a result

   Write results. The results are taken from a reconstructed
   Core::LinAlg::Vector<double>. In many cases this vector will contain just one
   variable (displacements) and thus is easy to write as a whole. At
   other times, however, there is more than one result (velocity,
   pressure) and we want to write just one part of it. So we have to
   specify which part.
   */
  virtual void write_result(
      const std::string groupname,  ///< name of the result group in the control file
      const std::string name,       ///< name of the result to be written
      const ResultType restype,     ///< type of the result to be written (nodal-/element-based)
      const int numdf,              ///< number of dofs per node to this result
      const int from = 0,           ///< start position of values in nodes
      const bool fillzeros = false  ///< zeros are filled to ensight file when no data is available
      ) = 0;

  /*!
   \brief write all time steps of a result in one time step

   Write results. The results are taken from a reconstructed
   Core::LinAlg::Vector<double>. In many cases this vector will contain just one
   variable (displacements) and thus is easy to write as a whole. At
   other times, however, there is more than one result (velocity,
   pressure) and we want to write just one part of it. So we have to
   specify which part. Currently file continuation is not supported
   because Paraview is not able to load it due to some weird wild card
   issue.
   */
  virtual void write_result_one_time_step(PostResult& result,  ///< result group in the control file
      const std::string groupname,  ///< name of the result group in the control file
      const std::string name,       ///< name of the result to be written
      const ResultType restype,     ///< type of the result to be written (nodal-/element-based)
      const int numdf,              ///< number of dofs per node to this result
      bool firststep,               ///< bool whether this is the first time step
      bool laststep,                ///< bool whether this is the last time step
      const int from = 0            ///< start position of values in nodes
      ) = 0;

  /*!
   \brief write a particular variable to file

   Write results. Some variables need interaction with the post filter,
   e.g. structural stresses that do some element computations before output.
   To allow for a generic interface, the calling site needs to supply a
   class derived from SpecialFieldInterface that knows which function to call.
   */
  virtual void write_special_field(SpecialFieldInterface& special,
      PostResult& result,  ///< result group in the control file
      const ResultType restype, const std::string& groupname,
      const std::vector<std::string>& fieldnames, const std::string& outinfo) = 0;

  /*!
   \brief Write one step of a nodal result
   */
  virtual void write_nodal_result_step(std::ofstream& file,
      const std::shared_ptr<Core::LinAlg::MultiVector<double>>& data,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string& groupname, const std::string& name, const int numdf) = 0;

  /*!
   \brief Write one step of an element result
   */
  virtual void write_element_result_step(std::ofstream& file,
      const std::shared_ptr<Core::LinAlg::MultiVector<double>>& data,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string& groupname, const std::string& name, const int numdf, const int from) = 0;

 protected:
  PostField* field_;
  std::string filename_;
  unsigned int myrank_;   ///< global processor id
  unsigned int numproc_;  ///< number of processors
};

FOUR_C_NAMESPACE_CLOSE

#endif
