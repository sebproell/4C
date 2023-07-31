/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief Write visualization output in vtk/vtp format at runtime

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/

/* headers */
#include "baci_io_runtime_vtp_writer.H"

#include "baci_io_vtp_writer.H"

#include "baci_utils_exceptions.H"


#include <iomanip>
#include <sstream>
#include <vector>

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
RuntimeVtpWriter::RuntimeVtpWriter() { vtp_writer_ = Teuchos::rcp(new VtpWriter()); }

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void RuntimeVtpWriter::Initialize(unsigned int myrank, unsigned int num_processors,
    unsigned int max_number_timesteps_to_be_written,
    const std::string& path_existing_output_directory, const std::string& simulation_name,
    const std::string& geometry_name, const std::string& restart_name, double restart_time,
    bool write_binary_output)
{
  vtp_writer_->Initialize(myrank, num_processors, max_number_timesteps_to_be_written,
      path_existing_output_directory, (simulation_name + "-vtp-files"), geometry_name, restart_name,
      restart_time, write_binary_output);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void RuntimeVtpWriter::SetupForNewTimeStepAndGeometry(
    double time, unsigned int timestep, const std::string& geometryname)
{
  vtp_writer_->ResetTimeAndTimeStep(time, timestep);

  vtp_writer_->ResetGeometryName(geometryname);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void RuntimeVtpWriter::ResetGeometry(const std::vector<double>& point_coordinates)
{
  point_coordinates_ = point_coordinates;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
// Todo template <typename T>, double, int,
void RuntimeVtpWriter::AppendVisualizationPointDataVector(const std::vector<double>& datavalues,
    const unsigned int num_components_per_point, const std::string& dataname)
{
  point_data_vectors_[dataname] = std::make_pair(datavalues, num_components_per_point);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
std::vector<double>& RuntimeVtpWriter::GetMutablePointCoordinateVector()
{
  return point_coordinates_;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void RuntimeVtpWriter::WriteFiles()
{
  vtp_writer_->InitializeVtkFileStreamsForNewGeometryAndOrTimeStep();

  vtp_writer_->WriteVtkHeaders();

  vtp_writer_->WriteVtkTimeAndOrCycle();

  vtp_writer_->WriteGeometryPolyData(point_coordinates_);

  WriteAllPointDataVectorsToFiles();

  vtp_writer_->WriteVtkFooters();
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void RuntimeVtpWriter::WriteAllPointDataVectorsToFiles()
{
  /* note: the implementation of the vtu writer expects to first write all point data and then
   *       write all cell data */

  // write all collected point data vectors to files
  typedef std::map<std::string, std::pair<std::vector<double>, unsigned int>> point_data_map;

  for (point_data_map::const_iterator point_data_iter = point_data_vectors_.begin();
       point_data_iter != point_data_vectors_.end(); ++point_data_iter)
  {
    vtp_writer_->WritePointDataVector(
        point_data_iter->second.first, point_data_iter->second.second, point_data_iter->first);
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void RuntimeVtpWriter::WriteCollectionFileOfAllWrittenFiles(const std::string& collectionfilename)
{
  vtp_writer_->WriteVtkCollectionFileForAllWrittenMasterFiles(collectionfilename);
}