/*----------------------------------------------------------------------*/
/*! \file

\brief Binning strategy for neighborhood search

\level 2

*----------------------------------------------------------------------*/

#include "4C_binstrategy.hpp"

#include "4C_binstrategy_meshfree_multibin.hpp"
#include "4C_binstrategy_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_dofset_independent.hpp"
#include "4C_fem_geometry_intersection_math.hpp"
#include "4C_fem_geometry_periodic_boundingbox.hpp"
#include "4C_inpar_binningstrategy.hpp"
#include "4C_io.hpp"
#include "4C_io_pstream.hpp"
#include "4C_linalg_utils_densematrix_communication.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_mortar_element.hpp"
#include "4C_mortar_node.hpp"
#include "4C_rebalance_graph_based.hpp"
#include "4C_rebalance_print.hpp"

#include <Teuchos_TimeMonitor.hpp>

#include <utility>

FOUR_C_NAMESPACE_OPEN


Core::Binstrategy::BinningStrategy::BinningStrategy(const Teuchos::ParameterList& binning_params,
    Teuchos::RCP<Core::IO::OutputControl> output_control, const Epetra_Comm& comm,
    const int my_rank,
    const std::function<Utils::SpecialElement(const Core::Elements::Element* element)>
        element_filter,
    const std::function<double(const Core::Elements::Element* element)> rigid_sphere_radius,
    const std::function<Core::Nodes::Node const*(Core::Nodes::Node const* node)>
        correct_beam_center_node,
    const std::vector<Teuchos::RCP<Core::FE::Discretization>>& discret,
    std::vector<Teuchos::RCP<const Epetra_Vector>> disnp)
    : bin_size_lower_bound_(binning_params.get<double>("BIN_SIZE_LOWER_BOUND")),
      deforming_simulation_domain_handler_(Teuchos::null),
      writebinstype_(
          Core::UTILS::IntegralValue<Inpar::BINSTRATEGY::Writebins>(binning_params, ("WRITEBINS"))),
      myrank_(my_rank),
      comm_(comm.Clone()),
      element_filter_(std::move(element_filter)),
      rigid_sphere_radius_(std::move(rigid_sphere_radius)),
      correct_beam_center_node_(std::move(correct_beam_center_node))
{
  // create binning discretization
  bindis_ = Teuchos::rcp(new Core::FE::Discretization("binning", comm_, 3));

  // create discretization writer
  auto spatial_approximation_type = Teuchos::getIntegralValue<Core::FE::ShapeFunctionType>(
      binning_params, "spatial_approximation_type");
  bindis_->SetWriter(Teuchos::rcp(
      new Core::IO::DiscretizationWriter(bindis_, output_control, spatial_approximation_type)));
  bindis_->fill_complete(false, false, false);

  visbindis_ = Teuchos::rcp(new Core::FE::Discretization("bins", comm_, 3));
  // create discretization writer
  visbindis_->SetWriter(Teuchos::rcp(
      new Core::IO::DiscretizationWriter(visbindis_, output_control, spatial_approximation_type)));

  // try to read valid input
  domain_bounding_box_corner_positions_.PutScalar(1.0e12);
  // get bounding box specified in the input file
  std::istringstream domain_bounding_box_stream(
      Teuchos::getNumericStringParameter(binning_params, "DOMAINBOUNDINGBOX"));
  for (int col = 0; col < 2; ++col)
  {
    for (int row = 0; row < 3; ++row)
    {
      std::string value;
      if (domain_bounding_box_stream >> value)
      {
        const double doublevalue = std::atof(value.c_str());
        domain_bounding_box_corner_positions_(row, col) = doublevalue;
      }
      else
        FOUR_C_THROW(
            "specify six values for bounding box in three dimensional problem. Fix input file");
    }
  }

  std::istringstream binstream(Teuchos::getNumericStringParameter(binning_params, "BIN_PER_DIR"));
  for (int idim = 0; idim < 3; ++idim)
  {
    std::string val;
    if (binstream >> val)
    {
      int intval = std::atoi(val.c_str());
      if (intval > 0) bin_per_dir_[idim] = intval;
    }
    else
    {
      FOUR_C_THROW(
          "You need to specify three figures for BIN_PER_DIR in input file for three dimensional "
          "problem. ");
    }
  }

  // check input
  bool feasiblebininput = true;

  for (int idim = 0; idim < 3; ++idim)
    if (bin_per_dir_[idim] < 1) feasiblebininput = false;

  // safety check
  if (feasiblebininput and bin_size_lower_bound_ > 0.0)
    FOUR_C_THROW("Choose either bin_size_lower_bound_ or binsperdir to specify binning domain.");

  boundaryrowbins_.clear();
  boundarycolbins_.clear();

  // init vectors for function calls
  if (disnp.size() == 0) disnp.resize(discret.size(), Teuchos::null);

  bool feasibleboxinput = true;
  for (int col = 0; col < 2; ++col)
    for (int row = 0; row < 3; ++row)
      if (domain_bounding_box_corner_positions_(row, col) > 1.0e11) feasibleboxinput = false;

  if (not feasibleboxinput)
  {
    FOUR_C_THROW_UNLESS(discret.size() != 0, "We need a discretization at this point.");
    compute_min_binning_domain_containing_all_elements_of_multiple_discrets(
        discret, disnp, domain_bounding_box_corner_positions_, bin_size_lower_bound_ < 0.0);
  }
  else if (bin_size_lower_bound_ < 0.0)
  {
    FOUR_C_THROW_UNLESS(discret.size() != 0, "We need a discretization at this point.");
    bin_size_lower_bound_ =
        compute_lower_bound_for_bin_size_as_max_edge_length_of_aabb_of_largest_ele(discret, disnp);
  }

  // init binning domain edge length
  for (int idim = 0; idim < 3; ++idim)
    edge_length_binning_domain_[idim] = domain_bounding_box_corner_positions_(idim, 1) -
                                        domain_bounding_box_corner_positions_(idim, 0);

  // create bins
  create_bins_based_on_bin_size_lower_bound_and_binning_domain_dimensions();

  // build periodic boundary condition
  BuildPeriodicBC(binning_params);
}

void Core::Binstrategy::BinningStrategy::GidsInijkRange(
    const int* ijk_range, std::set<int>& binIds, bool checkexistence) const
{
  if (checkexistence and bindis_ == Teuchos::null)
    FOUR_C_THROW("particle discretization is not set up correctly");

  for (int i = ijk_range[0]; i <= ijk_range[1]; ++i)
  {
    for (int j = ijk_range[2]; j <= ijk_range[3]; ++j)
    {
      for (int k = ijk_range[4]; k <= ijk_range[5]; ++k)
      {
        std::array ijk = {i, j, k};

        const int gid = ConvertijkToGid(ijk.data());
        if (gid != -1)
        {
          if (checkexistence)
          {
            if (bindis_->HaveGlobalElement(gid)) binIds.insert(gid);
          }
          else
          {
            binIds.insert(gid);
          }
        }
      }  // end for int k
    }    // end for int j
  }      // end for int i
}

void Core::Binstrategy::BinningStrategy::GidsInijkRange(
    const int* ijk_range, std::vector<int>& binIds, bool checkexistence) const
{
  if (checkexistence and bindis_ == Teuchos::null)
    FOUR_C_THROW("particle discretization is not set up correctly");

  for (int i = ijk_range[0]; i <= ijk_range[1]; ++i)
  {
    for (int j = ijk_range[2]; j <= ijk_range[3]; ++j)
    {
      for (int k = ijk_range[4]; k <= ijk_range[5]; ++k)
      {
        std::array ijk = {i, j, k};

        const int gid = ConvertijkToGid(ijk.data());
        if (gid != -1)
        {
          if (checkexistence)
          {
            if (bindis_->HaveGlobalElement(gid)) binIds.push_back(gid);
          }
          else
          {
            binIds.push_back(gid);
          }
        }
      }  // end for int k
    }    // end for int j
  }      // end for int i
}

int Core::Binstrategy::BinningStrategy::get_number_of_bins_inijk_range(int const ijk_range[6]) const
{
  return ((ijk_range[1] - ijk_range[0] + 1) * (ijk_range[3] - ijk_range[2] + 1) *
          (ijk_range[5] - ijk_range[4] + 1));
}

int Core::Binstrategy::BinningStrategy::ConvertijkToGid(int* ijk) const
{
  // might need to modify ijk connectivity in the presence of periodic boundary conditions
  if (havepbc_)
  {
    for (unsigned idim = 0; idim < 3; ++idim)
    {
      if (pbconoff_[idim])
      {
        if (ijk[idim] == -1)
          ijk[idim] = bin_per_dir_[idim] - 1;
        else if (ijk[idim] == bin_per_dir_[idim])
          ijk[idim] = 0;
      }
    }
  }

  // given ijk is outside of XAABB
  if (ijk[0] < 0 || ijk[1] < 0 || ijk[2] < 0 || ijk[0] >= bin_per_dir_[0] ||
      ijk[1] >= bin_per_dir_[1] || ijk[2] >= bin_per_dir_[2])
    return -1;

  return ijk[0] + ijk[1] * id_calc_bin_per_dir_[0] +
         ijk[2] * id_calc_bin_per_dir_[0] * id_calc_bin_per_dir_[1];
}

void Core::Binstrategy::BinningStrategy::ConvertGidToijk(const int gid, int* ijk) const
{
  // in order to efficiently compute the ijk triple from a given bin id,
  // use of the shift operator is made (right shift by one equals division by two)
  // therefore it is necessary that the number of bins per direction is
  // divisible by 2
  // (shift operation costs one cycle vs division (or modulo) costs 20--40 cycles on cpu)
  // Hence, two different number of bins per direction are needed
  // one for the used domain and the other one for converting gid <-> ijk

  // example: 2^n = bin_per_dir
  // example: gid >> n = (int)gid/bin_per_dir

  ijk[2] = gid >> (id_calc_exp_bin_per_dir_[0] + id_calc_exp_bin_per_dir_[1]);

  const int tmp = gid - ijk[2] * id_calc_bin_per_dir_[0] * id_calc_bin_per_dir_[1];

  ijk[1] = tmp >> id_calc_exp_bin_per_dir_[0];

  ijk[0] = tmp - ijk[1] * id_calc_bin_per_dir_[0];

  // alternative method - more expensive but only based on integer operations:
  //  {
  //    const int tmp1 = gid % (id_calc_bin_per_dir_[0] * id_calc_bin_per_dir_[1]);
  //    // compute i
  //    ijk[0] = tmp1 % id_calc_bin_per_dir_[0];
  //    // compute j
  //    ijk[1] = (tmp1 - ijk[0]) / id_calc_bin_per_dir_[0];
  //    // compute k
  //    ijk[2] = (gid - ijk[0] - ijk[1]*id_calc_bin_per_dir_[0]) / (id_calc_bin_per_dir_[0] *
  //    id_calc_bin_per_dir_[1]);
  //  }

  // found ijk is outside of XAABB
  if (ijk[0] < 0 || ijk[1] < 0 || ijk[2] < 0 || ijk[0] >= bin_per_dir_[0] ||
      ijk[1] >= bin_per_dir_[1] || ijk[2] >= bin_per_dir_[2])
    FOUR_C_THROW("ijk (%d %d %d) for given gid: %d is outside of range (bin per dir: %d %d %d)",
        ijk[0], ijk[1], ijk[2], gid, bin_per_dir_[0], bin_per_dir_[1], bin_per_dir_[2]);
}


int Core::Binstrategy::BinningStrategy::ConvertPosToGid(const double* pos) const
{
  int ijk[3];
  double pos_ud[3];
  if (deforming_simulation_domain_handler_ != Teuchos::null)
  {
    deforming_simulation_domain_handler_->transform_from_global_to_undeformed_bounding_box_system(
        pos, pos_ud);
    for (int dim = 0; dim < 3; ++dim)
      ijk[dim] = static_cast<int>(std::floor(
          (pos_ud[dim] - domain_bounding_box_corner_positions_(dim, 0)) * inv_bin_size_[dim]));
  }
  else
  {
    for (int dim = 0; dim < 3; ++dim)
      ijk[dim] = static_cast<int>(std::floor(
          (pos[dim] - domain_bounding_box_corner_positions_(dim, 0)) * inv_bin_size_[dim]));
  }

  return ConvertijkToGid(ijk);
}

void Core::Binstrategy::BinningStrategy::ConvertPosToijk(const double* pos, int* ijk) const
{
  double pos_ud[3];
  if (deforming_simulation_domain_handler_ != Teuchos::null)
  {
    deforming_simulation_domain_handler_->transform_from_global_to_undeformed_bounding_box_system(
        pos, pos_ud);
    for (int dim = 0; dim < 3; ++dim)
      ijk[dim] = static_cast<int>(std::floor(
          (pos_ud[dim] - domain_bounding_box_corner_positions_(dim, 0)) * inv_bin_size_[dim]));
  }
  else
  {
    for (int dim = 0; dim < 3; ++dim)
      ijk[dim] = static_cast<int>(std::floor(
          (pos[dim] - domain_bounding_box_corner_positions_(dim, 0)) * inv_bin_size_[dim]));
  }
}

void Core::Binstrategy::BinningStrategy::ConvertPosToijk(
    const Core::LinAlg::Matrix<3, 1>& pos, int* ijk) const
{
  Core::LinAlg::Matrix<3, 1> pos_ud;
  if (deforming_simulation_domain_handler_ != Teuchos::null)
  {
    deforming_simulation_domain_handler_->transform_from_global_to_undeformed_bounding_box_system(
        pos, pos_ud);
    for (int dim = 0; dim < 3; ++dim)
      ijk[dim] = static_cast<int>(std::floor(
          (pos_ud(dim) - domain_bounding_box_corner_positions_(dim, 0)) * inv_bin_size_[dim]));
  }
  else
  {
    for (int dim = 0; dim < 3; ++dim)
      ijk[dim] = static_cast<int>(std::floor(
          (pos(dim) - domain_bounding_box_corner_positions_(dim, 0)) * inv_bin_size_[dim]));
  }
}

int Core::Binstrategy::BinningStrategy::ConvertPosToGid(const Core::LinAlg::Matrix<3, 1>& pos) const
{
  int ijk[3];
  Core::LinAlg::Matrix<3, 1> pos_ud;
  if (deforming_simulation_domain_handler_ != Teuchos::null)
  {
    deforming_simulation_domain_handler_->transform_from_global_to_undeformed_bounding_box_system(
        pos, pos_ud);
    for (int dim = 0; dim < 3; ++dim)
      ijk[dim] = static_cast<int>(std::floor(
          (pos_ud(dim) - domain_bounding_box_corner_positions_(dim, 0)) * inv_bin_size_[dim]));
  }
  else
  {
    for (int dim = 0; dim < 3; ++dim)
      ijk[dim] = static_cast<int>(std::floor(
          (pos(dim) - domain_bounding_box_corner_positions_(dim, 0)) * inv_bin_size_[dim]));
  }

  return ConvertijkToGid(ijk);
}

void Core::Binstrategy::BinningStrategy::GetNeighborBinIds(
    const int binId, std::vector<int>& binIds) const
{
  int ijk_base[3];
  ConvertGidToijk(binId, ijk_base);

  for (int i = ijk_base[0] - 1; i <= ijk_base[0] + 1; ++i)
  {
    for (int j = ijk_base[1] - 1; j <= ijk_base[1] + 1; ++j)
    {
      for (int k = ijk_base[2] - 1; k <= ijk_base[2] + 1; ++k)
      {
        std::array ijk = {i, j, k};
        const int gid = ConvertijkToGid(ijk.data());
        if (gid != -1 and gid != binId)
        {
          binIds.push_back(gid);
        }
      }
    }
  }
}

void Core::Binstrategy::BinningStrategy::get_neighbor_and_own_bin_ids(
    const int binId, std::vector<int>& binIds) const
{
  // get neighbors
  GetNeighborBinIds(binId, binIds);

  // add myself
  binIds.push_back(binId);

  // in case of less than two bins in pbc direction, this is needed
  // to avoid double contact evaluation
  if (havepbc_)
  {
    std::sort(binIds.begin(), binIds.end());
    binIds.erase(std::unique(binIds.begin(), binIds.end()), binIds.end());
  }
}

void Core::Binstrategy::BinningStrategy::GetBinCorners(
    const int binId, std::vector<Core::LinAlg::Matrix<3, 1>>& bincorners) const
{
  bincorners.clear();
  bincorners.reserve(8);
  int ijk_base[3];
  ConvertGidToijk(binId, ijk_base);

  // order in bincorners is identical to ordering of i,j and k
  for (int k = ijk_base[2]; k < (ijk_base[2] + 2); ++k)
  {
    for (int j = ijk_base[1]; j < (ijk_base[1] + 2); ++j)
    {
      for (int i = ijk_base[0]; i < (ijk_base[0] + 2); ++i)
      {
        const std::array<int, 3> ijk_curr = {i, j, k};
        Core::LinAlg::Matrix<3, 1> curr_corner;
        for (int dim = 0; dim < 3; ++dim)
        {
          curr_corner(dim) =
              domain_bounding_box_corner_positions_(dim, 0) + bin_size_[dim] * ijk_curr[dim];
        }
        bincorners.push_back(curr_corner);

      }  // end for int k
    }    // end for int j
  }      // end for int i

  // change entries to get node numbering according to 4C convention
  std::swap(bincorners[2], bincorners[3]);
  std::swap(bincorners[6], bincorners[7]);
}

void Core::Binstrategy::BinningStrategy::GetAllBinCenters(
    Teuchos::RCP<Epetra_Map>& binrowmap, Teuchos::RCP<Epetra_MultiVector>& bincenters) const
{
  // loop over row bins and get center coordinates
  for (int i = 0; i < binrowmap->NumMyElements(); ++i)
  {
    // get global id of bin
    const int gidofbin = binrowmap->GID(i);

    // get coordinates of bin center
    Core::LinAlg::Matrix<3, 1> center = GetBinCentroid(gidofbin);

    for (int dim = 0; dim < 3; ++dim) bincenters->ReplaceMyValue(i, dim, center(dim));
  }
}

Core::LinAlg::Matrix<3, 1> Core::Binstrategy::BinningStrategy::GetBinCentroid(const int binId) const
{
  int ijk[3];
  ConvertGidToijk(binId, ijk);
  if (ijk[0] == -1)
    FOUR_C_THROW("given bin id is outside of bins; centroid of bin is does not make sense");

  Core::LinAlg::Matrix<3, 1> centroid;
  for (int dim = 0; dim < 3; ++dim)
    centroid(dim) =
        domain_bounding_box_corner_positions_(dim, 0) + bin_size_[dim] * (ijk[dim] + 0.5);

  return centroid;
}

double Core::Binstrategy::BinningStrategy::GetMinBinSize() const
{
  return std::min(bin_size_[0], std::min(bin_size_[1], bin_size_[2]));
}

double Core::Binstrategy::BinningStrategy::GetMaxBinSize() const
{
  return std::max(bin_size_[0], std::max(bin_size_[1], bin_size_[2]));
}

void Core::Binstrategy::BinningStrategy::BuildPeriodicBC(
    const Teuchos::ParameterList& binning_params)
{
  std::istringstream periodicbc(
      Teuchos::getNumericStringParameter(binning_params, "PERIODICONOFF"));

  // loop over all spatial directions
  for (int dim = 0; dim < 3; ++dim)
  {
    std::string val;
    if (periodicbc >> val)
    {
      int intval = std::atoi(val.c_str());
      if (intval)
      {
        // output pbc bounds based on XAABB of bins
        if (myrank_ == 0)
        {
          Core::IO::cout(Core::IO::verbose)
              << "INFO: PBC bounds for particles is computed automatically for direction " << dim
              << " based on XAABB of bins (left: " << domain_bounding_box_corner_positions_(dim, 0)
              << " , right: " << domain_bounding_box_corner_positions_(dim, 1) << " )"
              << Core::IO::endl;
        }

        // set flag
        pbconoff_[dim] = true;

        // set global flag
        havepbc_ = true;
      }
    }
    else
    {
      FOUR_C_THROW(
          "Enter three values to specify each direction as periodic or non periodic. Fix input "
          "file ...");
    }
  }
}

void Core::Binstrategy::BinningStrategy::determine_boundary_row_bins()
{
  // clear old content
  boundaryrowbins_.clear();

  // fill discret if necessary to obtain bin row map
  if (!bindis_->Filled()) bindis_->fill_complete(false, false, false);

  // determine maximal possible number of neighbors
  size_t nummaxneighbors = 26;

  // loop over row bins and decide whether they are located at the boundary
  const Epetra_Map* binrowmap = bindis_->ElementRowMap();
  for (int lid = 0; lid < binrowmap->NumMyElements(); ++lid)
  {
    Core::Elements::Element* currbin = bindis_->lRowElement(lid);
    std::vector<int> binvec;
    binvec.reserve(nummaxneighbors);
    // get neighboring bins
    GetNeighborBinIds(currbin->Id(), binvec);

    // a bin with less than 26 (or 8 in 2D) neighbors is a boundary bin
    if (binvec.size() < nummaxneighbors)
    {
      boundaryrowbins_.push_back(currbin);
      continue;
    }

    // a bin with less than 26 (or 8 in 2D) row neighbors is a boundary bin between processors
    std::vector<int> rowbinvec;
    rowbinvec.reserve(nummaxneighbors);
    for (int& it : binvec)
    {
      if (binrowmap->LID(it) != -1) rowbinvec.push_back(it);
    }

    if (rowbinvec.size() < nummaxneighbors) boundaryrowbins_.push_back(currbin);
  }
}

void Core::Binstrategy::BinningStrategy::determine_boundary_col_bins()
{
  boundarycolbins_.clear();

  if (boundaryrowbins_.size() == 0) determine_boundary_row_bins();

  // loop over boundary row bins and add neighbors
  std::list<Core::Elements::Element*>::const_iterator it;
  for (it = boundaryrowbins_.begin(); it != boundaryrowbins_.end(); ++it)
  {
    std::vector<int> binvec;
    binvec.reserve(26);
    // get neighboring bins
    GetNeighborBinIds((*it)->Id(), binvec);
    boundarycolbins_.insert(binvec.begin(), binvec.end());
  }
}

Teuchos::RCP<Epetra_Map> Core::Binstrategy::BinningStrategy::create_linear_map_for_numbin(
    const Epetra_Comm& comm) const
{
  // initial dummy distribution using a linear map
  const int numproc = comm.NumProc();
  const int numbin = bin_per_dir_[0] * bin_per_dir_[1] * bin_per_dir_[2];
  const int start = numbin / numproc * myrank_;
  int end;
  // special treatment for last proc
  if (myrank_ != numproc - 1)
    end = (int)(numbin / numproc * (myrank_ + 1));
  else
    end = numbin;

  std::vector<int> linearmap;
  linearmap.reserve(end - start);
  for (int k = 0; k < bin_per_dir_[2]; ++k)
  {
    for (int j = 0; j < bin_per_dir_[1]; ++j)
    {
      for (int i = 0; i < bin_per_dir_[0]; ++i)
      {
        int curr = i + j * bin_per_dir_[0] + k * bin_per_dir_[0] * bin_per_dir_[1];
        if (start <= curr and curr < end)
          linearmap.push_back(i + j * id_calc_bin_per_dir_[0] +
                              k * id_calc_bin_per_dir_[0] * id_calc_bin_per_dir_[1]);
      }
    }
  }

  return Teuchos::rcp(new Epetra_Map(numbin, linearmap.size(), linearmap.data(), 0, comm));
}

void Core::Binstrategy::BinningStrategy::WriteBinOutput(int const step, double const time)
{
  // no bin output
  if (writebinstype_ == Inpar::BINSTRATEGY::none) return;

  if (myrank_ == 0)
    Core::IO::cout(Core::IO::verbose) << "\nBinning discretization output (step " << step
                                      << ", time " << time << ") written." << Core::IO::endl;

  // -------------------------------------------------------------------------
  // note: this is a debug feature only (as very expensive)
  // -------------------------------------------------------------------------

  // -------------------------------------------------------------------------
  // create new discretization containing the bins as elements
  // -------------------------------------------------------------------------
  Teuchos::RCP<Epetra_Comm> com = Teuchos::rcp(bindis_->Comm().Clone());

  // store gids of ghosted elements
  std::map<int, std::vector<Core::LinAlg::Matrix<3, 1>>> ghostcorners;
  // add elements and nodes
  for (int i = 0; i < bindis_->NumMyColElements(); ++i)
  {
    Core::Elements::Element* ele = bindis_->lColElement(i);
    if (!ele) FOUR_C_THROW("Cannot find element with lid %", i);

    // get corner position as node positions
    const int numcorner = 8;
    int bingid = ele->Id();
    std::vector<Core::LinAlg::Matrix<3, 1>> bincorners;
    GetBinCorners(bingid, bincorners);

    // if element is a ghost
    if (ele->Owner() != myrank_)
    {
      ghostcorners[ele->Id()] = bincorners;
      continue;
    }

    // add new node
    std::vector<double> cornerpos(3, 0.0);
    std::vector<int> nids(8, 0);
    for (int corner_i = 0; corner_i < numcorner; ++corner_i)
    {
      for (int dim = 0; dim < 3; ++dim) cornerpos[dim] = bincorners[corner_i](dim);

      nids[corner_i] = (bingid * numcorner) + corner_i;
      Teuchos::RCP<Core::Nodes::Node> newnode =
          Teuchos::rcp(new Core::Nodes::Node(nids[corner_i], cornerpos, myrank_));
      visbindis_->AddNode(newnode);
    }

    // assign nodes to elements
    Teuchos::RCP<Core::Elements::Element> newele =
        Core::Communication::Factory("VELE3", "Polynomial", ele->Id(), myrank_);
    newele->SetNodeIds(nids.size(), nids.data());
    visbindis_->add_element(newele);
  }

  // get max gid before adding elements
  int maxgid = bindis_->ElementRowMap()->MaxAllGID() + 1;
  if (writebinstype_ == Inpar::BINSTRATEGY::cols)
  {
    // gather all numbers of ghosted bins that are going to be row eles
    int nummycol;
    nummycol = static_cast<int>(ghostcorners.size());
    // initialize std::vector for communication
    std::vector<int> numcol(com->NumProc(), 0);
    // communicate
    com->GatherAll(&nummycol, numcol.data(), 1);
    com->Barrier();

    // calculate starting index on myrank
    int startnewgid = 0;
    for (int i = 0; i < myrank_; ++i) startnewgid += numcol[i];

    // loop over all ghosted bins
    std::map<int, std::vector<Core::LinAlg::Matrix<3, 1>>>::const_iterator iter;
    int counter = 0;
    for (iter = ghostcorners.begin(); iter != ghostcorners.end(); ++iter)
    {
      // new elegid (unique over all procs)
      int newelegid = maxgid + startnewgid + counter;
      counter++;

      // add new node
      // get corner position as node positions
      const int numcorner = 8;
      std::vector<double> cornerpos(3, 0.0);
      std::vector<int> nids(8, 0);
      for (int corner_i = 0; corner_i < numcorner; ++corner_i)
      {
        for (int dim = 0; dim < 3; ++dim) cornerpos[dim] = iter->second[corner_i](dim);

        nids[corner_i] = (newelegid * numcorner) + corner_i;
        Teuchos::RCP<Core::Nodes::Node> newnode =
            Teuchos::rcp(new Core::Nodes::Node(nids[corner_i], cornerpos, myrank_));
        visbindis_->AddNode(newnode);
      }

      // assign nodes to elements
      Teuchos::RCP<Core::Elements::Element> newele =
          Core::Communication::Factory("VELE3", "Polynomial", newelegid, myrank_);
      newele->SetNodeIds(nids.size(), nids.data());
      visbindis_->add_element(newele);
    }
  }

  // complete new dis
  Teuchos::RCP<Core::DOFSets::IndependentDofSet> independentdofset =
      Teuchos::rcp(new Core::DOFSets::IndependentDofSet(true));
  visbindis_->ReplaceDofSet(independentdofset);
  visbindis_->fill_complete(true, true, false);

  // create vector that shows ghosting
  Teuchos::RCP<Epetra_Vector> ownedghostsvec =
      Core::LinAlg::CreateVector(*visbindis_->ElementRowMap(), true);
  for (int i = 0; i < visbindis_->NumMyRowElements(); ++i)
  {
    Core::Elements::Element* ele = visbindis_->lRowElement(i);

    if (ele->Id() < maxgid)
      (*ownedghostsvec)[i] = 0;  // owned
    else
      (*ownedghostsvec)[i] = 1;  // ghost
  }

  // write output
  visbindis_->Writer()->write_mesh(step, time);
  visbindis_->Writer()->new_step(step, time);
  visbindis_->Writer()->write_vector("owner0ghost1", ownedghostsvec, Core::IO::elementvector);
  visbindis_->Writer()->write_element_data(true);

  visbindis_->ClearDiscret();
}

void Core::Binstrategy::BinningStrategy::distribute_bins_recurs_coord_bisection(
    Teuchos::RCP<Epetra_Map>& binrowmap, Teuchos::RCP<Epetra_MultiVector>& bincenters,
    Teuchos::RCP<Epetra_MultiVector>& binweights) const
{
  // create a parameter list for partitioner
  Teuchos::ParameterList params;
  params.set("Partitioning Method", "RCB");

  // set low-level partitioning parameters (see Zoltan Users' Guide:
  // http://www.cs.sandia.gov/zoltan)
  Teuchos::ParameterList& sublist = params.sublist("Zoltan");

  // debug level (see http://www.cs.sandia.gov/zoltan/ug_html/ug_param.html)
  sublist.set("DEBUG_LEVEL", "0");

  // recursive coordinate bisection (see http://www.cs.sandia.gov/zoltan/ug_html/ug_alg_rcb.html)
  sublist.set("RCB_OUTPUT_LEVEL", "0");
  sublist.set("RCB_RECTILINEAR_BLOCKS", "1");

  std::tie(bincenters, binweights) =
      Core::Rebalance::RebalanceCoordinates(*bincenters.getConst(), params, *binweights.getConst());

  // create bin row map
  binrowmap = Teuchos::rcp(new Epetra_Map(-1, bincenters->Map().NumMyElements(),
      bincenters->Map().MyGlobalElements(), 0, BinDiscret()->Comm()));
}

void Core::Binstrategy::BinningStrategy::fill_bins_into_bin_discretization(
    Teuchos::RCP<Epetra_Map> const& rowbins)
{
  // fill bins into bindis_
  for (int i = 0; i < rowbins->NumMyElements(); ++i)
  {
    const int gid = rowbins->GID(i);
    Teuchos::RCP<Core::Elements::Element> bin =
        Core::Communication::Factory("MESHFREEMULTIBIN", "dummy", gid, myrank_);
    bindis_->add_element(bin);
  }
}

void Core::Binstrategy::BinningStrategy::addijk_to_axis_alignedijk_range_of_element(
    int const ijk[3], int ijk_range[6]) const
{
  for (int dim = 0; dim < 3; ++dim)
  {
    if (ijk[dim] < ijk_range[dim * 2])
    {
      ijk_range[dim * 2] = ijk[dim];
    }

    if (ijk[dim] > ijk_range[dim * 2 + 1])
    {
      ijk_range[dim * 2 + 1] = ijk[dim];
    }
  }
}

void Core::Binstrategy::BinningStrategy::addijk_to_axis_alignedijk_range_of_beam_element(
    int const ijk[3], int ijk_range[6]) const
{
  // this should be large enough
  int const cutcheckfac = 5;

  for (int dim = 0; dim < 3; ++dim)
  {
    if (ijk[dim] < ijk_range[dim * 2])
    {
      // this is needed if your element is cut by a periodic boundary and you don't
      // want a bounding box that is as big as the the whole binning domain
      if ((ijk[dim] == 0) && (abs(ijk[dim] - ijk_range[dim * 2]) > cutcheckfac))
      {
        ijk_range[dim * 2 + 1] = bin_per_dir_[dim];
        continue;
      }
      else
      {
        ijk_range[dim * 2] = ijk[dim];
      }
    }
    if (ijk[dim] > ijk_range[dim * 2 + 1])
    {
      // cut check
      if ((ijk[dim] == bin_per_dir_[dim] - 1) &&
          (abs(ijk[dim] - ijk_range[dim * 2 + 1]) > cutcheckfac))
      {
        ijk_range[dim * 2] = -1;
      }
      else
      {
        ijk_range[dim * 2 + 1] = ijk[dim];
      }
    }
  }
}

void Core::Binstrategy::BinningStrategy::build_axis_alignedijk_range_for_rigid_sphere(
    Core::Elements::Element const* const sphereele, double currpos[3], int ijk[3], int ijk_range[6],
    const double radius) const
{
  for (int j = 0; j < 3; ++j)
  {
    double* coords = currpos;
    coords[j] += radius;
    ConvertPosToijk(coords, ijk);
    addijk_to_axis_alignedijk_range_of_element(ijk, ijk_range);
    coords[j] -= (2.0 * radius);
    ConvertPosToijk(coords, ijk);
    addijk_to_axis_alignedijk_range_of_element(ijk, ijk_range);
    coords[j] += radius;
  }
}

void Core::Binstrategy::BinningStrategy::distribute_eles_to_bins(
    const Core::FE::Discretization& mortardis, std::map<int, std::set<int>>& binelemap,
    bool isslave) const
{
  // exploit bounding box idea for elements and bins
  for (int lid = 0; lid < mortardis.NumMyColElements(); ++lid)
  {
    Core::Elements::Element* ele = mortardis.lColElement(lid);
    if (dynamic_cast<Mortar::Element*>(ele)->IsSlave() == isslave)
    {
      Core::Nodes::Node** nodes = ele->Nodes();
      const int numnode = ele->num_node();

      // initialize ijk_range with ijk of first node of element
      int ijk[3];
      {
        Core::Nodes::Node* node = nodes[0];
        const double* coords = dynamic_cast<Mortar::Node*>(node)->xspatial();
        ConvertPosToijk(coords, ijk);
      }

      // ijk_range contains: i_min i_max j_min j_max k_min k_max
      int ijk_range[] = {ijk[0], ijk[0], ijk[1], ijk[1], ijk[2], ijk[2]};

      // fill in remaining nodes
      for (int j = 1; j < numnode; ++j)
      {
        Core::Nodes::Node* node = nodes[j];
        const double* coords = dynamic_cast<Mortar::Node*>(node)->xspatial();
        int ijk[3];
        ConvertPosToijk(coords, ijk);

        for (int dim = 0; dim < 3; ++dim)
        {
          if (ijk[dim] < ijk_range[dim * 2]) ijk_range[dim * 2] = ijk[dim];
          if (ijk[dim] > ijk_range[dim * 2 + 1]) ijk_range[dim * 2 + 1] = ijk[dim];
        }
      }

      // get corresponding bin ids in ijk range
      std::vector<int> binIds;
      binIds.reserve(get_number_of_bins_inijk_range(ijk_range));
      GidsInijkRange(ijk_range, binIds, false);

      // assign element to bins
      for (int binId : binIds) binelemap[binId].insert(ele->Id());
    }
  }
}

void Core::Binstrategy::BinningStrategy::distribute_row_elements_to_bins_using_ele_aabb(
    Teuchos::RCP<Core::FE::Discretization> const& discret,
    std::map<int, std::set<int>>& bintorowelemap, Teuchos::RCP<const Epetra_Vector> disnp) const
{
  bintorowelemap.clear();

  // exploit bounding box idea for elements in underlying discretization and bins
  // loop over all row elements
  for (int lid = 0; lid < discret->NumMyRowElements(); ++lid)
  {
    Core::Elements::Element* eleptr = discret->lRowElement(lid);
    // get corresponding bin ids in ijk range
    std::vector<int> binIds;
    distribute_single_element_to_bins_using_ele_aabb(discret, eleptr, binIds, disnp);

    // assign element to bins
    std::vector<int>::const_iterator biniter;
    for (biniter = binIds.begin(); biniter != binIds.end(); ++biniter)
      bintorowelemap[*biniter].insert(eleptr->Id());
  }
}

void Core::Binstrategy::BinningStrategy::distribute_col_elements_to_bins_using_ele_aabb(
    Teuchos::RCP<Core::FE::Discretization> const& discret,
    std::map<int, std::set<int>>& bintocolelemap, Teuchos::RCP<const Epetra_Vector> disnp) const
{
  bintocolelemap.clear();

  // exploit bounding box idea for elements in underlying discretization and bins
  // loop over all row elements
  for (int lid = 0; lid < discret->NumMyColElements(); ++lid)
  {
    Core::Elements::Element* eleptr = discret->lColElement(lid);
    // get corresponding bin ids in ijk range
    std::vector<int> binIds;
    distribute_single_element_to_bins_using_ele_aabb(discret, eleptr, binIds, disnp);

    // assign element to bins
    std::vector<int>::const_iterator biniter;
    for (biniter = binIds.begin(); biniter != binIds.end(); ++biniter)
      bintocolelemap[*biniter].insert(eleptr->Id());
  }
}

void Core::Binstrategy::BinningStrategy::distribute_single_element_to_bins_using_ele_aabb(
    Teuchos::RCP<Core::FE::Discretization> const& discret, Core::Elements::Element* eleptr,
    std::vector<int>& binIds, Teuchos::RCP<const Epetra_Vector> const& disnp) const
{
  binIds.clear();
  Core::Nodes::Node** nodes = eleptr->Nodes();

  // initialize ijk_range with ijk of first node of element
  int ijk[3];
  Core::Nodes::Node const* const node = nodes[0];
  getijk_of_single_node_in_current_position(discret, node, disnp, ijk);

  // ijk_range contains: i_min i_max j_min j_max k_min k_max
  int ijk_range[] = {ijk[0], ijk[0], ijk[1], ijk[1], ijk[2], ijk[2]};

  // bounding box idea for rigid sphere element with just one node needs
  // some special treatment
  switch (element_filter_(eleptr))
  {
    case Utils::SpecialElement::rigid_sphere:
    {
      double currpos[3] = {0.0, 0.0, 0.0};
      Utils::GetCurrentNodePos(
          discret, eleptr->Nodes()[0], correct_beam_center_node_, disnp, currpos);
      build_axis_alignedijk_range_for_rigid_sphere(
          eleptr, currpos, ijk, ijk_range, rigid_sphere_radius_(eleptr));
      break;
    }
    case Utils::SpecialElement::beam:
    {
      // fill in remaining nodes
      for (int j = 1; j < eleptr->num_node(); ++j)
      {
        Core::Nodes::Node const* const node = nodes[j];
        getijk_of_single_node_in_current_position(discret, node, disnp, ijk);
        addijk_to_axis_alignedijk_range_of_beam_element(ijk, ijk_range);
      }
      break;
    }
    default:
    {
      // fill in remaining nodes
      for (int j = 1; j < eleptr->num_node(); ++j)
      {
        Core::Nodes::Node const* const node = nodes[j];
        getijk_of_single_node_in_current_position(discret, node, disnp, ijk);
        addijk_to_axis_alignedijk_range_of_element(ijk, ijk_range);
      }
      break;
    }
  }
  // get corresponding bin ids in ijk range
  binIds.reserve(get_number_of_bins_inijk_range(ijk_range));
  GidsInijkRange(ijk_range, binIds, false);
}

void Core::Binstrategy::BinningStrategy::AssignElesToBins(
    Teuchos::RCP<Core::FE::Discretization> discret,
    std::map<int, std::set<int>> const& extended_bin_to_row_ele_map,
    const std::function<Utils::BinContentType(const Core::Elements::Element* element)>&
        ele_to_bin_type) const
{
  // loop over bins
  for (const auto& [bin_gid, ele_gids] : extended_bin_to_row_ele_map)
  {
    // extract bins from discretization after checking on existence
    const int bin_lid = bindis_->ElementColMap()->LID(bin_gid);
    if (bin_lid < 0) continue;

    // get current bin
    auto* currbin = dynamic_cast<Core::FE::MeshFree::MeshfreeMultiBin*>(bindis_->gElement(bin_gid));

    // loop over ele content of this bin
    for (const auto& ele_gid : ele_gids)
    {
      auto* ele = discret->gElement(ele_gid);
      // add eleid and elepointer to current bin
      currbin->AddAssociatedEle(ele_to_bin_type(ele), ele);
    }
  }
}

void Core::Binstrategy::BinningStrategy::GetBinContent(std::set<Core::Elements::Element*>& eles,
    const std::vector<Core::Binstrategy::Utils::BinContentType>& bincontent,
    std::vector<int>& binIds, bool roweles) const
{
  // loop over all bins
  std::vector<int>::const_iterator biniter;
  for (biniter = binIds.begin(); biniter != binIds.end(); ++biniter)
  {
    // extract bins from discretization after checking on existence
    const int lid = bindis_->ElementColMap()->LID(*biniter);
    if (lid < 0) continue;

    // get content of current bin
    auto* bin = static_cast<Core::FE::MeshFree::MeshfreeMultiBin*>(bindis_->lColElement(lid));

    // loop over bincontent you want to get
    for (const auto& bc_i : bincontent)
    {
      // gather elements of with specific bincontent type
      auto elements = bin->AssociatedEles(bc_i);
      for (const auto& ele : elements)
      {
        if (roweles && ele->Owner() != myrank_) continue;
        eles.insert(ele);
      }
    }
  }
}

void Core::Binstrategy::BinningStrategy::remove_all_eles_from_bins()
{
  // loop over all bins and remove assigned elements
  const int numcolbins = bindis_->NumMyColElements();
  for (int binlid = 0; binlid < numcolbins; ++binlid)
  {
    Core::Elements::Element* currentbin = bindis_->lColElement(binlid);
    dynamic_cast<Core::FE::MeshFree::MeshfreeMultiBin*>(currentbin)->remove_all_associated_eles();
  }
}

void Core::Binstrategy::BinningStrategy::getijk_of_single_node_in_current_position(
    Teuchos::RCP<Core::FE::Discretization> const& discret, Core::Nodes::Node const* const node,
    Teuchos::RCP<const Epetra_Vector> const& disnp, int ijk[3]) const
{
  double currpos[3] = {0.0, 0.0, 0.0};
  Utils::GetCurrentNodePos(discret, node, correct_beam_center_node_, disnp, currpos);
  double const* coords = currpos;
  ConvertPosToijk(coords, ijk);
}

void Core::Binstrategy::BinningStrategy::distribute_row_nodes_to_bins(
    Teuchos::RCP<Core::FE::Discretization> discret,
    std::map<int, std::vector<int>>& bin_to_rownodes_map,
    Teuchos::RCP<const Epetra_Vector> disnp) const
{
  // current position of nodes
  double currpos[3] = {0.0, 0.0, 0.0};

  // loop over row nodes
  for (int lid = 0; lid < discret->NumMyRowNodes(); ++lid)
  {
    Core::Nodes::Node* node = discret->lRowNode(lid);
    Utils::GetCurrentNodePos(discret, node, correct_beam_center_node_, disnp, currpos);

    const double* coords = currpos;
    int ijk[3];
    ConvertPosToijk(coords, ijk);
    const int binid = ConvertijkToGid(ijk);

    FOUR_C_THROW_UNLESS(binid != -1,
        "Node %i in your discretization resides outside the binning \n"
        "domain, this does not work at this point.",
        node->Id());

    // assign node to bin
    bin_to_rownodes_map[binid].push_back(node->Id());
  }
}

Teuchos::RCP<Epetra_Map> Core::Binstrategy::BinningStrategy::
    do_weighted_partitioning_of_bins_and_extend_ghosting_of_discret_to_one_bin_layer(
        std::vector<Teuchos::RCP<Core::FE::Discretization>> discret,
        std::vector<Teuchos::RCP<Epetra_Map>>& stdelecolmap,
        std::vector<Teuchos::RCP<Epetra_Map>>& stdnodecolmap)
{
  // initialize dummys
  std::vector<std::map<int, std::set<int>>> dummy1(discret.size());
  std::vector<Teuchos::RCP<const Epetra_Vector>> dummy2(discret.size());
  std::vector<Teuchos::RCP<Epetra_Vector>> mutabledummy2(discret.size());

  // ------------------------------------------------------------------------
  // create bins, weight them according to number of nodes (of discrets) they
  // contain, account for bin connectivity. Then an optimal distribution of
  // bins to procs can be obtained
  // ------------------------------------------------------------------------
  // nodes, that are owned by a proc, are distributed to the bins of this proc
  std::vector<std::map<int, std::vector<int>>> nodesinbin(discret.size());
  // default weight 10.0
  double const weight = 1.0;
  Teuchos::RCP<Epetra_Map> newrowbins =
      weighted_distribution_of_bins_to_procs(discret, dummy2, nodesinbin, weight);

  stdelecolmap.resize(discret.size());
  stdnodecolmap.resize(discret.size());

  // ------------------------------------------------------------------------
  // now we have an optimal distribution of bins (with respect to their content
  // and connectivity). Now we have to apply it by rebuilding the input discret,
  // i.e. we need to change the ownership of the nodes/elements according to
  // the bin they belong to (each proc then owns the nodes/eles laying in its
  // bins.
  // ------------------------------------------------------------------------
  // rebuild discretizations including extended ghosting
  for (size_t i = 0; i < discret.size(); ++i)
  {
    // ----------------------------------------------------------------------
    // start with standard ghosting
    // ----------------------------------------------------------------------
    standard_discretization_ghosting(
        discret[i], newrowbins, mutabledummy2[i], stdelecolmap[i], stdnodecolmap[i]);

    // ----------------------------------------------------------------------
    // extended ghosting
    // ----------------------------------------------------------------------
    // ----------------------------------------------------------------------
    // start with extended ghosting. Here this means the following: Each proc
    // ghosts all elements whose XAABB cuts a bin that is next to a bin that is
    // owned by a proc an not empty. All associated nodes are ghosted as well
    // ----------------------------------------------------------------------
    // here each proc assignes his owned elements in the means of a XAABB to
    // the global binids that do not need be owned by this proc.
    // binelemap on each proc than contains all bins (not neccesarily owned by
    // this proc) that are cut by the procs row elements
    std::map<int, std::set<int>> bintoelemap;
    distribute_row_elements_to_bins_using_ele_aabb(discret[i], bintoelemap, dummy2[i]);

    // ghosting is extended to one layer (two layer ghosting is excluded as it
    // is not needed, this case is covered by other procs then) around bins that
    // actually contain elements.
    // extbintoelemap[i] than contains all bins and its corresponding elements
    // that need to be owned or ghosted to ensure correct interaction handling
    // of the elements in the range of one layer
    Teuchos::RCP<Epetra_Map> extendedelecolmap = ExtendElementColMap(bintoelemap, bintoelemap,
        dummy1[i], Teuchos::null, newrowbins, discret[i]->ElementColMap());

    // adapt layout to extended ghosting in discret
    // first export the elements according to the processor local element column maps
    discret[i]->export_column_elements(*extendedelecolmap);

    // get the node ids of the elements that are to be ghosted
    // and create a proper node column map for their export
    std::set<int> nodes;
    for (int lid = 0; lid < extendedelecolmap->NumMyElements(); ++lid)
    {
      Core::Elements::Element* ele = discret[i]->gElement(extendedelecolmap->GID(lid));
      const int* nodeids = ele->NodeIds();
      for (int inode = 0; inode < ele->num_node(); ++inode) nodes.insert(nodeids[inode]);
    }

    std::vector<int> colnodes(nodes.begin(), nodes.end());
    Teuchos::RCP<Epetra_Map> nodecolmap = Teuchos::rcp(
        new Epetra_Map(-1, (int)colnodes.size(), colnodes.data(), 0, discret[i]->Comm()));

    // now ghost the nodes
    discret[i]->ExportColumnNodes(*nodecolmap);

    // fillcomplete discret with extended ghosting
    discret[i]->fill_complete();
    if (myrank_ == 0) std::cout << "parallel distribution with extended ghosting\n";
    Core::Rebalance::UTILS::print_parallel_distribution(*discret[i]);
  }

  return newrowbins;
}

Teuchos::RCP<Epetra_Map> Core::Binstrategy::BinningStrategy::weighted_distribution_of_bins_to_procs(
    std::vector<Teuchos::RCP<Core::FE::Discretization>>& discret,
    std::vector<Teuchos::RCP<const Epetra_Vector>>& disnp,
    std::vector<std::map<int, std::vector<int>>>& row_nodes_to_bin_map, double const& weight,
    bool repartition) const
{
  // calculate total number of bins
  const int numbin = bin_per_dir_[0] * bin_per_dir_[1] * bin_per_dir_[2];

  // some safety checks to ensure efficiency
  {
    if (numbin < discret[0]->Comm().NumProc() && myrank_ == 0)
    {
      FOUR_C_THROW(
          "ERROR:NumProc %i > NumBin %i. Too many processors to "
          "distribute your bins properly!!!",
          numbin, numbin < discret[0]->Comm().NumProc());
    }

    if (numbin < 8 * discret[0]->Comm().NumProc() && myrank_ == 0)
    {
      std::cout << "\n\nWARNING: partitioning not useful, choose less procs. "
                   " Owner distribution may be inefficient!\n\n";
    }
  }

  // row bin distribution
  Teuchos::RCP<Epetra_Map> rowbins = Teuchos::null;
  Teuchos::RCP<Epetra_CrsGraph> bingraph;
  if (repartition)
  {
    // use old bin distribution
    rowbins = Teuchos::rcp(const_cast<Epetra_Map*>(bindis_->ElementRowMap()));
    const Epetra_Map* oldrowmap = bindis_->ElementRowMap();

    const int maxband = 26;
    bingraph = Teuchos::rcp(new Epetra_CrsGraph(Copy, *oldrowmap, maxband, false));

    // fill all local entries into the graph
    for (int lid = 0; lid < oldrowmap->NumMyElements(); ++lid)
    {
      const int binId = oldrowmap->GID(lid);

      std::vector<int> neighbors;
      GetNeighborBinIds(binId, neighbors);

      int err = bingraph->InsertGlobalIndices(binId, (int)neighbors.size(), neighbors.data());
      if (err < 0)
        FOUR_C_THROW(
            "Epetra_CrsGraph::InsertGlobalIndices returned %d for global row %d", err, binId);
    }
  }
  else
  {
    // dummy row bin distribution (equally distributed over all procs as no
    // weighting done so far)
    rowbins = create_linear_map_for_numbin(discret[0]->Comm());
    // create nodal graph
    bingraph = Teuchos::rcp(new Epetra_CrsGraph(Copy, *rowbins, 108, false));
  }

  // Now we're going to create a Epetra_Vector with vertex/node weights to be
  // used for the partitioning operation (weights must be at least one for zoltan)
  Teuchos::RCP<Epetra_Vector> vweights = Core::LinAlg::CreateVector(*rowbins, true);

  // set weights of bins related to the number of nodes of discrets that are contained
  // empty bins have weight of 1
  vweights->PutScalar(1.0);

  // determine which node is in which bin and weight each bin according to
  // "weight" times the number of nodes it contains
  // assign all node gids to their corresponding bin
  for (int i = 0; i < static_cast<int>(discret.size()); ++i)
  {
    // distribute nodes, that are owned by a proc, to the bins of this proc
    distribute_row_nodes_to_bins(discret[i], row_nodes_to_bin_map[i], disnp[i]);

    std::map<int, std::vector<int>> nodesinmybins;
    // gather information of bin content from other procs (bin is owned by this
    // proc and there are some nodes on other procs which are located in this bin)
    // mynodesinbin then contains all node gids (vector) that reside in a owned bin (gid is map key)
    collect_information_about_content_of_bins_from_other_procs_via_round_robin(
        rowbins, row_nodes_to_bin_map[i], nodesinmybins);

    // weight each bin with 10 times the number of node it contains
    // empty bins remain with weight one
    std::map<int, std::vector<int>>::const_iterator biniter;
    for (biniter = nodesinmybins.begin(); biniter != nodesinmybins.end(); ++biniter)
    {
      int lid = rowbins->LID(biniter->first);
      // safety check
      if (lid < 0)
        FOUR_C_THROW("Proc %d: Cannot find gid=%d in Epetra_Vector", discret[i]->Comm().MyPID(),
            biniter->first);

      // weighting
      (*vweights)[lid] += weight * static_cast<double>(biniter->second.size());
    }
  }

  // fill bin connectivity into bin graph
  for (int lid = 0; lid < rowbins->NumMyElements(); ++lid)
  {
    int rowbinid = rowbins->GID(lid);
    // insert 26 (one level) neighboring bins to graph
    // (if active, periodic boundary conditions are considered here)
    std::vector<int> neighbors;
    GetNeighborBinIds(rowbinid, neighbors);

    int err = bingraph->InsertGlobalIndices(
        rowbinid, static_cast<int>(neighbors.size()), neighbors.data());
    if (err < 0)
      FOUR_C_THROW(
          "Epetra_CrsGraph::InsertGlobalIndices returned %d for global row %d", err, rowbinid);
  }

  // complete graph
  int err = bingraph->FillComplete();
  if (err) FOUR_C_THROW("graph->FillComplete() returned err=%d", err);
  err = bingraph->OptimizeStorage();
  if (err) FOUR_C_THROW("graph->OptimizeStorage() returned err=%d", err);

  Teuchos::ParameterList paramlist;
  paramlist.set("PARTITIONING METHOD", "GRAPH");
  Teuchos::ParameterList& sublist = paramlist.sublist("Zoltan");
  if (repartition)
    sublist.set("LB_APPROACH", "REPARTITION");
  else
    sublist.set("LB_APPROACH", "PARTITION");

  Teuchos::RCP<Epetra_CrsGraph> balanced_bingraph =
      Core::Rebalance::RebalanceGraph(*bingraph.getConst(), paramlist, vweights);

  // extract repartitioned bin row map
  const Epetra_BlockMap& rbinstmp = balanced_bingraph->RowMap();
  Teuchos::RCP<Epetra_Map> newrowbins = Teuchos::rcp(new Epetra_Map(
      -1, rbinstmp.NumMyElements(), rbinstmp.MyGlobalElements(), 0, discret[0]->Comm()));

  return newrowbins;
}

Teuchos::RCP<Epetra_Map> Core::Binstrategy::BinningStrategy::ExtendElementColMap(
    std::map<int, std::set<int>> const& bin_to_row_ele_map,
    std::map<int, std::set<int>>& bin_to_row_ele_map_to_lookup_requests,
    std::map<int, std::set<int>>& ext_bin_to_ele_map, Teuchos::RCP<Epetra_Map> bin_colmap,
    Teuchos::RCP<Epetra_Map> bin_rowmap, const Epetra_Map* ele_colmap_from_standardghosting) const
{
  // do communication to gather all elements for extended ghosting
  const int numproc = comm_->NumProc();
  for (int iproc = 0; iproc < numproc; ++iproc)
  {
    // gather set of column bins for each proc
    std::set<int> bins;
    if (iproc == myrank_)
    {
      // either use given column layout of bins ...
      if (bin_colmap != Teuchos::null)
      {
        int nummyeles = bin_colmap->NumMyElements();
        int* entries = bin_colmap->MyGlobalElements();
        bins.insert(entries, entries + nummyeles);
      }
      else  // ... or add an extra layer to the given bin distribution
      {
        std::map<int, std::set<int>>::const_iterator iter;
        for (iter = bin_to_row_ele_map.begin(); iter != bin_to_row_ele_map.end(); ++iter)
        {
          int binId = iter->first;
          // avoid getting two layer ghosting as this is not needed
          if (bin_rowmap != Teuchos::null)
          {
            const int lid = bin_rowmap->LID(binId);
            if (lid < 0) continue;
          }
          std::vector<int> binvec;
          // get neighboring bins
          get_neighbor_and_own_bin_ids(binId, binvec);
          bins.insert(binvec.begin(), binvec.end());
        }
      }
    }
    // copy set to vector in order to broadcast data
    std::vector<int> binids(bins.begin(), bins.end());

    // first: proc i tells all procs how many bins it has
    int numbin = binids.size();
    comm_->Broadcast(&numbin, 1, iproc);
    // second: proc i tells all procs which bins it has
    binids.resize(numbin);

    comm_->Broadcast(binids.data(), numbin, iproc);

    // loop over all own bins and find requested ones, fill in elements in these bins
    std::map<int, std::set<int>> sdata;
    std::map<int, std::set<int>> rdata;

    for (int i = 0; i < numbin; ++i)
    {
      if (bin_to_row_ele_map_to_lookup_requests.find(binids[i]) !=
          bin_to_row_ele_map_to_lookup_requests.end())
        sdata[binids[i]].insert(bin_to_row_ele_map_to_lookup_requests[binids[i]].begin(),
            bin_to_row_ele_map_to_lookup_requests[binids[i]].end());
    }

    Core::LinAlg::Gather<int>(sdata, rdata, 1, &iproc, *comm_);

    // proc i has to store the received data
    if (iproc == myrank_)
    {
      ext_bin_to_ele_map = rdata;
    }
  }

  // reduce map of sets to one set and copy to a vector to create extended elecolmap
  std::set<int> coleleset;
  std::map<int, std::set<int>>::iterator iter;
  for (iter = ext_bin_to_ele_map.begin(); iter != ext_bin_to_ele_map.end(); ++iter)
    coleleset.insert(iter->second.begin(), iter->second.end());

  // insert standard ghosting
  if (ele_colmap_from_standardghosting != nullptr)
    for (int lid = 0; lid < ele_colmap_from_standardghosting->NumMyElements(); ++lid)
      coleleset.insert(ele_colmap_from_standardghosting->GID(lid));

  std::vector<int> colgids(coleleset.begin(), coleleset.end());

  // return extended elecolmap
  return Teuchos::rcp(new Epetra_Map(-1, (int)colgids.size(), colgids.data(), 0, *comm_));
}

void Core::Binstrategy::BinningStrategy::extend_ghosting_of_binning_discretization(
    Teuchos::RCP<Epetra_Map> rowbins, std::set<int> const& colbins, bool assigndegreesoffreedom)
{
  // gather set of bins that need to be ghosted by myrank
  std::set<int> bins(colbins.begin(), colbins.end());

  // insert row bins
  for (int i = 0; i < rowbins->NumMyElements(); ++i) bins.insert(rowbins->GID(i));

  std::vector<int> bincolmapvec(bins.begin(), bins.end());
  Teuchos::RCP<Epetra_Map> bincolmap = Teuchos::rcp(new Epetra_Map(
      -1, static_cast<int>(bincolmapvec.size()), bincolmapvec.data(), 0, bindis_->Comm()));

  if (bincolmap->NumGlobalElements() == 1 && bindis_->Comm().NumProc() > 1)
    FOUR_C_THROW("one bin cannot be run in parallel -> reduce BIN_SIZE_LOWER_BOUND");

  Utils::ExtendDiscretizationGhosting(bindis_, bincolmap, assigndegreesoffreedom, false, true);

#ifdef FOUR_C_ENABLE_ASSERTIONS
  // check whether each proc has only particles that are within bins on this proc
  for (int k = 0; k < bindis_->NumMyColElements(); ++k)
  {
    int binid = bindis_->lColElement(k)->Id();
    Core::Nodes::Node** particles = bindis_->lColElement(k)->Nodes();

    for (int iparticle = 0; iparticle < bindis_->lColElement(k)->num_node(); ++iparticle)
    {
      double const* pos = particles[iparticle]->X().data();
      std::array ijk = {-1, -1, -1};
      ConvertPosToijk(pos, ijk.data());

      int gidofbin = ConvertijkToGid(ijk.data());
      if (gidofbin != binid)
        FOUR_C_THROW(
            "after ghosting: particle which should be in bin no. %i is in %i", gidofbin, binid);
    }
  }
#endif
}

void Core::Binstrategy::BinningStrategy::standard_discretization_ghosting(
    Teuchos::RCP<Core::FE::Discretization>& discret, Teuchos::RCP<Epetra_Map> const& rowbins,
    Teuchos::RCP<Epetra_Vector>& disnp, Teuchos::RCP<Epetra_Map>& stdelecolmap,
    Teuchos::RCP<Epetra_Map>& stdnodecolmap) const
{
  // each owner of a bin gets owner of the nodes this bin contains
  // all other nodes of elements, of which proc is owner of at least one
  // node, are ghosted
  Teuchos::RCP<Epetra_CrsGraph> initgraph = discret->BuildNodeGraph();

  // Todo introduced this export to column map due to special handling of
  //      beam nodes without own position DoFs in Utils::GetCurrentNodePos()
  Teuchos::RCP<Epetra_Vector> disnp_col = Teuchos::null;
  if (discret->HaveDofs() and disnp != Teuchos::null)
  {
    disnp_col = Teuchos::rcp(new Epetra_Vector(*discret->DofColMap()));
    Core::LinAlg::Export(*disnp, *disnp_col);
  }

  // distribute nodes, that are owned by a proc, to the bins of this proc
  std::map<int, std::vector<int>> nodesinbin;
  distribute_row_nodes_to_bins(discret, nodesinbin, disnp_col);

  std::map<int, std::vector<int>> nodesinmybins;
  // gather information of bin content from other procs (bin is owned by
  // this proc and there are some nodes on other procs which are located
  // in this bin here)
  collect_information_about_content_of_bins_from_other_procs_via_round_robin(
      rowbins, nodesinbin, nodesinmybins);

  // build new node row map
  std::vector<int> mynewrownodes;
  std::map<int, std::vector<int>>::const_iterator biniter;
  for (biniter = nodesinmybins.begin(); biniter != nodesinmybins.end(); ++biniter)
  {
    std::vector<int>::const_iterator nodeiter;
    for (nodeiter = biniter->second.begin(); nodeiter != biniter->second.end(); ++nodeiter)
    {
      mynewrownodes.push_back(*nodeiter);
    }
  }
  nodesinmybins.clear();

  Teuchos::RCP<Epetra_Map> newnoderowmap = Teuchos::rcp(
      new Epetra_Map(-1, mynewrownodes.size(), mynewrownodes.data(), 0, discret->Comm()));

  // create the new graph and export to it
  Teuchos::RCP<Epetra_CrsGraph> newnodegraph;

  newnodegraph = Teuchos::rcp(new Epetra_CrsGraph(Copy, *newnoderowmap, 108, false));
  Epetra_Export exporter(initgraph->RowMap(), *newnoderowmap);
  int err = newnodegraph->Export(*initgraph, exporter, Add);
  if (err < 0) FOUR_C_THROW("Graph export returned err=%d", err);
  newnodegraph->FillComplete();
  newnodegraph->OptimizeStorage();

  // the column map will become the new ghosted distribution of nodes (standard ghosting)
  const Epetra_BlockMap cntmp = newnodegraph->ColMap();
  stdnodecolmap = Teuchos::rcp(
      new Epetra_Map(-1, cntmp.NumMyElements(), cntmp.MyGlobalElements(), 0, discret->Comm()));

  // rebuild of the discretizations with new maps for standard ghosting
  Teuchos::RCP<Epetra_Map> roweles;
  std::tie(roweles, stdelecolmap) =
      discret->build_element_row_column(*newnoderowmap, *stdnodecolmap);
  discret->ExportRowNodes(*newnoderowmap);
  discret->ExportRowElements(*roweles);
  discret->ExportColumnNodes(*stdnodecolmap);
  discret->export_column_elements(*stdelecolmap);
  // in case we have a state vector, we need to build the dof map to enable its rebuild
  if (disnp == Teuchos::null)
  {
    discret->fill_complete(false, false, false);
  }
  else
  {
    discret->fill_complete(true, false, false);
    Teuchos::RCP<Epetra_Vector> old;
    old = disnp;
    disnp = Core::LinAlg::CreateVector(*discret->dof_row_map(), true);
    Core::LinAlg::Export(*old, *disnp);
  }

#ifdef FOUR_C_ENABLE_ASSERTIONS
  // print distribution after standard ghosting
  // some output after standard ghosting
  if (myrank_ == 0) std::cout << "parallel distribution with standard ghosting" << std::endl;
  Core::Rebalance::UTILS::print_parallel_distribution(*discret);
#endif
}

void Core::Binstrategy::BinningStrategy::
    collect_information_about_content_of_bins_from_other_procs_via_round_robin(
        Teuchos::RCP<Epetra_Map> rowbins, std::map<int, std::vector<int>>& mynodesinbins,
        std::map<int, std::vector<int>>& allnodesinmybins) const
{
  // do communication to gather all nodes
  const int numproc = rowbins->Comm().NumProc();
  for (int iproc = 0; iproc < numproc; ++iproc)
  {
    // vector with row bins on this proc
    std::vector<int> binids;
    int numbin;
    if (iproc == myrank_)
    {
      int* myrowbinsdata = rowbins->MyGlobalElements();
      numbin = rowbins->NumMyElements();
      binids.insert(binids.begin(), myrowbinsdata, myrowbinsdata + numbin);
    }

    // first: proc i tells all procs how many bins it has
    rowbins->Comm().Broadcast(&numbin, 1, iproc);
    binids.resize(numbin);
    // second: proc i tells all procs which bins it has, now each proc contains
    // rowbingids of iproc in vector binids
    rowbins->Comm().Broadcast(binids.data(), numbin, iproc);

    // loop over all own bins and find requested ones, fill in master elements in these bins
    // (map key is bin gid owned by iproc, vector contains all node gids of all procs in this bin)
    std::map<int, std::vector<int>> sdata;
    std::map<int, std::vector<int>> rdata;

    for (int i = 0; i < numbin; ++i)
    {
      // now each procs checks if row nodes lie in bins of iproc ...
      if (mynodesinbins.find(binids[i]) != mynodesinbins.end())
      {
        // ... if so, each proc assignes its node gids to iprocs bins
        sdata[binids[i]].insert(sdata[binids[i]].begin(), mynodesinbins[binids[i]].begin(),
            mynodesinbins[binids[i]].end());
      }
    }

    // iprocs gathers all this information from other procs
    Core::LinAlg::Gather<int>(sdata, rdata, 1, &iproc, rowbins->Comm());

    // iproc has to store the received data
    if (iproc == myrank_)
    {
      // clear data and refill
      allnodesinmybins.clear();
      allnodesinmybins.insert(rdata.begin(), rdata.end());
    }
  }
}

void Core::Binstrategy::BinningStrategy::revert_extended_ghosting(
    std::vector<Teuchos::RCP<Core::FE::Discretization>> dis,
    std::vector<Teuchos::RCP<Epetra_Map>>& stdelecolmap,
    std::vector<Teuchos::RCP<Epetra_Map>>& stdnodecolmap) const
{
  for (size_t i = 0; i < dis.size(); ++i)
  {
    //----------------------------
    // revert extended ghosting
    //----------------------------

    // adapt layout to standard ghosting in discret
    // first export the elements according to the processor local element column maps
    dis[i]->export_column_elements(*(stdelecolmap[i]));

    // now ghost the nodes
    dis[i]->ExportColumnNodes(*(stdnodecolmap[i]));

    // fillcomplete discret with standard ghosting
    dis[i]->fill_complete();
    if (myrank_ == 0) std::cout << "parallel distribution with reverted ghosting\n";
    Core::Rebalance::UTILS::print_parallel_distribution(*dis[i]);
  }
}

void Core::Binstrategy::BinningStrategy::
    compute_min_binning_domain_containing_all_elements_of_multiple_discrets(
        std::vector<Teuchos::RCP<Core::FE::Discretization>> discret,
        std::vector<Teuchos::RCP<const Epetra_Vector>> disnp,
        Core::LinAlg::Matrix<3, 2>& domain_bounding_box_corner_positions,
        bool set_bin_size_lower_bound_)
{
  // reset lower bound for bin size
  if (set_bin_size_lower_bound_) bin_size_lower_bound_ = 0.0;

  // safety check
  FOUR_C_THROW_UNLESS(discret[0]->NodeRowMap()->NumMyElements() > 0,
      "At least one proc does not even own at least one element, this leads to problems."
      " Choose less procs or change parallel distribution");

  // initialize XAABB_ as rectangle around the first node of first discret
  const Core::Nodes::Node* node = discret[0]->lRowNode(0);

  // calculate current position of this node
  double currpos[3] = {0.0, 0.0, 0.0};
  Utils::GetCurrentNodePos(discret[0], node, correct_beam_center_node_, disnp[0], currpos);

  for (int dim = 0; dim < 3; ++dim)
  {
    domain_bounding_box_corner_positions(dim, 0) = currpos[dim] - Core::Geo::TOL7;
    domain_bounding_box_corner_positions(dim, 1) = currpos[dim] + Core::Geo::TOL7;
  }

  // build XAABB_ from XAABB of all discrets and determine maximal element extension
  // to use as new set_bin_size_lower_bound_
  for (size_t i = 0; i < discret.size(); ++i)
  {
    Core::LinAlg::Matrix<3, 2> locXAABB;
    compute_min_binning_domain_containing_all_elements_of_single_discret(
        discret[i], locXAABB, disnp[i], set_bin_size_lower_bound_);

    // set XAABB_ considering all input discrets
    for (int dim = 0; dim < 3; ++dim)
    {
      domain_bounding_box_corner_positions(dim, 0) =
          std::min(domain_bounding_box_corner_positions(dim, 0), locXAABB(dim, 0));
      domain_bounding_box_corner_positions(dim, 1) =
          std::max(domain_bounding_box_corner_positions(dim, 1), locXAABB(dim, 1));
    }
  }

  // enlarge lower bound for bin size a little bit for safety reasons
  if (set_bin_size_lower_bound_) bin_size_lower_bound_ += Core::Geo::TOL7;
}

double Core::Binstrategy::BinningStrategy::
    compute_lower_bound_for_bin_size_as_max_edge_length_of_aabb_of_largest_ele(
        std::vector<Teuchos::RCP<Core::FE::Discretization>> discret,
        std::vector<Teuchos::RCP<const Epetra_Vector>> disnp)
{
  double bin_size_lower_bound = 0.0;

  // loop over all input discrets
  for (size_t ndis = 0; ndis < discret.size(); ++ndis)
  {
    // lower bound for bin size as largest element in discret
    double loc_max_bin_size_lower_bound = 0.0;
    double currpos[3] = {0.0, 0.0, 0.0};

    // loop over row elements of each proc
    for (int i = 0; i < discret[ndis]->NumMyRowElements(); ++i)
    {
      Core::Elements::Element* ele = discret[ndis]->lRowElement(i);

      // eleXAABB for each row element
      Core::LinAlg::Matrix<3, 2> eleXAABB(false);

      // initialize eleXAABB as rectangle around the first node of ele
      Utils::GetCurrentNodePos(
          discret[ndis], ele->Nodes()[0], correct_beam_center_node_, disnp[ndis], currpos);
      for (int dim = 0; dim < 3; ++dim)
      {
        eleXAABB(dim, 0) = currpos[dim] - Core::Geo::TOL7;
        eleXAABB(dim, 1) = currpos[dim] + Core::Geo::TOL7;
      }

      // rigid sphere elements needs to consider its radius
      if (element_filter_(ele) == Utils::SpecialElement::rigid_sphere)
      {
        for (int dim = 0; dim < 3; ++dim)
        {
          eleXAABB(dim, 0) = std::min(
              eleXAABB(dim, 0), eleXAABB(dim, 0) - rigid_sphere_radius_(ele) - Core::Geo::TOL7);
          eleXAABB(dim, 1) = std::max(
              eleXAABB(dim, 1), eleXAABB(dim, 0) + rigid_sphere_radius_(ele) + Core::Geo::TOL7);
        }
      }
      else
      {
        // loop over remaining nodes of current rowele
        for (int lid = 1; lid < ele->num_node(); ++lid)
        {
          const Core::Nodes::Node* node = ele->Nodes()[lid];
          Utils::GetCurrentNodePos(
              discret[ndis], node, correct_beam_center_node_, disnp[ndis], currpos);

          //  merge eleXAABB of all nodes of this element
          for (int dim = 0; dim < 3; ++dim)
          {
            eleXAABB(dim, 0) = std::min(eleXAABB(dim, 0), currpos[dim] - Core::Geo::TOL7);
            eleXAABB(dim, 1) = std::max(eleXAABB(dim, 1), currpos[dim] + Core::Geo::TOL7);
          }
        }
      }
      // compute lower bound for bin size as largest element in discret
      for (int dim = 0; dim < 3; ++dim)
        loc_max_bin_size_lower_bound =
            std::max(loc_max_bin_size_lower_bound, eleXAABB(dim, 1) - eleXAABB(dim, 0));
    }

    double globmax_bin_size_lower_bound = 0.0;
    discret[ndis]->Comm().MaxAll(&loc_max_bin_size_lower_bound, &globmax_bin_size_lower_bound, 1);
    // this is necessary if more than one discret is relevant
    bin_size_lower_bound = std::max(globmax_bin_size_lower_bound, bin_size_lower_bound);
  }

  return bin_size_lower_bound;
}

void Core::Binstrategy::BinningStrategy::
    create_bins_based_on_bin_size_lower_bound_and_binning_domain_dimensions(
        Teuchos::RCP<Core::FE::Discretization> dis)
{
  // create XAABB for discretization
  if (dis != Teuchos::null)
    compute_min_binning_domain_containing_all_elements_of_single_discret(
        dis, domain_bounding_box_corner_positions_);

  // divide global bounding box into bins
  for (int dim = 0; dim < 3; ++dim)
  {
    // determine number of bins per direction for prescribed bin_size_lower_bound_
    // std::floor leads to bins that are at least of size bin_size_lower_bound_
    if (bin_size_lower_bound_ > 0.0)
    {
      bin_per_dir_[dim] = std::max(1, (int)((domain_bounding_box_corner_positions_(dim, 1) -
                                                domain_bounding_box_corner_positions_(dim, 0)) /
                                            bin_size_lower_bound_));
    }

    // for detailed description of the difference between bin_per_dir
    // and id_calc_bin_per_dir_ see BinningStrategy::ConvertGidToijk;
    int n = 0;
    do
    {
      id_calc_bin_per_dir_[dim] = std::pow(2, n);
      id_calc_exp_bin_per_dir_[dim] = n;
      ++n;
    } while (id_calc_bin_per_dir_[dim] < bin_per_dir_[dim]);

    // calculate size of bins in each direction
    bin_size_[dim] = (domain_bounding_box_corner_positions_(dim, 1) -
                         domain_bounding_box_corner_positions_(dim, 0)) /
                     bin_per_dir_[dim];
    // calculate inverse of size of bins in each direction
    inv_bin_size_[dim] = 1.0 / bin_size_[dim];
  }

  FOUR_C_THROW_UNLESS(id_calc_bin_per_dir_[0] * id_calc_bin_per_dir_[1] * id_calc_bin_per_dir_[2] <
                          std::numeric_limits<int>::max(),
      "number of bins is larger than an integer can hold! Reduce number of bins by increasing "
      "the bin_size_lower_bound_");

  // determine lower bound for bin size if number of bins per dir was prescribed
  // was prescribed in input file
  if (bin_size_lower_bound_ <= 0.0)
    bin_size_lower_bound_ = std::min(bin_size_[0], std::min(bin_size_[1], bin_size_[2]));
}

void Core::Binstrategy::BinningStrategy::
    compute_min_binning_domain_containing_all_elements_of_single_discret(
        Teuchos::RCP<Core::FE::Discretization> discret, Core::LinAlg::Matrix<3, 2>& XAABB,
        Teuchos::RCP<const Epetra_Vector> disnp, bool set_bin_size_lower_bound_)
{
  // set_bin_size_lower_bound_ as largest element in discret on each proc
  double locmax_set_bin_size_lower_bound = 0.0;
  double currpos[3] = {0.0, 0.0, 0.0};
  // initialize XAABB of discret as rectangle around the first node of
  // discret on each proc
  Utils::GetCurrentNodePos(
      discret, discret->lRowNode(0), correct_beam_center_node_, disnp, currpos);
  for (int dim = 0; dim < 3; ++dim)
  {
    XAABB(dim, 0) = currpos[dim] - Core::Geo::TOL7;
    XAABB(dim, 1) = currpos[dim] + Core::Geo::TOL7;
  }

  // loop over row elements of each proc
  for (int i = 0; i < discret->NumMyRowElements(); ++i)
  {
    Core::Elements::Element* ele = discret->lRowElement(i);

    // eleXAABB for each row element
    Core::LinAlg::Matrix<3, 2> eleXAABB(false);

    // initialize eleXAABB as rectangle around the first node of ele
    Utils::GetCurrentNodePos(discret, ele->Nodes()[0], correct_beam_center_node_, disnp, currpos);
    for (int dim = 0; dim < 3; ++dim)
    {
      eleXAABB(dim, 0) = currpos[dim] - Core::Geo::TOL7;
      eleXAABB(dim, 1) = currpos[dim] + Core::Geo::TOL7;
    }

    // loop over remaining nodes of current rowele
    for (int lid = 1; lid < ele->num_node(); ++lid)
    {
      const Core::Nodes::Node* node = ele->Nodes()[lid];
      Utils::GetCurrentNodePos(discret, node, correct_beam_center_node_, disnp, currpos);

      //  merge eleXAABB of all nodes of this element
      for (int dim = 0; dim < 3; dim++)
      {
        eleXAABB(dim, 0) = std::min(eleXAABB(dim, 0), currpos[dim] - Core::Geo::TOL7);
        eleXAABB(dim, 1) = std::max(eleXAABB(dim, 1), currpos[dim] + Core::Geo::TOL7);
      }
    }

    // compute lower bound for bin size as largest element in discret
    if (set_bin_size_lower_bound_)
    {
      for (int dim = 0; dim < 3; ++dim)
        locmax_set_bin_size_lower_bound =
            std::max(locmax_set_bin_size_lower_bound, eleXAABB(dim, 1) - eleXAABB(dim, 0));
    }

    // merge XAABB of all roweles
    for (int dim = 0; dim < 3; dim++)
    {
      XAABB(dim, 0) = std::min(XAABB(dim, 0), eleXAABB(dim, 0));
      XAABB(dim, 1) = std::max(XAABB(dim, 1), eleXAABB(dim, 1));
    }
  }

  // local bounding box on each proc
  double locmin[3] = {XAABB(0, 0), XAABB(1, 0), XAABB(2, 0)};
  double locmax[3] = {XAABB(0, 1), XAABB(1, 1), XAABB(2, 1)};
  // global bounding box over all procs
  double globmin[3];
  double globmax[3];
  // do the necessary communication
  discret->Comm().MinAll(locmin, globmin, 3);
  discret->Comm().MaxAll(locmax, globmax, 3);

  // set global XAABB for discret
  for (int dim = 0; dim < 3; ++dim)
  {
    XAABB(dim, 0) = globmin[dim];
    XAABB(dim, 1) = globmax[dim];
  }

  // maxall of lower bound for cutoff
  if (set_bin_size_lower_bound_)
  {
    double globmax_set_bin_size_lower_bound = 0.0;
    discret->Comm().MaxAll(&locmax_set_bin_size_lower_bound, &globmax_set_bin_size_lower_bound, 1);
    // this is necessary if more than one discret is relevant
    bin_size_lower_bound_ = std::max(globmax_set_bin_size_lower_bound, bin_size_lower_bound_);
  }
}

void Core::Binstrategy::BinningStrategy::transfer_nodes_and_elements(
    Teuchos::RCP<Core::FE::Discretization>& discret, Teuchos::RCP<const Epetra_Vector> disnp,
    std::map<int, std::set<int>>& bintorowelemap)
{
  TEUCHOS_FUNC_TIME_MONITOR("Core::Binstrategy::BinningStrategy::transfer_nodes_and_elements");

  // clear map before setting up new one
  bintorowelemap.clear();

  // store elements that need treatment
  std::set<Core::Elements::Element*> elestoupdate;

  double currpos[3] = {0.0, 0.0, 0.0};
  // loop over all column nodes and check ownership
  for (int i = 0; i < discret->NodeColMap()->NumMyElements(); ++i)
  {
    // get current node and position
    Core::Nodes::Node* currnode = discret->lColNode(i);
    Utils::GetCurrentNodePos(discret, currnode, correct_beam_center_node_, disnp, currpos);

    int const gidofbin = ConvertPosToGid(currpos);

    if (bindis_->HaveGlobalElement(gidofbin))
    {
      int const hostbinowner = bindis_->gElement(gidofbin)->Owner();
      if (currnode->Owner() != hostbinowner)
      {
        // set new owner of node
        currnode->SetOwner(hostbinowner);
        // in case myrank is owner of associated element, add it to set
        Core::Elements::Element** curreles = currnode->Elements();
        for (int j = 0; j < currnode->NumElement(); ++j)
          if (curreles[j]->Owner() == myrank_) elestoupdate.insert(curreles[j]);
      }
    }
    /*else: in this case myrank was not owner of node and a corresponding element had and
    will only have ghost nodes on myrank, therefore we can leave the old owner because
    during the built up of the node col map all ghost nodes get deleted anyway  */
  }

  // store elements that need to be communicated
  std::map<int, std::vector<Core::Elements::Element*>> toranktosendeles;
  std::map<int, std::vector<std::pair<int, std::vector<int>>>> toranktosendbinids;

  // loop over row elements whose ownership may need an update
  std::set<Core::Elements::Element*>::const_iterator eleiter;
  for (eleiter = elestoupdate.begin(); eleiter != elestoupdate.end(); ++eleiter)
  {
    Core::Elements::Element* currele = *eleiter;
    Core::Nodes::Node** nodes = currele->Nodes();
    std::map<int, int> owner;
    for (int inode = 0; inode < currele->num_node(); ++inode) owner[nodes[inode]->Owner()] += 1;

    // check if any proc owns more nodes than myrank (for same number myrank_ stays owner)
    int newowner = myrank_;
    int numowned = (owner.find(myrank_) != owner.end()) ? owner[myrank_] : -1;
    std::map<int, int>::const_iterator i;
    for (i = owner.begin(); i != owner.end(); ++i)
      if (i->second > numowned) newowner = i->first;

    if (newowner != myrank_)
    {
      currele->SetOwner(newowner);
      toranktosendeles[newowner].push_back(currele);
      std::vector<int> binids;
      distribute_single_element_to_bins_using_ele_aabb(discret, currele, binids, disnp);
      std::pair<int, std::vector<int>> dummy(currele->Id(), binids);
      toranktosendbinids[newowner].push_back(dummy);
    }
  }

  // exploit bounding box idea for elements in underlying discretization and bins
  // loop over all row elements
  for (int lid = 0; lid < discret->NumMyRowElements(); ++lid)
  {
    Core::Elements::Element* eleptr = discret->lRowElement(lid);

    // check if owner did change
    if (eleptr->Owner() != myrank_) continue;

    // get corresponding bin ids in ijk range
    std::vector<int> binIds;
    distribute_single_element_to_bins_using_ele_aabb(discret, eleptr, binIds, disnp);

    // assign element to bins
    std::vector<int>::const_iterator biniter;
    for (biniter = binIds.begin(); biniter != binIds.end(); ++biniter)
      bintorowelemap[*biniter].insert(eleptr->Id());
  }

  // todo: send in one package
  // send and receive elements
  Utils::CommunicateElements(discret, toranktosendeles);
  // send and receive new elements to bin relation, like this no fillcomplete call necessary here
  Utils::CommunicateDistributionOfTransferredElementsToBins(
      discret, toranktosendbinids, bintorowelemap);
}

FOUR_C_NAMESPACE_CLOSE