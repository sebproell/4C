// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_porofluid_pressure_based_elast_scatra_artery_coupling_surfbased.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_utils_densematrix_communication.hpp"
#include "4C_porofluid_pressure_based_elast_scatra_artery_coupling_pair.hpp"


FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
PoroPressureBased::PorofluidElastScatraArteryCouplingSurfaceBasedAlgorithm::
    PorofluidElastScatraArteryCouplingSurfaceBasedAlgorithm(
        std::shared_ptr<Core::FE::Discretization> arterydis,
        std::shared_ptr<Core::FE::Discretization> contdis,
        const Teuchos::ParameterList& couplingparams, const std::string& condname,
        const std::string& artcoupleddofname, const std::string& contcoupleddofname)
    : PorofluidElastScatraArteryCouplingNonConformingAlgorithm(
          arterydis, contdis, couplingparams, condname, artcoupleddofname, contcoupleddofname)
{
  // user info
  if (my_mpi_rank_ == 0)
  {
    std::cout << "<                                                  >" << std::endl;
    print_coupling_method();
    std::cout << "<                                                  >" << std::endl;
    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    std::cout << "\n";
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroPressureBased::PorofluidElastScatraArteryCouplingSurfaceBasedAlgorithm::
    pre_evaluate_coupling_pairs()
{
  const int numpatch_axi = Global::Problem::instance()
                               ->poro_fluid_multi_phase_dynamic_params()
                               .sublist("ARTERY COUPLING")
                               .get<int>("NUMPATCH_AXI");
  const int numpatch_rad = Global::Problem::instance()
                               ->poro_fluid_multi_phase_dynamic_params()
                               .sublist("ARTERY COUPLING")
                               .get<int>("NUMPATCH_RAD");
  const int numartele = artery_dis_->num_global_elements();
  const int numgp_per_artele = numpatch_axi * numpatch_rad * 25;
  const int numgp_desired = numgp_per_artele * numartele;

  // this vector keeps track of evaluation of GPs
  std::shared_ptr<Core::LinAlg::MultiVector<double>> gp_vector =
      std::make_shared<Core::LinAlg::MultiVector<double>>(
          *artery_dis_->element_col_map(), numgp_per_artele);

  // pre-evaluate
  for (unsigned i = 0; i < coupled_ele_pairs_.size(); i++)
    coupled_ele_pairs_[i]->pre_evaluate(gp_vector);

  // delete the inactive pairs
  coupled_ele_pairs_.erase(std::remove_if(coupled_ele_pairs_.begin(), coupled_ele_pairs_.end(),
                               [](auto& pair) { return !pair->is_active(); }),
      coupled_ele_pairs_.end());

  // the following takes care of a very special case, namely, if a GP on the lateral surface lies
  // exactly in between two or more 3D elements owned by different procs.
  // In that case, the GP is duplicated across all owning procs.
  // We detect such cases and communicate them to all procs. below by adapting the
  // "gp_vector", to contain the multiplicity of the respective GP. Later the GP weight inside the
  // coupling pair is scaled by the inverse of the multiplicity

  int duplicates = 0;
  if (Core::Communication::num_mpi_ranks(get_comm()) > 1)
  {
    std::vector<int> mygpvec(numgp_per_artele, 0);
    std::vector<int> sumgpvec(numgp_per_artele, 0);
    // loop over all GIDs
    for (int gid = gp_vector->get_map().min_all_gid(); gid <= gp_vector->get_map().max_all_gid();
        gid++)
    {
      // reset
      std::fill(sumgpvec.data(), sumgpvec.data() + numgp_per_artele, 0);

      const int mylid = gp_vector->get_map().lid(gid);
      // if not owned or ghosted fill with zeros
      if (mylid < 0) std::fill(mygpvec.data(), mygpvec.data() + numgp_per_artele, 0);
      // else get the GP vector
      else
        for (int igp = 0; igp < numgp_per_artele; igp++)
          mygpvec[igp] = static_cast<int>(((*gp_vector)(igp))[mylid]);

      // communicate to all via summation
      Core::Communication::sum_all(mygpvec.data(), sumgpvec.data(), numgp_per_artele, get_comm());

      // this is ok for now, either the GID does not exist or the entire element protrudes.
      // Inform user and continue
      if (*std::max_element(sumgpvec.data(), sumgpvec.data() + numgp_per_artele) < 1)
      {
        std::cout << "WARNING! No GP of element  " << gid + 1 << " could be projected!"
                  << std::endl;
        continue;
      }

      // if one entry is equal to zero, this GP could not be projected
      if (*std::min_element(sumgpvec.data(), sumgpvec.data() + numgp_per_artele) < 1)
        FOUR_C_THROW("It seems as if one GP could not be projected");

      // find number of duplicates
      int sum = 0;
      sum = std::accumulate(sumgpvec.data(), sumgpvec.data() + numgp_per_artele, sum);
      duplicates += sum - numgp_per_artele;

      // if owned or ghosted by this proc. and if duplicates have been detected, replace entry in
      // gp_vector
      if (mylid >= 0 && sum > numgp_per_artele)
      {
        for (int igp = 0; igp < numgp_per_artele; igp++)
        {
          int err = gp_vector->ReplaceMyValue(mylid, igp, static_cast<double>(sumgpvec[igp]));
          if (err != 0) FOUR_C_THROW("ReplaceMyValue failed with error code {}!", err);
        }
      }
    }
  }

  for (unsigned i = 0; i < coupled_ele_pairs_.size(); i++)
    coupled_ele_pairs_[i]->delete_unnecessary_gps(gp_vector);

  int total_num_gp = 0;
  int numgp = 0;

  for (unsigned i = 0; i < coupled_ele_pairs_.size(); i++)
  {
    // segment ID not needed in this case, just set to zero
    coupled_ele_pairs_[i]->set_segment_id(0);
    numgp = numgp + coupled_ele_pairs_[i]->num_gp();
  }
  // safety check
  Core::Communication::sum_all(&numgp, &total_num_gp, 1, get_comm());
  if (numgp_desired != total_num_gp - duplicates)
    FOUR_C_THROW("It seems as if some GPs could not be projected");

  // output
  int total_numactive_pairs = 0;
  int numactive_pairs = static_cast<int>(coupled_ele_pairs_.size());
  Core::Communication::sum_all(&numactive_pairs, &total_numactive_pairs, 1, get_comm());
  if (homogenized_dis_->name() == "porofluid" && my_mpi_rank_ == 0)
    std::cout << "Only " << total_numactive_pairs
              << " Artery-to-PoroMultiphaseScatra coupling pairs are active" << std::endl;

  // print out summary of pairs
  if (homogenized_dis_->name() == "porofluid" &&
      coupling_params_.get<bool>("PRINT_OUT_SUMMARY_PAIRS"))
  {
    if (my_mpi_rank_ == 0)
      std::cout << "In total " << numgp_desired << " GPs (" << numgp_per_artele
                << " per artery element) required for lateral surface coupling" << std::endl;
    std::cout << "Proc. " << my_mpi_rank_ << " evaluates " << numgp << " GPs " << "("
              << (double)(numgp) / (double)(total_num_gp) * 100.0 << "% of all GPs)" << std::endl;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroPressureBased::PorofluidElastScatraArteryCouplingSurfaceBasedAlgorithm::setup()
{
  // call base class
  PoroPressureBased::PorofluidElastScatraArteryCouplingNonConformingAlgorithm::setup();

  // error-checks
  if (has_variable_diameter_)
    FOUR_C_THROW("Varying diameter not yet possible for surface-based coupling");
  if (!evaluate_in_ref_config_)
    FOUR_C_THROW("Evaluation in current configuration not yet possible for surface-based coupling");

  is_setup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroPressureBased::PorofluidElastScatraArteryCouplingSurfaceBasedAlgorithm::evaluate(
    std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase> sysmat,
    std::shared_ptr<Core::LinAlg::Vector<double>> rhs)
{
  if (!is_setup_) FOUR_C_THROW("setup() has not been called");

  if (!porofluid_managers_initialized_)
  {
    // pre-evaluate the pairs
    // --> has to be done here since the radius inside the material is required
    pre_evaluate_coupling_pairs();
  }

  // call base class
  PorofluidElastScatraArteryCouplingNonConformingAlgorithm::evaluate(sysmat, rhs);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroPressureBased::PorofluidElastScatraArteryCouplingSurfaceBasedAlgorithm::setup_system(
    std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase> sysmat,
    std::shared_ptr<Core::LinAlg::Vector<double>> rhs,
    std::shared_ptr<Core::LinAlg::SparseMatrix> sysmat_cont,
    std::shared_ptr<Core::LinAlg::SparseMatrix> sysmat_art,
    std::shared_ptr<const Core::LinAlg::Vector<double>> rhs_cont,
    std::shared_ptr<const Core::LinAlg::Vector<double>> rhs_art,
    std::shared_ptr<const Core::LinAlg::MapExtractor> dbcmap_cont,
    std::shared_ptr<const Core::LinAlg::MapExtractor> dbcmap_art)
{
  // call base class
  PorofluidElastScatraArteryCouplingNonConformingAlgorithm::setup_system(*sysmat, rhs, *sysmat_cont,
      *sysmat_art, rhs_cont, rhs_art, *dbcmap_cont, *dbcmap_art->cond_map(),
      *dbcmap_art->cond_map());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroPressureBased::PorofluidElastScatraArteryCouplingSurfaceBasedAlgorithm::
    apply_mesh_movement()
{
  if (!evaluate_in_ref_config_)
    FOUR_C_THROW("Evaluation in current configuration not possible for surface-based coupling");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::shared_ptr<const Core::LinAlg::Vector<double>> PoroPressureBased::
    PorofluidElastScatraArteryCouplingSurfaceBasedAlgorithm::blood_vessel_volume_fraction()
{
  FOUR_C_THROW("Output of vessel volume fraction not possible for surface-based coupling");
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroPressureBased::PorofluidElastScatraArteryCouplingSurfaceBasedAlgorithm::
    print_coupling_method() const
{
  std::cout << "<   surface-based formulation                      >" << '\n';
  PorofluidElastScatraArteryCouplingNonConformingAlgorithm::print_coupling_method();
}

FOUR_C_NAMESPACE_CLOSE
