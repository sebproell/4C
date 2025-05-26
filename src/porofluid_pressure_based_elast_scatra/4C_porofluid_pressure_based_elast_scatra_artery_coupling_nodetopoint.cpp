// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_porofluid_pressure_based_elast_scatra_artery_coupling_nodetopoint.hpp"

#include "4C_fem_condition_selector.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_utils_densematrix_communication.hpp"
#include "4C_porofluid_pressure_based_elast_scatra_artery_coupling_pair.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
PoroPressureBased::PorofluidElastScatraArteryCouplingNodeToPointAlgorithm::
    PorofluidElastScatraArteryCouplingNodeToPointAlgorithm(
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
void PoroPressureBased::PorofluidElastScatraArteryCouplingNodeToPointAlgorithm::setup()
{
  // call base class
  PorofluidElastScatraArteryCouplingNonConformingAlgorithm::setup();

  // pre-evaluate coupling pairs
  pre_evaluate_coupling_pairs();

  // print out summary of pairs
  if (homogenized_dis_->name() == "porofluid" &&
      coupling_params_.get<bool>("PRINT_OUT_SUMMARY_PAIRS"))
    output_coupling_pairs();

  // error-checks
  if (has_variable_diameter_)
    FOUR_C_THROW("Varying diameter not yet possible for node-to-point coupling");
  if (!evaluate_in_ref_config_)
    FOUR_C_THROW("Evaluation in current configuration not yet possible for node-to-point coupling");

  is_setup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroPressureBased::PorofluidElastScatraArteryCouplingNodeToPointAlgorithm::
    pre_evaluate_coupling_pairs()
{
  // pre-evaluate
  for (const auto& coupled_elepair : coupled_elepairs_) coupled_elepair->pre_evaluate(nullptr);

  // delete the inactive pairs
  coupled_elepairs_.erase(
      std::remove_if(coupled_elepairs_.begin(), coupled_elepairs_.end(),
          [](const std::shared_ptr<PoroMultiPhaseScatraArteryCouplingPairBase>& coupling_pair)
          { return not coupling_pair->is_active(); }),
      coupled_elepairs_.end());

  // output
  int total_num_active_pairs = 0;
  int num_active_pairs = static_cast<int>(coupled_elepairs_.size());
  Core::Communication::sum_all(&num_active_pairs, &total_num_active_pairs, 1, get_comm());
  if (my_mpi_rank_ == 0)
  {
    std::cout << total_num_active_pairs
              << " Artery-to-PoroMultiphaseScatra coupling pairs are active" << '\n';
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroPressureBased::PorofluidElastScatraArteryCouplingNodeToPointAlgorithm::evaluate(
    std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase> sysmat,
    std::shared_ptr<Core::LinAlg::Vector<double>> rhs)
{
  if (!is_setup_) FOUR_C_THROW("setup() has not been called");

  // call base class
  PorofluidElastScatraArteryCouplingNonConformingAlgorithm::evaluate(sysmat, rhs);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroPressureBased::PorofluidElastScatraArteryCouplingNodeToPointAlgorithm::setup_system(
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
void PoroPressureBased::PorofluidElastScatraArteryCouplingNodeToPointAlgorithm::
    apply_mesh_movement()
{
  if (!evaluate_in_ref_config_)
    FOUR_C_THROW("Evaluation in current configuration not possible for node-to-point coupling");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::shared_ptr<const Core::LinAlg::Vector<double>> PoroPressureBased::
    PorofluidElastScatraArteryCouplingNodeToPointAlgorithm::blood_vessel_volume_fraction()
{
  FOUR_C_THROW("Output of vessel volume fraction not possible for node-to-point coupling");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroPressureBased::PorofluidElastScatraArteryCouplingNodeToPointAlgorithm::
    print_coupling_method() const
{
  std::cout << "<Coupling-Method: 1D node to coincident point in 3D>" << '\n';
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroPressureBased::PorofluidElastScatraArteryCouplingNodeToPointAlgorithm::
    output_coupling_pairs() const
{
  if (my_mpi_rank_ == 0)
  {
    std::cout << "\nSummary of coupling pairs (segments):" << '\n';
    std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << '\n';
  }
  Core::Communication::barrier(get_comm());
  for (const auto& coupled_elepair : coupled_elepairs_)
  {
    std::cout << "Proc " << std::right << std::setw(2) << my_mpi_rank_ << ": Artery-ele "
              << std::right << std::setw(5) << coupled_elepair->ele1_gid()
              << ": <---> continuous-ele " << std::right << std::setw(7)
              << coupled_elepair->ele2_gid() << std::endl;
  }
  Core::Communication::barrier(get_comm());
  if (my_mpi_rank_ == 0) std::cout << "\n";
}

FOUR_C_NAMESPACE_CLOSE
