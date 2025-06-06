// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_porofluid_pressure_based_elast_scatra_artery_coupling_base.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_utils_parameter_list.hpp"

#include <Teuchos_ParameterListExceptions.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
PoroPressureBased::PoroMultiPhaseScaTraArtCouplBase::PoroMultiPhaseScaTraArtCouplBase(
    std::shared_ptr<Core::FE::Discretization> arterydis,
    std::shared_ptr<Core::FE::Discretization> contdis, const Teuchos::ParameterList& couplingparams,
    const std::string& condname, const std::string& artcoupleddofname,
    const std::string& contcoupleddofname)
    : arterydis_(arterydis),
      contdis_(contdis),
      myrank_(Core::Communication::my_mpi_rank(arterydis->get_comm())),
      evaluate_in_ref_config_(Global::Problem::instance()
              ->poro_fluid_multi_phase_dynamic_params()
              .sublist("ARTERY COUPLING")
              .get<bool>("EVALUATE_IN_REF_CONFIG")),
      comm_(arterydis->get_comm())
{
  // safety check
  if (arterydis_->num_global_nodes() == 0)
    FOUR_C_THROW("artery discretization does not seem to have any nodes");

  // get the actual coupled DOFs
  // 1) 1D artery discretization
  int word1;
  int dummy = 0;
  std::istringstream coupled_art_dof_stream(
      Teuchos::getNumericStringParameter(couplingparams, artcoupleddofname));
  while (coupled_art_dof_stream >> word1)
  {
    // check ascending order
    if (dummy > 0)
      if ((int)(word1 - 1) <= coupleddofs_art_[dummy - 1])
        FOUR_C_THROW("DOFs have to be ordered in ascending order");
    coupleddofs_art_.push_back((int)(word1 - 1));
    dummy++;
  }

  // 2) 2D, 3D continuous field discretization
  dummy = 0;
  std::istringstream coupled_poro_dof_stream(
      Teuchos::getNumericStringParameter(couplingparams, contcoupleddofname));
  while (coupled_poro_dof_stream >> word1)
  {
    // check ascending order
    if (dummy > 0)
      if ((int)(word1 - 1) <= coupleddofs_cont_[dummy - 1])
        FOUR_C_THROW("DOFs have to be ordered in ascending order");
    coupleddofs_cont_.push_back((int)(word1 - 1));
    dummy++;
  }

  // no coupling selected by user
  if (coupleddofs_cont_.size() == 1 and coupleddofs_art_.size() == 1 and
      coupleddofs_cont_[0] < 0 and coupleddofs_art_[0] < 0)
  {
    coupleddofs_cont_.resize(0);
    coupleddofs_art_.resize(0);
  }

  if (coupleddofs_cont_.size() != coupleddofs_art_.size())
    FOUR_C_THROW("size mismatch between COUPLEDDOFS_ART and COUPLEDDOFS_PORO");

  num_coupled_dofs_ = coupleddofs_cont_.size();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroPressureBased::PoroMultiPhaseScaTraArtCouplBase::recompute_coupled_do_fs_for_ntp(
    std::vector<const Core::Conditions::Condition*> coupcond, unsigned int couplingnode)
{
  coupleddofs_art_ =
      (coupcond[couplingnode]->parameters().get<std::vector<int>>("COUPLEDDOF_REDUCED"));
  coupleddofs_cont_ =
      (coupcond[couplingnode]->parameters().get<std::vector<int>>("COUPLEDDOF_PORO"));

  // decrease the value of all elements by 1, because we start counting from 0 here and in the input
  // file we start from 1
  std::for_each(coupleddofs_art_.begin(), coupleddofs_art_.end(), [](int& value) { value--; });
  std::for_each(coupleddofs_cont_.begin(), coupleddofs_cont_.end(), [](int& value) { value--; });
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const std::shared_ptr<const Core::LinAlg::Map>&
PoroPressureBased::PoroMultiPhaseScaTraArtCouplBase::full_map() const
{
  return globalex_->full_map();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const std::shared_ptr<Core::LinAlg::MultiMapExtractor>&
PoroPressureBased::PoroMultiPhaseScaTraArtCouplBase::global_extractor() const
{
  return globalex_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroPressureBased::PoroMultiPhaseScaTraArtCouplBase::set_solution_vectors(
    std::shared_ptr<const Core::LinAlg::Vector<double>> phinp_cont,
    std::shared_ptr<const Core::LinAlg::Vector<double>> phin_cont,
    std::shared_ptr<const Core::LinAlg::Vector<double>> phinp_art)
{
  // do nothing
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroPressureBased::PoroMultiPhaseScaTraArtCouplBase::set_nearby_ele_pairs(
    const std::map<int, std::set<int>>* nearbyelepairs)
{
  // do nothing
}

FOUR_C_NAMESPACE_CLOSE
