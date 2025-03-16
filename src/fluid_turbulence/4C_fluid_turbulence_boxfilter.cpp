// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fluid_turbulence_boxfilter.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_fluid_ele_action.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_mat_newtonianfluid.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_scatra_ele_action.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  Constructor (public)                                     krank 09/13|
 *----------------------------------------------------------------------*/
FLD::Boxfilter::Boxfilter(
    std::shared_ptr<Core::FE::Discretization> actdis, Teuchos::ParameterList& params)
    :  // call constructor for "nontrivial" objects
      discret_(actdis),
      params_(params),
      physicaltype_(
          Teuchos::getIntegralValue<Inpar::FLUID::PhysicalType>(params_, "Physical Type")),
      //  available control settings
      apply_dynamic_smagorinsky_(false),
      vreman_dynamic_(false),
      apply_box_filter_(false),
      loma_(false),
      incomp_(false),
      velocity_(false),
      reynoldsstress_(false),
      modeled_subgrid_stress_(false),
      expression_(false),
      strainrate_(false),
      alphaij_(false),
      alpha2_(false),
      finescale_velocity_(false),
      densvelocity_(false),
      densstrainrate_(false),
      density_(false),
      phi_(false),
      phi2_(false),
      phiexpression_(false),
      alphaijsc_(false)
{
  Teuchos::ParameterList* modelparams = &(params_.sublist("TURBULENCE MODEL"));

  if (modelparams->get<std::string>("TURBULENCE_APPROACH", "DNS_OR_RESVMM_LES") == "CLASSICAL_LES")
  {
    if (modelparams->get<std::string>("PHYSICAL_MODEL", "no_model") == "Dynamic_Smagorinsky")
      apply_dynamic_smagorinsky_ = true;

    if (modelparams->get<std::string>("PHYSICAL_MODEL", "no_model") ==
        "Multifractal_Subgrid_Scales")
      apply_box_filter_ = true;

    if (physicaltype_ == Inpar::FLUID::loma) loma_ = true;

    if (physicaltype_ == Inpar::FLUID::incompressible) incomp_ = true;

    if (modelparams->get<std::string>("PHYSICAL_MODEL", "no_model") == "Dynamic_Vreman")
      vreman_dynamic_ = true;
  }

  dynsmag_loma_on_ = (loma_ and apply_dynamic_smagorinsky_);

  if (apply_dynamic_smagorinsky_)
  {
    velocity_ = true;
    reynoldsstress_ = true;
    modeled_subgrid_stress_ = true;
  }
  if (apply_box_filter_)
  {
    velocity_ = true;
    reynoldsstress_ = true;
    finescale_velocity_ = true;
  }
  if (dynsmag_loma_on_)
  {
    densvelocity_ = true;
    densstrainrate_ = true;
    density_ = true;
    velocity_ = true;
    reynoldsstress_ = true;
    modeled_subgrid_stress_ = true;
  }
  if (loma_ and vreman_dynamic_) FOUR_C_THROW("Dynamic Vreman model not implemented for loma!");

  return;
}



/*----------------------------------------------------------------------*
 | add some scatra specific parameters                  rasthofer 08/12 |
 * ---------------------------------------------------------------------*/
void FLD::Boxfilter::add_scatra(std::shared_ptr<Core::FE::Discretization> scatradis)
{
  scatradiscret_ = scatradis;

  return;
}

void FLD::Boxfilter::initialize_vreman()
{
  strainrate_ = true;
  expression_ = true;
  alphaij_ = true;
  alpha2_ = true;

  return;
}

void FLD::Boxfilter::initialize_vreman_scatra(std::shared_ptr<Core::FE::Discretization> scatradis)
{
  scatradiscret_ = scatradis;

  phi_ = true;
  phi2_ = true;
  phiexpression_ = true;
  alphaijsc_ = true;
  return;
}


/*---------------------------------------------------------------------*
 | Perform box filter operation                                        |
 *---------------------------------------------------------------------*/
void FLD::Boxfilter::apply_filter(
    const std::shared_ptr<const Core::LinAlg::Vector<double>> velocity,
    const std::shared_ptr<const Core::LinAlg::Vector<double>> scalar, const double thermpress,
    const std::shared_ptr<const Core::LinAlg::Vector<double>> dirichtoggle)
{
  // perform filtering depending on the LES model
  apply_box_filter(velocity, scalar, thermpress, *dirichtoggle);

  return;
}

void FLD::Boxfilter::apply_filter_scatra(
    const std::shared_ptr<const Core::LinAlg::Vector<double>> scalar, const double thermpress,
    const std::shared_ptr<const Core::LinAlg::Vector<double>> dirichtoggle, const int ndsvel)
{
  // perform filtering depending on the LES model
  apply_box_filter_scatra(scalar, thermpress, *dirichtoggle, ndsvel);

  return;
}

/*----------------------------------------------------------------------*
 | perform box filtering                                      (private) |
 |                                                            rasthofer |
 *----------------------------------------------------------------------*/
void FLD::Boxfilter::apply_box_filter(
    const std::shared_ptr<const Core::LinAlg::Vector<double>> velocity,
    const std::shared_ptr<const Core::LinAlg::Vector<double>> scalar, const double thermpress,
    const Core::LinAlg::Vector<double>& dirichtoggle)
{
  TEUCHOS_FUNC_TIME_MONITOR("apply_filter_for_dynamic_computation_of_cs");

  // LES turbulence modeling is only valid for 3 dimensions
  const int numdim = 3;

  // generate a parameterlist for communication and control
  Teuchos::ParameterList filterparams;
  // action for elements
  filterparams.set<FLD::Action>("action", FLD::calc_fluid_box_filter);
  filterparams.set("thermpress", thermpress);

  // set state vector to pass distributed vector to the element
  discret_->clear_state();
  discret_->set_state("u and p (trial)", velocity);
  discret_->set_state("T (trial)", scalar);

  // dummies
  Core::LinAlg::SerialDenseMatrix emat1;
  Core::LinAlg::SerialDenseMatrix emat2;
  Core::LinAlg::SerialDenseVector evec1;
  Core::LinAlg::SerialDenseVector evec2;
  Core::LinAlg::SerialDenseVector evec3;

  // ---------------------------------------------------------------
  // get a vector layout from the discretization to construct
  const Epetra_Map* noderowmap = discret_->node_row_map();

  // alloc an additional vector to store/add up the patch volume
  Core::LinAlg::Vector<double> patchvol(*noderowmap, true);

  // free mem and reallocate to zero out vecs
  if (velocity_) filtered_vel_ = nullptr;
  if (reynoldsstress_) filtered_reynoldsstress_ = nullptr;
  if (modeled_subgrid_stress_) filtered_modeled_subgrid_stress_ = nullptr;
  if (densvelocity_) filtered_dens_vel_ = nullptr;
  if (density_) filtered_dens_ = nullptr;
  if (densstrainrate_) filtered_dens_strainrate_ = nullptr;
  if (finescale_velocity_) fs_vel_ = nullptr;
  if (strainrate_) filtered_strainrate_ = nullptr;
  if (expression_) filtered_expression_ = nullptr;
  if (alphaij_) filtered_alphaij_ = nullptr;
  if (alpha2_) filtered_alpha2_ = nullptr;

  if (velocity_)
    filtered_vel_ = std::make_shared<Core::LinAlg::MultiVector<double>>(*noderowmap, numdim, true);
  if (reynoldsstress_)
    filtered_reynoldsstress_ =
        std::make_shared<Core::LinAlg::MultiVector<double>>(*noderowmap, numdim * numdim, true);
  if (modeled_subgrid_stress_)
    filtered_modeled_subgrid_stress_ =
        std::make_shared<Core::LinAlg::MultiVector<double>>(*noderowmap, numdim * numdim, true);
  if (densvelocity_)
    filtered_dens_vel_ =
        std::make_shared<Core::LinAlg::MultiVector<double>>(*noderowmap, numdim, true);
  if (density_) filtered_dens_ = std::make_shared<Core::LinAlg::Vector<double>>(*noderowmap, true);
  if (densstrainrate_)
    filtered_dens_strainrate_ = std::make_shared<Core::LinAlg::Vector<double>>(*noderowmap, true);
  if (strainrate_)
    filtered_strainrate_ =
        std::make_shared<Core::LinAlg::MultiVector<double>>(*noderowmap, numdim * numdim, true);
  if (expression_)
    filtered_expression_ = std::make_shared<Core::LinAlg::Vector<double>>(*noderowmap, true);
  if (alphaij_)
    filtered_alphaij_ =
        std::make_shared<Core::LinAlg::MultiVector<double>>(*noderowmap, numdim * numdim, true);
  if (alpha2_) filtered_alpha2_ = std::make_shared<Core::LinAlg::Vector<double>>(*noderowmap, true);

  if (finescale_velocity_)
    fs_vel_ = std::make_shared<Core::LinAlg::MultiVector<double>>(*noderowmap, numdim, true);

  // ---------------------------------------------------------------
  // do the integration of the (not normalized) box filter function
  // on the element

  // loop all elements on this proc (including ghosted ones)
  for (int nele = 0; nele < discret_->num_my_col_elements(); ++nele)
  {
    // get the element
    Core::Elements::Element* ele = discret_->l_col_element(nele);

    // provide vectors for filtered quantities
    std::shared_ptr<std::vector<double>> vel_hat =
        std::make_shared<std::vector<double>>((numdim), 0.0);
    std::shared_ptr<std::vector<std::vector<double>>> reynoldsstress_hat =
        std::make_shared<std::vector<std::vector<double>>>();
    std::shared_ptr<std::vector<std::vector<double>>> modeled_subgrid_stress =
        std::make_shared<std::vector<std::vector<double>>>();
    // set to dimensions
    if (reynoldsstress_) (*reynoldsstress_hat).resize(numdim);
    if (modeled_subgrid_stress_) (*modeled_subgrid_stress).resize(numdim);
    for (int rr = 0; rr < numdim; rr++)
    {
      if (reynoldsstress_) ((*reynoldsstress_hat)[rr]).resize(numdim);
      if (modeled_subgrid_stress_) ((*modeled_subgrid_stress)[rr]).resize(numdim);
    }
    // initialize with zeros
    for (int rr = 0; rr < numdim; rr++)
    {
      for (int ss = 0; ss < numdim; ss++)
      {
        if (reynoldsstress_) (*reynoldsstress_hat)[rr][ss] = 0.0;
        if (modeled_subgrid_stress_) (*modeled_subgrid_stress)[rr][ss] = 0.0;
      }
    }
    std::shared_ptr<std::vector<double>> densvel_hat =
        std::make_shared<std::vector<double>>((numdim), 0.0);
    // and set them in parameter list
    filterparams.set<std::shared_ptr<std::vector<double>>>("densvel_hat", densvel_hat);

    filterparams.set<std::shared_ptr<std::vector<double>>>("vel_hat", vel_hat);
    filterparams.set<std::shared_ptr<std::vector<std::vector<double>>>>(
        "reynoldsstress_hat", reynoldsstress_hat);
    filterparams.set<std::shared_ptr<std::vector<std::vector<double>>>>(
        "modeled_subgrid_stress", modeled_subgrid_stress);

    // Vreman_initialization
    std::shared_ptr<std::vector<std::vector<double>>> strainrate_hat =
        std::make_shared<std::vector<std::vector<double>>>();
    std::shared_ptr<std::vector<std::vector<double>>> alphaij_hat =
        std::make_shared<std::vector<std::vector<double>>>();
    if (strainrate_) (*strainrate_hat).resize(numdim);
    if (alphaij_) (*alphaij_hat).resize(numdim);
    for (int rr = 0; rr < numdim; rr++)
    {
      if (strainrate_) ((*strainrate_hat)[rr]).resize(numdim);
      if (alphaij_) ((*alphaij_hat)[rr]).resize(numdim);
    }
    // initialize with zeros
    for (int rr = 0; rr < numdim; rr++)
    {
      for (int ss = 0; ss < numdim; ss++)
      {
        if (strainrate_) (*strainrate_hat)[rr][ss] = 0.0;
        if (alphaij_) (*alphaij_hat)[rr][ss] = 0.0;
      }
    }

    // if(strainrate_)
    filterparams.set<std::shared_ptr<std::vector<std::vector<double>>>>(
        "strainrate_hat", strainrate_hat);
    // if(alphaij_)
    filterparams.set<std::shared_ptr<std::vector<std::vector<double>>>>("alphaij_hat", alphaij_hat);

    // get element location vector, dirichlet flags and ownerships
    std::vector<int> lm;
    std::vector<int> lmowner;
    std::vector<int> lmstride;
    ele->location_vector(*discret_, lm, lmowner, lmstride);

    // call the element evaluate method to integrate functions
    // against heaviside function element
    int err = ele->evaluate(filterparams, *discret_, lm, emat1, emat2, evec1, evec2, evec2);
    if (err)
      FOUR_C_THROW("Proc {}: Element {} returned err={}",
          Core::Communication::my_mpi_rank(discret_->get_comm()), ele->id(), err);

    // get contribution to patch volume of this element. Add it up.
    double volume_contribution = filterparams.get<double>("volume_contribution");
    // if(density_)
    double dens_hat = filterparams.get<double>("dens_hat");
    // if(densstrainrate_)
    double dens_strainrate_hat = filterparams.get<double>("dens_strainrate_hat");

    double expression_hat = filterparams.get<double>("expression_hat");
    double alpha2_hat = filterparams.get<double>("alpha2_hat");

    // loop all nodes of this element, add values to the global vectors
    Core::Nodes::Node** elenodes = ele->nodes();
    for (int nn = 0; nn < ele->num_node(); ++nn)
    {
      Core::Nodes::Node* node = (elenodes[nn]);

      // we are interested only in  row nodes
      if (node->owner() == Core::Communication::my_mpi_rank(discret_->get_comm()))
      {
        // now assemble the computed values into the global vector
        int id = (node->id());

        patchvol.sum_into_global_values(1, &volume_contribution, &id);


        if (density_) filtered_dens_->sum_into_global_values(1, &dens_hat, &id);
        if (densstrainrate_)
          filtered_dens_strainrate_->sum_into_global_values(1, &dens_strainrate_hat, &id);
        if (expression_) filtered_expression_->sum_into_global_values(1, &expression_hat, &id);
        if (alpha2_) filtered_alpha2_->sum_into_global_values(1, &alpha2_hat, &id);

        for (int idim = 0; idim < numdim; ++idim)
        {
          if (velocity_)
          {
            double val = (*vel_hat)[idim];
            ((*filtered_vel_)(idim)).sum_into_global_values(1, &val, &id);
          }
          if (densvelocity_)
          {
            double val = (*densvel_hat)[idim];
            ((*filtered_dens_vel_)(idim)).sum_into_global_values(1, &val, &id);
          }


          for (int jdim = 0; jdim < numdim; ++jdim)
          {
            const int ij = numdim * idim + jdim;
            if (reynoldsstress_)
            {
              double val = (*reynoldsstress_hat)[idim][jdim];
              ((*filtered_reynoldsstress_)(ij)).sum_into_global_values(1, &val, &id);
            }
            if (modeled_subgrid_stress_)
            {
              double val = (*modeled_subgrid_stress)[idim][jdim];
              ((*filtered_modeled_subgrid_stress_)(ij)).sum_into_global_values(1, &val, &id);
            }
            if (strainrate_)
            {
              double val = (*strainrate_hat)[idim][jdim];
              ((*filtered_strainrate_)(ij)).sum_into_global_values(1, &val, &id);
            }
            if (alphaij_)
            {
              double val = (*alphaij_hat)[idim][jdim];
              ((*filtered_alphaij_)(ij)).sum_into_global_values(1, &val, &id);
            }
          }
        }
      }
    }
  }  // end elementloop

  // ---------------------------------------------------------------
  // send add values from masters and slaves
  {
    double val;
    // if(velocity_)
    std::vector<double> vel_val(3);
    std::vector<std::vector<double>> reystress_val;
    if (reynoldsstress_)
    {
      reystress_val.resize(3);
      for (int rr = 0; rr < 3; rr++) (reystress_val[rr]).resize(3);
    }
    std::vector<std::vector<double>> modeled_subgrid_stress_val;
    if (modeled_subgrid_stress_)
    {
      modeled_subgrid_stress_val.resize(3);
      for (int rr = 0; rr < 3; rr++) (modeled_subgrid_stress_val[rr]).resize(3);
    }
    std::vector<std::vector<double>> strainrate_val;
    if (strainrate_)
    {
      strainrate_val.resize(3);
      for (int rr = 0; rr < 3; rr++) (strainrate_val[rr]).resize(3);
    }
    std::vector<std::vector<double>> alphaij_val;
    if (alphaij_)
    {
      alphaij_val.resize(3);
      for (int rr = 0; rr < 3; rr++) (alphaij_val[rr]).resize(3);
    }
    // loma specific quantities
    std::vector<double> dens_vel_val(3);
    double dens_val;
    double dens_strainrate_val;
    double expression_val;
    double alpha2_val;

    std::map<int, std::vector<int>>* pbcmapmastertoslave =
        discret_->get_all_pbc_coupled_col_nodes();
    // loop all master nodes on this proc
    if (pbcmapmastertoslave)
    {
      for (const auto& [master_gid, slave_gids] : *pbcmapmastertoslave)
      {
        // loop only owned nodes
        if ((discret_->g_node(master_gid))->owner() !=
            Core::Communication::my_mpi_rank(discret_->get_comm()))
          continue;

        int lid = noderowmap->LID(master_gid);
        if (lid < 0) FOUR_C_THROW("nodelid < 0 ?");

        val = (patchvol)[lid];

        if (density_) dens_val = (*filtered_dens_)[lid];
        if (densstrainrate_) dens_strainrate_val = (*filtered_dens_strainrate_)[lid];
        if (expression_) expression_val = (*filtered_expression_)[lid];
        if (alpha2_) alpha2_val = (*filtered_alpha2_)[lid];

        for (int idim = 0; idim < numdim; ++idim)
        {
          if (velocity_) vel_val[idim] = ((((*filtered_vel_)(idim)))[lid]);

          if (densvelocity_) dens_vel_val[idim] = ((((*filtered_dens_vel_)(idim)))[lid]);

          for (int jdim = 0; jdim < numdim; ++jdim)
          {
            const int ij = numdim * idim + jdim;
            if (reynoldsstress_)
              reystress_val[idim][jdim] = (((*filtered_reynoldsstress_)(ij)))[lid];
            if (modeled_subgrid_stress_)
              modeled_subgrid_stress_val[idim][jdim] =
                  (((*filtered_modeled_subgrid_stress_)(ij)))[lid];
            if (strainrate_) strainrate_val[idim][jdim] = (((*filtered_strainrate_)(ij)))[lid];
            if (alphaij_) alphaij_val[idim][jdim] = (((*filtered_alphaij_)(ij)))[lid];
          }
        }

        // loop all this masters slaves
        for (auto slave_gid : slave_gids)
        {
          lid = noderowmap->LID(slave_gid);
          val += (patchvol)[lid];

          if (density_) dens_val += (*filtered_dens_)[lid];
          if (densstrainrate_) dens_strainrate_val += (*filtered_dens_strainrate_)[lid];
          if (expression_) expression_val += (*filtered_expression_)[lid];
          if (alpha2_) alpha2_val += (*filtered_alpha2_)[lid];

          for (int idim = 0; idim < numdim; ++idim)
          {
            if (velocity_) vel_val[idim] += ((((*filtered_vel_)(idim)))[lid]);

            if (densvelocity_) dens_vel_val[idim] += ((((*filtered_dens_vel_)(idim)))[lid]);

            for (int jdim = 0; jdim < numdim; ++jdim)
            {
              const int ij = numdim * idim + jdim;
              if (reynoldsstress_)
                reystress_val[idim][jdim] += (((*filtered_reynoldsstress_)(ij)))[lid];
              if (modeled_subgrid_stress_)
                modeled_subgrid_stress_val[idim][jdim] +=
                    (((*filtered_modeled_subgrid_stress_)(ij)))[lid];
              if (strainrate_) strainrate_val[idim][jdim] += (((*filtered_strainrate_)(ij)))[lid];
              if (alphaij_) alphaij_val[idim][jdim] += (((*filtered_alphaij_)(ij)))[lid];
            }  // end loop jdim
          }  // end loop idim
        }  // end loop slaves

        // replace value by sum
        lid = noderowmap->LID(master_gid);
        int error = patchvol.replace_local_values(1, &val, &lid);
        if (error != 0) FOUR_C_THROW("dof not on proc");

        {
          int err = 0;
          if (density_) err += filtered_dens_->replace_local_values(1, &dens_val, &lid);
          if (densstrainrate_)
            err += filtered_dens_strainrate_->replace_local_values(1, &dens_strainrate_val, &lid);
          if (expression_)
            err += filtered_expression_->replace_local_values(1, &expression_val, &lid);
          if (alpha2_) err += filtered_alpha2_->replace_local_values(1, &alpha2_val, &lid);
          if (err != 0) FOUR_C_THROW("dof not on proc");
        }

        for (int idim = 0; idim < numdim; ++idim)
        {
          int err = 0;
          if (velocity_)
            err += ((*filtered_vel_)(idim)).replace_local_values(1, &(vel_val[idim]), &lid);

          if (densvelocity_)
            err +=
                ((*filtered_dens_vel_)(idim)).replace_local_values(1, &(dens_vel_val[idim]), &lid);

          for (int jdim = 0; jdim < numdim; ++jdim)
          {
            const int ij = numdim * idim + jdim;
            if (reynoldsstress_)
              err += ((*filtered_reynoldsstress_)(ij))
                         .replace_local_values(1, &(reystress_val[idim][jdim]), &lid);
            if (modeled_subgrid_stress_)
              err += ((*filtered_modeled_subgrid_stress_)(ij))
                         .replace_local_values(1, &(modeled_subgrid_stress_val[idim][jdim]), &lid);
            if (strainrate_)
              err += ((*filtered_strainrate_)(ij))
                         .replace_local_values(1, &(strainrate_val[idim][jdim]), &lid);
            if (alphaij_)
              err += ((*filtered_alphaij_)(ij))
                         .replace_local_values(1, &(alphaij_val[idim][jdim]), &lid);
          }  // end loop jdim
          if (err != 0) FOUR_C_THROW("dof not on proc");
        }  // end loop idim

        // loop all this masters slaves
        for (auto slave_gid : slave_gids)
        {
          int err = 0;
          lid = noderowmap->LID(slave_gid);
          err += patchvol.replace_local_values(1, &val, &lid);

          {
            int err = 0;
            if (density_) err += filtered_dens_->replace_local_values(1, &dens_val, &lid);
            if (densstrainrate_)
              err += filtered_dens_strainrate_->replace_local_values(1, &dens_strainrate_val, &lid);
            if (expression_)
              err += filtered_expression_->replace_local_values(1, &expression_val, &lid);
            if (alpha2_) err += filtered_alpha2_->replace_local_values(1, &alpha2_val, &lid);
            if (err != 0) FOUR_C_THROW("dof not on proc");
          }

          for (int idim = 0; idim < numdim; ++idim)
          {
            if (velocity_)
              err += ((*filtered_vel_)(idim)).replace_local_values(1, &(vel_val[idim]), &lid);

            if (densvelocity_)
              err += ((*filtered_dens_vel_)(idim))
                         .replace_local_values(1, &(dens_vel_val[idim]), &lid);

            for (int jdim = 0; jdim < numdim; ++jdim)
            {
              const int ij = numdim * idim + jdim;
              if (reynoldsstress_)
                err += ((*filtered_reynoldsstress_)(ij))
                           .replace_local_values(1, &(reystress_val[idim][jdim]), &lid);
              if (modeled_subgrid_stress_)
                err +=
                    ((*filtered_modeled_subgrid_stress_)(ij))
                        .replace_local_values(1, &(modeled_subgrid_stress_val[idim][jdim]), &lid);
              if (strainrate_)
                err += ((*filtered_strainrate_)(ij))
                           .replace_local_values(1, &(strainrate_val[idim][jdim]), &lid);
              if (alphaij_)
                err += ((*filtered_alphaij_)(ij))
                           .replace_local_values(1, &(alphaij_val[idim][jdim]), &lid);
            }  // end loop jdim
          }  // end loop idim
          if (err != 0) FOUR_C_THROW("dof not on proc");
        }  // end loop slaves
      }  // end loop masters
    }
  }

  // ---------------------------------------------------------------
  // replace values at dirichlet nodes

  {
    // get a rowmap for the dofs
    const Epetra_Map* dofrowmap = discret_->dof_row_map();

    // loop all nodes on the processor
    for (int lnodeid = 0; lnodeid < discret_->num_my_row_nodes(); ++lnodeid)
    {
      // get the processor local node
      Core::Nodes::Node* lnode = discret_->l_row_node(lnodeid);

      // the set of degrees of freedom associated with the node
      std::vector<int> nodedofset = discret_->dof(lnode);

      // check whether the node is on a wall, i.e. all velocity dofs
      // are Dirichlet constrained
      int is_dirichlet_node = 0;
      int is_no_slip_node = 0;
      for (int index = 0; index < numdim; ++index)
      {
        int gid = nodedofset[index];
        int lid = dofrowmap->LID(gid);

        if ((dirichtoggle)[lid] == 1)  // this is a dirichlet node
        {
          is_dirichlet_node++;
          double vel_i = (*velocity)[lid];
          if (abs(vel_i) < 1e-12)
          {
            is_no_slip_node++;
          }
        }
      }

      // this node is on a dirichlet boundary
      if (is_dirichlet_node == numdim)
      {
        int err = 0;

        // determine volume
        double thisvol = (patchvol)[lnodeid];

        // determine density
        double dens = 1.0;
        if (physicaltype_ == Inpar::FLUID::incompressible)  // this is important to have here,
        {  //  since, for the pure box filter application,
           // get fluid viscosity from material definition                                //  we do
           // not want to multiply the reynolds stress by density
          int id =
              Global::Problem::instance()->materials()->first_id_by_type(Core::Materials::m_fluid);
          if (id == -1)
            FOUR_C_THROW("Could not find Newtonian fluid material");
          else
          {
            const Core::Mat::PAR::Parameter* mat =
                Global::Problem::instance()->materials()->parameter_by_id(id);
            const Mat::PAR::NewtonianFluid* actmat =
                static_cast<const Mat::PAR::NewtonianFluid*>(mat);
            // we need the kinematic viscosity here
            dens = actmat->density_;
          }
        }
        if (density_) dens = (*filtered_dens_)[lnodeid] / thisvol;


        // set density (only required for loma)
        // set value to mean value
        // we already divide by the corresponding volume of all contributing elements,
        // since we set the volume to 1.0 in the next step in order not to modify the dirichlet
        // values
        if (density_) err += filtered_dens_->replace_local_values(1, &dens, &lnodeid);

        // this node is on a wall
        if (is_no_slip_node == numdim)
        {
          // Peter style
          double val = 0.0;
          if (densstrainrate_)
            err += filtered_dens_strainrate_->replace_local_values(1, &val, &lnodeid);
          if (expression_) err += filtered_expression_->replace_local_values(1, &val, &lnodeid);
          if (alpha2_) err += filtered_alpha2_->replace_local_values(1, &val, &lnodeid);
        }
        else
        {
          if (densstrainrate_)
          {
            double val = (*filtered_dens_strainrate_)[lnodeid] / thisvol;
            err += filtered_dens_strainrate_->replace_local_values(1, &val, &lnodeid);
          }
          if (expression_)
          {
            double val = (*filtered_expression_)[lnodeid] / thisvol;
            err += filtered_expression_->replace_local_values(1, &val, &lnodeid);
          }
          if (alpha2_)
          {
            double val = (*filtered_alpha2_)[lnodeid] / thisvol;
            err += filtered_alpha2_->replace_local_values(1, &val, &lnodeid);
          }
        }



        for (int idim = 0; idim < numdim; ++idim)
        {
          int gid_i = nodedofset[idim];
          int lid_i = dofrowmap->LID(gid_i);
          double valvel_i = (*velocity)[lid_i];
          if (velocity_)
          {
            err += ((*filtered_vel_)(idim)).replace_local_values(1, &valvel_i, &lnodeid);
          }
          // dens*reynoldsstress not in parameter list until now?
          if (densvelocity_)  //=loma
          {
            // note: for incompressible flow, this vector is rebuild in calculation of Lij and Mij
            double valdensvel_i = dens * valvel_i;
            err += ((*filtered_dens_vel_)(idim)).replace_local_values(1, &valdensvel_i, &lnodeid);
          }

          for (int jdim = 0; jdim < numdim; ++jdim)
          {
            const int ij = numdim * idim + jdim;

            if (reynoldsstress_)
            {
              int gid_j = nodedofset[jdim];
              int lid_j = dofrowmap->LID(gid_j);

              double valvel_j = (*velocity)[lid_j];
              double valvel_ij = dens * valvel_i * valvel_j;
              // remember: density = 1.0 for pure box filter application
              err +=
                  ((*filtered_reynoldsstress_)(ij)).replace_local_values(1, &valvel_ij, &lnodeid);
            }

            if (is_no_slip_node == numdim)
            {
              // set value to zero (original Peter style)
              double val = 0.0;
              if (modeled_subgrid_stress_)
                err += ((*filtered_modeled_subgrid_stress_)(ij))
                           .replace_local_values(1, &val, &lnodeid);
              // remark: setting the modeled stresses equal to zero improves the estimated friction
              // Reynolds number!
              if (strainrate_)
                err += ((*filtered_strainrate_)(ij)).replace_local_values(1, &val, &lnodeid);
              if (alphaij_)
                err += ((*filtered_alphaij_)(ij)).replace_local_values(1, &val, &lnodeid);
            }
            else
            {
              // set value to mean value
              // we already divide by the corresponding volume of all contributing elements,
              // since we set the volume to 1.0 in the next step in order not to modify the
              // dirichlet values
              if (modeled_subgrid_stress_)
              {
                double val = ((((*filtered_modeled_subgrid_stress_)(ij)))[lnodeid]) / thisvol;
                err += ((*filtered_modeled_subgrid_stress_)(ij))
                           .replace_local_values(1, &val, &lnodeid);
              }
              if (strainrate_)
              {
                double val = ((((*filtered_strainrate_)(ij)))[lnodeid]) / thisvol;
                err += ((*filtered_strainrate_)(ij)).replace_local_values(1, &val, &lnodeid);
              }
              if (alphaij_)
              {
                double val = ((((*filtered_alphaij_)(ij)))[lnodeid]) / thisvol;
                err += ((*filtered_alphaij_)(ij)).replace_local_values(1, &val, &lnodeid);
              }
            }
          }  // end loop jdim
        }  // end loop idim

        double volval = 1.0;
        err += patchvol.replace_local_values(1, &volval, &lnodeid);
        if (err != 0) FOUR_C_THROW("dof/node not on proc");
      }  // is dirichlet node
    }  // end loop all nodes
  }


  // ---------------------------------------------------------------
  // scale vectors by element patch sizes --- this corresponds to
  // the normalization of the box filter function

  // loop all nodes on the processor
  for (int lnodeid = 0; lnodeid < discret_->num_my_row_nodes(); ++lnodeid)
  {
    double thisvol = (patchvol)[lnodeid];

    int err = 0;
    double val = 0.0;
    if (density_)
    {
      val = (*filtered_dens_)[lnodeid] / thisvol;
      err += filtered_dens_->replace_local_values(1, &val, &lnodeid);
    }
    if (densstrainrate_)
    {
      val = (*filtered_dens_strainrate_)[lnodeid] / thisvol;
      err += filtered_dens_strainrate_->replace_local_values(1, &val, &lnodeid);
    }
    if (expression_)
    {
      val = (*filtered_expression_)[lnodeid] / thisvol;
      err += filtered_expression_->replace_local_values(1, &val, &lnodeid);
    }
    if (alpha2_)
    {
      val = (*filtered_alpha2_)[lnodeid] / thisvol;
      err += filtered_alpha2_->replace_local_values(1, &val, &lnodeid);
    }

    for (int idim = 0; idim < 3; ++idim)
    {
      if (velocity_)
      {
        val = ((((*filtered_vel_)(idim)))[lnodeid]) / thisvol;
        err += ((*filtered_vel_)(idim)).replace_local_values(1, &val, &lnodeid);
      }

      if (densvelocity_)
      {
        val = ((((*filtered_dens_vel_)(idim)))[lnodeid]) / thisvol;
        err += ((*filtered_dens_vel_)(idim)).replace_local_values(1, &val, &lnodeid);
      }

      for (int jdim = 0; jdim < 3; ++jdim)
      {
        const int ij = numdim * idim + jdim;

        if (reynoldsstress_)
        {
          val = ((((*filtered_reynoldsstress_)(ij)))[lnodeid]) / thisvol;
          err += ((*filtered_reynoldsstress_)(ij)).replace_local_values(1, &val, &lnodeid);
        }
        if (modeled_subgrid_stress_)
        {
          val = ((((*filtered_modeled_subgrid_stress_)(ij)))[lnodeid]) / thisvol;
          err += ((*filtered_modeled_subgrid_stress_)(ij)).replace_local_values(1, &val, &lnodeid);
        }
        if (strainrate_)
        {
          val = ((((*filtered_strainrate_)(ij)))[lnodeid]) / thisvol;
          err += ((*filtered_strainrate_)(ij)).replace_local_values(1, &val, &lnodeid);
        }
        if (alphaij_)
        {
          val = ((((*filtered_alphaij_)(ij)))[lnodeid]) / thisvol;
          err += ((*filtered_alphaij_)(ij)).replace_local_values(1, &val, &lnodeid);
        }
      }  // end loop jdim
      if (err != 0) FOUR_C_THROW("dof not on proc");
    }  // end loop idim
  }  // end loop nodes

  // clean up
  discret_->clear_state();

  // calculate fine scale velocities
  if (finescale_velocity_)
  {
    // fine scale veocity requires filtered velocity
    if (not velocity_)
      FOUR_C_THROW(
          "filtered velocity is required in the box filter to calculate the fine scale velocity");
    // loop all elements on this proc
    for (int nid = 0; nid < discret_->num_my_row_nodes(); ++nid)
    {
      // get the node
      Core::Nodes::Node* node = discret_->l_row_node(nid);
      // get global ids of all dofs of the node
      std::vector<int> dofs = discret_->dof(node);
      // we only loop over all velocity dofs
      for (int d = 0; d < discret_->num_dof(node) - 1; ++d)
      {
        // get global id of the dof
        int gid = dofs[d];
        // get local dof id corresponding to the global id
        int lid = discret_->dof_row_map()->LID(gid);
        // filtered velocity and all scale velocity
        double filteredvel = (((*filtered_vel_)(d)))[nid];
        double vel = (*velocity)[lid];
        // calculate fine scale velocity
        double val = vel - filteredvel;
        // calculate fine scale velocity
        int err = ((*fs_vel_)(d)).replace_local_values(1, &val, &nid);
        if (err != 0) FOUR_C_THROW("dof not on proc");
      }
    }
  }

  // ----------------------------------------------------------
  // the communication part: Export from row to column map

  // get the column map in order to communicate the result to all ghosted nodes
  const Epetra_Map* nodecolmap = discret_->node_col_map();

  // allocate distributed vectors in col map format to have the filtered
  // quantities available on ghosted nodes
  if (velocity_)
    col_filtered_vel_ = std::make_shared<Core::LinAlg::MultiVector<double>>(*nodecolmap, 3, true);
  if (reynoldsstress_)
    col_filtered_reynoldsstress_ =
        std::make_shared<Core::LinAlg::MultiVector<double>>(*nodecolmap, 9, true);
  if (modeled_subgrid_stress_)
    col_filtered_modeled_subgrid_stress_ =
        std::make_shared<Core::LinAlg::MultiVector<double>>(*nodecolmap, 9, true);
  if (finescale_velocity_)
    col_fs_vel_ = std::make_shared<Core::LinAlg::MultiVector<double>>(*nodecolmap, 3, true);
  if (densvelocity_)
    col_filtered_dens_vel_ =
        std::make_shared<Core::LinAlg::MultiVector<double>>(*nodecolmap, 3, true);
  if (density_)
    col_filtered_dens_ = std::make_shared<Core::LinAlg::Vector<double>>(*nodecolmap, true);
  if (densstrainrate_)
    col_filtered_dens_strainrate_ =
        std::make_shared<Core::LinAlg::Vector<double>>(*nodecolmap, true);
  if (strainrate_)
    col_filtered_strainrate_ =
        std::make_shared<Core::LinAlg::MultiVector<double>>(*nodecolmap, 9, true);
  if (alphaij_)
    col_filtered_alphaij_ =
        std::make_shared<Core::LinAlg::MultiVector<double>>(*nodecolmap, 9, true);
  if (expression_)
    col_filtered_expression_ = std::make_shared<Core::LinAlg::Vector<double>>(*nodecolmap, true);
  if (alpha2_)
    col_filtered_alpha2_ = std::make_shared<Core::LinAlg::Vector<double>>(*nodecolmap, true);

  // export filtered vectors in rowmap to columnmap format
  if (velocity_) Core::LinAlg::export_to(*filtered_vel_, *col_filtered_vel_);
  if (reynoldsstress_)
    Core::LinAlg::export_to(*filtered_reynoldsstress_, *col_filtered_reynoldsstress_);
  if (modeled_subgrid_stress_)
    Core::LinAlg::export_to(
        *filtered_modeled_subgrid_stress_, *col_filtered_modeled_subgrid_stress_);
  if (finescale_velocity_) Core::LinAlg::export_to(*fs_vel_, *col_fs_vel_);
  if (densvelocity_) Core::LinAlg::export_to(*filtered_dens_vel_, *col_filtered_dens_vel_);
  if (density_) Core::LinAlg::export_to(*filtered_dens_, *col_filtered_dens_);
  if (densstrainrate_)
    Core::LinAlg::export_to(*filtered_dens_strainrate_, *col_filtered_dens_strainrate_);
  if (strainrate_) Core::LinAlg::export_to(*filtered_strainrate_, *col_filtered_strainrate_);
  if (alphaij_) Core::LinAlg::export_to(*filtered_alphaij_, *col_filtered_alphaij_);
  if (expression_) Core::LinAlg::export_to(*filtered_expression_, *col_filtered_expression_);
  if (alpha2_) Core::LinAlg::export_to(*filtered_alpha2_, *col_filtered_alpha2_);
  return;
}



/*----------------------------------------------------------------------*
 | perform box filtering                                      (private) |
 |                                                      rasthofer 08/12 |
 *----------------------------------------------------------------------*/
void FLD::Boxfilter::apply_box_filter_scatra(
    const std::shared_ptr<const Core::LinAlg::Vector<double>> scalar, const double thermpress,
    const Core::LinAlg::Vector<double>& dirichtoggle, const int ndsvel)
{
  TEUCHOS_FUNC_TIME_MONITOR("apply_filter_for_dynamic_computation_of_prt");
  if (apply_box_filter_ == true) FOUR_C_THROW("not yet considered");
  // LES turbulence modeling is only valid for 3 dimensions
  const int numdim = 3;

  // generate a parameterlist for communication and control
  Teuchos::ParameterList filterparams;
  // action for elements
  Core::Utils::add_enum_class_to_parameter_list<ScaTra::Action>(
      "action", ScaTra::Action::calc_scatra_box_filter, filterparams);

  filterparams.set("thermpress", thermpress);

  // set state vector to pass distributed vector to the element
  scatradiscret_->clear_state();
  scatradiscret_->set_state("scalar", scalar);

  // dummies
  Core::LinAlg::SerialDenseMatrix emat1;
  Core::LinAlg::SerialDenseMatrix emat2;
  Core::LinAlg::SerialDenseVector evec1;
  Core::LinAlg::SerialDenseVector evec2;
  Core::LinAlg::SerialDenseVector evec3;

  // ---------------------------------------------------------------
  // get a vector layout from the discretization to construct
  const Epetra_Map* noderowmap = scatradiscret_->node_row_map();

  // alloc an additional vector to store/add up the patch volume
  Core::LinAlg::Vector<double> patchvol(*noderowmap, true);

  // free mem and reallocate to zero out vecs
  filtered_dens_vel_temp_ = nullptr;
  filtered_dens_rateofstrain_temp_ = nullptr;
  filtered_vel_ = nullptr;
  filtered_dens_vel_ = nullptr;
  filtered_temp_ = nullptr;
  filtered_dens_temp_ = nullptr;
  filtered_dens_ = nullptr;
  if (phi_) filtered_phi_ = nullptr;
  if (phi2_) filtered_phi2_ = nullptr;
  if (phiexpression_) filtered_phiexpression_ = nullptr;
  if (alphaijsc_) filtered_alphaijsc_ = nullptr;

  filtered_dens_vel_temp_ =
      std::make_shared<Core::LinAlg::MultiVector<double>>(*noderowmap, numdim, true);
  filtered_dens_rateofstrain_temp_ =
      std::make_shared<Core::LinAlg::MultiVector<double>>(*noderowmap, numdim, true);
  filtered_vel_ = std::make_shared<Core::LinAlg::MultiVector<double>>(*noderowmap, numdim, true);
  filtered_dens_vel_ =
      std::make_shared<Core::LinAlg::MultiVector<double>>(*noderowmap, numdim, true);
  filtered_temp_ = std::make_shared<Core::LinAlg::Vector<double>>(*noderowmap, true);
  filtered_dens_temp_ = std::make_shared<Core::LinAlg::Vector<double>>(*noderowmap, true);
  filtered_dens_ = std::make_shared<Core::LinAlg::Vector<double>>(*noderowmap, true);
  if (phi_)
    filtered_phi_ = std::make_shared<Core::LinAlg::MultiVector<double>>(*noderowmap, numdim, true);
  if (phi2_) filtered_phi2_ = std::make_shared<Core::LinAlg::Vector<double>>(*noderowmap, true);
  if (phiexpression_)
    filtered_phiexpression_ = std::make_shared<Core::LinAlg::Vector<double>>(*noderowmap, true);
  if (alphaijsc_)
    filtered_alphaijsc_ =
        std::make_shared<Core::LinAlg::MultiVector<double>>(*noderowmap, numdim * numdim, true);

  // ---------------------------------------------------------------
  // do the integration of the (not normalized) box filter function
  // on the element

  // loop all elements on this proc (including ghosted ones)
  for (int nele = 0; nele < scatradiscret_->num_my_col_elements(); ++nele)
  {
    // get the element
    Core::Elements::Element* ele = scatradiscret_->l_col_element(nele);

    // provide vectors for filtered quantities //declaration necessary even if not used
    std::shared_ptr<std::vector<double>> vel_hat =
        std::make_shared<std::vector<double>>((numdim), 0.0);
    std::shared_ptr<std::vector<double>> densvel_hat =
        std::make_shared<std::vector<double>>((numdim), 0.0);
    std::shared_ptr<std::vector<double>> densveltemp_hat =
        std::make_shared<std::vector<double>>((numdim), 0.0);
    std::shared_ptr<std::vector<double>> densstraintemp_hat =
        std::make_shared<std::vector<double>>((numdim), 0.0);
    std::shared_ptr<std::vector<double>> phi_hat =
        std::make_shared<std::vector<double>>((numdim), 0.0);
    std::shared_ptr<std::vector<std::vector<double>>> alphaijsc_hat =
        std::make_shared<std::vector<std::vector<double>>>();
    if (alphaijsc_)
    {
      (*alphaijsc_hat).resize(numdim);
      for (int rr = 0; rr < numdim; rr++) ((*alphaijsc_hat)[rr]).resize(numdim);
      // initialize with zeros
      for (int rr = 0; rr < numdim; rr++)
      {
        for (int ss = 0; ss < numdim; ss++) (*alphaijsc_hat)[rr][ss] = 0.0;
      }
    }
    // and set them in parameter list
    filterparams.set<std::shared_ptr<std::vector<double>>>("vel_hat", vel_hat);
    filterparams.set<std::shared_ptr<std::vector<double>>>("densvel_hat", densvel_hat);
    filterparams.set<std::shared_ptr<std::vector<double>>>("densveltemp_hat", densveltemp_hat);
    filterparams.set<std::shared_ptr<std::vector<double>>>(
        "densstraintemp_hat", densstraintemp_hat);
    filterparams.set<std::shared_ptr<std::vector<double>>>("phi_hat", phi_hat);
    filterparams.set<std::shared_ptr<std::vector<std::vector<double>>>>(
        "alphaijsc_hat", alphaijsc_hat);

    // initialize variables for filtered scalar quantities
    double dens_hat = 0.0;
    double temp_hat = 0.0;
    double dens_temp_hat = 0.0;
    double phi2_hat = 0.0;
    double phiexpression_hat = 0.0;

    // initialize volume contribution
    double volume_contribution = 0.0;

    // get element location vector, dirichlet flags and ownerships
    Core::Elements::LocationArray la(scatradiscret_->num_dof_sets());
    ele->location_vector(*scatradiscret_, la, false);

    // call the element evaluate method to integrate functions
    // against heaviside function element
    int err = ele->evaluate(filterparams, *scatradiscret_, la, emat1, emat2, evec1, evec2, evec2);
    if (err)
      FOUR_C_THROW("Proc {}: Element {} returned err={}",
          Core::Communication::my_mpi_rank(scatradiscret_->get_comm()), ele->id(), err);

    // get contribution to patch volume of this element. Add it up.
    // double volume_contribution = filterparams.get<double>("volume_contribution");
    volume_contribution = filterparams.get<double>("volume_contribution");

    // filtered scalar quantities
    dens_hat = filterparams.get<double>("dens_hat");
    temp_hat = filterparams.get<double>("temp_hat");
    dens_temp_hat = filterparams.get<double>("dens_temp_hat");
    if (phi2_) phi2_hat = filterparams.get<double>("phi2_hat");
    if (phiexpression_) phiexpression_hat = filterparams.get<double>("phiexpression_hat");
    // loop all nodes of this element, add values to the global vectors
    Core::Nodes::Node** elenodes = ele->nodes();
    for (int nn = 0; nn < ele->num_node(); ++nn)
    {
      Core::Nodes::Node* node = (elenodes[nn]);

      // we are interested only in  row nodes
      if (node->owner() == Core::Communication::my_mpi_rank(scatradiscret_->get_comm()))
      {
        // now assemble the computed values into the global vector
        int id = (node->id());

        patchvol.sum_into_global_values(1, &volume_contribution, &id);
        filtered_dens_->sum_into_global_values(1, &dens_hat, &id);
        filtered_temp_->sum_into_global_values(1, &temp_hat, &id);
        filtered_dens_temp_->sum_into_global_values(1, &dens_temp_hat, &id);
        if (phi2_) filtered_phi2_->sum_into_global_values(1, &phi2_hat, &id);
        if (phiexpression_)
          filtered_phiexpression_->sum_into_global_values(1, &phiexpression_hat, &id);
        for (int idim = 0; idim < numdim; ++idim)
        {
          double val = (*vel_hat)[idim];
          ((*filtered_vel_)(idim)).sum_into_global_values(1, &val, &id);
          val = (*densveltemp_hat)[idim];
          ((*filtered_dens_vel_temp_)(idim)).sum_into_global_values(1, &val, &id);
          val = (*densstraintemp_hat)[idim];
          ((*filtered_dens_rateofstrain_temp_)(idim)).sum_into_global_values(1, &val, &id);
          val = (*densvel_hat)[idim];
          ((*filtered_dens_vel_)(idim)).sum_into_global_values(1, &val, &id);
          if (phi_)
          {
            val = (*phi_hat)[idim];
            ((*filtered_phi_)(idim)).sum_into_global_values(1, &val, &id);
          }
          if (alphaijsc_)
          {
            for (int jdim = 0; jdim < numdim; ++jdim)
            {
              const int ij = numdim * idim + jdim;
              double val = (*alphaijsc_hat)[idim][jdim];
              ((*filtered_alphaijsc_)(ij)).sum_into_global_values(1, &val, &id);
            }
          }
        }
      }
    }
  }  // end elementloop

  // ---------------------------------------------------------------
  // send add values from masters and slaves
  {
    double val = 0.0;
    std::vector<double> vel_val(3);
    std::vector<double> dens_vel_val(3);
    std::vector<double> dens_vel_temp_val(3);
    std::vector<double> dens_strain_temp_val(3);
    std::vector<double> phi_val(3);
    std::vector<std::vector<double>> alphaijsc_val;
    if (alphaijsc_)
    {
      alphaijsc_val.resize(3);
      for (int rr = 0; rr < 3; rr++) (alphaijsc_val[rr]).resize(3);
    }
    double temp_val = 0.0;
    double dens_val = 0.0;
    double dens_temp_val = 0.0;
    double phi2_val = 0.0;
    double phiexpression_val = 0.0;

    // loop all master nodes on this proc
    std::map<int, std::vector<int>>* pbcmapmastertoslave =
        scatradiscret_->get_all_pbc_coupled_col_nodes();
    // loop all master nodes on this proc
    if (pbcmapmastertoslave)
    {
      for (const auto& [master_gid, slave_gids] : *pbcmapmastertoslave)
      {
        // loop only owned nodes
        if ((scatradiscret_->g_node(master_gid))->owner() !=
            Core::Communication::my_mpi_rank(scatradiscret_->get_comm()))
          continue;

        int lid = noderowmap->LID(master_gid);
        if (lid < 0) FOUR_C_THROW("nodelid < 0 ?");

        val = (patchvol)[lid];

        dens_val = (*filtered_dens_)[lid];
        dens_temp_val = (*filtered_dens_temp_)[lid];
        temp_val = (*filtered_temp_)[lid];
        if (phi2_) phi2_val = (*filtered_phi2_)[lid];
        if (phiexpression_) phiexpression_val = (*filtered_phiexpression_)[lid];

        for (int idim = 0; idim < numdim; ++idim)
        {
          vel_val[idim] = ((((*filtered_vel_)(idim)))[lid]);
          dens_vel_val[idim] = ((((*filtered_dens_vel_)(idim)))[lid]);
          dens_vel_temp_val[idim] = ((((*filtered_dens_vel_temp_)(idim)))[lid]);
          dens_strain_temp_val[idim] = ((((*filtered_dens_rateofstrain_temp_)(idim)))[lid]);
          if (phi_) phi_val[idim] = ((((*filtered_phi_)(idim)))[lid]);
          if (alphaijsc_)
          {
            for (int jdim = 0; jdim < numdim; ++jdim)
            {
              const int ij = numdim * idim + jdim;
              alphaijsc_val[idim][jdim] = (((*filtered_alphaijsc_)(ij)))[lid];
            }
          }
        }

        // loop all this masters slaves
        for (auto slave_gid : slave_gids)
        {
          lid = noderowmap->LID(slave_gid);
          val += (patchvol)[lid];

          dens_val += (*filtered_dens_)[lid];
          dens_temp_val += (*filtered_dens_temp_)[lid];
          temp_val += (*filtered_temp_)[lid];
          if (phi2_) phi2_val += (*filtered_phi2_)[lid];
          if (phiexpression_) phiexpression_val += (*filtered_phiexpression_)[lid];

          for (int idim = 0; idim < numdim; ++idim)
          {
            vel_val[idim] += ((((*filtered_vel_)(idim)))[lid]);
            dens_vel_val[idim] += ((((*filtered_dens_vel_)(idim)))[lid]);
            dens_vel_temp_val[idim] += ((((*filtered_dens_vel_temp_)(idim)))[lid]);
            dens_strain_temp_val[idim] += ((((*filtered_dens_rateofstrain_temp_)(idim)))[lid]);
            if (phi_) phi_val[idim] += ((((*filtered_phi_)(idim)))[lid]);
            if (alphaijsc_)
            {
              for (int jdim = 0; jdim < numdim; ++jdim)
              {
                const int ij = numdim * idim + jdim;
                alphaijsc_val[idim][jdim] += (((*filtered_alphaijsc_)(ij)))[lid];
              }
            }
          }
        }  // end loop slaves

        // replace value by sum
        lid = noderowmap->LID(master_gid);
        int error = patchvol.replace_local_values(1, &val, &lid);
        if (error != 0) FOUR_C_THROW("dof not on proc");

        int e = 0;
        e += filtered_dens_->replace_local_values(1, &dens_val, &lid);
        e += filtered_dens_temp_->replace_local_values(1, &dens_temp_val, &lid);
        e += filtered_temp_->replace_local_values(1, &temp_val, &lid);
        if (phi2_) e += filtered_phi2_->replace_local_values(1, &phi2_val, &lid);
        if (phiexpression_)
          e += filtered_phiexpression_->replace_local_values(1, &phiexpression_val, &lid);
        if (e != 0) FOUR_C_THROW("dof not on proc");

        for (int idim = 0; idim < numdim; ++idim)
        {
          int err = 0;
          err += ((*filtered_vel_)(idim)).replace_local_values(1, &(vel_val[idim]), &lid);
          err += ((*filtered_dens_vel_)(idim)).replace_local_values(1, &(dens_vel_val[idim]), &lid);
          err += ((*filtered_dens_vel_temp_)(idim))
                     .replace_local_values(1, &(dens_vel_temp_val[idim]), &lid);
          err += ((*filtered_dens_rateofstrain_temp_)(idim))
                     .replace_local_values(1, &(dens_strain_temp_val[idim]), &lid);
          if (phi_) err += ((*filtered_phi_)(idim)).replace_local_values(1, &(phi_val[idim]), &lid);
          if (alphaijsc_)
          {
            for (int jdim = 0; jdim < numdim; ++jdim)
            {
              const int ij = numdim * idim + jdim;
              err += ((*filtered_alphaijsc_)(ij))
                         .replace_local_values(1, &(alphaijsc_val[idim][jdim]), &lid);
            }  // end loop jdim
          }
          if (err != 0) FOUR_C_THROW("dof not on proc");
        }

        // loop all this masters slaves
        for (auto slave_gid : slave_gids)
        {
          int err = 0;
          lid = noderowmap->LID(slave_gid);
          err += patchvol.replace_local_values(1, &val, &lid);

          err += filtered_dens_->replace_local_values(1, &dens_val, &lid);
          err += filtered_dens_temp_->replace_local_values(1, &dens_temp_val, &lid);
          err += filtered_temp_->replace_local_values(1, &temp_val, &lid);
          if (phi2_) err += filtered_phi2_->replace_local_values(1, &phi2_val, &lid);
          if (phiexpression_)
            err += filtered_phiexpression_->replace_local_values(1, &phiexpression_val, &lid);

          for (int idim = 0; idim < numdim; ++idim)
          {
            err += ((*filtered_vel_)(idim)).replace_local_values(1, &(vel_val[idim]), &lid);
            err +=
                ((*filtered_dens_vel_)(idim)).replace_local_values(1, &(dens_vel_val[idim]), &lid);
            err += ((*filtered_dens_vel_temp_)(idim))
                       .replace_local_values(1, &(dens_vel_temp_val[idim]), &lid);
            err += ((*filtered_dens_rateofstrain_temp_)(idim))
                       .replace_local_values(1, &(dens_strain_temp_val[idim]), &lid);
            if (phi_)
              err += ((*filtered_phi_)(idim)).replace_local_values(1, &(phi_val[idim]), &lid);
            if (alphaijsc_)
            {
              for (int jdim = 0; jdim < numdim; ++jdim)
              {
                const int ij = numdim * idim + jdim;
                err += ((*filtered_alphaijsc_)(ij))
                           .replace_local_values(1, &(alphaijsc_val[idim][jdim]), &lid);
              }  // end loop jdim
            }
          }

          if (err != 0) FOUR_C_THROW("dof not on proc");
        }  // end loop slaves
      }  // end loop masters
    }
  }

  // ---------------------------------------------------------------
  // extract convective velocity from scatra discretization
  std::shared_ptr<const Core::LinAlg::Vector<double>> convel =
      scatradiscret_->get_state(ndsvel, "convective velocity field");
  if (convel == nullptr) FOUR_C_THROW("Cannot extract convective velocity field");

  // replace values at dirichlet nodes
  {
    // get a rowmap for the dofs
    const Epetra_Map* dofrowmap = scatradiscret_->dof_row_map();

    // as we want to identify nodes at walls,
    // we have to be sure that fluid and scatra are still matching
    if (not scatradiscret_->node_row_map()->SameAs(*(discret_->node_row_map())))
      FOUR_C_THROW("Fluid and ScaTra noderowmaps are NOT identical.");

    // loop all nodes on the processor
    for (int lnodeid = 0; lnodeid < scatradiscret_->num_my_row_nodes(); ++lnodeid)
    {
      // get the processor local node
      Core::Nodes::Node* lnode = scatradiscret_->l_row_node(lnodeid);
      std::vector<int> nodedofs = scatradiscret_->dof(ndsvel, lnode);
      // get the corresponding processor local fluid node
      Core::Nodes::Node* fluidlnode = discret_->l_row_node(lnodeid);

      // do we have a dirichlet boundary conditions in the fluid
      std::vector<Core::Conditions::Condition*> dbccond;
      fluidlnode->get_condition("Dirichlet", dbccond);

      // yes, we have a dirichlet boundary condition
      if (dbccond.size() > 0)
      {
#ifdef FOUR_C_ENABLE_ASSERTIONS
        if ((lnode->x()[0] != fluidlnode->x()[0]) or (lnode->x()[1] != fluidlnode->x()[1]) or
            (lnode->x()[2] != fluidlnode->x()[2]))
          FOUR_C_THROW("Nodes do not match.");
#endif
        // we only want to modify nodes at the wall, as the model should vanish there
        // check, whether we have a no-slip node
        int no_slip_node = 0;
        for (int idim = 0; idim < numdim; idim++)
        {
          // get global and local dof IDs
          const int gid = nodedofs[idim];
          const int lid = convel->get_map().LID(gid);
          if (lid < 0) FOUR_C_THROW("Local ID not found in map for given global ID!");

          double vel_i = (*convel)[lid];
          if (abs(vel_i) < 1e-12) no_slip_node++;
        }

        // yes, we have a no-slip node
        if (no_slip_node == numdim)
        {
          // do we also have a temperature dirichlet boundary condition
          // get the set of temperature degrees of freedom associated with the node
          std::vector<int> nodedofset = scatradiscret_->dof(0, lnode);
          if (nodedofset.size() > 1)
            FOUR_C_THROW(
                "Dynamic Smagorinsky or dynamic Vreman currently only implemented for one scalar "
                "field!");

          // check whether the dofs are Dirichlet constrained
          bool is_dirichlet_node = false;
          int gid = nodedofset[0];
          int lid = dofrowmap->LID(gid);

          // this is a dirichlet node
          if ((dirichtoggle)[lid] == 1) is_dirichlet_node = true;

          // get volume
          double thisvol = (patchvol)[lnodeid];
          // and density
          double dens = (*filtered_dens_)[lnodeid] / thisvol;
          int err = 0;
          err += filtered_dens_->replace_local_values(1, &dens, &lnodeid);

          double temp = 0.0;
          if (is_dirichlet_node)
          {
            temp = (*scalar)[lid];
            err += filtered_temp_->replace_local_values(1, &temp, &lnodeid);
            double val = dens * temp;
            err += filtered_dens_temp_->replace_local_values(1, &val, &lnodeid);
            if (phi2_)
            {
              double val = 0.0;
              err += filtered_phi2_->replace_local_values(1, &val, &lnodeid);
            }
            if (phiexpression_)
            {
              double val = 0.0;
              err += filtered_phiexpression_->replace_local_values(1, &val, &lnodeid);
            }
          }
          else
          {
            temp = (*filtered_temp_)[lnodeid] / thisvol;
            err += filtered_temp_->replace_local_values(1, &temp, &lnodeid);
            double val = (*filtered_dens_temp_)[lnodeid] / thisvol;
            err += filtered_dens_temp_->replace_local_values(1, &val, &lnodeid);
          }

          for (int idim = 0; idim < numdim; idim++)
          {
            // get global and local dof IDs
            const int gid = nodedofs[idim];
            const int lid = convel->get_map().LID(gid);
            if (lid < 0) FOUR_C_THROW("Local ID not found in map for given global ID!");

            double valvel_i = (*convel)[lid];
            err += ((*filtered_vel_)(idim)).replace_local_values(1, &valvel_i, &lnodeid);

            double valdensvel_i = dens * valvel_i;
            err += ((*filtered_dens_vel_)(idim)).replace_local_values(1, &valdensvel_i, &lnodeid);

            double dvtval = dens * temp * valvel_i;
            err += ((*filtered_dens_vel_temp_)(idim)).replace_local_values(1, &dvtval, &lnodeid);

            // Peter style
            double drtval = 0.0;
            err += ((*filtered_dens_rateofstrain_temp_)(idim))
                       .replace_local_values(1, &drtval, &lnodeid);

            if (phi_)
            {
              double drtval = 0.0;
              err += ((*filtered_phi_)(idim)).replace_local_values(1, &drtval, &lnodeid);
            }
            if (alphaijsc_)
            {
              for (int jdim = 0; jdim < numdim; ++jdim)
              {
                const int ij = numdim * idim + jdim;
                if (no_slip_node == numdim)
                {
                  // set value to zero (original Peter style)
                  double val = 0.0;
                  err += ((*filtered_alphaijsc_)(ij)).replace_local_values(1, &val, &lnodeid);
                }
                else
                {
                  // set value to mean value
                  // we already divide by the corresponding volume of all contributing elements,
                  // since we set the volume to 1.0 in the next step in order not to modify the
                  // dirichlet values
                  double val = ((((*filtered_alphaijsc_)(ij)))[lnodeid]) / thisvol;
                  err += ((*filtered_alphaijsc_)(ij)).replace_local_values(1, &val, &lnodeid);
                }
              }  // end loop jdim
            }

            // alternative: see comment in apply_box_filter() for velocity field
            // double drtval = ((*((*filtered_dens_rateofstrain_temp_)(idim)))[lnodeid])/thisvol;
            // err +=
            // ((*filtered_dens_rateofstrain_temp_)(idim)).ReplaceMyValues(1,&drtval,&lnodeid);
          }

          double volval = 1.0;
          err += patchvol.replace_local_values(1, &volval, &lnodeid);
          if (err != 0) FOUR_C_THROW("dof/node not on proc");
        }
      }
    }
  }

  // ---------------------------------------------------------------
  // scale vectors by element patch sizes --- this corresponds to
  // the normalization of the box filter function

  // loop all nodes on the processor
  for (int lnodeid = 0; lnodeid < scatradiscret_->num_my_row_nodes(); ++lnodeid)
  {
    double thisvol = (patchvol)[lnodeid];

    int err = 0;
    double val = 0.0;

    val = (*filtered_temp_)[lnodeid] / thisvol;
    err += filtered_temp_->replace_local_values(1, &val, &lnodeid);
    val = (*filtered_dens_)[lnodeid] / thisvol;
    err += filtered_dens_->replace_local_values(1, &val, &lnodeid);
    val = (*filtered_dens_temp_)[lnodeid] / thisvol;
    err += filtered_dens_temp_->replace_local_values(1, &val, &lnodeid);
    if (phi2_)
    {
      val = (*filtered_phi2_)[lnodeid] / thisvol;
      err += filtered_phi2_->replace_local_values(1, &val, &lnodeid);
    }
    if (phiexpression_)
    {
      val = (*filtered_phiexpression_)[lnodeid] / thisvol;
      err += filtered_phiexpression_->replace_local_values(1, &val, &lnodeid);
    }
    for (int idim = 0; idim < 3; ++idim)
    {
      val = ((((*filtered_vel_)(idim)))[lnodeid]) / thisvol;
      err += ((*filtered_vel_)(idim)).replace_local_values(1, &val, &lnodeid);
      val = ((((*filtered_dens_vel_)(idim)))[lnodeid]) / thisvol;
      err += ((*filtered_dens_vel_)(idim)).replace_local_values(1, &val, &lnodeid);
      val = ((((*filtered_dens_vel_temp_)(idim)))[lnodeid]) / thisvol;
      err += ((*filtered_dens_vel_temp_)(idim)).replace_local_values(1, &val, &lnodeid);
      val = ((((*filtered_dens_rateofstrain_temp_)(idim)))[lnodeid]) / thisvol;
      err += ((*filtered_dens_rateofstrain_temp_)(idim)).replace_local_values(1, &val, &lnodeid);
      if (phi_)
      {
        val = ((((*filtered_phi_)(idim)))[lnodeid]) / thisvol;
        err += ((*filtered_phi_)(idim)).replace_local_values(1, &val, &lnodeid);
      }
      if (alphaijsc_)
      {
        for (int jdim = 0; jdim < 3; ++jdim)
        {
          const int ij = numdim * idim + jdim;
          val = ((((*filtered_alphaijsc_)(ij)))[lnodeid]) / thisvol;
          err += ((*filtered_alphaijsc_)(ij)).replace_local_values(1, &val, &lnodeid);
        }  // end loop jdim
      }
    }  // end loop idim
    if (err != 0) FOUR_C_THROW("dof not on proc");
  }  // end loop nodes

  // clean up
  scatradiscret_->clear_state();

  // ----------------------------------------------------------
  // the communication part: Export from row to column map

  // get the column map in order to communicate the result to all ghosted nodes
  const Epetra_Map* nodecolmap = scatradiscret_->node_col_map();

  // allocate distributed vectors in col map format to have the filtered
  // quantities available on ghosted nodes
  col_filtered_vel_ = std::make_shared<Core::LinAlg::MultiVector<double>>(*nodecolmap, 3, true);
  col_filtered_dens_vel_ =
      std::make_shared<Core::LinAlg::MultiVector<double>>(*nodecolmap, 3, true);
  col_filtered_dens_vel_temp_ =
      std::make_shared<Core::LinAlg::MultiVector<double>>(*nodecolmap, 3, true);
  col_filtered_dens_rateofstrain_temp_ =
      std::make_shared<Core::LinAlg::MultiVector<double>>(*nodecolmap, 3, true);
  col_filtered_temp_ = std::make_shared<Core::LinAlg::Vector<double>>(*nodecolmap, true);
  col_filtered_dens_ = std::make_shared<Core::LinAlg::Vector<double>>(*nodecolmap, true);
  col_filtered_dens_temp_ = std::make_shared<Core::LinAlg::Vector<double>>(*nodecolmap, true);
  if (phi_)
    col_filtered_phi_ = std::make_shared<Core::LinAlg::MultiVector<double>>(*nodecolmap, 3, true);
  if (phi2_) col_filtered_phi2_ = std::make_shared<Core::LinAlg::Vector<double>>(*nodecolmap, true);
  if (phiexpression_)
    col_filtered_phiexpression_ = std::make_shared<Core::LinAlg::Vector<double>>(*nodecolmap, true);
  if (alphaijsc_)
    col_filtered_alphaijsc_ =
        std::make_shared<Core::LinAlg::MultiVector<double>>(*nodecolmap, 9, true);

  // export filtered vectors in rowmap to columnmap format
  Core::LinAlg::export_to(*filtered_vel_, *col_filtered_vel_);
  Core::LinAlg::export_to(*filtered_dens_vel_, *col_filtered_dens_vel_);
  Core::LinAlg::export_to(*filtered_dens_vel_temp_, *col_filtered_dens_vel_temp_);
  Core::LinAlg::export_to(*filtered_dens_rateofstrain_temp_, *col_filtered_dens_rateofstrain_temp_);
  Core::LinAlg::export_to(*filtered_temp_, *col_filtered_temp_);
  Core::LinAlg::export_to(*filtered_dens_, *col_filtered_dens_);
  Core::LinAlg::export_to(*filtered_dens_temp_, *col_filtered_dens_temp_);
  if (phi_) Core::LinAlg::export_to(*filtered_phi_, *col_filtered_phi_);
  if (phi2_) Core::LinAlg::export_to(*filtered_phi2_, *col_filtered_phi2_);
  if (phiexpression_)
    Core::LinAlg::export_to(*filtered_phiexpression_, *col_filtered_phiexpression_);
  if (alphaijsc_) Core::LinAlg::export_to(*filtered_alphaijsc_, *col_filtered_alphaijsc_);

  return;
}

FOUR_C_NAMESPACE_CLOSE
