// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_adapter_str_fsiwrapper.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fsi_str_model_evaluator_partitioned.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_structure_aux.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

namespace
{
  bool prestress_is_active(const double currentTime)
  {
    Inpar::Solid::PreStress pstype = Teuchos::getIntegralValue<Inpar::Solid::PreStress>(
        Global::Problem::instance()->structural_dynamic_params(), "PRESTRESS");
    const double pstime =
        Global::Problem::instance()->structural_dynamic_params().get<double>("PRESTRESSTIME");
    return pstype != Inpar::Solid::PreStress::none && currentTime <= pstime + 1.0e-15;
  }
}  // namespace

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Adapter::FSIStructureWrapper::FSIStructureWrapper(std::shared_ptr<Structure> structure)
    : StructureWrapper(structure)
{
  // set-up FSI interface
  interface_ = std::make_shared<Solid::MapExtractor>();

  if (Global::Problem::instance()->get_problem_type() != Core::ProblemType::fpsi)
    interface_->setup(*discretization(), *discretization()->dof_row_map());
  else
    interface_->setup(*discretization(), *discretization()->dof_row_map(),
        true);  // create overlapping maps for fpsi problem

  const Teuchos::ParameterList& fsidyn = Global::Problem::instance()->fsi_dynamic_params();
  const Teuchos::ParameterList& fsipart = fsidyn.sublist("PARTITIONED SOLVER");
  predictor_ = fsipart.get<std::string>("PREDICTOR");
}

/*------------------------------------------------------------------------------------*
 * Rebuild FSI interface on structure side                              sudhakar 09/13
 * This is necessary if elements are added/deleted from interface
 *------------------------------------------------------------------------------------*/
void Adapter::FSIStructureWrapper::rebuild_interface()
{
  interface_ = std::make_shared<Solid::MapExtractor>();
  interface_->setup(*discretization(), *discretization()->dof_row_map());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FSIStructureWrapper::use_block_matrix()
{
  StructureWrapper::use_block_matrix(interface_, interface_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::Vector<double>> Adapter::FSIStructureWrapper::relaxation_solve(
    std::shared_ptr<Core::LinAlg::Vector<double>> iforce)
{
  apply_interface_forces(iforce);
  fsi_model_evaluator()->set_is_relaxation_solve(true);
  std::shared_ptr<const Core::LinAlg::Vector<double>> idisi =
      fsi_model_evaluator()->solve_relaxation_linear(structure_);
  fsi_model_evaluator()->set_is_relaxation_solve(false);

  // we are just interested in the incremental interface displacements
  return interface_->extract_fsi_cond_vector(*idisi);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::Vector<double>>
Adapter::FSIStructureWrapper::predict_interface_dispnp()
{
  // prestressing business
  std::shared_ptr<Core::LinAlg::Vector<double>> idis;

  if (predictor_ == "d(n)")
  {
    // respect Dirichlet conditions at the interface (required for pseudo-rigid body)
    if (prestress_is_active(time()))
    {
      idis = std::make_shared<Core::LinAlg::Vector<double>>(*interface_->fsi_cond_map(), true);
    }
    else
    {
      idis = interface_->extract_fsi_cond_vector(*dispn());
    }
  }
  else if (predictor_ == "d(n)+dt*(1.5*v(n)-0.5*v(n-1))")
  {
    FOUR_C_THROW("interface velocity v(n-1) not available");
  }
  else if (predictor_ == "d(n)+dt*v(n)")
  {
    if (prestress_is_active(time()))
      FOUR_C_THROW("only constant interface predictor useful for prestressing");

    double current_dt = dt();

    idis = interface_->extract_fsi_cond_vector(*dispn());
    std::shared_ptr<Core::LinAlg::Vector<double>> ivel =
        interface_->extract_fsi_cond_vector(*veln());

    idis->update(current_dt, *ivel, 1.0);
  }
  else if (predictor_ == "d(n)+dt*v(n)+0.5*dt^2*a(n)")
  {
    if (prestress_is_active(time()))
      FOUR_C_THROW("only constant interface predictor useful for prestressing");

    double current_dt = dt();

    idis = interface_->extract_fsi_cond_vector(*dispn());
    std::shared_ptr<Core::LinAlg::Vector<double>> ivel =
        interface_->extract_fsi_cond_vector(*veln());
    std::shared_ptr<Core::LinAlg::Vector<double>> iacc =
        interface_->extract_fsi_cond_vector(*accn());

    idis->update(current_dt, *ivel, 0.5 * current_dt * current_dt, *iacc, 1.0);
  }
  else
  {
    FOUR_C_THROW("unknown interface displacement predictor '{}'", Global::Problem::instance()
                                                                      ->fsi_dynamic_params()
                                                                      .sublist("PARTITIONED SOLVER")
                                                                      .get<std::string>("PREDICTOR")
                                                                      .c_str());
  }

  return idis;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::Vector<double>>
Adapter::FSIStructureWrapper::extract_interface_dispn()
{
  FOUR_C_ASSERT(interface_->full_map()->PointSameAs(dispn()->get_map()),
      "Full map of map extractor and Dispn() do not match.");

  // prestressing business
  if (prestress_is_active(time()))
  {
    return std::make_shared<Core::LinAlg::Vector<double>>(*interface_->fsi_cond_map(), true);
  }
  else
  {
    return interface_->extract_fsi_cond_vector(*dispn());
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::Vector<double>>
Adapter::FSIStructureWrapper::extract_interface_dispnp()
{
  FOUR_C_ASSERT(interface_->full_map()->PointSameAs(dispnp()->get_map()),
      "Full map of map extractor and Dispnp() do not match.");

  // prestressing business
  if (prestress_is_active(time()))
  {
    if (Core::Communication::my_mpi_rank(discretization()->get_comm()) == 0)
      std::cout << "Applying no displacements to the fluid since we do prestressing" << std::endl;

    return std::make_shared<Core::LinAlg::Vector<double>>(*interface_->fsi_cond_map(), true);
  }
  else
  {
    return interface_->extract_fsi_cond_vector(*dispnp());
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// Apply interface forces
void Adapter::FSIStructureWrapper::apply_interface_forces(
    std::shared_ptr<Core::LinAlg::Vector<double>> iforce)
{
  fsi_model_evaluator()->get_interface_force_np_ptr()->put_scalar(0.0);
  interface_->add_fsi_cond_vector(*iforce, *fsi_model_evaluator()->get_interface_force_np_ptr());
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// Apply interface forces deprecated version ! Remove as soon as possible!
void Adapter::FSIStructureWrapper::apply_interface_forces_temporary_deprecated(
    std::shared_ptr<Core::LinAlg::Vector<double>> iforce)
{
  std::shared_ptr<Core::LinAlg::MultiVector<double>> fifc =
      Core::LinAlg::create_multi_vector(*dof_row_map(), 1, true);

  interface_->add_fsi_cond_vector(*iforce, (*fifc)(0));

  set_force_interface(fifc);

  prepare_partition_step();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Solid::ModelEvaluator::PartitionedFSI>
Adapter::FSIStructureWrapper::fsi_model_evaluator()
{
  return fsi_model_evaluator_;
};

FOUR_C_NAMESPACE_CLOSE
