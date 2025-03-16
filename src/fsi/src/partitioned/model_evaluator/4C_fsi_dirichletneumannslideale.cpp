// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fsi_dirichletneumannslideale.hpp"

#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_coupling_adapter_mortar.hpp"
#include "4C_fem_geometry_searchtree.hpp"
#include "4C_fsi_debugwriter.hpp"
#include "4C_fsi_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_fsi.hpp"
#include "4C_mortar_interface.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::DirichletNeumannSlideale::DirichletNeumannSlideale(MPI_Comm comm) : DirichletNeumann(comm)
{
  // empty constructor
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::DirichletNeumannSlideale::setup()
{
  // call setup of base class
  FSI::DirichletNeumann::setup();

  const Teuchos::ParameterList& fsidyn = Global::Problem::instance()->fsi_dynamic_params();
  const Teuchos::ParameterList& fsipart = fsidyn.sublist("PARTITIONED SOLVER");
  set_kinematic_coupling(Teuchos::getIntegralValue<Inpar::FSI::CoupVarPart>(
                             fsipart, "COUPVARIABLE") == Inpar::FSI::CoupVarPart::disp);

  auto aletype = Teuchos::getIntegralValue<Inpar::FSI::SlideALEProj>(
      Global::Problem::instance()->fsi_dynamic_params(), "SLIDEALEPROJ");

  slideale_ = std::make_shared<FSI::Utils::SlideAleUtils>(structure_field()->discretization(),
      mb_fluid_field()->discretization(), structure_fluid_coupling_mortar(), true, aletype);

  islave_ = std::make_shared<Core::LinAlg::Vector<double>>(
      *structure_fluid_coupling_mortar().slave_dof_map(), true);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::DirichletNeumannSlideale::remeshing()
{
  // dispn and dispnp of structure, used for surface integral and velocity of the fluid in the
  // interface
  std::shared_ptr<Core::LinAlg::Vector<double>> idisptotal =
      structure_field()->extract_interface_dispnp();

  slideale_->remeshing(*structure_field(), *mb_fluid_field()->discretization(), *idisptotal,
      *islave_, structure_fluid_coupling_mortar(), get_comm());

  // Evaluate solid/fluid Mortar coupling
  slideale_->evaluate_mortar(
      *structure_field()->extract_interface_dispnp(), *islave_, structure_fluid_coupling_mortar());
  // Evaluate solid/ale Mortar coupling
  slideale_->evaluate_fluid_mortar(idisptotal, islave_);

  std::shared_ptr<Core::LinAlg::Vector<double>> unew =
      slideale_->interpolate_fluid(*mb_fluid_field()->extract_interface_velnp());
  mb_fluid_field()->apply_interface_values(islave_, unew);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::Vector<double>> FSI::DirichletNeumannSlideale::fluid_op(
    std::shared_ptr<Core::LinAlg::Vector<double>> idispcurr, const FillType fillFlag)
{
  FSI::Partitioned::fluid_op(idispcurr, fillFlag);

  if (fillFlag == User)
  {
    FOUR_C_THROW("not implemented");
    // SD relaxation calculation
    return fluid_to_struct(mb_fluid_field()->relaxation_solve(struct_to_fluid(idispcurr), dt()));
  }
  else
  {
    // normal fluid solve

    // the displacement -> velocity conversion at the interface
    const std::shared_ptr<Core::LinAlg::Vector<double>> ivel = interface_velocity(*idispcurr);

    // A rather simple hack. We need something better!
    const int itemax = mb_fluid_field()->itemax();
    if (fillFlag == MF_Res and mfresitemax_ > 0) mb_fluid_field()->set_itemax(mfresitemax_ + 1);

    // new Core::LinAlg::Vector<double> for aledisp in interface
    std::shared_ptr<Core::LinAlg::Vector<double>> iale =
        std::make_shared<Core::LinAlg::Vector<double>>(
            *(structure_fluid_coupling_mortar().master_dof_map()), true);

    std::shared_ptr<Core::LinAlg::Vector<double>> idispn =
        structure_field()->extract_interface_dispn();

    iale->update(1.0, *idispcurr, 0.0);

    // iale reduced by old displacement dispn and instead added the real last displacements
    iale->update(1.0, *ft_stemp_, -1.0, *idispn, 1.0);

    mb_fluid_field()->nonlinear_solve(struct_to_fluid(iale), struct_to_fluid(ivel));

    mb_fluid_field()->set_itemax(itemax);

    return fluid_to_struct(mb_fluid_field()->extract_interface_forces());
  }
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::Vector<double>> FSI::DirichletNeumannSlideale::struct_op(
    std::shared_ptr<Core::LinAlg::Vector<double>> iforce, const FillType fillFlag)
{
  FSI::Partitioned::struct_op(iforce, fillFlag);

  if (fillFlag == User)
  {
    // SD relaxation calculation
    return structure_field()->relaxation_solve(iforce);
  }
  else
  {
    // normal structure solve
    structure_field()->apply_interface_forces(iforce);
    structure_field()->solve();
    return structure_field()->extract_interface_dispnp();
  }
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::Vector<double>> FSI::DirichletNeumannSlideale::initial_guess()
{
  if (get_kinematic_coupling())
  {
    // real displacement of slave side at time step begin on master side --> for calculation of
    // FluidOp
    ft_stemp_ = fluid_to_struct(islave_);
    // predict displacement
    return structure_field()->predict_interface_dispnp();
  }
  else
  {
    const Teuchos::ParameterList& fsidyn = Global::Problem::instance()->fsi_dynamic_params();
    const Teuchos::ParameterList& fsipart = fsidyn.sublist("PARTITIONED SOLVER");
    if (fsipart.get<std::string>("PREDICTOR") != "d(n)")
    {
      FOUR_C_THROW(
          "unknown interface force predictor '{}'", fsipart.get<std::string>("PREDICTOR").c_str());
    }
    return interface_force();
  }
}

FOUR_C_NAMESPACE_CLOSE
