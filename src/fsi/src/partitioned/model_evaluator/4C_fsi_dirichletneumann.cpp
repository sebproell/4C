// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fsi_dirichletneumann.hpp"

#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_global_data.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::DirichletNeumann::DirichletNeumann(MPI_Comm comm)
    : Partitioned(comm), kinematiccoupling_(false)
{
  // empty constructor
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::DirichletNeumann::setup()
{
  /// call setup of base class
  FSI::Partitioned::setup();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::DirichletNeumann::fsi_op(
    const Core::LinAlg::Vector<double>& x, Core::LinAlg::Vector<double>& F, const FillType fillFlag)
{
  if (kinematiccoupling_)  // coupling variable: interface displacements/velocity
  {
    const std::shared_ptr<Core::LinAlg::Vector<double>> icoupn =
        std::make_shared<Core::LinAlg::Vector<double>>(x);

    const std::shared_ptr<Core::LinAlg::Vector<double>> iforce = fluid_op(icoupn, fillFlag);

    const std::shared_ptr<Core::LinAlg::Vector<double>> icoupnp = struct_op(iforce, fillFlag);

    F.update(1.0, *icoupnp, -1.0, *icoupn, 0.0);
  }
  else  // coupling variable: interface forces
  {
    const std::shared_ptr<Core::LinAlg::Vector<double>> iforcen =
        std::make_shared<Core::LinAlg::Vector<double>>(x);

    const std::shared_ptr<Core::LinAlg::Vector<double>> icoupn = struct_op(iforcen, fillFlag);

    const std::shared_ptr<Core::LinAlg::Vector<double>> iforcenp = fluid_op(icoupn, fillFlag);

    F.update(1.0, *iforcenp, -1.0, *iforcen, 0.0);
  }
}

FOUR_C_NAMESPACE_CLOSE
