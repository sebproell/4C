// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fbi_calc_utils.hpp"

#include "4C_beam3_base.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_linalg_blocksparsematrix.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/

void FBI::Utils::get_fbi_element_centerline_dof_indices(Core::FE::Discretization const& discret,
    const Core::Elements::Element* ele, std::vector<unsigned int>& ele_centerline_dof_indices,
    unsigned int& num_dof)
{
  // Todo implement method in Core::Elements::Element or find alternative way of doing this
  // find out the elements' number of Dofs (=dimension of element vector/matrices)
  std::vector<int> lmrow;
  std::vector<int> dummy1, dummy2;

  ele->location_vector(discret, lmrow, dummy1, dummy2);
  num_dof = lmrow.size();

  const Discret::Elements::Beam3Base* beamele =
      dynamic_cast<const Discret::Elements::Beam3Base*>(ele);

  if (beamele != nullptr)
  {
    beamele->centerline_dof_indices_of_element(ele_centerline_dof_indices);
  }
  else
  {
    ele_centerline_dof_indices.resize(num_dof * 3 / 4);
    int j = 0;
    for (unsigned int i = 0; i < num_dof; ++i)
    {
      if ((i + 1) % 4 != 0)
      {
        ele_centerline_dof_indices[j] = i;
        j++;
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void FBI::Utils::assemble_centerline_dof_force_stiff_into_fbi_element_force_stiff(
    const Core::FE::Discretization& discretization1,
    const Core::FE::Discretization& discretization2, std::vector<int> const& elegid,
    std::vector<Core::LinAlg::SerialDenseVector> const& eleforce_centerlineDOFs,
    std::vector<std::vector<Core::LinAlg::SerialDenseMatrix>> const& elestiff_centerlineDOFs,
    std::vector<Core::LinAlg::SerialDenseVector>* eleforce,
    std::vector<std::vector<Core::LinAlg::SerialDenseMatrix>>* elestiff)
{
  std::vector<unsigned int> numdof_ele(2);
  std::vector<std::vector<unsigned int>> ele_centerlinedofindices(2);

  // Get DOFs for beam element
  Core::Elements::Element* ele = discretization1.g_element(elegid[0]);
  get_fbi_element_centerline_dof_indices(
      discretization1, ele, ele_centerlinedofindices[0], numdof_ele[0]);

  // Get DOFs for fluid element
  ele = discretization2.g_element(elegid[1]);
  get_fbi_element_centerline_dof_indices(
      discretization2, ele, ele_centerlinedofindices[1], numdof_ele[1]);


  // assemble centerline DOF values correctly into element DOFvec vectors/matrices
  if (eleforce != nullptr)
  {
    for (unsigned int iele = 0; iele < 2; ++iele)
    {
      // resize and clear variable
      ((*eleforce)[iele]).size(numdof_ele[iele]);

      // safety check: dimensions
      if ((unsigned int)eleforce_centerlineDOFs[iele].numRows() !=
          ele_centerlinedofindices[iele].size())
        FOUR_C_THROW(
            "size mismatch! need to assemble {} values of centerline-Dof based "
            "force vector into element vector but only got {} element-local Dof indices",
            eleforce_centerlineDOFs[iele].numRows(), ele_centerlinedofindices[iele].size());

      // Todo maybe use a more general 'SerialDenseAssemble' method here
      for (unsigned int idof = 0; idof < ele_centerlinedofindices[iele].size(); ++idof)
        ((*eleforce)[iele])(ele_centerlinedofindices[iele][idof]) =
            eleforce_centerlineDOFs[iele](idof);
    }
  }

  if (elestiff != nullptr)
  {
    for (unsigned int iele = 0; iele < 2; ++iele)
    {
      for (unsigned int jele = 0; jele < 2; ++jele)
      {
        // resize and clear variable
        ((*elestiff)[iele][jele]).shape(numdof_ele[iele], numdof_ele[jele]);

        // safety check: dimensions
        if ((unsigned int)elestiff_centerlineDOFs[iele][jele].numRows() !=
            ele_centerlinedofindices[iele].size())
          FOUR_C_THROW(
              "size mismatch! need to assemble {} row values of centerline-Dof based "
              "stiffness matrix into element matrix but only got {} element-local Dof indices",
              elestiff_centerlineDOFs[iele][jele].numRows(), ele_centerlinedofindices[iele].size());

        if ((unsigned int)elestiff_centerlineDOFs[iele][jele].numCols() !=
            ele_centerlinedofindices[jele].size())
          FOUR_C_THROW(
              "size mismatch! need to assemble {} column values of centerline-Dof based "
              "stiffness matrix into element matrix but only got {} element-local Dof indices",
              elestiff_centerlineDOFs[iele][jele].numCols(), ele_centerlinedofindices[jele].size());

        for (unsigned int idof = 0; idof < ele_centerlinedofindices[iele].size(); ++idof)
          for (unsigned int jdof = 0; jdof < ele_centerlinedofindices[jele].size(); ++jdof)
            ((*elestiff)[iele][jele])(
                ele_centerlinedofindices[iele][idof], ele_centerlinedofindices[jele][jdof]) =
                elestiff_centerlineDOFs[iele][jele](idof, jdof);
      }
    }
  }
}

FOUR_C_NAMESPACE_CLOSE
