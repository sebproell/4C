// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_linear_solver_method_parameters.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_discretization_nullspace.hpp"
#include "4C_fem_general_elementtype.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_utils_exceptions.hpp"

#include <Xpetra_EpetraIntMultiVector.hpp>

FOUR_C_NAMESPACE_OPEN

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void Core::LinearSolver::Parameters::compute_solver_parameters(
    Core::FE::Discretization& dis, Teuchos::ParameterList& solverlist)
{
  std::shared_ptr<Epetra_Map> nullspaceMap =
      solverlist.get<std::shared_ptr<Epetra_Map>>("null space: map", nullptr);

  int numdf = 1;
  int dimns = 1;
  int nv = 0;
  int np = 0;

  // set parameter information for solver
  {
    if (nullspaceMap == nullptr and dis.num_my_row_nodes() > 0)
    {
      // no map given, just grab the block information on the first element that appears
      Core::Elements::Element* dwele = dis.l_row_element(0);
      dwele->element_type().nodal_block_information(dwele, numdf, dimns, nv, np);
    }
    else
    {
      // if a map is given, grab the block information of the first element in that map
      for (int i = 0; i < dis.num_my_row_nodes(); ++i)
      {
        Core::Nodes::Node* actnode = dis.l_row_node(i);
        std::vector<int> dofs = dis.dof(0, actnode);

        const int localIndex = nullspaceMap->LID(dofs[0]);

        if (localIndex == -1) continue;

        Core::Elements::Element* dwele = dis.l_row_element(localIndex);
        actnode->elements()[0]->element_type().nodal_block_information(dwele, numdf, dimns, nv, np);
        break;
      }
    }

    // communicate data to procs without row element
    std::array<int, 4> ldata{numdf, dimns, nv, np};
    std::array<int, 4> gdata{0, 0, 0, 0};
    Core::Communication::max_all(ldata.data(), gdata.data(), 4, dis.get_comm());
    numdf = gdata[0];
    dimns = gdata[1];
    nv = gdata[2];
    np = gdata[3];

    // store nullspace information in solver list
    solverlist.set("PDE equations", numdf);
    solverlist.set("null space: dimension", dimns);
    solverlist.set("null space: type", "pre-computed");
    solverlist.set("null space: add default vectors", false);
  }

  // set coordinate information
  {
    std::shared_ptr<Core::LinAlg::MultiVector<double>> coordinates;
    if (nullspaceMap == nullptr)
      coordinates = dis.build_node_coordinates();
    else
      coordinates = dis.build_node_coordinates(nullspaceMap);

    solverlist.set<std::shared_ptr<Core::LinAlg::MultiVector<double>>>("Coordinates", coordinates);
  }

  // set nullspace information
  {
    if (nullspaceMap == nullptr)
    {
      // if no map is given, we calculate the nullspace on the map describing the
      // whole discretization
      nullspaceMap = std::make_shared<Epetra_Map>(*dis.dof_row_map());
    }

    auto nullspace = Core::FE::compute_null_space(dis, numdf, dimns, *nullspaceMap);

    solverlist.set<std::shared_ptr<Core::LinAlg::MultiVector<double>>>("nullspace", nullspace);
    solverlist.set("null space: vectors", nullspace->Values());
    solverlist.set<bool>("ML validate parameter list", false);
  }
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void Core::LinearSolver::Parameters::fix_null_space(std::string field, const Epetra_Map& oldmap,
    const Epetra_Map& newmap, Teuchos::ParameterList& solveparams)
{
  if (!Core::Communication::my_mpi_rank(Core::Communication::unpack_epetra_comm(oldmap.Comm())))
    printf("Fixing %s Nullspace\n", field.c_str());

  // find the ML or MueLu list
  Teuchos::ParameterList* params_ptr = nullptr;
  if (solveparams.isSublist("ML Parameters"))
    params_ptr = &(solveparams.sublist("ML Parameters"));
  else if (solveparams.isSublist("MueLu Parameters"))
    params_ptr = &(solveparams.sublist("MueLu Parameters"));
  else
    params_ptr = &(solveparams);
  Teuchos::ParameterList& params = *params_ptr;

  const int ndim = params.get("null space: dimension", -1);
  if (ndim == -1) FOUR_C_THROW("List does not contain nullspace dimension");

  std::shared_ptr<Core::LinAlg::MultiVector<double>> nullspace =
      params.get<std::shared_ptr<Core::LinAlg::MultiVector<double>>>("nullspace", nullptr);
  if (nullspace == nullptr) FOUR_C_THROW("List does not contain nullspace");

  const int nullspaceLength = nullspace->MyLength();
  const int newmapLength = newmap.NumMyElements();

  if (nullspaceLength == newmapLength) return;
  if (nullspaceLength != oldmap.NumMyElements())
    FOUR_C_THROW("Nullspace map of length {} does not match old map length of {}", nullspaceLength,
        oldmap.NumMyElements());
  if (newmapLength > nullspaceLength)
    FOUR_C_THROW("New problem size larger than old - full rebuild of nullspace necessary");

  std::shared_ptr<Core::LinAlg::MultiVector<double>> nullspaceNew =
      std::make_shared<Core::LinAlg::MultiVector<double>>(newmap, ndim, true);

  for (int i = 0; i < ndim; i++)
  {
    auto& nullspaceData = (*nullspace)(i);
    auto& nullspaceDataNew = (*nullspaceNew)(i);
    const int myLength = nullspaceDataNew.local_length();

    for (int j = 0; j < myLength; j++)
    {
      int gid = newmap.GID(j);
      int olid = oldmap.LID(gid);
      if (olid == -1) continue;
      nullspaceDataNew[j] = nullspaceData[olid];
    }
  }

  params.set<std::shared_ptr<Core::LinAlg::MultiVector<double>>>("nullspace", nullspaceNew);
  params.set("null space: vectors", nullspaceNew->Values());
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
Core::LinearSolver::Parameters::extract_nullspace_from_parameterlist(
    const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& row_map,
    const Teuchos::ParameterList& list)
{
  if (!list.isParameter("null space: dimension"))
    FOUR_C_THROW("Multigrid parameter 'null space: dimension' missing  in solver parameter list.");

  const int nullspace_dimension = list.get<int>("null space: dimension");
  if (nullspace_dimension < 1)
    FOUR_C_THROW("Multigrid parameter 'null space: dimension' wrong. It has to be > 0.");

  auto nullspace_data = list.get<std::shared_ptr<Core::LinAlg::MultiVector<double>>>("nullspace");
  if (!nullspace_data) FOUR_C_THROW("Nullspace data is null.");

  Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> nullspace =
      Teuchos::make_rcp<Xpetra::EpetraMultiVectorT<GlobalOrdinal, Node>>(
          Teuchos::rcpFromRef(*nullspace_data->get_ptr_of_Epetra_MultiVector()));

  nullspace->replaceMap(Teuchos::rcpFromRef(row_map));

  return nullspace;
}

FOUR_C_NAMESPACE_CLOSE
