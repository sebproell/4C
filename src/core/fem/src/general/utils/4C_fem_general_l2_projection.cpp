// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fem_general_l2_projection.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linear_solver_method.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::MultiVector<double>> Core::FE::evaluate_and_solve_nodal_l2_projection(
    Core::FE::Discretization& dis, const Epetra_Map& noderowmap, const std::string& statename,
    const int& numvec, Teuchos::ParameterList& params, const Teuchos::ParameterList& solverparams,
    const std::function<const Teuchos::ParameterList&(int)> get_solver_params,
    const Epetra_Map& fullnoderowmap, const std::map<int, int>& slavetomastercolnodesmap)
{
  // create empty matrix
  Core::LinAlg::SparseMatrix massmatrix(noderowmap, 108, false, true);
  // create empty right hand side
  Core::LinAlg::MultiVector<double> rhs(noderowmap, numvec);

  std::vector<int> lm;
  std::vector<int> lmowner;
  std::vector<int> lmstride;
  Core::Elements::LocationArray la(dis.num_dof_sets());

  // define element matrices and vectors
  Core::LinAlg::SerialDenseMatrix elematrix1;
  Core::LinAlg::SerialDenseMatrix elematrix2;
  Core::LinAlg::SerialDenseVector elevector1;
  Core::LinAlg::SerialDenseVector elevector2;
  Core::LinAlg::SerialDenseVector elevector3;

  // loop column elements

  for (auto* actele : dis.my_col_element_range())
  {
    const int numnode = actele->num_node();

    actele->location_vector(dis, la, false);
    lmowner = la[0].lmowner_;
    lmstride = la[0].stride_;
    lm = la[0].lm_;

    // Reshape element matrices and vectors and initialize to zero
    elevector1.size(numnode);
    elematrix1.shape(numnode, numnode);
    elematrix2.shape(numnode, numvec);

    // call the element specific evaluate method (elemat1 = mass matrix, elemat2 = rhs)
    int err = actele->evaluate(
        params, dis, la, elematrix1, elematrix2, elevector1, elevector2, elevector3);
    if (err) FOUR_C_THROW("Element {} returned err={}", actele->id(), err);


    // get element location vector for nodes
    lm.resize(numnode);
    lmowner.resize(numnode);

    Core::Nodes::Node** nodes = actele->nodes();
    for (int n = 0; n < numnode; ++n)
    {
      const int nodeid = nodes[n]->id();
      if (!slavetomastercolnodesmap.empty())
      {
        auto slavemasterpair = slavetomastercolnodesmap.find(nodeid);
        if (slavemasterpair != slavetomastercolnodesmap.end())
          lm[n] = slavemasterpair->second;
        else
          lm[n] = nodeid;
      }
      else
        lm[n] = nodeid;

      // owner of pbc master and slave nodes are identical
      lmowner[n] = nodes[n]->owner();
    }

    // mass matrix assembling into node map
    massmatrix.assemble(actele->id(), elematrix1, lm, lmowner);

    // assemble numvec entries sequentially
    for (int n = 0; n < numvec; ++n)
    {
      // copy results into Serial_DenseVector for assembling
      for (int inode = 0; inode < numnode; ++inode) elevector1(inode) = elematrix2(inode, n);
      // assemble into nth vector of MultiVector
      Core::LinAlg::assemble(rhs, n, elevector1, lm, lmowner);
    }
  }  // end element loop

  // finalize the matrix
  massmatrix.complete();

  return solve_nodal_l2_projection(massmatrix, rhs, dis.get_comm(), numvec, solverparams,
      get_solver_params, noderowmap, fullnoderowmap, slavetomastercolnodesmap);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::MultiVector<double>> Core::FE::compute_nodal_l2_projection(
    Core::FE::Discretization& dis, const std::string& statename, const int& numvec,
    Teuchos::ParameterList& params, const Teuchos::ParameterList& solverparams,
    const std::function<const Teuchos::ParameterList&(int)> get_solver_params)
{
  // check if the statename has been set
  if (!dis.has_state(statename))
  {
    FOUR_C_THROW(
        "The discretization does not know about this statename. Please "
        "review how you call this function.");
  }

  // check whether action type is set
  if (params.getEntryRCP("action") == Teuchos::null)
    FOUR_C_THROW("action type for element is missing");

  // handle pbcs if existing
  // build inverse map from slave to master nodes
  std::map<int, int> slavetomastercolnodesmap;

  std::map<int, std::vector<int>>* allcoupledcolnodes = dis.get_all_pbc_coupled_col_nodes();
  if (allcoupledcolnodes)
  {
    for (auto [master_gid, slave_gids] : *allcoupledcolnodes)
    {
      for (const auto slave_gid : slave_gids)
      {
        slavetomastercolnodesmap[slave_gid] = master_gid;
      }
    }
  }

  // get reduced node row map of fluid field --> will be used for setting up linear system
  const auto* fullnoderowmap = dis.node_row_map();
  // remove pbc slave nodes from full noderowmap
  std::vector<int> reducednoderowmap;
  // a little more memory than necessary is possibly reserved here
  reducednoderowmap.reserve(fullnoderowmap->NumMyElements());
  for (int i = 0; i < fullnoderowmap->NumMyElements(); ++i)
  {
    const int nodeid = fullnoderowmap->GID(i);
    // do not add slave pbc nodes here
    if (slavetomastercolnodesmap.empty() or slavetomastercolnodesmap.count(nodeid) == 0)
    {
      reducednoderowmap.push_back(nodeid);
    }
  }

  // build node row map which does not include slave pbc nodes
  Epetra_Map noderowmap(-1, static_cast<int>(reducednoderowmap.size()), reducednoderowmap.data(), 0,
      fullnoderowmap->Comm());

  auto nodevec = evaluate_and_solve_nodal_l2_projection(dis, noderowmap, statename, numvec, params,
      solverparams, get_solver_params, *fullnoderowmap, slavetomastercolnodesmap);

  // if no pbc are involved leave here
  if (slavetomastercolnodesmap.empty() or noderowmap.PointSameAs(*fullnoderowmap)) return nodevec;

  // solution vector based on full row map in which the solution of the master node is inserted into
  // slave nodes
  auto fullnodevec = std::make_shared<Core::LinAlg::MultiVector<double>>(*fullnoderowmap, numvec);

  for (int i = 0; i < fullnoderowmap->NumMyElements(); ++i)
  {
    const int nodeid = fullnoderowmap->GID(i);

    auto slavemasterpair = slavetomastercolnodesmap.find(nodeid);
    if (slavemasterpair != slavetomastercolnodesmap.end())
    {
      const int mastergid = slavemasterpair->second;
      const int masterlid = noderowmap.LID(mastergid);
      for (int j = 0; j < numvec; ++j) fullnodevec->ReplaceMyValue(i, j, (*nodevec)(j)[masterlid]);
    }
    else
    {
      const int lid = noderowmap.LID(nodeid);
      for (int j = 0; j < numvec; ++j) fullnodevec->ReplaceMyValue(i, j, (*nodevec)(j)[lid]);
    }
  }

  return fullnodevec;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::MultiVector<double>> Core::FE::solve_nodal_l2_projection(
    Core::LinAlg::SparseMatrix& massmatrix, Core::LinAlg::MultiVector<double>& rhs, MPI_Comm comm,
    const int& numvec, const Teuchos::ParameterList& solverparams,
    const std::function<const Teuchos::ParameterList&(int)> get_solver_params,
    const Epetra_Map& noderowmap, const Epetra_Map& fullnoderowmap,
    const std::map<int, int>& slavetomastercolnodesmap)
{
  // get solver parameter list of linear solver
  const auto solvertype =
      Teuchos::getIntegralValue<Core::LinearSolver::SolverType>(solverparams, "SOLVER");

  Core::LinAlg::Solver solver(
      solverparams, comm, get_solver_params, Core::IO::Verbositylevel::standard);

  // skip setup of preconditioner in case of a direct solver
  if (solvertype != Core::LinearSolver::SolverType::umfpack and
      solvertype != Core::LinearSolver::SolverType::superlu)
  {
    const auto prectype =
        Teuchos::getIntegralValue<Core::LinearSolver::PreconditionerType>(solverparams, "AZPREC");
    switch (prectype)
    {
      case Core::LinearSolver::PreconditionerType::multigrid_muelu:
      {
        Teuchos::ParameterList* preclist_ptr = nullptr;
        // Parameter for MueLu
        if (prectype == Core::LinearSolver::PreconditionerType::multigrid_muelu)
          preclist_ptr = &((solver.params()).sublist("MueLu Parameters"));
        else
          FOUR_C_THROW("please add correct parameter list");

        Teuchos::ParameterList& preclist = *preclist_ptr;
        preclist.set("PDE equations", 1);
        preclist.set("null space: dimension", 1);
        preclist.set("null space: type", "pre-computed");
        preclist.set("null space: add default vectors", false);

        std::shared_ptr<Core::LinAlg::MultiVector<double>> nullspace =
            std::make_shared<Core::LinAlg::MultiVector<double>>(noderowmap, 1, true);
        nullspace->PutScalar(1.0);

        preclist.set<std::shared_ptr<Core::LinAlg::MultiVector<double>>>("nullspace", nullspace);
        preclist.set("null space: vectors", nullspace->Values());
        preclist.set("ML validate parameter list", false);
      }
      break;
      case Core::LinearSolver::PreconditionerType::ilu:
        // do nothing
        break;
      default:
        FOUR_C_THROW("You have to choose ML, MueLu or ILU preconditioning");
        break;
    }
  }

  // solution vector based on reduced node row map
  auto nodevec = std::make_shared<Core::LinAlg::MultiVector<double>>(noderowmap, numvec);

  switch (solvertype)
  {
    case Core::LinearSolver::SolverType::belos:
    {
      // solve for numvec rhs at the same time using Belos solver
      Core::LinAlg::SolverParams solver_params;
      solver_params.refactor = true;
      solver_params.reset = true;
      solver.solve_with_multi_vector(massmatrix.epetra_operator(), nodevec,
          Core::Utils::shared_ptr_from_ref(rhs), solver_params);
      break;
    }
    default:
    {
      if (numvec != 1 and Core::Communication::my_mpi_rank(comm) == 0)
      {
        std::cout << "Think about using a Belos solver which can handle several rhs vectors at the "
                     "same time\n";
      }

      // solve for numvec rhs iteratively
      for (int i = 0; i < numvec; i++)
      {
        Core::LinAlg::SolverParams solver_params;
        solver_params.refactor = true;
        solver_params.reset = true;
        solver.solve_with_multi_vector(massmatrix.epetra_operator(),
            (*nodevec)(i).get_ptr_of_multi_vector(), rhs(i).get_ptr_of_multi_vector(),
            solver_params);
      }
      break;
    }
  }

  return nodevec;
}

FOUR_C_NAMESPACE_CLOSE