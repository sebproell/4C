// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_inpar_solver.hpp"

#include "4C_linear_solver_method.hpp"
#include "4C_utils_parameter_list.hpp"

#include <BelosTypes.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Inpar::SOLVER
{
  void set_valid_solver_parameters(Core::IO::InputSpec& spec)
  {
    using namespace Core::IO::InputSpecBuilders;
    Core::Utils::SectionSpecs list{"dummy"};
    // Solver options
    {
      list.specs.emplace_back(deprecated_selection<Core::LinearSolver::SolverType>("SOLVER",
          {
              {"UMFPACK", Core::LinearSolver::SolverType::umfpack},
              {"Superlu", Core::LinearSolver::SolverType::superlu},
              {"Belos", Core::LinearSolver::SolverType::belos},
              {"undefined", Core::LinearSolver::SolverType::undefined},
          },
          {.description = "The solver to attack the system of linear equations arising of FE "
                          "approach with.",
              .default_value = Core::LinearSolver::SolverType::undefined}));
    }

    // Iterative solver options
    {
      list.specs.emplace_back(
          deprecated_selection<Core::LinearSolver::IterativeSolverType>("AZSOLVE",
              {
                  {"CG", Core::LinearSolver::IterativeSolverType::cg},
                  {"GMRES", Core::LinearSolver::IterativeSolverType::gmres},
                  {"BiCGSTAB", Core::LinearSolver::IterativeSolverType::bicgstab},
              },
              {.description = "Type of linear solver algorithm to use.",
                  .default_value = Core::LinearSolver::IterativeSolverType::gmres}));
    }

    // Preconditioner options
    {
      list.specs.emplace_back(deprecated_selection<Core::LinearSolver::PreconditionerType>("AZPREC",
          {
              {"ILU", Core::LinearSolver::PreconditionerType::ilu},
              {"MueLu", Core::LinearSolver::PreconditionerType::multigrid_muelu},
              {"AMGnxn", Core::LinearSolver::PreconditionerType::multigrid_nxn},
              {"Teko", Core::LinearSolver::PreconditionerType::block_teko},
          },
          {.description = "Type of internal preconditioner to use.\nNote! this preconditioner will "
                          "only be used if the input operator\nsupports the Epetra_RowMatrix "
                          "interface and the client does not pass\nin an external preconditioner!",
              .default_value = Core::LinearSolver::PreconditionerType::ilu}));
    }

    // Ifpack options
    {
      list.specs.emplace_back(parameter<int>("IFPACKOVERLAP",
          {.description = "The amount of overlap used for the ifpack \"ilu\" preconditioner.",
              .default_value = 0}));

      list.specs.emplace_back(parameter<int>("IFPACKGFILL",
          {.description = "The amount of fill allowed for an internal \"ilu\" preconditioner.",
              .default_value = 0}));

      std::vector<std::string> ifpack_combine_valid_input = {"Add", "Insert", "Zero"};
      list.specs.emplace_back(
          deprecated_selection<std::string>("IFPACKCOMBINE", ifpack_combine_valid_input,
              {.description = "Combine mode for Ifpack Additive Schwarz", .default_value = "Add"}));
    }

    // Iterative solver options
    {
      list.specs.emplace_back(parameter<int>(
          "AZITER", {.description = "The maximum number of iterations the underlying iterative "
                                    "solver is allowed to perform",
                        .default_value = 1000}));

      list.specs.emplace_back(parameter<double>("AZTOL",
          {.description =
                  "The level the residual norms must reach to decide about successful convergence",
              .default_value = 1e-8}));

      list.specs.emplace_back(deprecated_selection<Belos::ScaleType>("AZCONV",
          {
              {"AZ_r0", Belos::ScaleType::NormOfInitRes},
              {"AZ_noscaled", Belos::ScaleType::None},
          },
          {.description = "The implicit residual norm scaling type to use for terminating the "
                          "iterative solver.",
              .default_value = Belos::ScaleType::NormOfInitRes}));

      list.specs.emplace_back(parameter<int>(
          "AZOUTPUT", {.description = "The number of iterations between each output of the "
                                      "solver's progress is written to "
                                      "screen",
                          .default_value = 0}));

      list.specs.emplace_back(parameter<int>("AZREUSE",
          {.description = "The number specifying how often to recompute some preconditioners",
              .default_value = 0}));

      list.specs.emplace_back(parameter<int>(
          "AZSUB", {.description = "The maximum size of the Krylov subspace used with \"GMRES\" "
                                   "before\n a restart is performed.",
                       .default_value = 50}));

      list.specs.emplace_back(
          Core::IO::InputSpecBuilders::parameter<std::optional<std::filesystem::path>>(
              "SOLVER_XML_FILE", {.description = "xml file defining any iterative solver"}));
    }

    // MueLu options
    {
      list.specs.emplace_back(
          Core::IO::InputSpecBuilders::parameter<std::optional<std::filesystem::path>>(
              "MUELU_XML_FILE", {.description = "xml file defining any MueLu preconditioner"}));
    }

    // Teko options
    {
      list.specs.emplace_back(
          Core::IO::InputSpecBuilders::parameter<std::optional<std::filesystem::path>>(
              "TEKO_XML_FILE", {.description = "xml file defining any Teko preconditioner"}));
    }

    // user-given name of solver block (just for beauty)
    list.specs.emplace_back(parameter<std::string>("NAME",
        {.description = "User specified name for solver block", .default_value = "No_name"}));

    // Parameters for AMGnxn Preconditioner
    {
      list.specs.emplace_back(parameter<std::string>(
          "AMGNXN_TYPE", {.description = "Name of the pre-built preconditioner to be used. If set "
                                         "to\"XML\" the preconditioner is defined using a xml file",
                             .default_value = "AMG(BGS)"}));
      list.specs.emplace_back(
          Core::IO::InputSpecBuilders::parameter<std::optional<std::filesystem::path>>(
              "AMGNXN_XML_FILE", {.description = "xml file defining the AMGnxn preconditioner"}));
    }

    spec = Core::IO::InputSpecBuilders::all_of(std::move(list.specs));
  }


  void set_valid_parameters(std::map<std::string, Core::IO::InputSpec>& list)
  {
    // set valid parameters for solver blocks

    // Note: the maximum number of solver blocks is hardwired here. If you change this,
    // don't forget to edit the corresponding parts in globalproblems.cpp, too.
    for (int i = 1; i < 10; i++)
    {
      std::stringstream ss;
      ss << "SOLVER " << i;
      std::stringstream ss_description;
      ss_description << "solver parameters for solver block " << i;
      set_valid_solver_parameters(list[ss.str()]);
    }
  }

}  // namespace Inpar::SOLVER

FOUR_C_NAMESPACE_CLOSE
