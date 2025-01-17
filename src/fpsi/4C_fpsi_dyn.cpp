// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fpsi_dyn.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fpsi.hpp"
#include "4C_fpsi_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_fpsi.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*------------------------------------------------------------------------------------------------*
 | main control routine for fluid-porous-structure-interaction problems                rauch 11/12 |
 *------------------------------------------------------------------------------------------------*/
void fpsi_drt()
{
  Global::Problem* problem = Global::Problem::instance();

  // 1.- Get Communicator
  MPI_Comm comm = problem->get_dis("structure")->get_comm();

  // print the chuck
  if (Core::Communication::my_mpi_rank(comm) == 0)
  {
    std::cout << "                            ``;@@@@@@@@@#@`               " << std::endl;
    std::cout << "                            .@#@#@@@@@@@@@@,+.`           " << std::endl;
    std::cout << "                          '##@@@@@@@@@@@@@@@+;`           " << std::endl;
    std::cout << "                         :@#@@@@@@@@@@@@@@@@'`#           " << std::endl;
    std::cout << "                         @@@#@@@@@@@@@@@@@@@#@#:          " << std::endl;
    std::cout << "                        #@@@@@@@@@@@@@@@@@@@@@@@          " << std::endl;
    std::cout << "                       ,##@@@@@@@@@@@@@@@@@@@@@@:         " << std::endl;
    std::cout << "                       ;###@@@@@@@@#@@@@@@@@@@@@#         " << std::endl;
    std::cout << "                        @@@@@@@@@@@@@#@#@@@@@@@@@         " << std::endl;
    std::cout << "                       `@#@@@@@@@@@@@@@@@@@@@@@@@         " << std::endl;
    std::cout << "                       ,@#@@@@@@@@@#':#@@##@@@@@@`        " << std::endl;
    std::cout << "                       +##@@@@@@@@@+: '#@@#@@@@@@@        " << std::endl;
    std::cout << "                       :@#@@@@@@@@@;    .`#@@@@@@@        " << std::endl;
    std::cout << "                       ,@@@@@@@@##;.     `@@#@@@##.       " << std::endl;
    std::cout << "                       .@##@@;.@@@        @@#@@@@#+       " << std::endl;
    std::cout << "                        @@@@#@#@#`         @##@@#@#       " << std::endl;
    std::cout << "                        @@@@@  `           #@#@@@##       " << std::endl;
    std::cout << "                        ,#@@````           +@@#@#@#       " << std::endl;
    std::cout << "                         @,@.   `..        +@#@,@##       " << std::endl;
    std::cout << "                         #@@,@# #@@@@#      @@ # @#       " << std::endl;
    std::cout << "                         '##@@@ ;@@@#@'     ` + :##       " << std::endl;
    std::cout << "                         `@@@@@  @ ``      ';   #@`       " << std::endl;
    std::cout << "                          @@@##            .@   @@        " << std::endl;
    std::cout << "                          #  @+            @:  :##        " << std::endl;
    std::cout << "                          '  @.            @,  @@,        " << std::endl;
    std::cout << "                          #  @            `@` '@@         " << std::endl;
    std::cout << "                          , `@            .#: ;#.         " << std::endl;
    std::cout << "                           +'@     #       #   ,          " << std::endl;
    std::cout << "                           @@@   +       :.`   .          " << std::endl;
    std::cout << "                           @##@@        ,@'@   @;         " << std::endl;
    std::cout << "                           @@@@@@#': #. @#@+   ;@         " << std::endl;
    std::cout << "                     @;+#@@@@####@@@@@@@@@@`   +#'        " << std::endl;
    std::cout << "                    @@#@@@#@@@@+;   ;@@#:#'    #@@        " << std::endl;
    std::cout << "                    ;#@@@#@@@@@'.   :#@:@@    ;#@@        " << std::endl;
    std::cout << "                     @#@@@@@@@@@@@@@@###@@   ,@@@@;       " << std::endl;
    std::cout << "                      @#@@@@@@@@#@@#@@@@@+   @@@#@#+      " << std::endl;
    std::cout << "                     :.#@@@@@@@#@#@##@###  `@#@@@@@@#     " << std::endl;
    std::cout << "                      _______ ______ _______  ______      " << std::endl;
    std::cout << "                        ||____ ||___|||_____    ||        " << std::endl;
    std::cout << "                        ||     ||     _____|| __||__      " << std::endl;
    std::cout << std::endl << std::endl;
  }

  // 2.- Parameter reading
  const Teuchos::ParameterList& fpsidynparams = problem->fpsi_dynamic_params();
  const Teuchos::ParameterList& poroelastdynparams = problem->poroelast_dynamic_params();

  FPSI::InterfaceUtils* FPSI_UTILS = FPSI::InterfaceUtils::instance();

  // 3.- Creation of Poroelastic + Fluid problem. (discretization called inside)
  std::shared_ptr<FPSI::FpsiBase> fpsi = nullptr;
  fpsi = FPSI_UTILS->setup_discretizations(comm, fpsidynparams, poroelastdynparams);

  // 3.1- Read restart if needed.
  const int restartstep = problem->restart();
  if (restartstep)
  {
    fpsi->read_restart(restartstep);
  }

  // 3.2.- redistribute the FPSI interface
  fpsi->redistribute_interface();

  //////////////////////////////////
  // 4.- Run of the actual problem.//
  //////////////////////////////////

  // 4.1.- Coupling and creation of combined dofmap
  fpsi->setup_system();
  // setup the linear solver, if necessary
  fpsi->setup_solver();
  // 4.2.- Solve the whole problem
  fpsi->timeloop();
  Teuchos::TimeMonitor::summarize();

  // 5. - perform the result test
  fpsi->test_results(comm);


  return;
}  // fpsi_drt()

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_CLOSE
