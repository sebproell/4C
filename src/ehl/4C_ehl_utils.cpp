// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_ehl_utils.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*                                                        AN,JH 09/2014 */
/* Function for checking that the different time steps are a
 multiplicative of each other                                           */

int EHL::Utils::check_time_stepping(double dt1, double dt2)
{
  double workdt1 = std::min(dt1, dt2);
  double workdt2 = std::max(dt1, dt2);
  double t1 = 0.0;
  int i = 0;

  while (true)
  {
    i++;
    t1 = i * workdt1;

    if (std::abs(t1 - workdt2) < 10E-10)
      break;

    else if (t1 > workdt2)
      FOUR_C_THROW("Chosen time steps {} and {} are not a multiplicative of each other", dt1, dt2);
  }
  return i;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*                                                        AN,JH 10/2014 */
// Modification of time parameter list for problem with different time step size

void EHL::Utils::change_time_parameter(MPI_Comm comm, Teuchos::ParameterList& ehlparams,
    Teuchos::ParameterList& lubricationdyn, Teuchos::ParameterList& sdyn)
{
  const bool difftimestep = ehlparams.get<bool>("DIFFTIMESTEPSIZE");

  if (difftimestep)  // Create subproblems with different time steps
  {
    // Check correct choice of time stepping for single fields
    double lubricationstep = lubricationdyn.get<double>("TIMESTEP");
    double solidstep = sdyn.get<double>("TIMESTEP");

    check_time_stepping(lubricationstep, solidstep);

    // modify global time step size
    ehlparams.set<double>("TIMESTEP", std::min(lubricationstep, solidstep));
  }
  else
  {
    // -------------------------------------------------------------------
    // overrule certain parameters for coupled problems
    // -------------------------------------------------------------------
    // the default time step size
    lubricationdyn.set<double>("TIMESTEP", ehlparams.get<double>("TIMESTEP"));
    sdyn.set<double>("TIMESTEP", ehlparams.get<double>("TIMESTEP"));
    // maximum simulation time
    lubricationdyn.set<double>("MAXTIME", ehlparams.get<double>("MAXTIME"));
    sdyn.set<double>("MAXTIME", ehlparams.get<double>("MAXTIME"));
    // maximum number of timesteps
    lubricationdyn.set<int>("NUMSTEP", ehlparams.get<int>("NUMSTEP"));
    sdyn.set<int>("NUMSTEP", ehlparams.get<int>("NUMSTEP"));
  }

  // Check correct input of restart. Code relies that both time value RESTARTEVERYTIME and
  // RESULTSEVERYTIME are given if restart from time is applied
  double restarttime = ehlparams.get<double>("RESTARTEVERYTIME");
  double updatetime = ehlparams.get<double>("RESULTSEVERYTIME");
  if ((updatetime > 0.0) or (restarttime > 0.0))
    if (!(updatetime > 0.0) and !(restarttime > 0.0))
      FOUR_C_THROW(
          "If time controlled output and restart is desired, both parameters RESTARTEVERYTIME and "
          "RESULTSEVERYTIME has to be set");

  // set restart params
  int lubricationrestart;
  int structurerestart;

  if (restarttime > 0.0)
  {
    lubricationrestart = check_time_stepping(lubricationdyn.get<double>("TIMESTEP"), restarttime);
    structurerestart = check_time_stepping(sdyn.get<double>("TIMESTEP"), restarttime);
  }
  else
  {
    int restart = ehlparams.get<int>("RESTARTEVERY");
    lubricationrestart = restart;
    structurerestart = restart;
  }

  // set output params
  int lubricationupres;
  int structureupres;

  if (updatetime > 0.0)
  {
    lubricationupres = check_time_stepping(lubricationdyn.get<double>("TIMESTEP"), updatetime);
    structureupres = check_time_stepping(sdyn.get<double>("TIMESTEP"), updatetime);
  }
  else
  {
    int update = ehlparams.get<int>("RESULTSEVERY");
    lubricationupres = update;
    structureupres = update;
  }

  // restart
  lubricationdyn.set<int>("RESTARTEVERY", lubricationrestart);
  sdyn.set<int>("RESTARTEVERY", structurerestart);
  // solution output
  lubricationdyn.set<int>("RESULTSEVERY", lubricationupres);
  sdyn.set<int>("RESULTSEVERY", structureupres);

  if (Core::Communication::my_mpi_rank(comm) == 0)
  {
    std::cout << "====================== Overview of chosen time stepping: "
                 "==============================\n"
              << "\t Timestep lubrication:           " << lubricationdyn.get<double>("TIMESTEP")
              << "\n"
              << "\t Timestep structure:        " << sdyn.get<double>("TIMESTEP") << "\n"
              << "\t Result step lubrication:        " << lubricationdyn.get<int>("RESULTSEVERY")
              << "\n"
              << "\t Result step structure:     " << sdyn.get<int>("RESULTSEVERY") << "\n"
              << "\t Restart step lubrication:       " << lubricationdyn.get<int>("RESTARTEVERY")
              << "\n"
              << "\t Restart step structure:    " << sdyn.get<int>("RESTARTEVERY") << "\n"
              << "================================================================================="
                 "=======\n \n";
  }
}

/*----------------------------------------------------------------------*
 | print EHL-logo                                            Faraji 05/19 |
 *----------------------------------------------------------------------*/
void EHL::printlogo()
{
  // more at http://www.ascii-art.de
  std::cout
      << "                                    ,/                                              \n"
      << "                                    #%                                              \n"
      << "                                   ,&&                                              \n"
      << "                                  ,&%%(                                             \n"
      << "                                  %&&&&%                                            \n"
      << "                                 /&&%%%%#                                           \n"
      << "                                ,&&%&&&&%#                                          \n"
      << "                                %&&&&&&&&&.                                         \n"
      << "                               .#&&&&&&&%%*                                         \n"
      << "                                %&&%&&&&%&                                          \n"
      << "                     *%%(        *&&%%%&/           .&&&%&                          \n"
      << "            .#&%,    &&&&%,                 /*(     *%&&&&,    %&&&%                \n"
      << "           %%%&&%%%&&&&&&&,  /&&&          %&&&&/  /%%&&&&&&%#&&&&&&                \n"
      << "            #%&&&&%%&%&&&&&#&%&&&          &&&&&%&&&&&&&&&&&&&&&&/&&                \n"
      << "           .#&&&%,   *#&&&&&&%&&/          %&&&&&&&&&&&&&&&&&&&&&*%&                \n"
      << "       &%&%&&&%%.         /&&&&%%*   *&&  ,%&&&&&&%&&&%###%&&%%&&&&&&%%&&&%&        \n"
      << "      *%&&%&&%%            .&&&&.   /&&&&%&&&&%&&&(           .%&&&&&&&&&&&/        \n"
      << "      (&&%&&&&*             %%&%%,  %&&&&&&&&&&&/                #&&&&&&&/          \n"
      << "         .%&&&*             %%&&&&&&. *%%&&&&&&                   (&&&&&&&          \n"
      << "          /%&%&            *&&&&&&&&   %&&&&&&.                    /%&&&&&&&&       \n"
      << "        ,&&&&&&&(        .%&&&&(,**,.*#&&&&&%%                      &&&&&&&&&&&.    \n"
      << "       .%%&%&&&&&%&%#(#%&&&&&&    #%&&&&&&&&%%                      &%&&&&&%&&%,    \n"
      << "         #&&( ,&&&&&&%%&&&&&&%(   /&&&&&&&&&%&                     .&&&&&%/         \n"
      << "                (%&&&%*  ,&&&.       /%&&&&%(                   ,%&&&&&&,           \n"
      << "                %&%%&,     *.           %&&&&&&(                 ,&&&&&&&&%&(.      \n"
      << "                                      .%&&&&&&&&&(             ,%%&&&&&&&&&%&/      \n"
      << "                                     (%&&&&&&&&&&&&&%/,    *(%&&&&&&&%&&&%&&,       \n"
      << "                                      (&&&%(%&&&&&&&&&%&&&&&&&&&&&&&                \n"
      << "                                              %&&&&&&&&&&&&&&&&&&&&&%&              \n"
      << "                                             .&&&&%&&%%&&&&&&&&(&&%&&&%.            \n"
      << "                                             #%&%%%    %%&&&&     %,                \n"
      << "                                               ,(*     *&&&&%                       \n"
      << "                                                        ...                         \n"
      << "                                                                                    \n"
      << "\n"
      << std::endl;

}  // printlogo()

/*----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_CLOSE
