// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_particle_wall_result_test.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_particle_wall_datastate.hpp"
#include "4C_particle_wall_interface.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEWALL::WallResultTest::WallResultTest() : Core::Utils::ResultTest("PARTICLEWALL")
{
  // empty constructor
}

void PARTICLEWALL::WallResultTest::init()
{
  // nothing to do
}

void PARTICLEWALL::WallResultTest::setup(
    const std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface)
{
  // set interface to particle wall handler
  particlewallinterface_ = particlewallinterface;

  // get wall discretization
  walldiscretization_ = particlewallinterface_->get_wall_discretization();
}

void PARTICLEWALL::WallResultTest::test_node(
    const Core::IO::InputParameterContainer& container, int& nerr, int& test_count)
{
  // extract and check discretization name
  std::string dis = container.get<std::string>("DIS");
  if (dis != walldiscretization_->name()) return;

  // extract node id
  int node = container.get<int>("NODE");
  node -= 1;

  int havenode(walldiscretization_->have_global_node(node));
  int havenodeonanyproc(0);
  Core::Communication::sum_all(&havenode, &havenodeonanyproc, 1, walldiscretization_->get_comm());

  // safety check
  if (not havenodeonanyproc)
    FOUR_C_THROW("node {} does not belong to discretization {}", node + 1,
        walldiscretization_->name().c_str());


  if (walldiscretization_->have_global_node(node))
  {
    const Core::Nodes::Node* actnode = walldiscretization_->g_node(node);

    // node not owned on this processor
    if (actnode->owner() != Core::Communication::my_mpi_rank(walldiscretization_->get_comm()))
      return;

    // get wall data state container
    std::shared_ptr<PARTICLEWALL::WallDataState> walldatastate =
        particlewallinterface_->get_wall_data_state();

    // extract test quantity
    std::string quantity = container.get<std::string>("QUANTITY");

    // init actual result
    double actresult = 0.0;

    // position
    if (quantity == "posx" or quantity == "posy" or quantity == "posz")
    {
      // get wall displacements
      std::shared_ptr<const Core::LinAlg::Vector<double>> disp = walldatastate->get_disp_col();

      int idx = -1;
      if (quantity == "posx")
        idx = 0;
      else if (quantity == "posy")
        idx = 1;
      else if (quantity == "posz")
        idx = 2;

      if (idx >= 0)
      {
        actresult = actnode->x()[idx];

        if (disp != nullptr)
        {
          const Epetra_BlockMap& disnpmap = disp->get_map();
          int lid = disnpmap.LID(walldiscretization_->dof(0, actnode, idx));
          if (lid < 0)
            FOUR_C_THROW("You tried to test {} on nonexistent dof {} on node {}", quantity.c_str(),
                idx, actnode->id());
          actresult += (*disp)[lid];
        }
      }
    }
    // displacement
    else if (quantity == "dispx" or quantity == "dispy" or quantity == "dispz")
    {
      // get wall displacements
      std::shared_ptr<const Core::LinAlg::Vector<double>> disp = walldatastate->get_disp_col();

      if (disp == nullptr) return;

      int idx = -1;
      if (quantity == "dispx")
        idx = 0;
      else if (quantity == "dispy")
        idx = 1;
      else if (quantity == "dispz")
        idx = 2;

      if (idx >= 0)
      {
        const Epetra_BlockMap& disnpmap = disp->get_map();
        int lid = disnpmap.LID(walldiscretization_->dof(0, actnode, idx));
        if (lid < 0)
          FOUR_C_THROW("You tried to test {} on nonexistent dof {} on node {}", quantity.c_str(),
              idx, actnode->id());
        actresult = (*disp)[lid];
      }
    }
    else
      FOUR_C_THROW("result check failed with unknown quantity '{}'!", quantity.c_str());

    // compare values
    const int err = compare_values(actresult, "NODE", container);
    nerr += err;
    test_count++;
  }
}

void PARTICLEWALL::WallResultTest::test_special(
    const Core::IO::InputParameterContainer& container, int& nerr, int& test_count)
{
  // check results only for processor 0
  if (Core::Communication::my_mpi_rank(walldiscretization_->get_comm()) != 0) return;

  // extract and check discretization name
  std::string dis = container.get<std::string>("DIS");
  if (dis != walldiscretization_->name()) return;

  // extract test quantity
  std::string quantity = container.get<std::string>("QUANTITY");

  // init actual result
  double actresult = 0.0;

  // number of total wall elements
  if (quantity == "nwalleles") actresult = walldiscretization_->num_global_elements();
  // number of total wall nodes
  else if (quantity == "nwallnodes")
    actresult = walldiscretization_->num_global_nodes();
  else
    FOUR_C_THROW("result check failed with unknown quantity '{}'!", quantity.c_str());

  // compare values
  const int err = compare_values(actresult, "SPECIAL", container);
  nerr += err;
  test_count++;
}

FOUR_C_NAMESPACE_CLOSE
