// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_cut_parallel.hpp"

#include "4C_comm_exporter.hpp"
#include "4C_comm_mpi_utils.hpp"
#include "4C_comm_pack_buffer.hpp"
#include "4C_cut_output.hpp"
#include "4C_cut_volumecell.hpp"
#include "4C_fem_discretization.hpp"

#include <Teuchos_Time.hpp>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN


/*------------------------------------------------------------------------------------------------*
 * basic CUT parallel constructor                                                    schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
Cut::Parallel::Parallel(const std::shared_ptr<Core::FE::Discretization>& discret, Cut::Mesh& mesh,
    Cut::ParentIntersection& parentintersection)
    : discret_(discret),
      myrank_(Core::Communication::my_mpi_rank(discret_->get_comm())),
      numproc_(Core::Communication::num_mpi_ranks(discret_->get_comm())),
      mesh_(mesh),
      parentintersection_(parentintersection)
{
  return;
}  // end constructor


/*------------------------------------------------------------------------------------------------*
 * communicate node positions                                                        schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void Cut::Parallel::communicate_node_positions()
{
  // wait for all processors before the Communication starts
  Core::Communication::barrier(discret_->get_comm());

  TEUCHOS_FUNC_TIME_MONITOR("Cut --- 5/6 --- cut_positions_dofsets --- communicate_node_positions");

  //  if(myrank_==0) std::cout << "\n\t ... communicate_node_positions" << std::flush;
  //
  //  const double t_start = Teuchos::Time::wallTime();

  //----------------------------------------------------------


  int counter = 0;  // loop counter to avoid infinite loops

  while (true)  // loop as long as not all procs finished with decided node positions for all their
                // nodes
  {
    counter += 1;

    // avoid infinite loops, after numproc rounds all procs should know their node positions
    if (counter > numproc_ + 10)
    {
      break;
      FOUR_C_THROW(
          "number of rounds for exchanging data in communicate_node_positions > numproc_+10");
    }

    // check if there are undecided node positions on this proc
    bool undecided_nodes = mesh_.check_for_undecided_node_positions(
        curr_undecided_node_pos_, curr_undecided_node_pos_shadow_);

    // no undecided node positions -> procDone=true
    // undecided node positions -> procDone=false
    bool procDone = !undecided_nodes;

    // export if current proc is done or not
    // explanation: export is not finished when another proc (especially previous proc) is not done
    // * after one round all procs have status procDone=false if at least one proc has undecided
    // node positions
    // * after one round all procs have status procDone=true  only if all procs have not undecided
    // node positions
    export_communication_finished(procDone);

    // check: procDone == 1 just if all procs have finished
    if (procDone)
      break;
    else
    {
      // perform a new Robin round to gather data from other procs
      // (send from current proc to next proc and receive info from proc before, numproc times)
      // fill the current maps with information (node positions) from myproc
      export_node_position_data();

      //-----------------------------------------------------------------------
      // ... now (back to the original proc) some of the ordered data could have been obtained from
      // other procs try do distribute the information (node positions)
      distribute_my_received_node_position_data();
    }
  }  // end while loop

#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (myrank_ == 0)
  {
    std::cout << "number of round Robin loops to check finished procs:\t" << counter << std::endl;
    std::cout << "number of round Robin loops to exchange node positions:\t" << counter - 1
              << std::endl;
  }
#endif

  Core::Communication::barrier(discret_->get_comm());

  //----------------------------------------------------------

  //  const double t_diff = Teuchos::Time::wallTime()-t_start;
  //  if ( myrank_ == 0 )
  //  {
  //    std::cout << " ... Success (" << t_diff  <<  " secs)";
  //  }


  return;
}



/*------------------------------------------------------------------------------------------------*
 * export data whether proc has finished to set node positions                       schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void Cut::Parallel::export_communication_finished(bool& procDone)
{
  // Initialization
  int dest = myrank_ + 1;  // destination proc (the "next" one)
  if (myrank_ == (numproc_ - 1)) dest = 0;

  int source = myrank_ - 1;  // source proc (the "last" one)
  if (myrank_ == 0) source = numproc_ - 1;


  /*-------------------------------------------*
   * first part: send procfinished in order to *
   * check whether all procs have finished     *
   *-------------------------------------------*/
  for (int iproc = 0; iproc < numproc_ - 1;
      iproc++)  // proc obtains information from numproc_-1 other processors
  {
    Core::Communication::PackBuffer dataSend;

    add_to_pack(dataSend, procDone);

    std::vector<char> dataRecv;
    send_data(dataSend, dest, source, dataRecv);

    // pointer to current position of group of cells in global std::string (counts bytes)
    bool allProcsDone = 0;

    Core::Communication::UnpackBuffer data_recv_buffer(dataRecv);
    // unpack received data
    extract_from_pack(data_recv_buffer, allProcsDone);

    // if the received information is allProcsDone==false, then set the current proc also to
    // procDone=false within the next round-iteration the next proc is also set to procDone=false
    if ((bool)allProcsDone == false) procDone = false;

    // processors wait for each other
    Core::Communication::barrier(discret_->get_comm());
  }


  return;
}  // end exportFinishedData



/*------------------------------------------------------------------------------------------------*
 * export position data to neighbor proc and receive data from previous proc         schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void Cut::Parallel::export_node_position_data()
{
  // destination proc (the "next" one)
  int dest = myrank_ + 1;
  if (myrank_ == (numproc_ - 1)) dest = 0;

  // source proc (the "last" one)
  int source = myrank_ - 1;
  if (myrank_ == 0) source = numproc_ - 1;


  // loop over processors
  for (int procid = 0; procid < numproc_; procid++)
  {
    //---------------------------------------------------------------------------------------------------------------
    //--------------------- send data to next proc and receive from previous proc
    //-----------------------------------
    //---------------------------------------------------------------------------------------------------------------
    // send current node position map to next proc and receive a new map from previous proc
    {
      // copy from sorted_vector to std::vector
      // sorting wont be changed during copying
      std::map<std::vector<int>, int> tmp_curr_undecidedNodePos_shadow;
      tmp_curr_undecidedNodePos_shadow.clear();

      //--------------------
      // copy from std::map<plain_int_set, int> -> std::map<std::vector<int>, int>
      for (std::map<plain_int_set, int>::iterator it = curr_undecided_node_pos_shadow_.begin();
          it != curr_undecided_node_pos_shadow_.end(); it++)
      {
        std::vector<int> tmp_vec;
        tmp_vec.clear();

        for (plain_int_set::const_iterator vec_it = (it->first).begin();
            vec_it != (it->first).end(); vec_it++)
        {
          tmp_vec.push_back(*vec_it);
        }

        tmp_curr_undecidedNodePos_shadow.insert(
            std::pair<std::vector<int>, int>(tmp_vec, it->second));
      }
      //--------------------


      Core::Communication::PackBuffer dataSend;  // data to be sent

      add_to_pack(dataSend, curr_undecided_node_pos_);
      add_to_pack(dataSend, tmp_curr_undecidedNodePos_shadow);

      std::vector<char> dataRecv;
      send_data(dataSend, dest, source, dataRecv);

      // clear vector that should be filled
      curr_undecided_node_pos_.clear();
      curr_undecided_node_pos_shadow_.clear();
      tmp_curr_undecidedNodePos_shadow.clear();

      // unpack received data
      Core::Communication::UnpackBuffer data_recv_buffer(dataRecv);
      while (!data_recv_buffer.at_end())
      {
        extract_from_pack(data_recv_buffer, curr_undecided_node_pos_);
        extract_from_pack(data_recv_buffer, tmp_curr_undecidedNodePos_shadow);

        //--------------------
        // copy from std::map<std::vector<int>, int> -> std::map<plain_int_set, int>
        for (std::map<std::vector<int>, int>::iterator it =
                 tmp_curr_undecidedNodePos_shadow.begin();
            it != tmp_curr_undecidedNodePos_shadow.end(); it++)
        {
          Cut::plain_int_set tmp_set;
          for (std::vector<int>::const_iterator vec_it = (it->first).begin();
              vec_it != (it->first).end(); vec_it++)
          {
            tmp_set.insert(*vec_it);
          }

          curr_undecided_node_pos_shadow_.insert(
              std::pair<plain_int_set, int>(tmp_set, it->second));
        }
        //--------------------
      }

      Core::Communication::barrier(discret_->get_comm());  // processors wait for each other
    }


    //---------------------------------------------------------------------------------------------------------------
    //---------------------------------- fill maps with data from current proc
    //-------------------------------------
    //---------------------------------------------------------------------------------------------------------------
    // find node positions for received nodes and set position if possible
    for (std::map<int, int>::iterator it = curr_undecided_node_pos_.begin();
        it != curr_undecided_node_pos_.end(); it++)
    {
      int nid = it->first;
      int pos = it->second;

      // find the node on current proc
      Node* n = mesh_.get_node(nid);

      if (n != nullptr)  // node on this proc found
      {
        Point* p = n->point();

        Point::PointPosition my_pos = p->position();


        if ((pos == Point::undecided) and (my_pos != Point::undecided))
        {
          // set the new position
          it->second = my_pos;
        }
        // both position are undecided
        else if ((pos == Point::undecided) and (my_pos == Point::undecided))
        {
          // no gain in information
        }
        else if ((pos != Point::undecided) and (my_pos != Point::undecided) and (pos != my_pos))
        {  // only an additional check for not unique node positions
          FOUR_C_THROW(
              "current position on myproc and received position for node {} are set, but not equal",
              nid);
        }
      }

    }  // end iteration map

    //---------------------------------------------------------------------------------------------------------------
    // do the same for shadow nodes in case of quadratic elements
    // find node positions for received shadow! nodes and set position if possible
    for (std::map<plain_int_set, int>::iterator it = curr_undecided_node_pos_shadow_.begin();
        it != curr_undecided_node_pos_shadow_.end(); it++)
    {
      const plain_int_set& nids = it->first;
      int pos = it->second;

      // find the shadow node on current proc (uniquely identified with other nodes of element
      // boundary or element itself)
      Node* n = mesh_.get_node(nids);

      if (n != nullptr)  // node on this proc found
      {
        Point* p = n->point();

        Point::PointPosition my_pos = p->position();


        if ((pos == Point::undecided) and (my_pos != Point::undecided))
        {
          // set the new position
          it->second = my_pos;
        }
        // both position are undecided
        else if ((pos == Point::undecided) and (my_pos == Point::undecided))
        {
          // no gain in information
        }
        else if ((pos != Point::undecided) and (my_pos != Point::undecided) and (pos != my_pos))
        {  // only an additional check for not unique node positions
          FOUR_C_THROW(
              "current position on myproc and received position for shadow node are set, but not "
              "equal");
        }
      }

    }  // end iteration map

    //---------------------------------------------------------------------------------------------------------------


  }  // end loop over procs


}  // end export_node_position_data



/*------------------------------------------------------------------------------------------------*
 * distribute received node positions on my processor                                schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void Cut::Parallel::distribute_my_received_node_position_data()
{
  // distribute the received data for myproc
  for (std::map<int, int>::iterator it = curr_undecided_node_pos_.begin();
      it != curr_undecided_node_pos_.end(); it++)
  {
    int nid = it->first;
    Point::PointPosition received_pos = (Point::PointPosition)it->second;

    // find the node on current proc (remark: this node has to be on this proc, it's the original
    // source proc)
    Node* n = mesh_.get_node(nid);

    // set the node position for the node and distribute it to facets, vcs ...
    if (n != nullptr)
    {
      set_position_for_node(n, received_pos);
    }
    else
      FOUR_C_THROW(
          " this node should be available on its original proc that asked other procs for "
          "nodepositions");
  }

  //--------------------------------------------------------------------------------------------------
  // do the same for shadow nodes in case of quadratic elements

  // distribute the received data for myproc
  for (std::map<plain_int_set, int>::iterator it = curr_undecided_node_pos_shadow_.begin();
      it != curr_undecided_node_pos_shadow_.end(); it++)
  {
    const plain_int_set& nids = it->first;
    Point::PointPosition received_pos = (Point::PointPosition)it->second;

    // find the shadow node on current proc (remark: this node has to be on this proc, it's the
    // original source proc)
    Node* n = mesh_.get_node(nids);

    // set the node position for the node and distribute it to facets, vcs ...
    if (n != nullptr)
    {
      set_position_for_node(n, received_pos);
    }
    else
      FOUR_C_THROW(
          " this shadow node should be available on its original proc that asked other procs for "
          "nodepositions");
  }

  return;
}

/*------------------------------------------------------------------------------------------------*
 * set received node positions for node and distribute it to facets, vcs ...         schott 05/14 *
 *------------------------------------------------------------------------------------------------*/
void Cut::Parallel::set_position_for_node(const Node* n, const Point::PointPosition& pos)
{
  Point* p = n->point();

  // set the new position for this point and distribute the information via facets and volumecells
  if (pos != Point::undecided)
  {
    // std::cout << "reset the position for node " << nid << std::endl;
    p->position(pos);
  }

  if (pos == Point::outside or pos == Point::inside)
  {
    // The nodal position is already known. Set it to my facets. If the
    // facets are already set, this will not have much effect anyway. But on
    // multiple cuts we avoid unset facets this way.
    const plain_facet_set& facets = p->facets();
    for (plain_facet_set::const_iterator i = facets.begin(); i != facets.end(); ++i)
    {
      Facet* f = *i;
      f->position(pos);
    }
  }

  return;
}


/*------------------------------------------------------------------------------------------------*
 * communicate the node dofset number for single volumecells                         schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void Cut::Parallel::communicate_node_dof_set_numbers(bool include_inner)
{
  // wait for all processors before the Communication starts
  Core::Communication::barrier(discret_->get_comm());

  TEUCHOS_FUNC_TIME_MONITOR(
      "Cut --- 5/6 --- cut_positions_dofsets --- communicate_node_dof_set_numbers");

  //  if(myrank_==0) std::cout << "\n\t ... communicate_node_dof_set_numbers" << std::flush;
  //
  //  const double t_start = Teuchos::Time::wallTime();

  dof_set_data_.clear();

  // check if there are missing dofset for volumecells on this proc
  parentintersection_.fill_parallel_dof_set_data(dof_set_data_, *discret_, include_inner);


  // perform just one Robin round to gather data from other procs
  // (send from current proc to next proc and receive info from proc before)
  // fill the current maps with information (dofset number for vc and the current row node) from
  // myproc
  export_dof_set_data(include_inner);

  //-----------------------------------------------------------------------
  // ... now (back to the original proc) all the ordered data should have been obtained from other
  // procs set the information
  distribute_dof_set_data();

  Core::Communication::barrier(discret_->get_comm());

  //  const double t_diff = Teuchos::Time::wallTime()-t_start;
  //  if ( myrank_ == 0 )
  //  {
  //    std::cout << " ... Success (" << t_diff  <<  " secs)\n";
  //  }


  return;
}


/*------------------------------------------------------------------------------------------------*
 * export dofset data to neighbor proc and receive data from previous proc           schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void Cut::Parallel::export_dof_set_data(bool include_inner)
{
  //  bool include_inner = false;

  // destination proc (the "next" one)
  int dest = myrank_ + 1;
  if (myrank_ == (numproc_ - 1)) dest = 0;

  // source proc (the "last" one)
  int source = myrank_ - 1;
  if (myrank_ == 0) source = numproc_ - 1;


  // loop over processors
  for (int procid = 0; procid < numproc_; procid++)
  {
    //---------------------------------------------------------------------------------------------------------------
    //--------------------- send data to next proc and receive from previous proc
    //-----------------------------------
    //---------------------------------------------------------------------------------------------------------------
    // send current DofSetData to next proc and receive a new map from previous proc
    {
      Core::Communication::PackBuffer dataSend;  // data to be sent
      for (std::vector<std::shared_ptr<MeshIntersection::DofSetData>>::iterator data =
               dof_set_data_.begin();
          data != dof_set_data_.end(); data++)
      {
        add_to_pack(dataSend, (*data)->set_index_);
        add_to_pack(dataSend, (*data)->inside_cell_);
        pack_points(dataSend, (*data)->cut_points_coords_);
        add_to_pack(dataSend, (*data)->peid_);
        add_to_pack(dataSend, (*data)->node_dofsetnumber_map_);
      }

      std::vector<char> dataRecv;
      send_data(dataSend, dest, source, dataRecv);

      // clear vector that should be filled
      dof_set_data_.clear();

      Core::Communication::UnpackBuffer data_recv_buffer(dataRecv);
      while (!data_recv_buffer.at_end())
      {
        // unpack volumecell
        int set_index = -1;                                         // set index for Volumecell
        bool inside_cell = false;                                   // inside or outside cell
        std::vector<Core::LinAlg::Matrix<3, 1>> cut_points_coords;  // coordinates of cut points
        int peid = -1;                             // parent element id for volume cell
        std::map<int, int> node_dofsetnumber_map;  // map <nid, current dofset number>

        // unpack volumecell data
        extract_from_pack(data_recv_buffer, set_index);
        extract_from_pack(data_recv_buffer, inside_cell);
        unpack_points(data_recv_buffer, cut_points_coords);
        extract_from_pack(data_recv_buffer, peid);
        extract_from_pack(data_recv_buffer, node_dofsetnumber_map);

        // create a new dofSetData object with unpacked data
        dof_set_data_.push_back(std::make_shared<Cut::MeshIntersection::DofSetData>(
            set_index, (bool)inside_cell, cut_points_coords, peid, node_dofsetnumber_map));
      }



      Core::Communication::barrier(discret_->get_comm());  // processors wait for each other
    }


    //---------------------------------------------------------------------------------------------------------------
    //---------------------------------- fill maps with data from current proc
    //-------------------------------------
    //---------------------------------------------------------------------------------------------------------------
    for (std::vector<std::shared_ptr<MeshIntersection::DofSetData>>::iterator vc_data =
             dof_set_data_.begin();
        vc_data != dof_set_data_.end(); vc_data++)
    {
      bool find_volumecell = false;  // do we have to identify the received volumecell on myrank?


      // safety check
      if (include_inner == false)
      {
        if ((*vc_data)->inside_cell_ == true)
          FOUR_C_THROW(
              "why did you communicate volumecells with inside position in element {}, where "
              "include_inner is set to false",
              (*vc_data)->peid_);
      }

      std::map<int, int>& node_dofsetnumber_map = (*vc_data)->node_dofsetnumber_map_;


      // check if it necessary to find the received volumecell on myrank
      // required, if at least one node in the received node_dofsetnumber_map is a row node on
      // myrank
      for (std::map<int, int>::iterator node_dofsetnumber_it = node_dofsetnumber_map.begin();
          node_dofsetnumber_it != node_dofsetnumber_map.end(); node_dofsetnumber_it++)
      {
        int nid = node_dofsetnumber_it->first;

        bool haveGlobalNode = discret_->have_global_node(nid);

        if (haveGlobalNode)
        {
          Core::Nodes::Node* node = discret_->g_node(nid);

          if (node->owner() == myrank_)
          {
            find_volumecell = true;
            // std::cout << "find the vc!!!" << std::endl;
            break;  // at least one node found as row node, we have to identify the received
                    // volumecell on myrank
          }
        }
      }

      //----------------------------------------------------------------------
      if (!find_volumecell) continue;  // check the next received dofSetData_
      //----------------------------------------------------------------------

      // parent element Id for current vc data
      int peid = (*vc_data)->peid_;

      // find the volumecell on myrank
      VolumeCell* my_vc = nullptr;
      double tol = PARALLEL_COORD_TOL;
      int step = 0;
      while (my_vc == nullptr && step < 10)
      {
        my_vc = find_volume_cell(*(*vc_data), tol);
        if (!my_vc)
        {
          std::cout << "==| Identification of Volumecells with tolerance = " << tol
                    << " was not successful |==!" << std::endl;
          tol *= 10;
        }
        ++step;
      }

      if (my_vc == nullptr)
        FOUR_C_THROW("no corresponding volumecell for vc in element {} found", peid);


      for (std::map<int, int>::iterator node_dofsetnumber_it = node_dofsetnumber_map.begin();
          node_dofsetnumber_it != node_dofsetnumber_map.end(); node_dofsetnumber_it++)
      {
        // std::cout << "loop the node_dofsetnumber_map" << std::endl;
        int nid = node_dofsetnumber_it->first;
        int curr_dofset_number = node_dofsetnumber_it->second;

        bool haveGlobalNode = discret_->have_global_node(nid);

        // decide if the current proc carries the required information
        if (haveGlobalNode)  // node on this proc available as row or col node
        {
          // std::cout << "in haveGlobalNode for node " << nid << std::endl;
          Core::Nodes::Node* node = discret_->g_node(nid);
          if (node->owner() == myrank_)
          {
            int new_dofset_number = -1;

            // find the local index of the current node w.r.t the element
            int index = -1;

            index = get_dof_set_vec_index(nid, peid);

            // std::cout << "index for node " << nid << " is: " << index << std::endl;

            // get the right dofsetnumber for the node
            if (my_vc != nullptr)
            {
              // std::cout << "in my_vc != nullptr" << std::endl;
              const std::vector<int> nds = my_vc->nodal_dof_set();

              if ((int)(nds.size()) == 0)
              {
                std::cout << "position of found vc is " << my_vc->position() << std::endl;

                std::stringstream str;
                str << "cut_element" << my_vc->parent_element()->id() << "_fail.pos";
                std::ofstream file(str.str().c_str());
                Cut::Output::gmsh_complete_cut_element(file, my_vc->parent_element());
                Cut::Output::gmsh_new_section(file, "MyVC");
                Cut::Output::gmsh_volumecell_dump(file, my_vc);
                Cut::Output::gmsh_new_section(file, "OtherPointsVC", true);
                for (std::vector<Core::LinAlg::Matrix<3, 1>>::iterator rec_it =
                         ((*vc_data)->cut_points_coords_).begin();
                    rec_it != ((*vc_data)->cut_points_coords_).end(); rec_it++)
                {
                  Cut::Output::gmsh_coord_dump(file, (*rec_it), 0);
                }
                Cut::Output::gmsh_end_section(file, true);
                std::cout << "Is the element cut = " << my_vc->parent_element()->is_cut() << "\n";
                std::cout << "VC_DATA is inside: Is the element cut = " << (*vc_data)->inside_cell_
                          << "\n";
                FOUR_C_THROW("the nds-vector of volume cell in element {} on proc {} has size {}",
                    my_vc->parent_element()->id(), myrank_, (int)(nds.size()));
              }
              if (index >= (int)(nds.size()))
              {
                FOUR_C_THROW(
                    " index {} exceeds the nds vector of my vc with size {} in element {} on proc "
                    "{}",
                    index, nds.size(), my_vc->parent_element()->id(), myrank_);
              }

              new_dofset_number = nds[index];

              if (new_dofset_number == -1)
                FOUR_C_THROW("the new dofset number for node {} is not valid", nid);
            }
            else
              FOUR_C_THROW("Volumecell Pointer is nullptr");



            // set the new dofset-number
            if (curr_dofset_number == -1)
            {
              node_dofsetnumber_it->second = new_dofset_number;
            }
            else
              FOUR_C_THROW(
                  "dofset for this node for this volumecell already set by another proc, this "
                  "should be done just by the row-node proc");


          }  // end if myrank_
        }  // have global node on this proc
      }  // end loop node_dofsetnumber_it ( some elemental nodes for the current received
         // volumecell)
    }  // end loop dofSetData_ (data for volumecells)


    //    std::cout << "replaced data: " << std::endl;
    //    print_dof_set_data();


    //---------------------------------------------------------------------------------------------------------------

  }  // end loop over procs

  return;
}  // end export_node_position_data



/*------------------------------------------------------------------------------------------------*
 * distribute received dofset number for nodes of volumecell on my processor         schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void Cut::Parallel::distribute_dof_set_data()
{
  // back to the original proc

  for (std::vector<std::shared_ptr<MeshIntersection::DofSetData>>::iterator data =
           dof_set_data_.begin();
      data != dof_set_data_.end(); data++)
  {
    // set data in first volumecell of set with setindex

    // get the Volumecell
    // find volumecell sets and non-row nodes for that dofset numbers has to be communicated
    // parallel the communication is done element wise for all its sets of volumecells when there is
    // a non-row node in this element

    int peid = (*data)->peid_;

    Cut::ElementHandle* e = parentintersection_.get_element(peid);

    std::vector<int> nds;

    if (e != nullptr)
    {
      if ((*data)->inside_cell_)
      {
        const std::vector<plain_volumecell_set>& ele_vc_sets_inside = e->get_vc_sets_inside();
        std::vector<std::vector<int>>& nodaldofset_vc_sets_inside =
            e->get_nodal_dof_set_vc_sets_inside();

        replace_nds_vectors(e, ele_vc_sets_inside, nodaldofset_vc_sets_inside, (*data)->set_index_,
            (*data)->node_dofsetnumber_map_);
      }
      else
      {
        const std::vector<plain_volumecell_set>& ele_vc_sets_outside = e->get_vc_sets_outside();
        std::vector<std::vector<int>>& nodaldofset_vc_sets_outside =
            e->get_nodal_dof_set_vc_sets_outside();

        replace_nds_vectors(e, ele_vc_sets_outside, nodaldofset_vc_sets_outside,
            (*data)->set_index_, (*data)->node_dofsetnumber_map_);
      }
    }
    else
    {
      FOUR_C_THROW(
          "there must be an elementhandle for element {} on my original proc {}", peid, myrank_);
    }


    //... and set data in elementhandle


    // check if all nodes are filled with a valid dofset number
    for (std::map<int, int>::iterator dofnumber_map = (*data)->node_dofsetnumber_map_.begin();
        dofnumber_map != (*data)->node_dofsetnumber_map_.end(); dofnumber_map++)
    {
      // safety check if setting dofset for data was successful
      if (dofnumber_map->second == -1)
      {
        FOUR_C_THROW("dofset number for node {} for vc in element {} could not be determined",
            dofnumber_map->first, (*data)->peid_);
      }
    }
  }


  return;
}

/*------------------------------------------------------------------------------------------------*
 * find the volumecell on myrank for which we received data stored in vc_data        schott 10/12 *
 *------------------------------------------------------------------------------------------------*/
Cut::VolumeCell* Cut::Parallel::find_volume_cell(
    MeshIntersection::DofSetData&
        vc_data,  ///< volumecell data which have to be identified on myrank
    double tol)
{
  bool vc_found = true;

  VolumeCell* my_vc = nullptr;

  ElementHandle* pele = parentintersection_.get_element(vc_data.peid_);

  if (pele == nullptr)
    FOUR_C_THROW("element with Id {} not found on proc {}", vc_data.peid_, myrank_);

  plain_volumecell_set my_vcs;

  pele->get_volume_cells(my_vcs);

  // find all the received volumecell's points
  for (plain_volumecell_set::iterator c = my_vcs.begin(); c != my_vcs.end(); c++)
  {
    VolumeCell* cell = *c;

    // assume that the corresponding vc on myrank_ has been found
    vc_found = true;

    // get points for my_vc
    plain_point_set my_cut_points;
    const plain_facet_set& facets = cell->facets();

    for (plain_facet_set::const_iterator i = facets.begin(); i != facets.end(); i++)
    {
      Facet* f = *i;
      const std::vector<Point*>& facetpoints = f->points();
      std::copy(facetpoints.begin(), facetpoints.end(),
          std::inserter(my_cut_points, my_cut_points.begin()));
    }


    //-----------------------------------------------------------------------------
    // find received cut_points in my_cut_points

    // compare the size of vecs
    //    if(my_cut_points.size() != (vc_data.cut_points_coords_).size())
    //    {
    //
    //      //std::cout << "my_cut_points.size() != (vc_data->cut_points_coords_).size()" <<
    //      std::endl;
    //      // no identical number of points found!
    //      // this volumecell is not a candidate
    //      vc_found = false;
    //    }
    //    else
    {
      // std::cout << "my_cut_points.size() == (vc_data->cut_points_coords_).size()" << std::endl;
      // identical number of points found! this vc is a candidate: check the points"

      // brute force search for identical point coords
      for (plain_point_set::iterator my_it = my_cut_points.begin(); my_it != my_cut_points.end();
          my_it++)
      {
        bool point_found = true;

        for (std::vector<Core::LinAlg::Matrix<3, 1>>::iterator rec_it =
                 (vc_data.cut_points_coords_).begin();
            rec_it != (vc_data.cut_points_coords_).end(); rec_it++)
        {
          point_found = true;

          // compare the coordinates
          for (int i = 0; i < 3; i++)
          {
            if (fabs((*my_it)->x()[i] - (*rec_it)(i)) > tol)
            {
              point_found = false;
              break;  // stop the loop over coordinates
            }
          }

          if (point_found == true)
            break;  // stop the inner! loop over received vc's points if the point is already found

        }  // inner loop over received vc's points

        if (point_found == false)
        {
          vc_found = false;
          break;  // do not search for the next points, it is the wrong volumecell
        }
      }  // end loop over my_points

      // now check the other way around
      if (vc_found)
      {
        // brute force search for identical point coords
        for (std::vector<Core::LinAlg::Matrix<3, 1>>::iterator rec_it =
                 (vc_data.cut_points_coords_).begin();
            rec_it != (vc_data.cut_points_coords_).end(); rec_it++)
        {
          bool point_found = true;
          for (plain_point_set::iterator my_it = my_cut_points.begin();
              my_it != my_cut_points.end(); my_it++)
          {
            point_found = true;

            // compare the coordinates
            for (int i = 0; i < 3; i++)
            {
              if (fabs((*my_it)->x()[i] - (*rec_it)(i)) > tol)
              {
                point_found = false;
                break;  // stop the loop over coordinates
              }
            }

            if (point_found == true)
              break;  // stop the inner! loop over received vc's points if the point is already
                      // found

          }  // inner loop over received vc's points

          if (point_found == false)
          {
            vc_found = false;
            break;  // do not search for the next points, it is the wrong volumecell
          }
        }  // end loop over my_points
      }
    }  // end else

    // additional check which we can make for sure (If there are no DOF Sets this cannot be the
    // correct Volumecell)
    if (vc_found == true)
      if (!cell->nodal_dof_set().size()) vc_found = false;


    if (vc_found == true)
    {
      my_vc = cell;
      // std::cout << "volumecell has been found!!! Happy :-)" << std::endl;
      break;  // stop the loop over volumecells
    }

    //-----------------------------------------------------------------------------
  }  // it over volumecells

  //  if (vc_found == false)  //Do a fallbackstrategy
  //  {
  //    //identify by center
  //    std::cout << "==| Do center idendification |== " << std::endl;
  //    Core::LinAlg::Matrix<3,1> vc_center(true);
  //    for(std::vector<Core::LinAlg::Matrix<3,1> >::iterator
  //    rec_it=(vc_data.cut_points_coords_).begin();
  //       rec_it != (vc_data.cut_points_coords_).end();
  //       rec_it++)
  //    {
  //     vc_center.update(1.0,(*rec_it),1.0);
  //    }
  //    vc_center.scale(1.0/(vc_data.cut_points_coords_).size());
  //
  //    for(plain_volumecell_set::iterator c = my_vcs.begin(); c!=my_vcs.end(); c++)
  //    {
  //     vc_found=true;
  //
  //     VolumeCell* cell = *c;
  //     // get points for my_vc
  //     plain_point_set my_cut_points;
  //     const plain_facet_set & facets = cell->Facets();
  //
  //     for(plain_facet_set::const_iterator i=facets.begin(); i!=facets.end(); i++)
  //     {
  //       Facet* f = *i;
  //       const std::vector<Point*> & facetpoints = f->Points();
  //       std::copy(facetpoints.begin(), facetpoints.end(), std::inserter(my_cut_points,
  //       my_cut_points.begin()));
  //     }
  //     Core::LinAlg::Matrix<3,1> my_center(true);
  //     for(plain_point_set::iterator my_it=my_cut_points.begin(); my_it!=my_cut_points.end();
  //     my_it++)
  //     {
  //       Core::LinAlg::Matrix<3,1> coord((*my_it)->X(),true);
  //       my_center.update(1.0,coord,1.0);
  //     }
  //     my_center.scale(1.0/(my_cut_points.size()));
  //
  //     for(int i=0; i<3; i++)
  //     {
  //       if(fabs(my_center(i) - vc_center(i)) > PARALLEL_COORD_TOL)
  //       {
  //         vc_found=false;
  //         break; // do not search for the next points, it is the wrong volumecell
  //       }
  //     }
  //     my_vc = cell;
  //    }
  //    if (vc_found == true)
  //     std::cout << "==| identified vc by center |== " << std::endl;
  //   }

  if (vc_found == false)
  {
    std::stringstream str;
    str << "NIVC_" << myrank_ << ".pos";
    std::ofstream file(str.str());

    Cut::Output::gmsh_new_section(file, "OtherPoints");
    for (std::vector<Core::LinAlg::Matrix<3, 1>>::iterator rec_it =
             (vc_data.cut_points_coords_).begin();
        rec_it != (vc_data.cut_points_coords_).end(); rec_it++)
    {
      Cut::Output::gmsh_coord_dump(file, (*rec_it), 0);
    }
    for (plain_volumecell_set::iterator c = my_vcs.begin(); c != my_vcs.end(); c++)
    {
      Cut::Output::gmsh_new_section(file, "MyPoints", true);
      VolumeCell* cell = *c;
      // get points for my_vc
      plain_point_set my_cut_points;
      const plain_facet_set& facets = cell->facets();

      for (plain_facet_set::const_iterator i = facets.begin(); i != facets.end(); i++)
      {
        Facet* f = *i;
        const std::vector<Point*>& facetpoints = f->points();
        std::copy(facetpoints.begin(), facetpoints.end(),
            std::inserter(my_cut_points, my_cut_points.begin()));
      }
      for (plain_point_set::iterator my_it = my_cut_points.begin(); my_it != my_cut_points.end();
          my_it++)
      {
        Core::LinAlg::Matrix<3, 1> coord((*my_it)->x(), true);
        Cut::Output::gmsh_coord_dump(file, coord, 0);
      }
    }

    for (plain_volumecell_set::iterator c = my_vcs.begin(); c != my_vcs.end(); c++)
    {
      VolumeCell* cell = *c;
      const plain_facet_set& facete = cell->facets();
      Cut::Output::gmsh_new_section(file, "VolumeCell", true);
      Cut::Output::gmsh_volumecell_dump(file, cell, "lines", true);
      for (plain_facet_set::const_iterator it = facete.begin(); it != facete.end(); ++it)
      {
        Cut::Output::gmsh_new_section(file, "Facet", true);
        Cut::Output::gmsh_facet_dump(file, *it, "sides", true);
      }
    }
    Cut::Output::gmsh_end_section(file, true);
  }

  return my_vc;
}



void Cut::Parallel::replace_nds_vectors(ElementHandle* e,
    const std::vector<plain_volumecell_set>& ele_vc_sets,
    std::vector<std::vector<int>>& nodaldofset_vc_sets, int set_index,
    std::map<int, int>& node_dofsetnumber_map)
{
  // get the original set
  plain_volumecell_set cells = ele_vc_sets[set_index];

  // get the old nds vector for all cells in current set, represented by the first Volumecell
  const std::vector<int> nds_old = (*cells.begin())->nodal_dof_set();

  // get vector of nids in current element
  const std::vector<Cut::Node*> nodes = e->nodes();

  // create the new dofset Vector
  std::vector<int> nds_new;

  // fill the new nds-vector
  for (int i = 0; i < (int)nds_old.size(); i++)
  {
    if (nds_old[i] == -1)  // entry that has to be replaced with received data
    {
      int nid = nodes[i]->id();  // key for received node_dofsetnumber_map
      nds_new.push_back(node_dofsetnumber_map.find(nid)->second);
    }
    else  // set original dofset number
    {
      nds_new.push_back(nds_old[i]);
    }
  }

  // set the new nds vector for all volumecells in the current set
  for (plain_volumecell_set::iterator c = cells.begin(); c != cells.end(); c++)
  {
    (*c)->set_nodal_dof_set(nds_new);
  }

  nodaldofset_vc_sets[set_index] = nds_new;

  return;
}


/*------------------------------------------------------------------------------------------------*
 * packing a point for parallel communication only with the basic point data                      *
 * without an underlying discretization fitting to the node's new prozessor          schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void Cut::Parallel::pack_points(Core::Communication::PackBuffer& dataSend,
    std::vector<Core::LinAlg::Matrix<3, 1>>& points_coords) const
{
  add_to_pack(dataSend, points_coords);
}



/*------------------------------------------------------------------------------------------------*
 * unpacking a point for parallel communication only with the basic point data                    *
 * without an underlying discretization fitting to the node's new prozessor          schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void Cut::Parallel::unpack_points(Core::Communication::UnpackBuffer& buffer,
    std::vector<Core::LinAlg::Matrix<3, 1>>& points_coords) const
{
  extract_from_pack(buffer, points_coords);
}  // end function unpackNodes



/*------------------------------------------------------------------------------------------------*
 * basic function sending data to dest and receiving data from source                schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void Cut::Parallel::send_data(Core::Communication::PackBuffer& dataSend, int& dest, int& source,
    std::vector<char>& dataRecv) const
{
  std::vector<int> lengthSend(1, 0);
  lengthSend[0] = dataSend().size();
  int size_one = 1;

#ifdef FOUR_C_ENABLE_ASSERTIONS
  std::cout << "--- sending " << lengthSend[0] << " bytes: from proc " << myrank_ << " to proc "
            << dest << std::endl;
#endif

  // exporter for sending
  Core::Communication::Exporter exporter(discret_->get_comm());

  // send length of the data to be received ...
  MPI_Request req_length_data;
  int length_tag = 0;
  exporter.i_send(myrank_, dest, lengthSend.data(), size_one, length_tag, req_length_data);

  // ... and receive length
  std::vector<int> lengthRecv(1, 0);
  exporter.receive(source, length_tag, lengthRecv, size_one);
  exporter.wait(req_length_data);

  //-----------------------------------------

  // send actual data ...
  int data_tag = 4;
  MPI_Request req_data;
  exporter.i_send(myrank_, dest, dataSend().data(), lengthSend[0], data_tag, req_data);

  // ... and receive data
  dataRecv.clear();
  dataRecv.resize(lengthRecv[0]);
  exporter.receive_any(source, data_tag, dataRecv, lengthRecv[0]);
  exporter.wait(req_data);

#ifdef FOUR_C_ENABLE_ASSERTIONS
  std::cout << "--- receiving " << lengthRecv[0] << " bytes: to proc " << myrank_ << " from proc "
            << source << std::endl;
#endif
}  // end send_data



/*------------------------------------------------------------------------------------------------*
 * print current stored dofSetData_                                                  schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void Cut::Parallel::print_dof_set_data()
{
  for (std::vector<std::shared_ptr<MeshIntersection::DofSetData>>::iterator i =
           dof_set_data_.begin();
      i != dof_set_data_.end(); i++)
  {
    (*i)->print();
  }
}

/*------------------------------------------------------------------------------------------------*
 * get the index of nid in the vector of elements (eid) node Ids                     schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
int Cut::Parallel::get_dof_set_vec_index(int nid, int eid)
{
  Core::Elements::Element* ele = discret_->g_element(eid);

  if (ele == nullptr) FOUR_C_THROW("element {} not available on proc {}", eid, myrank_);

  int numnode = ele->num_node();

  const int* ele_n_ids = ele->node_ids();

  int index = 0;

  // find the local node Index for received nid w.r.t the element
  for (; index < numnode; ++index)
  {
    if (ele_n_ids[index] == nid) return index;
  }

  if (index > numnode)
  {
    FOUR_C_THROW("nid {} not found in nodes for element {}!", nid, eid);
    return -1;
  }

  return -1;
}

FOUR_C_NAMESPACE_CLOSE
