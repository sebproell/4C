// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_comm_exporter.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_fluid_turbulence_statistics_bfs.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

#include <fstream>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*!
  \brief Standard Constructor (public)

    o Create sets for lines

  o Allocate distributed vector for squares

*/
/*----------------------------------------------------------------------*/
FLD::TurbulenceStatisticsBfs::TurbulenceStatisticsBfs(
    std::shared_ptr<Core::FE::Discretization> actdis, Teuchos::ParameterList& params,
    const std::string& statistics_outfilename, const std::string& geotype)
    : discret_(actdis),
      params_(params),
      geotype_(TurbulenceStatisticsBfs::none),
      inflowchannel_(params_.sublist("TURBULENT INFLOW").get<bool>("TURBULENTINFLOW")),
      inflowmax_(params_.sublist("TURBULENT INFLOW").get<double>("INFLOW_CHA_SIDE", 0.0)),
      statistics_outfilename_(statistics_outfilename)
{
  if (Core::Communication::my_mpi_rank(discret_->get_comm()) == 0)
  {
    std::cout << "This is the turbulence statistics manager of backward-facing step problems:"
              << std::endl;
    std::cout << "It can deal with the following two geometries:" << std::endl;
    std::cout << "- geometry of DNS by Le,Moin,Kim (expansion ratio 1.2) and" << std::endl;
    std::cout << "- geometry of experiment by Kasagi,Matsunaga (expansion ratio 1.5) and"
              << std::endl;
    std::cout << "- geometry of experiment by Vogel, Eaton (expansion ratio 1.25) at Re 28,000."
              << std::endl;
    std::cout << "If additional output in front of the step is required, it has to be set manually "
                 "(look for numx1supplocations_)."
              << std::endl;
  }

  //----------------------------------------------------------------------
  // plausibility check
  int numdim = params_.get<int>("number of velocity degrees of freedom");
  if (numdim != 3) FOUR_C_THROW("Evaluation of turbulence statistics only for 3d flow problems!");

  // type of fluid flow solver: incompressible, Boussinesq approximation, varying density, loma
  const auto physicaltype =
      Teuchos::getIntegralValue<Inpar::FLUID::PhysicalType>(params_, "Physical Type");

  // geometry of bfs
  convert_string_to_geo_type(geotype);

  //----------------------------------------------------------------------
  // allocate some (toggle) vectors
  const Epetra_Map* dofrowmap = discret_->dof_row_map();

  squaredvelnp_ = Core::LinAlg::create_vector(*dofrowmap, true);
  squaredscanp_ = Core::LinAlg::create_vector(*dofrowmap, true);
  invscanp_ = Core::LinAlg::create_vector(*dofrowmap, true);
  squaredinvscanp_ = Core::LinAlg::create_vector(*dofrowmap, true);

  toggleu_ = Core::LinAlg::create_vector(*dofrowmap, true);
  togglev_ = Core::LinAlg::create_vector(*dofrowmap, true);
  togglew_ = Core::LinAlg::create_vector(*dofrowmap, true);
  togglep_ = Core::LinAlg::create_vector(*dofrowmap, true);

  // bounds for extension of flow domain in x2-direction
  x2min_ = +10e+19;
  x2max_ = -10e+19;
  // bounds for extension of flow domain in x3-direction
  x3min_ = +10e+19;
  x3max_ = -10e+19;

  //----------------------------------------------------------------------
  // create sets of coordinates
  //----------------------------------------------------------------------
  x1coordinates_ = std::make_shared<std::vector<double>>();
  x2coordinates_ = std::make_shared<std::vector<double>>();

  // the criterion allows differences in coordinates by 1e-9
  std::set<double, LineSortCriterion> x1avcoords;
  std::set<double, LineSortCriterion> x2avcoords;

  // loop nodes and build sets of lines in x1- and x2-direction
  // accessible on this proc
  // For x1-direction: consider horizontal line at x2=0
  // and assume no change in discretization behind the step
  // For x2-direction: consider vertical line at x1=0
  // and assume no change in discretization behind the step
  for (int i = 0; i < discret_->num_my_row_nodes(); ++i)
  {
    Core::Nodes::Node* node = discret_->l_row_node(i);

    if (inflowchannel_ == true)
    {
      // store also x-coordinate of outflow of inflow channel
      if ((node->x()[1] < 2e-9 && node->x()[1] > -2e-9) && (node->x()[0] > (inflowmax_ - 2e-9)))
        x1avcoords.insert(node->x()[0]);
    }
    else
    {
      if (node->x()[1] < 2e-9 && node->x()[1] > -2e-9) x1avcoords.insert(node->x()[0]);
    }
    if (node->x()[0] < 2e-9 && node->x()[0] > -2e-9) x2avcoords.insert(node->x()[1]);

    // find mins and maxs
    // we do not look for x1min and x1max as they depend
    // on the inflow generation technique
    if (x2min_ > node->x()[1]) x2min_ = node->x()[1];
    if (x2max_ < node->x()[1]) x2max_ = node->x()[1];

    if (x3min_ > node->x()[2]) x3min_ = node->x()[2];
    if (x3max_ < node->x()[2]) x3max_ = node->x()[2];
  }

  // communicate x2mins and x2maxs
  double min2;
  Core::Communication::min_all(&x2min_, &min2, 1, discret_->get_comm());
  x2min_ = min2;

  double max2;
  Core::Communication::max_all(&x2max_, &max2, 1, discret_->get_comm());
  x2max_ = max2;

  // communicate x3mins and x3maxs
  double min3;
  Core::Communication::min_all(&x3min_, &min3, 1, discret_->get_comm());
  x3min_ = min3;

  double max3;
  Core::Communication::max_all(&x3max_, &max3, 1, discret_->get_comm());
  x3max_ = max3;

  //--------------------------------------------------------------------
  // round robin loop to communicate coordinates to all procs
  //--------------------------------------------------------------------
  {
    int myrank = Core::Communication::my_mpi_rank(discret_->get_comm());
    int numprocs = Core::Communication::num_mpi_ranks(discret_->get_comm());

    std::vector<char> sblock;
    std::vector<char> rblock;

    // create an exporter for point to point communication
    Core::Communication::Exporter exporter(discret_->get_comm());

    // first, communicate coordinates in x1-direction
    for (int np = 0; np < numprocs; ++np)
    {
      Core::Communication::PackBuffer data;

      for (std::set<double, LineSortCriterion>::iterator x1line = x1avcoords.begin();
          x1line != x1avcoords.end(); ++x1line)
      {
        add_to_pack(data, *x1line);
      }
      std::swap(sblock, data());

      MPI_Request request;
      int tag = myrank;

      int frompid = myrank;
      int topid = (myrank + 1) % numprocs;

      int length = sblock.size();

      exporter.i_send(frompid, topid, sblock.data(), sblock.size(), tag, request);

      rblock.clear();

      // receive from predecessor
      frompid = (myrank + numprocs - 1) % numprocs;
      exporter.receive_any(frompid, tag, rblock, length);

      if (tag != (myrank + numprocs - 1) % numprocs)
      {
        FOUR_C_THROW("received wrong message (ReceiveAny)");
      }

      exporter.wait(request);

      {
        // for safety
        Core::Communication::barrier(exporter.get_comm());
      }

      //--------------------------------------------------
      // Unpack received block into set of all planes.
      {
        std::vector<double> coordsvec;

        coordsvec.clear();

        Core::Communication::UnpackBuffer buffer(rblock);
        while (!buffer.at_end())
        {
          double onecoord;
          extract_from_pack(buffer, onecoord);
          x1avcoords.insert(onecoord);
        }
      }
    }

    // second, communicate coordinates in x2-direction
    for (int np = 0; np < numprocs; ++np)
    {
      Core::Communication::PackBuffer data;
      for (std::set<double, LineSortCriterion>::iterator x2line = x2avcoords.begin();
          x2line != x2avcoords.end(); ++x2line)
      {
        add_to_pack(data, *x2line);
      }
      std::swap(sblock, data());

      MPI_Request request;
      int tag = myrank;

      int frompid = myrank;
      int topid = (myrank + 1) % numprocs;

      int length = sblock.size();

      exporter.i_send(frompid, topid, sblock.data(), sblock.size(), tag, request);

      rblock.clear();

      // receive from predecessor
      frompid = (myrank + numprocs - 1) % numprocs;
      exporter.receive_any(frompid, tag, rblock, length);

      if (tag != (myrank + numprocs - 1) % numprocs)
      {
        FOUR_C_THROW("received wrong message (ReceiveAny)");
      }

      exporter.wait(request);

      {
        // for safety
        Core::Communication::barrier(exporter.get_comm());
      }

      //--------------------------------------------------
      // Unpack received block into set of all planes.
      {
        std::vector<double> coordsvec;

        coordsvec.clear();

        Core::Communication::UnpackBuffer buffer(rblock);
        while (!buffer.at_end())
        {
          double onecoord;
          extract_from_pack(buffer, onecoord);
          x2avcoords.insert(onecoord);
        }
      }
    }
  }

  //----------------------------------------------------------------------
  // push coordinates in vectors
  //----------------------------------------------------------------------
  {
    x1coordinates_ = std::make_shared<std::vector<double>>();
    x2coordinates_ = std::make_shared<std::vector<double>>();

    for (std::set<double, LineSortCriterion>::iterator coord1 = x1avcoords.begin();
        coord1 != x1avcoords.end(); ++coord1)
    {
      x1coordinates_->push_back(*coord1);
      // std::cout << *coord1 << std::endl;
    }

    for (std::set<double, LineSortCriterion>::iterator coord2 = x2avcoords.begin();
        coord2 != x2avcoords.end(); ++coord2)
    {
      x2coordinates_->push_back(*coord2);
      // std::cout << *coord2 << std::endl;
    }
  }

  //----------------------------------------------------------------------
  // number of coordinates in x1- and x2-direction
  //----------------------------------------------------------------------
  numx1coor_ = x1coordinates_->size();
  numx2coor_ = x2coordinates_->size();

  //----------------------------------------------------------------------
  // number of locations in x1-direction for statistical evaluation
  //----------------------------------------------------------------------
  numx1statlocations_ = 21;
  numx1supplocations_ = 0;

  //----------------------------------------------------------------------
  // define locations in x1-direction for statistical evaluation
  //----------------------------------------------------------------------

  // compute step height
  double h = 0.0;
  if (geotype_ == TurbulenceStatisticsBfs::geometry_LES_flow_with_heating)
  {
    h = (x2max_ - x2min_) / 3.0;
  }
  else if (geotype_ == TurbulenceStatisticsBfs::geometry_DNS_incomp_flow)
  {
    h = (x2max_ - x2min_) / 6.0;
  }
  else if (geotype_ == TurbulenceStatisticsBfs::geometry_EXP_vogel_eaton)
  {
    h = (x2max_ - x2min_) / 5.0;
    numx1statlocations_ = 10;
  }

  if (geotype_ == TurbulenceStatisticsBfs::geometry_EXP_vogel_eaton)
  {
    // locactions given by Vogel&Eaton
    std::array<double, 10> givenpos = {2.2, 3.0, 3.73, 4.47, 5.2, 5.93, 6.67, 7.4, 8.13, 8.87};
    if (numx1statlocations_ != 10) FOUR_C_THROW("wrong number of numx1statlocations_");
    // find locations
    for (int rr = 0; rr < numx1statlocations_; rr++)
    {
      double actpos = givenpos[rr];
      double dist = 30.0 * h;
      double mindist = 30.0 * h;
      int pos = 0;

      for (int ll = 0; ll < numx1coor_; ll++)
      {
        dist = abs(actpos * h - x1coordinates_->at(ll));
        if (dist < mindist)
        {
          mindist = dist;
          pos = ll;
        }
      }

      x1statlocations_(rr) = x1coordinates_->at(pos);
    }
  }
  else
  {
    // find locations x/h=0 ... x/h=20
    for (int rr = 0; rr < numx1statlocations_; rr++)
    {
      int actpos = rr;
      double dist = 30 * h;
      double mindist = 30.0 * h;
      int pos = 0;

      for (int ll = 0; ll < numx1coor_; ll++)
      {
        dist = abs(actpos * h - x1coordinates_->at(ll));
        if (dist < mindist)
        {
          mindist = dist;
          pos = ll;
        }
      }

      x1statlocations_(rr) = x1coordinates_->at(pos);
    }
  }

  // find supplementary locations x/h=-2 and -1
  // remark1: number of supplementary location depends on length of inlet section
  // remark2: as separate channel may be located in front of the step, we have to
  //         use this rather complicated way
  numx1supplocations_ = 2;
  numx1statlocations_ += numx1supplocations_;
  for (int rr = 0; rr < numx1supplocations_; rr++)
  {
    int actpos = rr - numx1supplocations_;
    double dist = 30 * h;
    double mindist = 30.0 * h;
    int pos = 0;

    for (int ll = 0; ll < numx1coor_; ll++)
    {
      dist = abs(actpos * h - x1coordinates_->at(ll));
      if (dist < mindist)
      {
        mindist = dist;
        pos = ll;
      }
    }

    x1supplocations_(rr) = x1coordinates_->at(pos);
    // std::cout << x1supplocations_(rr) << std::endl;
  }

  //----------------------------------------------------------------------
  // define locations in x2-direction for statistical evaluation
  // (lower and upper wall)
  //----------------------------------------------------------------------
  if (geotype_ == TurbulenceStatisticsBfs::geometry_LES_flow_with_heating ||
      geotype_ == TurbulenceStatisticsBfs::geometry_EXP_vogel_eaton)
  {
    // num2statlocations_ also defines number of supplocations
    numx2statlocations_ = 2;

    x2statlocations_(0) = x2min_;
    x2statlocations_(1) = x2max_;
    //----------------------------------------------------------------------
    // define supplementary locations in x2-direction for statistical
    // evaluation of velocity derivative at wall
    // (first nodes off lower and upper wall, respectively)
    //----------------------------------------------------------------------
    x2supplocations_(0) = x2coordinates_->at(1);
    x2supplocations_(1) = x2coordinates_->at(x2coordinates_->size() - 2);
  }
  else if (geotype_ == TurbulenceStatisticsBfs::geometry_DNS_incomp_flow)
  {
    // num2statlocations_ also defines number of supplocations
    numx2statlocations_ = 1;

    x2statlocations_(0) = x2min_;
    x2statlocations_(1) = x2max_;  // not needed here, upper wall is slip wall
    //----------------------------------------------------------------------
    // define supplementary locations in x2-direction for statistical
    // evaluation of velocity derivative at wall
    // (first nodes off lower and upper wall, respectively)
    //----------------------------------------------------------------------
    x2supplocations_(0) = x2coordinates_->at(1);
    x2supplocations_(1) =
        x2coordinates_->at(x2coordinates_->size() - 2);  // not needed here, upper wall is slip wall
  }
  else
    FOUR_C_THROW("Unknown geometry of backward facing step!");

  //----------------------------------------------------------------------
  // allocate arrays for sums of mean values
  //----------------------------------------------------------------------
  // x1-direction
  x1sump_ = std::make_shared<Core::LinAlg::SerialDenseMatrix>();
  x1sump_->reshape(numx2statlocations_, numx1coor_);

  x1sumu_ = std::make_shared<Core::LinAlg::SerialDenseMatrix>();
  x1sumu_->reshape(numx2statlocations_, numx1coor_);

  // the following vectors are only necessary for low-Mach-number flow
  x1sumrho_ = std::make_shared<Core::LinAlg::SerialDenseMatrix>();
  x1sumrho_->reshape(numx2statlocations_, numx1coor_);

  x1sum_t_ = std::make_shared<Core::LinAlg::SerialDenseMatrix>();
  x1sum_t_->reshape(numx2statlocations_, numx1coor_);

  x1sumtauw_ = std::make_shared<Core::LinAlg::SerialDenseMatrix>();
  x1sumtauw_->reshape(numx2statlocations_, numx1coor_);

  // x2-direction
  // first-order moments
  x2sumu_ = std::make_shared<Core::LinAlg::SerialDenseMatrix>();
  x2sumu_->reshape(numx1statlocations_, numx2coor_);

  x2sumv_ = std::make_shared<Core::LinAlg::SerialDenseMatrix>();
  x2sumv_->reshape(numx1statlocations_, numx2coor_);

  x2sumw_ = std::make_shared<Core::LinAlg::SerialDenseMatrix>();
  x2sumw_->reshape(numx1statlocations_, numx2coor_);

  x2sump_ = std::make_shared<Core::LinAlg::SerialDenseMatrix>();
  x2sump_->reshape(numx1statlocations_, numx2coor_);

  // second-order moments
  x2sumsqu_ = std::make_shared<Core::LinAlg::SerialDenseMatrix>();
  x2sumsqu_->reshape(numx1statlocations_, numx2coor_);

  x2sumsqv_ = std::make_shared<Core::LinAlg::SerialDenseMatrix>();
  x2sumsqv_->reshape(numx1statlocations_, numx2coor_);

  x2sumsqw_ = std::make_shared<Core::LinAlg::SerialDenseMatrix>();
  x2sumsqw_->reshape(numx1statlocations_, numx2coor_);

  x2sumsqp_ = std::make_shared<Core::LinAlg::SerialDenseMatrix>();
  x2sumsqp_->reshape(numx1statlocations_, numx2coor_);

  x2sumuv_ = std::make_shared<Core::LinAlg::SerialDenseMatrix>();
  x2sumuv_->reshape(numx1statlocations_, numx2coor_);

  x2sumuw_ = std::make_shared<Core::LinAlg::SerialDenseMatrix>();
  x2sumuw_->reshape(numx1statlocations_, numx2coor_);

  x2sumvw_ = std::make_shared<Core::LinAlg::SerialDenseMatrix>();
  x2sumvw_->reshape(numx1statlocations_, numx2coor_);

  // the following vectors are only necessary for low-Mach-number flow
  // first-order moments
  x2sumrho_ = std::make_shared<Core::LinAlg::SerialDenseMatrix>();
  x2sumrho_->reshape(numx1statlocations_, numx2coor_);

  x2sum_t_ = std::make_shared<Core::LinAlg::SerialDenseMatrix>();
  x2sum_t_->reshape(numx1statlocations_, numx2coor_);

  // second-order moments
  x2sumsqrho_ = std::make_shared<Core::LinAlg::SerialDenseMatrix>();
  x2sumsqrho_->reshape(numx1statlocations_, numx2coor_);

  x2sumsq_t_ = std::make_shared<Core::LinAlg::SerialDenseMatrix>();
  x2sumsq_t_->reshape(numx1statlocations_, numx2coor_);

  x2sumrhou_ = std::make_shared<Core::LinAlg::SerialDenseMatrix>();
  x2sumrhou_->reshape(numx1statlocations_, numx2coor_);

  x2sumu_t_ = std::make_shared<Core::LinAlg::SerialDenseMatrix>();
  x2sumu_t_->reshape(numx1statlocations_, numx2coor_);

  x2sumrhov_ = std::make_shared<Core::LinAlg::SerialDenseMatrix>();
  x2sumrhov_->reshape(numx1statlocations_, numx2coor_);

  x2sumv_t_ = std::make_shared<Core::LinAlg::SerialDenseMatrix>();
  x2sumv_t_->reshape(numx1statlocations_, numx2coor_);

  // set number of samples to zero
  numsamp_ = 0;

  //----------------------------------------------------------------------
  // define homogeneous direction to compute averages of Smagorinsky constant

  Teuchos::ParameterList* modelparams = &(params_.sublist("TURBULENCE MODEL"));
  // check if we want to compute averages of Smagorinsky constant
  if (modelparams->get<std::string>("PHYSICAL_MODEL", "no_model") == "Dynamic_Smagorinsky")
  {
    // store them in parameterlist for access on the element
    modelparams->set<std::shared_ptr<std::vector<double>>>("dir1coords_", x1coordinates_);
    modelparams->set<std::shared_ptr<std::vector<double>>>("dir2coords_", x2coordinates_);
  }

  //----------------------------------------------------------------------
  // initialize output and initially open respective statistics output file

  std::shared_ptr<std::ofstream> log;

  if (Core::Communication::my_mpi_rank(discret_->get_comm()) == 0)
  {
    std::string s(statistics_outfilename_);

    if (physicaltype == Inpar::FLUID::loma)
    {
      s.append(".loma_statistics");

      log = std::make_shared<std::ofstream>(s.c_str(), std::ios::out);
      (*log) << "# Statistics for turbulent variable-density flow over a backward-facing step at "
                "low Mach number (first- and second-order moments)\n\n";
    }
    else
    {
      s.append(".flow_statistics");

      log = std::make_shared<std::ofstream>(s.c_str(), std::ios::out);
      (*log) << "# Statistics for turbulent incompressible flow over a backward-facing step "
                "(first- and second-order moments)\n\n";
    }

    log->flush();
  }

  return;
}  // TurbulenceStatisticsBfs::TurbulenceStatisticsBfs


//----------------------------------------------------------------------
// sampling of velocity/pressure values
//----------------------------------------------------------------------
void FLD::TurbulenceStatisticsBfs::do_time_sample(
    Core::LinAlg::Vector<double>& velnp, Core::LinAlg::Vector<double>& stresses)
{
  // compute squared values of velocity
  squaredvelnp_->multiply(1.0, velnp, velnp, 0.0);

  //----------------------------------------------------------------------
  // increase sample counter
  //----------------------------------------------------------------------
  numsamp_++;

  int x1nodnum = -1;
  //----------------------------------------------------------------------
  // values at lower and upper wall
  //----------------------------------------------------------------------
  for (std::vector<double>::iterator x1line = x1coordinates_->begin();
      x1line != x1coordinates_->end(); ++x1line)
  {
    x1nodnum++;

    for (int x2nodnum = 0; x2nodnum < numx2statlocations_; ++x2nodnum)
    {
      // current x2-coordinate of respective wall
      double x2cwall = x2statlocations_(x2nodnum);

      // current x2-coordinate of supplementary location to respective wall
      double x2csupp = x2supplocations_(x2nodnum);

      // toggle vectors are one in the position of a dof of this node,
      // else 0
      toggleu_->put_scalar(0.0);
      togglep_->put_scalar(0.0);
      // misuse togglev for stresses in u direction
      togglev_->put_scalar(0.0);

      // count the number of nodes in x3-direction contributing to this nodal value
      int countnodes = 0;

      for (int nn = 0; nn < discret_->num_my_row_nodes(); ++nn)
      {
        Core::Nodes::Node* node = discret_->l_row_node(nn);

        // this is the wall node
        if ((node->x()[0] < (*x1line + 2e-9) and node->x()[0] > (*x1line - 2e-9)) and
            (node->x()[1] < (x2cwall + 2e-5) and node->x()[1] > (x2cwall - 2e-5)))
        {
          std::vector<int> dof = discret_->dof(node);
          double one = 1.0;

          togglep_->replace_global_values(1, &one, dof.data() + 3);
          // stresses
          togglev_->replace_global_values(1, &one, dof.data());

          countnodes++;
        }
        // this is the supplementary node
        else if ((node->x()[0] < (*x1line + 2e-9) and node->x()[0] > (*x1line - 2e-9)) and
                 (node->x()[1] < (x2csupp + 2e-5) and node->x()[1] > (x2csupp - 2e-5)))
        {
          std::vector<int> dof = discret_->dof(node);
          double one = 1.0;

          toggleu_->replace_global_values(1, &one, dof.data());
        }
      }

      int countnodesonallprocs = 0;

      Core::Communication::sum_all(&countnodes, &countnodesonallprocs, 1, discret_->get_comm());

      // reduce by 1 due to periodic boundary condition
      countnodesonallprocs -= 1;

      if (countnodesonallprocs)
      {
        //----------------------------------------------------------------------
        // get values for velocity derivative and pressure
        //----------------------------------------------------------------------
        double u;
        velnp.dot(*toggleu_, &u);
        double p;
        velnp.dot(*togglep_, &p);
        // stresses
        double tauw;
        stresses.dot(*togglev_, &tauw);

        //----------------------------------------------------------------------
        // calculate spatial means
        //----------------------------------------------------------------------
        double usm = u / countnodesonallprocs;
        double psm = p / countnodesonallprocs;
        double tauwsm = tauw / countnodesonallprocs;

        //----------------------------------------------------------------------
        // add spatial mean values to statistical sample
        //----------------------------------------------------------------------
        (*x1sumu_)(x2nodnum, x1nodnum) += usm;
        (*x1sump_)(x2nodnum, x1nodnum) += psm;
        (*x1sumtauw_)(x2nodnum, x1nodnum) += tauwsm;
      }
    }
  }

  //----------------------------------------------------------------------
  // loop locations for statistical evaluation in x1-direction
  //----------------------------------------------------------------------
  for (int x1nodnum = 0; x1nodnum < numx1statlocations_; ++x1nodnum)
  {
    // current x1-coordinate
    // caution: if there are supplementary locations in x1-direction, we loop
    //          them first (only DNS geometry)
    double x1c = 1.0e20;
    if (x1nodnum < numx1supplocations_)
    {
      x1c = x1supplocations_(x1nodnum);
    }
    else
    {
      x1c = x1statlocations_(x1nodnum - numx1supplocations_);
    }

    int x2nodnum = -1;
    //----------------------------------------------------------------------
    // loop nodes in x2-direction
    //----------------------------------------------------------------------
    for (std::vector<double>::iterator x2line = x2coordinates_->begin();
        x2line != x2coordinates_->end(); ++x2line)
    {
      x2nodnum++;

      // toggle vectors are one in the position of a dof of this node,
      // else 0
      toggleu_->put_scalar(0.0);
      togglev_->put_scalar(0.0);
      togglew_->put_scalar(0.0);
      togglep_->put_scalar(0.0);

      // count the number of nodes in x3-direction contributing to this nodal value
      int countnodes = 0;

      for (int nn = 0; nn < discret_->num_my_row_nodes(); ++nn)
      {
        Core::Nodes::Node* node = discret_->l_row_node(nn);

        // this is the node
        if ((node->x()[0] < (x1c + 2e-5) and node->x()[0] > (x1c - 2e-5)) and
            (node->x()[1] < (*x2line + 2e-9) and node->x()[1] > (*x2line - 2e-9)))
        {
          std::vector<int> dof = discret_->dof(node);
          double one = 1.0;

          toggleu_->replace_global_values(1, &one, &(dof[0]));
          togglev_->replace_global_values(1, &one, &(dof[1]));
          togglew_->replace_global_values(1, &one, &(dof[2]));
          togglep_->replace_global_values(1, &one, &(dof[3]));

          countnodes++;
        }
      }

      int countnodesonallprocs = 0;

      Core::Communication::sum_all(&countnodes, &countnodesonallprocs, 1, discret_->get_comm());

      // reduce by 1 due to periodic boundary condition
      countnodesonallprocs -= 1;

      if (countnodesonallprocs)
      {
        //----------------------------------------------------------------------
        // get values for velocity and pressure on this line
        //----------------------------------------------------------------------
        double u;
        double v;
        double w;
        double p;
        velnp.dot(*toggleu_, &u);
        velnp.dot(*togglev_, &v);
        velnp.dot(*togglew_, &w);
        velnp.dot(*togglep_, &p);

        double uu;
        double vv;
        double ww;
        double pp;
        squaredvelnp_->dot(*toggleu_, &uu);
        squaredvelnp_->dot(*togglev_, &vv);
        squaredvelnp_->dot(*togglew_, &ww);
        squaredvelnp_->dot(*togglep_, &pp);

        double uv;
        double uw;
        double vw;
        double locuv = 0.0;
        double locuw = 0.0;
        double locvw = 0.0;
        for (int rr = 1; rr < velnp.local_length(); ++rr)
        {
          locuv += ((velnp)[rr - 1] * (*toggleu_)[rr - 1]) * ((velnp)[rr] * (*togglev_)[rr]);
        }
        Core::Communication::sum_all(&locuv, &uv, 1, discret_->get_comm());
        for (int rr = 2; rr < velnp.local_length(); ++rr)
        {
          locuw += ((velnp)[rr - 2] * (*toggleu_)[rr - 2]) * ((velnp)[rr] * (*togglew_)[rr]);
        }
        Core::Communication::sum_all(&locuw, &uw, 1, discret_->get_comm());
        for (int rr = 2; rr < velnp.local_length(); ++rr)
        {
          locvw += ((velnp)[rr - 1] * (*togglev_)[rr - 1]) * ((velnp)[rr] * (*togglew_)[rr]);
        }
        Core::Communication::sum_all(&locvw, &vw, 1, discret_->get_comm());

        //----------------------------------------------------------------------
        // calculate spatial means on this line
        // add spatial mean values to statistical sample
        //----------------------------------------------------------------------
        (*x2sumu_)(x1nodnum, x2nodnum) += u / countnodesonallprocs;
        (*x2sumv_)(x1nodnum, x2nodnum) += v / countnodesonallprocs;
        (*x2sumw_)(x1nodnum, x2nodnum) += w / countnodesonallprocs;
        (*x2sump_)(x1nodnum, x2nodnum) += p / countnodesonallprocs;

        (*x2sumsqu_)(x1nodnum, x2nodnum) += uu / countnodesonallprocs;
        (*x2sumsqv_)(x1nodnum, x2nodnum) += vv / countnodesonallprocs;
        (*x2sumsqw_)(x1nodnum, x2nodnum) += ww / countnodesonallprocs;
        (*x2sumsqp_)(x1nodnum, x2nodnum) += pp / countnodesonallprocs;

        (*x2sumuv_)(x1nodnum, x2nodnum) += uv / countnodesonallprocs;
        (*x2sumuw_)(x1nodnum, x2nodnum) += uw / countnodesonallprocs;
        (*x2sumvw_)(x1nodnum, x2nodnum) += vw / countnodesonallprocs;
      }
    }
  }

  return;
}  // TurbulenceStatisticsBfs::DoTimeSample


//----------------------------------------------------------------------
// sampling of velocity, pressure and temperature values
//----------------------------------------------------------------------
void FLD::TurbulenceStatisticsBfs::do_loma_time_sample(
    Core::LinAlg::Vector<double>& velnp, Core::LinAlg::Vector<double>& scanp, const double eosfac)
{
  // compute squared values of velocity
  squaredvelnp_->multiply(1.0, velnp, velnp, 0.0);
  squaredscanp_->multiply(1.0, scanp, scanp, 0.0);
  // compute 1/T and (1/T)^2
  for (int rr = 0; rr < squaredscanp_->local_length(); ++rr)
  {
    if ((scanp)[rr] > 0)  // temperature in kelvin is always > 0
      (*invscanp_)[rr] = 1 / ((scanp)[rr]);
  }
  squaredinvscanp_->multiply(1.0, *invscanp_, *invscanp_, 0.0);

  //----------------------------------------------------------------------
  // increase sample counter
  //----------------------------------------------------------------------
  numsamp_++;

  int x1nodnum = -1;
  //----------------------------------------------------------------------
  // values at lower and upper wall
  //----------------------------------------------------------------------
  for (std::vector<double>::iterator x1line = x1coordinates_->begin();
      x1line != x1coordinates_->end(); ++x1line)
  {
    x1nodnum++;

    for (int x2nodnum = 0; x2nodnum < numx2statlocations_; ++x2nodnum)
    {
      // current x2-coordinate of respective wall
      double x2cwall = x2statlocations_(x2nodnum);

      // current x2-coordinate of supplementary location to respective wall
      double x2csupp = x2supplocations_(x2nodnum);

      // toggle vectors are one in the position of a dof of this node,
      // else 0
      toggleu_->put_scalar(0.0);
      togglep_->put_scalar(0.0);

      // count the number of nodes in x3-direction contributing to this nodal value
      int countnodes = 0;

      for (int nn = 0; nn < discret_->num_my_row_nodes(); ++nn)
      {
        Core::Nodes::Node* node = discret_->l_row_node(nn);

        // this is the wall node
        if ((node->x()[0] < (*x1line + 2e-9) and node->x()[0] > (*x1line - 2e-9)) and
            (node->x()[1] < (x2cwall + 2e-5) and node->x()[1] > (x2cwall - 2e-5)))
        {
          std::vector<int> dof = discret_->dof(node);
          double one = 1.0;

          togglep_->replace_global_values(1, &one, &(dof[3]));

          countnodes++;
        }
        // this is the supplementary node
        else if ((node->x()[0] < (*x1line + 2e-9) and node->x()[0] > (*x1line - 2e-9)) and
                 (node->x()[1] < (x2csupp + 2e-5) and node->x()[1] > (x2csupp - 2e-5)))
        {
          std::vector<int> dof = discret_->dof(node);
          double one = 1.0;

          toggleu_->replace_global_values(1, &one, dof.data());
        }
      }

      int countnodesonallprocs = 0;

      Core::Communication::sum_all(&countnodes, &countnodesonallprocs, 1, discret_->get_comm());

      // reduce by 1 due to periodic boundary condition
      countnodesonallprocs -= 1;

      if (countnodesonallprocs)
      {
        //----------------------------------------------------------------------
        // get values for velocity derivative, pressure and temperature
        //----------------------------------------------------------------------
        double u;
        velnp.dot(*toggleu_, &u);
        double p;
        velnp.dot(*togglep_, &p);
        double T;
        scanp.dot(*togglep_, &T);

        double rho;
        invscanp_->dot(*togglep_, &rho);
        // compute density: rho = eosfac/T
        rho *= eosfac;

        //----------------------------------------------------------------------
        // calculate spatial means
        //----------------------------------------------------------------------
        double usm = u / countnodesonallprocs;
        double psm = p / countnodesonallprocs;
        double Tsm = T / countnodesonallprocs;
        double rhosm = rho / countnodesonallprocs;

        //----------------------------------------------------------------------
        // add spatial mean values to statistical sample
        //----------------------------------------------------------------------
        (*x1sumu_)(x2nodnum, x1nodnum) += usm;
        (*x1sump_)(x2nodnum, x1nodnum) += psm;
        (*x1sumrho_)(x2nodnum, x1nodnum) += rhosm;
        (*x1sum_t_)(x2nodnum, x1nodnum) += Tsm;
      }
    }
  }

  //----------------------------------------------------------------------
  // loop locations for statistical evaluation in x1-direction
  //----------------------------------------------------------------------
  for (int x1nodnum = 0; x1nodnum < numx1statlocations_; ++x1nodnum)
  {
    // current x1-coordinate
    // caution: if there are supplementary locations in x1-direction, we loop
    //          them first (only DNS geometry)
    double x1c = 1.0e20;
    if (x1nodnum < numx1supplocations_)
    {
      x1c = x1supplocations_(x1nodnum);
    }
    else
    {
      x1c = x1statlocations_(x1nodnum - numx1supplocations_);
    }

    int x2nodnum = -1;
    //----------------------------------------------------------------------
    // loop nodes in x2-direction and calculate pointwise means
    //----------------------------------------------------------------------
    for (std::vector<double>::iterator x2line = x2coordinates_->begin();
        x2line != x2coordinates_->end(); ++x2line)
    {
      x2nodnum++;

      // toggle vectors are one in the position of a dof of this node,
      // else 0
      toggleu_->put_scalar(0.0);
      togglev_->put_scalar(0.0);
      togglew_->put_scalar(0.0);
      togglep_->put_scalar(0.0);

      // count the number of nodes in x3-direction contributing to this nodal value
      int countnodes = 0;

      for (int nn = 0; nn < discret_->num_my_row_nodes(); ++nn)
      {
        Core::Nodes::Node* node = discret_->l_row_node(nn);

        // this is the node
        if ((node->x()[0] < (x1c + 2e-5) and node->x()[0] > (x1c - 2e-5)) and
            (node->x()[1] < (*x2line + 2e-9) and node->x()[1] > (*x2line - 2e-9)))
        {
          std::vector<int> dof = discret_->dof(node);
          double one = 1.0;

          toggleu_->replace_global_values(1, &one, &(dof[0]));
          togglev_->replace_global_values(1, &one, &(dof[1]));
          togglew_->replace_global_values(1, &one, &(dof[2]));
          togglep_->replace_global_values(1, &one, &(dof[3]));

          countnodes++;
        }
      }

      int countnodesonallprocs = 0;

      Core::Communication::sum_all(&countnodes, &countnodesonallprocs, 1, discret_->get_comm());

      // reduce by 1 due to periodic boundary condition
      countnodesonallprocs -= 1;

      if (countnodesonallprocs)
      {
        //----------------------------------------------------------------------
        // get values for velocity, pressure, density, temperature, and
        // subgrid viscosity on this line
        //----------------------------------------------------------------------
        double u;
        double v;
        double w;
        double p;
        velnp.dot(*toggleu_, &u);
        velnp.dot(*togglev_, &v);
        velnp.dot(*togglew_, &w);
        velnp.dot(*togglep_, &p);

        double T;
        scanp.dot(*togglep_, &T);

        double uu;
        double vv;
        double ww;
        double pp;
        squaredvelnp_->dot(*toggleu_, &uu);
        squaredvelnp_->dot(*togglev_, &vv);
        squaredvelnp_->dot(*togglew_, &ww);
        squaredvelnp_->dot(*togglep_, &pp);

        double uv;
        double uw;
        double vw;
        double locuv = 0.0;
        double locuw = 0.0;
        double locvw = 0.0;
        for (int rr = 1; rr < velnp.local_length(); ++rr)
        {
          locuv += ((velnp)[rr - 1] * (*toggleu_)[rr - 1]) * ((velnp)[rr] * (*togglev_)[rr]);
        }
        Core::Communication::sum_all(&locuv, &uv, 1, discret_->get_comm());
        for (int rr = 2; rr < velnp.local_length(); ++rr)
        {
          locuw += ((velnp)[rr - 2] * (*toggleu_)[rr - 2]) * ((velnp)[rr] * (*togglew_)[rr]);
        }
        Core::Communication::sum_all(&locuw, &uw, 1, discret_->get_comm());
        for (int rr = 2; rr < velnp.local_length(); ++rr)
        {
          locvw += ((velnp)[rr - 1] * (*togglev_)[rr - 1]) * ((velnp)[rr] * (*togglew_)[rr]);
        }
        Core::Communication::sum_all(&locvw, &vw, 1, discret_->get_comm());

        double TT;
        squaredscanp_->dot(*togglep_, &TT);

        double uT;
        double vT;
        double locuT = 0.0;
        double locvT = 0.0;
        for (int rr = 3; rr < velnp.local_length(); ++rr)
        {
          locuT += ((velnp)[rr - 3] * (*toggleu_)[rr - 3]) * ((scanp)[rr] * (*togglep_)[rr]);
        }
        Core::Communication::sum_all(&locuT, &uT, 1, discret_->get_comm());
        for (int rr = 3; rr < velnp.local_length(); ++rr)
        {
          locvT += ((velnp)[rr - 2] * (*togglev_)[rr - 2]) * ((scanp)[rr] * (*togglep_)[rr]);
        }
        Core::Communication::sum_all(&locvT, &vT, 1, discret_->get_comm());

        double rho;
        invscanp_->dot(*togglep_, &rho);
        // compute density: rho = eosfac/T
        rho *= eosfac;
        double rhorho;
        squaredinvscanp_->dot(*togglep_, &rhorho);
        rhorho *= eosfac * eosfac;

        double rhou;
        double rhov;
        double locrhou = 0.0;
        double locrhov = 0.0;
        for (int rr = 3; rr < velnp.local_length(); ++rr)
        {
          locrhou += (eosfac * ((*invscanp_)[rr] * (*togglep_)[rr])) *
                     ((velnp)[rr - 3] * (*toggleu_)[rr - 3]);
        }
        Core::Communication::sum_all(&locrhou, &rhou, 1, discret_->get_comm());
        for (int rr = 3; rr < velnp.local_length(); ++rr)
        {
          locrhov += (eosfac * ((*invscanp_)[rr] * (*togglep_)[rr])) *
                     ((velnp)[rr - 2] * (*togglev_)[rr - 2]);
        }
        Core::Communication::sum_all(&locrhov, &rhov, 1, discret_->get_comm());


        //----------------------------------------------------------------------
        // calculate spatial means on this line
        // add spatial mean values to statistical sample
        //----------------------------------------------------------------------
        (*x2sumu_)(x1nodnum, x2nodnum) += u / countnodesonallprocs;
        (*x2sumv_)(x1nodnum, x2nodnum) += v / countnodesonallprocs;
        (*x2sumw_)(x1nodnum, x2nodnum) += w / countnodesonallprocs;
        (*x2sump_)(x1nodnum, x2nodnum) += p / countnodesonallprocs;

        (*x2sum_t_)(x1nodnum, x2nodnum) += T / countnodesonallprocs;
        (*x2sumrho_)(x1nodnum, x2nodnum) += rho / countnodesonallprocs;

        (*x2sumsqu_)(x1nodnum, x2nodnum) += uu / countnodesonallprocs;
        (*x2sumsqv_)(x1nodnum, x2nodnum) += vv / countnodesonallprocs;
        (*x2sumsqw_)(x1nodnum, x2nodnum) += ww / countnodesonallprocs;
        (*x2sumsqp_)(x1nodnum, x2nodnum) += pp / countnodesonallprocs;

        (*x2sumsq_t_)(x1nodnum, x2nodnum) += TT / countnodesonallprocs;
        (*x2sumsqrho_)(x1nodnum, x2nodnum) += rhorho / countnodesonallprocs;

        (*x2sumuv_)(x1nodnum, x2nodnum) += uv / countnodesonallprocs;
        (*x2sumuw_)(x1nodnum, x2nodnum) += uw / countnodesonallprocs;
        (*x2sumvw_)(x1nodnum, x2nodnum) += vw / countnodesonallprocs;

        (*x2sumrhou_)(x1nodnum, x2nodnum) += rhou / countnodesonallprocs;
        (*x2sumu_t_)(x1nodnum, x2nodnum) += uT / countnodesonallprocs;
        (*x2sumrhov_)(x1nodnum, x2nodnum) += rhov / countnodesonallprocs;
        (*x2sumv_t_)(x1nodnum, x2nodnum) += vT / countnodesonallprocs;
      }
    }
  }

  return;
}  // TurbulenceStatisticsBfc::DoLomaTimeSample


//----------------------------------------------------------------------
// sampling of velocity, pressure and scalar values
//----------------------------------------------------------------------
void FLD::TurbulenceStatisticsBfs::do_scatra_time_sample(
    Core::LinAlg::Vector<double>& velnp, Core::LinAlg::Vector<double>& scanp)
{
  // compute squared values of velocity
  squaredvelnp_->multiply(1.0, velnp, velnp, 0.0);
  squaredscanp_->multiply(1.0, scanp, scanp, 0.0);

  //----------------------------------------------------------------------
  // increase sample counter
  //----------------------------------------------------------------------
  numsamp_++;

  int x1nodnum = -1;
  //----------------------------------------------------------------------
  // values at lower and upper wall
  //----------------------------------------------------------------------
  for (std::vector<double>::iterator x1line = x1coordinates_->begin();
      x1line != x1coordinates_->end(); ++x1line)
  {
    x1nodnum++;

    for (int x2nodnum = 0; x2nodnum < numx2statlocations_; ++x2nodnum)
    {
      // current x2-coordinate of respective wall
      double x2cwall = x2statlocations_(x2nodnum);

      // current x2-coordinate of supplementary location to respective wall
      double x2csupp = x2supplocations_(x2nodnum);

      // toggle vectors are one in the position of a dof of this node,
      // else 0
      toggleu_->put_scalar(0.0);
      togglep_->put_scalar(0.0);

      // count the number of nodes in x3-direction contributing to this nodal value
      int countnodes = 0;

      for (int nn = 0; nn < discret_->num_my_row_nodes(); ++nn)
      {
        Core::Nodes::Node* node = discret_->l_row_node(nn);

        // this is the wall node
        if ((node->x()[0] < (*x1line + 2e-9) and node->x()[0] > (*x1line - 2e-9)) and
            (node->x()[1] < (x2cwall + 2e-5) and node->x()[1] > (x2cwall - 2e-5)))
        {
          std::vector<int> dof = discret_->dof(node);
          double one = 1.0;

          togglep_->replace_global_values(1, &one, &(dof[3]));

          countnodes++;
        }
        // this is the supplementary node
        else if ((node->x()[0] < (*x1line + 2e-9) and node->x()[0] > (*x1line - 2e-9)) and
                 (node->x()[1] < (x2csupp + 2e-5) and node->x()[1] > (x2csupp - 2e-5)))
        {
          std::vector<int> dof = discret_->dof(node);
          double one = 1.0;

          toggleu_->replace_global_values(1, &one, dof.data());
        }
      }

      int countnodesonallprocs = 0;

      Core::Communication::sum_all(&countnodes, &countnodesonallprocs, 1, discret_->get_comm());

      // reduce by 1 due to periodic boundary condition
      countnodesonallprocs -= 1;

      if (countnodesonallprocs)
      {
        //----------------------------------------------------------------------
        // get values for velocity derivative, pressure and temperature
        //----------------------------------------------------------------------
        double u;
        velnp.dot(*toggleu_, &u);
        double p;
        velnp.dot(*togglep_, &p);
        double T;
        scanp.dot(*togglep_, &T);

        //----------------------------------------------------------------------
        // calculate spatial means
        //----------------------------------------------------------------------
        double usm = u / countnodesonallprocs;
        double psm = p / countnodesonallprocs;
        double Tsm = T / countnodesonallprocs;

        //----------------------------------------------------------------------
        // add spatial mean values to statistical sample
        //----------------------------------------------------------------------
        (*x1sumu_)(x2nodnum, x1nodnum) += usm;
        (*x1sump_)(x2nodnum, x1nodnum) += psm;
        (*x1sum_t_)(x2nodnum, x1nodnum) += Tsm;
      }
    }
  }

  //----------------------------------------------------------------------
  // loop locations for statistical evaluation in x1-direction
  //----------------------------------------------------------------------
  for (int x1nodnum = 0; x1nodnum < numx1statlocations_; ++x1nodnum)
  {
    // current x1-coordinate
    // caution: if there are supplementary locations in x1-direction, we loop
    //          them first (only DNS geometry)
    double x1c = 1.0e20;
    if (x1nodnum < numx1supplocations_)
    {
      x1c = x1supplocations_(x1nodnum);
    }
    else
    {
      x1c = x1statlocations_(x1nodnum - numx1supplocations_);
    }

    int x2nodnum = -1;
    //----------------------------------------------------------------------
    // loop nodes in x2-direction and calculate pointwise means
    //----------------------------------------------------------------------
    for (std::vector<double>::iterator x2line = x2coordinates_->begin();
        x2line != x2coordinates_->end(); ++x2line)
    {
      x2nodnum++;

      // toggle vectors are one in the position of a dof of this node,
      // else 0
      toggleu_->put_scalar(0.0);
      togglev_->put_scalar(0.0);
      togglew_->put_scalar(0.0);
      togglep_->put_scalar(0.0);

      // count the number of nodes in x3-direction contributing to this nodal value
      int countnodes = 0;

      for (int nn = 0; nn < discret_->num_my_row_nodes(); ++nn)
      {
        Core::Nodes::Node* node = discret_->l_row_node(nn);

        // this is the node
        if ((node->x()[0] < (x1c + 2e-5) and node->x()[0] > (x1c - 2e-5)) and
            (node->x()[1] < (*x2line + 2e-9) and node->x()[1] > (*x2line - 2e-9)))
        {
          std::vector<int> dof = discret_->dof(node);
          double one = 1.0;

          toggleu_->replace_global_values(1, &one, &(dof[0]));
          togglev_->replace_global_values(1, &one, &(dof[1]));
          togglew_->replace_global_values(1, &one, &(dof[2]));
          togglep_->replace_global_values(1, &one, &(dof[3]));

          countnodes++;
        }
      }

      int countnodesonallprocs = 0;

      Core::Communication::sum_all(&countnodes, &countnodesonallprocs, 1, discret_->get_comm());

      // reduce by 1 due to periodic boundary condition
      countnodesonallprocs -= 1;

      if (countnodesonallprocs)
      {
        //----------------------------------------------------------------------
        // get values for velocity, pressure, density, temperature, and
        // subgrid viscosity on this line
        //----------------------------------------------------------------------
        double u;
        double v;
        double w;
        double p;
        velnp.dot(*toggleu_, &u);
        velnp.dot(*togglev_, &v);
        velnp.dot(*togglew_, &w);
        velnp.dot(*togglep_, &p);

        double T;
        scanp.dot(*togglep_, &T);

        double uu;
        double vv;
        double ww;
        double pp;
        squaredvelnp_->dot(*toggleu_, &uu);
        squaredvelnp_->dot(*togglev_, &vv);
        squaredvelnp_->dot(*togglew_, &ww);
        squaredvelnp_->dot(*togglep_, &pp);

        double uv;
        double uw;
        double vw;
        double locuv = 0.0;
        double locuw = 0.0;
        double locvw = 0.0;
        for (int rr = 1; rr < velnp.local_length(); ++rr)
        {
          locuv += ((velnp)[rr - 1] * (*toggleu_)[rr - 1]) * ((velnp)[rr] * (*togglev_)[rr]);
        }
        Core::Communication::sum_all(&locuv, &uv, 1, discret_->get_comm());
        for (int rr = 2; rr < velnp.local_length(); ++rr)
        {
          locuw += ((velnp)[rr - 2] * (*toggleu_)[rr - 2]) * ((velnp)[rr] * (*togglew_)[rr]);
        }
        Core::Communication::sum_all(&locuw, &uw, 1, discret_->get_comm());
        for (int rr = 2; rr < velnp.local_length(); ++rr)
        {
          locvw += ((velnp)[rr - 1] * (*togglev_)[rr - 1]) * ((velnp)[rr] * (*togglew_)[rr]);
        }
        Core::Communication::sum_all(&locvw, &vw, 1, discret_->get_comm());

        double TT;
        squaredscanp_->dot(*togglep_, &TT);

        double uT;
        double vT;
        double locuT = 0.0;
        double locvT = 0.0;
        for (int rr = 3; rr < velnp.local_length(); ++rr)
        {
          locuT += ((velnp)[rr - 3] * (*toggleu_)[rr - 3]) * ((scanp)[rr] * (*togglep_)[rr]);
        }
        Core::Communication::sum_all(&locuT, &uT, 1, discret_->get_comm());
        for (int rr = 3; rr < velnp.local_length(); ++rr)
        {
          locvT += ((velnp)[rr - 2] * (*togglev_)[rr - 2]) * ((scanp)[rr] * (*togglep_)[rr]);
        }
        Core::Communication::sum_all(&locvT, &vT, 1, discret_->get_comm());

        //----------------------------------------------------------------------
        // calculate spatial means on this line
        // add spatial mean values to statistical sample
        //----------------------------------------------------------------------
        (*x2sumu_)(x1nodnum, x2nodnum) += u / countnodesonallprocs;
        (*x2sumv_)(x1nodnum, x2nodnum) += v / countnodesonallprocs;
        (*x2sumw_)(x1nodnum, x2nodnum) += w / countnodesonallprocs;
        (*x2sump_)(x1nodnum, x2nodnum) += p / countnodesonallprocs;

        (*x2sum_t_)(x1nodnum, x2nodnum) += T / countnodesonallprocs;

        (*x2sumsqu_)(x1nodnum, x2nodnum) += uu / countnodesonallprocs;
        (*x2sumsqv_)(x1nodnum, x2nodnum) += vv / countnodesonallprocs;
        (*x2sumsqw_)(x1nodnum, x2nodnum) += ww / countnodesonallprocs;
        (*x2sumsqp_)(x1nodnum, x2nodnum) += pp / countnodesonallprocs;

        (*x2sumsq_t_)(x1nodnum, x2nodnum) += TT / countnodesonallprocs;

        (*x2sumuv_)(x1nodnum, x2nodnum) += uv / countnodesonallprocs;
        (*x2sumuw_)(x1nodnum, x2nodnum) += uw / countnodesonallprocs;
        (*x2sumvw_)(x1nodnum, x2nodnum) += vw / countnodesonallprocs;

        (*x2sumu_t_)(x1nodnum, x2nodnum) += uT / countnodesonallprocs;
        (*x2sumv_t_)(x1nodnum, x2nodnum) += vT / countnodesonallprocs;
      }
    }
  }

  return;
}  // TurbulenceStatisticsBfc::DoScatraTimeSample


/*----------------------------------------------------------------------*
 *
 *----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsBfs::dump_statistics(int step)
{
  //----------------------------------------------------------------------
  // output to log-file
  std::shared_ptr<std::ofstream> log;
  if (Core::Communication::my_mpi_rank(discret_->get_comm()) == 0)
  {
    std::string s(statistics_outfilename_);
    s.append(".flow_statistics");

    log = std::make_shared<std::ofstream>(s.c_str(), std::ios::out);
    (*log) << "# Statistics for turbulent incompressible flow over a backward-facing step (first- "
              "and second-order moments)";
    (*log) << "\n\n";
    (*log) << "# Statistics record ";
    (*log) << " (Steps " << step - numsamp_ + 1 << "--" << step << ")\n\n\n";
    (*log) << std::scientific;

    (*log) << "\n\n\n";
    (*log) << "# lower wall behind step\n";
    (*log) << "#     x1";
    (*log) << "           duxdy         pmean         tauw\n";

    // distance from wall to first node off wall
    double dist = x2supplocations_(0) - x2statlocations_(0);

    for (unsigned i = 0; i < x1coordinates_->size(); ++i)
    {
      if ((*x1coordinates_)[i] > -2e-9)
      {
        double lwx1u = (*x1sumu_)(0, i) / numsamp_;
        double lwx1duxdy = lwx1u / dist;
        double lwx1p = (*x1sump_)(0, i) / numsamp_;
        double lwx1tauw = (*x1sumtauw_)(0, i) / numsamp_;

        (*log) << " " << std::setw(11) << std::setprecision(4) << (*x1coordinates_)[i];
        (*log) << "   " << std::setw(11) << std::setprecision(4) << lwx1duxdy;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << lwx1p;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << lwx1tauw;
        (*log) << "\n";
      }
    }

    if (geotype_ == TurbulenceStatisticsBfs::geometry_LES_flow_with_heating ||
        geotype_ == TurbulenceStatisticsBfs::geometry_EXP_vogel_eaton)
    {
      (*log) << "\n\n\n";
      (*log) << "# upper wall\n";
      (*log) << "#     x1";
      (*log) << "           duxdy         pmean\n";

      // distance from wall to first node off wall
      dist = x2statlocations_(1) - x2supplocations_(1);

      for (unsigned i = 0; i < x1coordinates_->size(); ++i)
      {
        double uwx1u = (*x1sumu_)(1, i) / numsamp_;
        double uwx1duxdy = uwx1u / dist;
        double uwx1p = (*x1sump_)(1, i) / numsamp_;

        (*log) << " " << std::setw(11) << std::setprecision(4) << (*x1coordinates_)[i];
        (*log) << "   " << std::setw(11) << std::setprecision(4) << uwx1duxdy;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << uwx1p;
        (*log) << "\n";
      }
    }

    for (int i = 0; i < numx1statlocations_; ++i)
    {
      // current x1-coordinate
      // caution: if there are supplementary locations in x1-direction, we loop
      //          them first
      double x1 = 1.0e20;
      if (i < numx1supplocations_)
      {
        x1 = x1supplocations_(i);
      }
      else
      {
        x1 = x1statlocations_(i - numx1supplocations_);
      }


      (*log) << "\n\n\n";
      (*log) << "# line in x2-direction at x1 = " << std::setw(11) << std::setprecision(4) << x1
             << "\n";
      (*log) << "#     x2";
      (*log) << "           umean         vmean         wmean         pmean";
      (*log) << "         urms          vrms          wrms          prms";
      (*log) << "          u'v'          u'w'          v'w'\n";

      for (unsigned j = 0; j < x2coordinates_->size(); ++j)
      {
        double x2u = (*x2sumu_)(i, j) / numsamp_;
        double x2v = (*x2sumv_)(i, j) / numsamp_;
        double x2w = (*x2sumw_)(i, j) / numsamp_;
        double x2p = (*x2sump_)(i, j) / numsamp_;

        double x2urms = 0.0;
        double x2vrms = 0.0;
        double x2wrms = 0.0;
        double x2prms = 0.0;

        if (((*x2sumsqu_)(i, j) / numsamp_ - x2u * x2u) > 0.0)
          x2urms = std::sqrt((*x2sumsqu_)(i, j) / numsamp_ - x2u * x2u);
        if (((*x2sumsqv_)(i, j) / numsamp_ - x2v * x2v) > 0.0)
          x2vrms = std::sqrt((*x2sumsqv_)(i, j) / numsamp_ - x2v * x2v);
        if (((*x2sumsqw_)(i, j) / numsamp_ - x2w * x2w) > 0.0)
          x2wrms = std::sqrt((*x2sumsqw_)(i, j) / numsamp_ - x2w * x2w);
        if (((*x2sumsqp_)(i, j) / numsamp_ - x2p * x2p) > 0.0)
          x2prms = std::sqrt((*x2sumsqp_)(i, j) / numsamp_ - x2p * x2p);

        double x2uv = (*x2sumuv_)(i, j) / numsamp_ - x2u * x2v;
        double x2uw = (*x2sumuw_)(i, j) / numsamp_ - x2u * x2w;
        double x2vw = (*x2sumvw_)(i, j) / numsamp_ - x2v * x2w;

        (*log) << " " << std::setw(11) << std::setprecision(4) << (*x2coordinates_)[j];
        (*log) << "   " << std::setw(11) << std::setprecision(4) << x2u;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << x2v;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << x2w;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << x2p;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << x2urms;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << x2vrms;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << x2wrms;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << x2prms;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << x2uv;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << x2uw;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << x2vw;
        (*log) << "\n";
      }
    }
    log->flush();
  }

  return;

}  // TurbulenceStatisticsBfs::DumpStatistics


/*----------------------------------------------------------------------*
 *
 *----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsBfs::dump_loma_statistics(int step)
{
  //----------------------------------------------------------------------
  // output to log-file
  std::shared_ptr<std::ofstream> log;
  if (Core::Communication::my_mpi_rank(discret_->get_comm()) == 0)
  {
    std::string s(statistics_outfilename_);
    s.append(".loma_statistics");

    log = std::make_shared<std::ofstream>(s.c_str(), std::ios::out);
    (*log) << "# Statistics for turbulent variable-density flow over a backward-facing step at low "
              "Mach number (first- and second-order moments)";
    (*log) << "\n\n";
    (*log) << "# Caution: The following statistics have to be used carefully:\n";
    (*log) << "#          rhoumean, uTmean, rhovmean, vTmean, rhou'T', rhov'T'\n";
    (*log) << "#          there are not any reference values for rhoumean, uTmean, rhovmean, "
              "vTmean and rhou'T'\n";
    (*log) << "#          it is not clear what is denoted by rhou'T' and rhov'T' in the paper of "
              "Avnacha and Pletcher(2002)\n\n";
    (*log) << "# Statistics record ";
    (*log) << " (Steps " << step - numsamp_ + 1 << "--" << step << ")\n\n\n";
    (*log) << std::scientific;

    (*log) << "\n\n\n";
    (*log) << "# lower wall behind step\n";
    (*log) << "#        x1";
    (*log) << "                 duxdy               pmean             rhomean              Tmean\n";

    // distance from wall to first node off wall
    double dist = x2supplocations_(0) - x2statlocations_(0);

    for (unsigned i = 0; i < x1coordinates_->size(); ++i)
    {
      if ((*x1coordinates_)[i] > -2e-9)
      {
        double lwx1u = (*x1sumu_)(0, i) / numsamp_;
        double lwx1duxdy = lwx1u / dist;
        double lwx1p = (*x1sump_)(0, i) / numsamp_;

        double lwx1rho = (*x1sumrho_)(0, i) / numsamp_;
        double lwx1T = (*x1sum_t_)(0, i) / numsamp_;

        (*log) << " " << std::setw(17) << std::setprecision(10) << (*x1coordinates_)[i];
        (*log) << "   " << std::setw(17) << std::setprecision(10) << lwx1duxdy;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << lwx1p;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << lwx1rho;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << lwx1T;
        (*log) << "\n";
      }
    }

    if (geotype_ == TurbulenceStatisticsBfs::geometry_LES_flow_with_heating)
    {
      (*log) << "\n\n\n";
      (*log) << "# upper wall\n";
      (*log) << "#        x1";
      (*log)
          << "                 duxdy               pmean             rhomean              Tmean\n";

      // distance from wall to first node off wall
      dist = x2statlocations_(1) - x2supplocations_(1);

      for (unsigned i = 0; i < x1coordinates_->size(); ++i)
      {
        double uwx1u = (*x1sumu_)(1, i) / numsamp_;
        double uwx1duxdy = uwx1u / dist;
        double uwx1p = (*x1sump_)(1, i) / numsamp_;

        double uwx1rho = (*x1sumrho_)(1, i) / numsamp_;
        double uwx1T = (*x1sum_t_)(1, i) / numsamp_;

        (*log) << " " << std::setw(17) << std::setprecision(10) << (*x1coordinates_)[i];
        (*log) << "   " << std::setw(17) << std::setprecision(10) << uwx1duxdy;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << uwx1p;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << uwx1rho;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << uwx1T;
        (*log) << "\n";
      }
    }
    else if (geotype_ == TurbulenceStatisticsBfs::geometry_EXP_vogel_eaton)
      FOUR_C_THROW("geometry not implemented for loma yet");

    for (int i = 0; i < numx1statlocations_; ++i)
    {
      // current x1-coordinate
      // caution: if there are supplementary locations in x1-direction, we loop
      //          them first
      double x1 = 1.0e20;
      if (i < numx1supplocations_)
      {
        x1 = x1supplocations_(i);
      }
      else
      {
        x1 = x1statlocations_(i - numx1supplocations_);
      }

      (*log) << "\n\n\n";
      (*log) << "# line in x2-direction at x1 = " << std::setw(11) << std::setprecision(10) << x1
             << "\n";
      (*log) << "#        x2";
      (*log) << "                 umean               vmean               wmean               "
                "pmean             rhomean               Tmean            rhoumean        "
                "rhouTmean            rhovmean        rhovTmean";
      (*log) << "               urms                vrms                wrms                prms   "
                "            rhorms                Trms";
      (*log) << "                u'v'                u'w'                v'w'             rhou'T'  "
                "           rhov'T'\n";

      for (unsigned j = 0; j < x2coordinates_->size(); ++j)
      {
        double x2u = (*x2sumu_)(i, j) / numsamp_;
        double x2v = (*x2sumv_)(i, j) / numsamp_;
        double x2w = (*x2sumw_)(i, j) / numsamp_;
        double x2p = (*x2sump_)(i, j) / numsamp_;

        double x2rho = (*x2sumrho_)(i, j) / numsamp_;
        double x2T = (*x2sum_t_)(i, j) / numsamp_;

        double x2urms = 0.0;
        double x2vrms = 0.0;
        double x2wrms = 0.0;
        double x2prms = 0.0;

        if (((*x2sumsqu_)(i, j) / numsamp_ - x2u * x2u) > 0.0)
          x2urms = std::sqrt((*x2sumsqu_)(i, j) / numsamp_ - x2u * x2u);
        if (((*x2sumsqv_)(i, j) / numsamp_ - x2v * x2v) > 0.0)
          x2vrms = std::sqrt((*x2sumsqv_)(i, j) / numsamp_ - x2v * x2v);
        if (((*x2sumsqw_)(i, j) / numsamp_ - x2w * x2w) > 0.0)
          x2wrms = std::sqrt((*x2sumsqw_)(i, j) / numsamp_ - x2w * x2w);
        if (((*x2sumsqp_)(i, j) / numsamp_ - x2p * x2p) > 0.0)
          x2prms = std::sqrt((*x2sumsqp_)(i, j) / numsamp_ - x2p * x2p);

        // as T and rho are constant in the inflow section
        // <T(rho)^2>-<T(rho)>*<T(rho)> should be zero
        // however, due to small errors, <T(rho)^2>-<T(rho)>*<T(rho)>
        // is only approximately equal zero
        // hence, zero negative values should be excluded
        // as they produce nans
        double x2rhorms = 0.0;
        double x2Trms = 0.0;
        if (std::abs((*x2sumsqrho_)(i, j) / numsamp_ - x2rho * x2rho) > 1e-9)
          x2rhorms = std::sqrt((*x2sumsqrho_)(i, j) / numsamp_ - x2rho * x2rho);
        if (std::abs((*x2sumsq_t_)(i, j) / numsamp_ - x2T * x2T) > 1e-9)
          x2Trms = std::sqrt((*x2sumsq_t_)(i, j) / numsamp_ - x2T * x2T);

        double x2uv = (*x2sumuv_)(i, j) / numsamp_ - x2u * x2v;
        double x2uw = (*x2sumuw_)(i, j) / numsamp_ - x2u * x2w;
        double x2vw = (*x2sumvw_)(i, j) / numsamp_ - x2v * x2w;

        double x2rhou = (*x2sumrhou_)(i, j) / numsamp_ - x2u * x2rho;
        double x2uT = (*x2sumu_t_)(i, j) / numsamp_ - x2u * x2T;
        double x2rhov = (*x2sumrhov_)(i, j) / numsamp_ - x2v * x2rho;
        double x2vT = (*x2sumv_t_)(i, j) / numsamp_ - x2v * x2T;

        double x2rhouppTpp = x2rho * (x2uT - x2u * x2T);
        double x2rhovppTpp = x2rho * (x2vT - x2v * x2T);

        (*log) << " " << std::setw(17) << std::setprecision(10) << (*x2coordinates_)[j];
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2u;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2v;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2w;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2p;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2rho;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2T;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2rhou;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2uT;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2rhov;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2vT;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2urms;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2vrms;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2wrms;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2prms;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2rhorms;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2Trms;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2uv;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2uw;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2vw;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2rhouppTpp;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2rhovppTpp;
        (*log) << "\n";
      }
    }
    log->flush();
  }

  return;

}  // TurbulenceStatisticsBfs::DumpLomaStatistics


/*----------------------------------------------------------------------*
 *
 *----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsBfs::dump_scatra_statistics(int step)
{
  //----------------------------------------------------------------------
  // output to log-file
  std::shared_ptr<std::ofstream> log;
  if (Core::Communication::my_mpi_rank(discret_->get_comm()) == 0)
  {
    std::string s(statistics_outfilename_);
    s.append(".flow_statistics");

    log = std::make_shared<std::ofstream>(s.c_str(), std::ios::out);
    (*log) << "# Statistics for turbulent flow with passive scalar over a backward-facing step "
              "(first- and second-order moments)";
    (*log) << "\n\n";
    (*log) << "# Statistics record ";
    (*log) << " (Steps " << step - numsamp_ + 1 << "--" << step << ")\n\n\n";
    (*log) << std::scientific;

    (*log) << "\n\n\n";
    (*log) << "# lower wall behind step\n";
    (*log) << "#        x1";
    (*log) << "                 duxdy               pmean              phimean\n";

    // distance from wall to first node off wall
    double dist = x2supplocations_(0) - x2statlocations_(0);

    for (unsigned i = 0; i < x1coordinates_->size(); ++i)
    {
      if ((*x1coordinates_)[i] > -2e-9)
      {
        double lwx1u = (*x1sumu_)(0, i) / numsamp_;
        double lwx1duxdy = lwx1u / dist;
        double lwx1p = (*x1sump_)(0, i) / numsamp_;

        double lwx1T = (*x1sum_t_)(0, i) / numsamp_;

        (*log) << " " << std::setw(17) << std::setprecision(10) << (*x1coordinates_)[i];
        (*log) << "   " << std::setw(17) << std::setprecision(10) << lwx1duxdy;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << lwx1p;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << lwx1T;
        (*log) << "\n";
      }
    }

    if (geotype_ == TurbulenceStatisticsBfs::geometry_LES_flow_with_heating)
    {
      (*log) << "\n\n\n";
      (*log) << "# upper wall\n";
      (*log) << "#        x1";
      (*log) << "                 duxdy               pmean              phimean\n";

      // distance from wall to first node off wall
      dist = x2statlocations_(1) - x2supplocations_(1);

      for (unsigned i = 0; i < x1coordinates_->size(); ++i)
      {
        double uwx1u = (*x1sumu_)(1, i) / numsamp_;
        double uwx1duxdy = uwx1u / dist;
        double uwx1p = (*x1sump_)(1, i) / numsamp_;

        double uwx1T = (*x1sum_t_)(1, i) / numsamp_;

        (*log) << " " << std::setw(17) << std::setprecision(10) << (*x1coordinates_)[i];
        (*log) << "   " << std::setw(17) << std::setprecision(10) << uwx1duxdy;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << uwx1p;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << uwx1T;
        (*log) << "\n";
      }
    }
    else if (geotype_ == TurbulenceStatisticsBfs::geometry_EXP_vogel_eaton)
      FOUR_C_THROW("geometry not implemented for scatra yet");

    for (int i = 0; i < numx1statlocations_; ++i)
    {
      // current x1-coordinate
      // caution: if there are supplementary locations in x1-direction, we loop
      //          them first
      double x1 = 1.0e20;
      if (i < numx1supplocations_)
      {
        x1 = x1supplocations_(i);
      }
      else
      {
        x1 = x1statlocations_(i - numx1supplocations_);
      }

      (*log) << "\n\n\n";
      (*log) << "# line in x2-direction at x1 = " << std::setw(11) << std::setprecision(10) << x1
             << "\n";
      (*log) << "#        x2";
      (*log) << "                 umean               vmean               wmean               "
                "pmean               phimean           uphimean           vphimean";
      (*log) << "               urms                vrms                wrms                prms   "
                "             phirms";
      (*log) << "                u'v'                u'w'                v'w'\n";

      for (unsigned j = 0; j < x2coordinates_->size(); ++j)
      {
        double x2u = (*x2sumu_)(i, j) / numsamp_;
        double x2v = (*x2sumv_)(i, j) / numsamp_;
        double x2w = (*x2sumw_)(i, j) / numsamp_;
        double x2p = (*x2sump_)(i, j) / numsamp_;

        double x2T = (*x2sum_t_)(i, j) / numsamp_;

        double x2urms = 0.0;
        double x2vrms = 0.0;
        double x2wrms = 0.0;
        double x2prms = 0.0;

        if (((*x2sumsqu_)(i, j) / numsamp_ - x2u * x2u) > 0.0)
          x2urms = std::sqrt((*x2sumsqu_)(i, j) / numsamp_ - x2u * x2u);
        if (((*x2sumsqv_)(i, j) / numsamp_ - x2v * x2v) > 0.0)
          x2vrms = std::sqrt((*x2sumsqv_)(i, j) / numsamp_ - x2v * x2v);
        if (((*x2sumsqw_)(i, j) / numsamp_ - x2w * x2w) > 0.0)
          x2wrms = std::sqrt((*x2sumsqw_)(i, j) / numsamp_ - x2w * x2w);
        if (((*x2sumsqp_)(i, j) / numsamp_ - x2p * x2p) > 0.0)
          x2prms = std::sqrt((*x2sumsqp_)(i, j) / numsamp_ - x2p * x2p);

        // as T is constant in the inflow section
        // <T^2>-<T>*<T> should be zero
        // however, due to small errors, <T^2>-<T>*<T>
        // is only approximately equal zero
        // hence, zero negative values should be excluded
        // as they produce nans
        double x2Trms = 0.0;
        if (std::abs((*x2sumsq_t_)(i, j) / numsamp_ - x2T * x2T) > 1e-9)
          x2Trms = std::sqrt((*x2sumsq_t_)(i, j) / numsamp_ - x2T * x2T);

        double x2uv = (*x2sumuv_)(i, j) / numsamp_ - x2u * x2v;
        double x2uw = (*x2sumuw_)(i, j) / numsamp_ - x2u * x2w;
        double x2vw = (*x2sumvw_)(i, j) / numsamp_ - x2v * x2w;

        double x2uT = (*x2sumu_t_)(i, j) / numsamp_ - x2u * x2T;
        double x2vT = (*x2sumv_t_)(i, j) / numsamp_ - x2v * x2T;


        (*log) << " " << std::setw(17) << std::setprecision(10) << (*x2coordinates_)[j];
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2u;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2v;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2w;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2p;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2T;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2uT;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2vT;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2urms;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2vrms;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2wrms;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2prms;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2Trms;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2uv;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2uw;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2vw;
        (*log) << "\n";
      }
    }
    log->flush();
  }

  return;

}  // TurbulenceStatisticsBfs::dump_scatra_statistics


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsBfs::convert_string_to_geo_type(const std::string& geotype)
{
  FOUR_C_ASSERT(geotype != "none", "No geometry supplied");

  geotype_ = TurbulenceStatisticsBfs::none;
  if (geotype == "geometry_DNS_incomp_flow")
    geotype_ = TurbulenceStatisticsBfs::geometry_DNS_incomp_flow;
  else if (geotype == "geometry_LES_flow_with_heating")
    geotype_ = TurbulenceStatisticsBfs::geometry_LES_flow_with_heating;
  else if (geotype == "geometry_EXP_vogel_eaton")
    geotype_ = TurbulenceStatisticsBfs::geometry_EXP_vogel_eaton;
  else
    FOUR_C_THROW("({}) geometry for backward facing step", geotype.c_str());
  return;
}

FOUR_C_NAMESPACE_CLOSE
