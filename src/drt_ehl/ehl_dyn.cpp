/*--------------------------------------------------------------------------*/
/*!
\file ehl_dyn.cpp

\brief dummy

<pre>
Maintainer: Andy Wirtz
            wirtz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15270
</pre>
*/
/*--------------------------------------------------------------------------*/

#include "../drt_ehl/ehl_dyn.H"

#include "../drt_adapter/adapter_lubrication.H"

#include "../drt_lib/drt_utils_createdis.H"

#include "../drt_lubrication/lubrication_timint_implicit.H"


/*----------------------------------------------------------------------*
 | Main control routine for EHL problems                    wirtz 11/15 |
 *----------------------------------------------------------------------*/
void ehl_dyn(int restart)
{
  // access the communicator
  const Epetra_Comm& comm = DRT::Problem::Instance()->GetDis("lubrication")->Comm();

  // print problem type
  if (comm.MyPID() == 0)
  {
    std::cout << "###################################################"
        << std::endl;
    std::cout << "# YOUR PROBLEM TYPE: "
        << DRT::Problem::Instance()->ProblemName() << std::endl;
    std::cout << "###################################################"
        << std::endl;
  }

  // access the problem-specific parameter list
  const Teuchos::ParameterList& lubricationdyn =
      DRT::Problem::Instance()->LubricationDynamicParams();

  // access the lubrication discretization
  Teuchos::RCP<DRT::Discretization> lubricationdis =
      DRT::Problem::Instance()->GetDis("lubrication");

  lubricationdis->FillComplete();

  // we directly use the elements from the Lubrication elements section
  if (lubricationdis->NumGlobalNodes() == 0)
    dserror("No elements in the ---LUBRICATION ELEMENTS section");

  // add proxy of velocity related degrees of freedom to Lubrication discretization
  if (lubricationdis->BuildDofSetAuxProxy(DRT::Problem::Instance()->NDim() + 1, 0, 0,
      true) != 1)
    dserror("Lubrication discretization has illegal number of dofsets!");

  // finalize discretization
  lubricationdis->FillComplete(true, false, false);

  // get linear solver id from LUBRICATION DYNAMIC
  const int linsolvernumber = lubricationdyn.get<int>("LINEAR_SOLVER");
  if (linsolvernumber == (-1))
    dserror(
        "no linear solver defined for LUBRICATION problem. Please set LINEAR_SOLVER in LUBRICATION DYNAMIC to a valid number!");

  // create instance of Lubrication basis algorithm
  Teuchos::RCP<ADAPTER::LubricationBaseAlgorithm> lubricationonly = Teuchos::rcp(
      new ADAPTER::LubricationBaseAlgorithm());

  // setup Lubrication basis algorithm
  lubricationonly->Setup(lubricationdyn, lubricationdyn,
          DRT::Problem::Instance()->SolverParams(linsolvernumber));

  // read the restart information, set vectors and variables
  if (restart)
    lubricationonly->LubricationField()->ReadRestart(restart);

  // enter time loop to solve problem
  (lubricationonly->LubricationField())->TimeLoop();

  // perform the result test if required
  DRT::Problem::Instance()->AddFieldTest(lubricationonly->CreateLubricationFieldTest());
  DRT::Problem::Instance()->TestAll(comm);

return;

} // end of ehl_dyn()