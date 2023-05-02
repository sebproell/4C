/*-----------------------------------------------------------*/
/*! \file

\brief Stationary problem driver for porous fluid


\level 2

*/
/*-----------------------------------------------------------*/

#include "fluid_timint_poro_stat.H"
#include "io.H"


FLD::TimIntPoroStat::TimIntPoroStat(const Teuchos::RCP<DRT::Discretization>& actdis,
    const Teuchos::RCP<LINALG::Solver>& solver, const Teuchos::RCP<Teuchos::ParameterList>& params,
    const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid)
    : FluidImplicitTimeInt(actdis, solver, params, output, alefluid),
      TimIntStationary(actdis, solver, params, output, alefluid),
      TimIntPoro(actdis, solver, params, output, alefluid)
{
}

void FLD::TimIntPoroStat::Init()
{
  // call Init()-functions of base classes
  // note: this order is important
  TimIntStationary::Init();
  TimIntPoro::Init();
}

void FLD::TimIntPoroStat::ReadRestart(int step)
{
  // call of base classes
  TimIntStationary::ReadRestart(step);
  TimIntPoro::ReadRestart(step);
}