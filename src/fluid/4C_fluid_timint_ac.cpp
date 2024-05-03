/*-----------------------------------------------------------*/
/*! \file

\brief Fluid time integrator for FS3I-AC problems


\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_fluid_timint_ac.hpp"

#include "4C_global_data.hpp"
#include "4C_io.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |  Constructor (public)                                     Thon 12/14 |
 *----------------------------------------------------------------------*/
FLD::TimIntAC::TimIntAC(const Teuchos::RCP<DRT::Discretization>& actdis,
    const Teuchos::RCP<CORE::LINALG::Solver>& solver,
    const Teuchos::RCP<Teuchos::ParameterList>& params,
    const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid /*= false*/)
    : FluidImplicitTimeInt(actdis, solver, params, output, alefluid)
{
  return;
}

/*----------------------------------------------------------------------*
 | output of solution vector to binio                        Thon 12/14 |
 *----------------------------------------------------------------------*/
void FLD::TimIntAC::ReadRestart(int step)
{
  const Teuchos::ParameterList& fs3idyn = GLOBAL::Problem::Instance()->FS3IDynamicParams();
  const bool restartfrompartfsi = CORE::UTILS::IntegralValue<int>(fs3idyn, "RESTART_FROM_PART_FSI");

  if (not restartfrompartfsi)  // standard restart
  {
    IO::DiscretizationReader reader(
        discret_, GLOBAL::Problem::Instance()->InputControlFile(), step);

    reader.ReadVector(trueresidual_, "trueresidual");
  }
}

/*----------------------------------------------------------------------*
 | output of solution vector to binio                        Thon 12/14 |
 *----------------------------------------------------------------------*/
void FLD::TimIntAC::Output()
{
  FluidImplicitTimeInt::Output();

  // output of solution
  if (uprestart_ > 0 and step_ % uprestart_ == 0)
  {
    output_->WriteVector("trueresidual", trueresidual_);
  }
  return;
}

FOUR_C_NAMESPACE_CLOSE