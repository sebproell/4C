/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of an isotropic growth strategy for the growth remodel mixture rule

\level 3
*/
/*----------------------------------------------------------------------*/

#include "4C_mixture_growth_strategy_isotropic.hpp"

#include "4C_mixture_growth_strategy.hpp"

FOUR_C_NAMESPACE_OPEN

MIXTURE::PAR::IsotropicGrowthStrategy::IsotropicGrowthStrategy(
    const Teuchos::RCP<CORE::MAT::PAR::Material>& matdata)
    : MIXTURE::PAR::MixtureGrowthStrategy(matdata)
{
}

std::unique_ptr<MIXTURE::MixtureGrowthStrategy>
MIXTURE::PAR::IsotropicGrowthStrategy::CreateGrowthStrategy()
{
  return std::make_unique<MIXTURE::IsotropicGrowthStrategy>();
}

void MIXTURE::IsotropicGrowthStrategy::EvaluateInverseGrowthDeformationGradient(
    CORE::LINALG::Matrix<3, 3>& iFgM, const MIXTURE::MixtureRule& mixtureRule,
    double currentReferenceGrowthScalar, int gp) const
{
  iFgM.Clear();

  for (int i = 0; i < 3; ++i)
  {
    iFgM(i, i) = std::pow(currentReferenceGrowthScalar, -1.0 / 3.0);
  }
}

void MIXTURE::IsotropicGrowthStrategy::EvaluateGrowthStressCmat(
    const MIXTURE::MixtureRule& mixtureRule, double currentReferenceGrowthScalar,
    const CORE::LINALG::Matrix<1, 6>& dCurrentReferenceGrowthScalarDC,
    const CORE::LINALG::Matrix<3, 3>& F, const CORE::LINALG::Matrix<6, 1>& E_strain,
    Teuchos::ParameterList& params, CORE::LINALG::Matrix<6, 1>& S_stress,
    CORE::LINALG::Matrix<6, 6>& cmat, const int gp, const int eleGID) const
{
  S_stress.Clear();
  cmat.Clear();
}
FOUR_C_NAMESPACE_CLOSE