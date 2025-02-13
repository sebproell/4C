// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_poromultiphase_scatra_function_parameters.hpp"

#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
PoroMultiPhaseScaTra::TumorGrowthLawHeavisideParameters::TumorGrowthLawHeavisideParameters(
    const std::map<std::string, double>& funct_params)
    : gamma_T_growth(funct_params.at("gamma_T_growth")),
      w_nl_crit(funct_params.at("w_nl_crit")),
      w_nl_env(funct_params.at("w_nl_env")),
      lambda(funct_params.at("lambda")),
      p_t_crit(funct_params.at("p_t_crit"))
{
  constexpr int NUMBER_OF_PARAMS = 5;
  FOUR_C_ASSERT_ALWAYS(funct_params.size() == NUMBER_OF_PARAMS, "Wrong size of funct_params");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
PoroMultiPhaseScaTra::NecrosisLawHeavisideParameters::NecrosisLawHeavisideParameters(
    const std::map<std::string, double>& funct_params)
    : gamma_t_necr(funct_params.at("gamma_t_necr")),
      w_nl_crit(funct_params.at("w_nl_crit")),
      w_nl_env(funct_params.at("w_nl_env")),
      delta_a_t(funct_params.at("delta_a_t")),
      p_t_crit(funct_params.at("p_t_crit"))
{
  constexpr int NUMBER_OF_PARAMS = 5;
  FOUR_C_ASSERT_ALWAYS(funct_params.size() == NUMBER_OF_PARAMS, "Wrong size of funct_params");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
PoroMultiPhaseScaTra::OxygenConsumptionLawHeavisideParameters::
    OxygenConsumptionLawHeavisideParameters(const std::map<std::string, double>& funct_params)
    : gamma_nl_growth(funct_params.at("gamma_nl_growth")),
      gamma_0_nl(funct_params.at("gamma_0_nl")),
      w_nl_crit(funct_params.at("w_nl_crit")),
      w_nl_env(funct_params.at("w_nl_env")),
      p_t_crit(funct_params.at("p_t_crit"))
{
  constexpr int NUMBER_OF_PARAMS = 5;
  FOUR_C_ASSERT_ALWAYS(funct_params.size() == NUMBER_OF_PARAMS, "Wrong size of funct_params");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
PoroMultiPhaseScaTra::TumorGrowthLawHeavisideNecroOxyParameters::
    TumorGrowthLawHeavisideNecroOxyParameters(const std::map<std::string, double>& funct_params)
    : gamma_T_growth(funct_params.at("gamma_T_growth")),
      w_nl_crit(funct_params.at("w_nl_crit")),
      w_nl_env(funct_params.at("w_nl_env")),
      lambda(funct_params.at("lambda")),
      p_t_crit(funct_params.at("p_t_crit"))
{
  constexpr int NUMBER_OF_PARAMS = 5;
  FOUR_C_ASSERT_ALWAYS(funct_params.size() == NUMBER_OF_PARAMS, "Wrong size of funct_params");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
PoroMultiPhaseScaTra::OxygenTransvascularExchangeLawContParameters::
    OxygenTransvascularExchangeLawContParameters(const std::map<std::string, double>& funct_params)
    : n(funct_params.at("n")),
      Pb50(funct_params.at("Pb50")),
      CaO2_max(funct_params.at("CaO2_max")),
      alpha_bl_eff(funct_params.at("alpha_bl_eff")),
      gammarhoSV(funct_params.at("gammarhoSV")),
      rho_oxy(funct_params.at("rho_oxy")),
      rho_if(funct_params.at("rho_if")),
      rho_bl(funct_params.at("rho_bl")),
      alpha_IF(funct_params.at("alpha_IF"))
{
  constexpr int NUMBER_OF_PARAMS = 9;
  FOUR_C_ASSERT_ALWAYS(funct_params.size() == NUMBER_OF_PARAMS, "Wrong size of funct_params");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
PoroMultiPhaseScaTra::OxygenTransvascularExchangeLawDiscParameters::
    OxygenTransvascularExchangeLawDiscParameters(const std::map<std::string, double>& funct_params)
    : n(funct_params.at("n")),
      Pb50(funct_params.at("Pb50")),
      CaO2_max(funct_params.at("CaO2_max")),
      alpha_bl_eff(funct_params.at("alpha_bl_eff")),
      gammarho(funct_params.at("gamma*rho")),
      rho_oxy(funct_params.at("rho_oxy")),
      rho_if(funct_params.at("rho_IF")),
      rho_bl(funct_params.at("rho_bl")),
      S2_max(funct_params.at("S2_max")),
      alpha_IF(funct_params.at("alpha_IF"))
{
  constexpr int NUMBER_OF_PARAMS = 10;
  FOUR_C_ASSERT_ALWAYS(funct_params.size() == NUMBER_OF_PARAMS, "Wrong size of funct_params");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
PoroMultiPhaseScaTra::LungOxygenExchangeLawParameters::LungOxygenExchangeLawParameters(
    const std::map<std::string, double>& funct_params)
    : rho_oxy(funct_params.at("rho_oxy")),
      DiffAdVTLC(funct_params.at("DiffAdVTLC")),
      alpha_oxy(funct_params.at("alpha_oxy")),
      rho_air(funct_params.at("rho_air")),
      rho_bl(funct_params.at("rho_bl")),
      n(funct_params.at("n")),
      P_oB50(funct_params.at("P_oB50")),
      NC_Hb(funct_params.at("NC_Hb")),
      P_atmospheric(funct_params.at("P_atmospheric")),
      volfrac_blood_ref(funct_params.at("volfrac_blood_ref"))
{
  constexpr int NUMBER_OF_PARAMS = 10;
  FOUR_C_ASSERT_ALWAYS(funct_params.size() == NUMBER_OF_PARAMS, "Wrong size of funct_params");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
PoroMultiPhaseScaTra::LungCarbonDioxideExchangeLawParameters::
    LungCarbonDioxideExchangeLawParameters(const std::map<std::string, double>& funct_params)
    : rho_CO2(funct_params.at("rho_CO2")),
      DiffsolAdVTLC(funct_params.at("DiffsolAdVTLC")),
      pH(funct_params.at("pH")),
      rho_air(funct_params.at("rho_air")),
      rho_bl(funct_params.at("rho_bl")),
      rho_oxy(funct_params.at("rho_oxy")),
      n(funct_params.at("n")),
      P_oB50(funct_params.at("P_oB50")),
      C_Hb(funct_params.at("C_Hb")),
      NC_Hb(funct_params.at("NC_Hb")),
      alpha_oxy(funct_params.at("alpha_oxy")),
      P_atmospheric(funct_params.at("P_atmospheric")),
      ScalingFormmHg(funct_params.at("ScalingFormmHg")),
      volfrac_blood_ref(funct_params.at("volfrac_blood_ref"))
{
  constexpr int NUMBER_OF_PARAMS = 14;
  FOUR_C_ASSERT_ALWAYS(funct_params.size() == NUMBER_OF_PARAMS, "Wrong size of funct_params");
}

FOUR_C_NAMESPACE_CLOSE
