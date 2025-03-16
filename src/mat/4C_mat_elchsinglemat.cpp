// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mat_elchsinglemat.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_utils_function_of_time.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Mat::PAR::ElchSingleMat::ElchSingleMat(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      diffusion_coefficient_concentration_dependence_funct_num_(
          matdata.parameters.get<int>("DIFF_COEF_CONC_DEP_FUNCT")),
      diffusion_coefficient_temperature_scaling_funct_num_(
          matdata.parameters.get<int>("DIFF_COEF_TEMP_SCALE_FUNCT")),
      number_diffusion_coefficient_params_(matdata.parameters.get<int>("DIFF_PARA_NUM")),
      diffusion_coefficient_params_(matdata.parameters.get<std::vector<double>>("DIFF_PARA")),
      number_diffusion_temp_scale_funct_params_(
          matdata.parameters.get<int>("DIFF_COEF_TEMP_SCALE_FUNCT_PARA_NUM")),
      diffusion_temp_scale_funct_params_(
          matdata.parameters.get<std::vector<double>>("DIFF_COEF_TEMP_SCALE_FUNCT_PARA")),
      conductivity_concentration_dependence_funct_num_(
          matdata.parameters.get<int>("COND_CONC_DEP_FUNCT")),
      conductivity_temperature_scaling_funct_num_(
          matdata.parameters.get<int>("COND_TEMP_SCALE_FUNCT")),
      number_conductivity_params_(matdata.parameters.get<int>("COND_PARA_NUM")),
      conductivity_params_(matdata.parameters.get<std::vector<double>>("COND_PARA")),
      number_conductivity_temp_scale_funct_params_(
          matdata.parameters.get<int>("COND_TEMP_SCALE_FUNCT_PARA_NUM")),
      conductivity_temp_scale_funct_params_(
          matdata.parameters.get<std::vector<double>>("COND_TEMP_SCALE_FUNCT_PARA")),
      R_(Global::Problem::instance()->elch_control_params().get<double>("GAS_CONSTANT"))
{
  // safety checks
  if (number_diffusion_coefficient_params_ !=
      static_cast<int>(diffusion_coefficient_params_.size()))
    FOUR_C_THROW("Mismatch in number of parameters for diffusion coefficient!");
  if (number_conductivity_params_ != static_cast<int>(conductivity_params_.size()))
    FOUR_C_THROW("Mismatch in number of parameters for conductivity!");
  if (number_diffusion_temp_scale_funct_params_ !=
      static_cast<int>(diffusion_temp_scale_funct_params_.size()))
    FOUR_C_THROW(
        "Mismatch in number of parameters for temp scale function for diffusion coefficient!");
  if (number_conductivity_temp_scale_funct_params_ !=
      static_cast<int>(conductivity_temp_scale_funct_params_.size()))
    FOUR_C_THROW("Mismatch in number of parameters for temp scale function for conductivity!");
  check_provided_params(
      diffusion_coefficient_concentration_dependence_funct_num_, diffusion_coefficient_params_);
  check_provided_params(conductivity_concentration_dependence_funct_num_, conductivity_params_);
  check_provided_params(
      diffusion_coefficient_temperature_scaling_funct_num_, diffusion_temp_scale_funct_params_);
  check_provided_params(
      conductivity_temperature_scaling_funct_num_, conductivity_temp_scale_funct_params_);
}


/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void Mat::PAR::ElchSingleMat::check_provided_params(
    const int functnr, const std::vector<double>& functparams)
{
  // name of specified curve
  std::string functionname;

  // expected number of parameters for specified curve
  unsigned int nfunctparams = 0;

  // check set of implemented functions with negative curve number
  if (functnr < 0)
  {
    switch (functnr)
    {
      case Mat::ElchSingleMat::CONSTANT_FUNCTION:
      {
        // constant value: functval=functparams[0];
        functionname = "'constant function'";
        nfunctparams = 1;
        break;
      }
      case Mat::ElchSingleMat::LINEAR_FUNCTION:
      {
        // linear function: functval=functparams[0]+functparams[1]*concentration;
        functionname = "'linear function'";
        nfunctparams = 2;
        break;
      }
      case Mat::ElchSingleMat::QUADRATIC_FUNCTION:
      {
        // quadratic function:
        // functval=functparams[0]+functparams[1]*concentration+functparams[2]*concentration*concentration;
        functionname = "'quadratic function'";
        nfunctparams = 3;
        break;
      }
      case Mat::ElchSingleMat::POWER_FUNCTION:
      {
        // power function: functval=functparams[0]*pow(concentration,functparams[1]);
        functionname = "'power function'";
        nfunctparams = 2;
        break;
      }
      case Mat::ElchSingleMat::CONDUCT:
      {
        // function 1 for conductivity;
        functionname = "'function 1 for conductivity'";
        nfunctparams = 4;
        break;
      }
      case Mat::ElchSingleMat::MOD_CUBIC_FUNCTION:
      {
        // a0*c + a1*c^1.5 + a2*c^3
        functionname = "'a0*c + a1*c^1.5 + a2*c^3'";
        nfunctparams = 3;
        break;
      }
      case Mat::ElchSingleMat::CUBIC_FUNCTION:
      {
        // a0 + a1*c + a2*c^2 + a3*c^3
        functionname = "'a0 + a1*c + a2*c^2 + a3*c^3'";
        nfunctparams = 4;
        break;
      }
      case Mat::ElchSingleMat::NYMAN:
      {
        // thermodynamic factor Nyman 2008
        functionname = "'function thermodynamic factor (Nyman 2008)'";
        nfunctparams = 7;
        break;
      }
      case Mat::ElchSingleMat::DEBYE_HUECKEL:
      {
        // linear thermodynamic factor including Debye-Hueckel theory
        functionname = "'function linear thermodynamic factor (including Debye Hueckel theory)'";
        nfunctparams = 2;
        break;
      }
      case Mat::ElchSingleMat::KOHLRAUSCH_SQUAREROOT:
      {
        // function 1 for conductivity
        functionname = "'function 1 for conductivity: own definition'";
        nfunctparams = 6;
        break;
      }
      case Mat::ElchSingleMat::GOLDIN:
      {
        // conductivity as a function of concentration according to Goldin, Colclasure, Wiedemann,
        // Kee (2012) kappa = a0*c*exp(a1*c^a2)
        functionname =
            "'conductivity as a function of concentration according to Goldin, Colclasure, "
            "Wiedemann, Kee (2012)'";
        nfunctparams = 3;
        break;
      }
      case Mat::ElchSingleMat::STEWART_NEWMAN:
      {
        // diffusion coefficient based on a function defined in
        // Stewart, S. G. & Newman, J. The Use of UV/vis Absorption to Measure Diffusion
        // Coefficients in LiPF6 Electrolytic Solutions Journal of The Electrochemical Society,
        // 2008, 155, F13-F16 diff = a0*exp(-a1*c^a2)
        functionname = "'diffusion coefficient as an exponential function: a1*exp(a2*c)'";
        nfunctparams = 2;
        break;
      }
      case Mat::ElchSingleMat::TDF:
      {
        // TDF based on a function defined in
        // J. Landesfeind, A. Ehrl, M. Graf, W.A. Wall, H.A. Gasteiger: Direct electrochemical
        // determination of activity coefficients in aprotic binary electrolytes TDF = 1.0 - 0.5 a1
        // sqrt(c)/(1+a2*sqrt(c)) + a2*c
        functionname =
            "'TDF as as a function of concentration according to Landesfeind, Ehrl, Graf, Wall, "
            "Gasteiger (2015)'";
        nfunctparams = 3;
        break;
      }
      case Mat::ElchSingleMat::ARRHENIUS:
      {
        // Arrhenius Ansatz for temperature dependent diffusion coefficient in solids D0 *
        // exp(-Q/(R*T)) Q: activation energy, R: universal gas constant, T: temperature, D0: max
        // diffusion coefficient D0 is provided by DIFF PARA and R is a constant which is already
        // defined
        functionname =
            "'Arrhenius Ansatz for temperature dependent diffusion coefficient in solids'";
        nfunctparams = 1;
        break;
      }
      case Mat::ElchSingleMat::INVERSE_LINEAR:
      {
        // Temperature dependent factor for electric conductivity (sigma)
        // electric conductivity is the inverse electric resistivity (rho)
        // for "small" temperature differences the temperature dependence of the specific electric
        // resistance follows rho = rho_0 * (1 + alpha(T - T_0))
        functionname =
            "'Linear approximation of specific electrical resistance/electric conductivity'";
        nfunctparams = 2;
        break;
      }
      default:
      {
        FOUR_C_THROW("Curve number {} is not implemented", functnr);
        break;
      }
    }

    // safety check
    if (functparams.size() != nfunctparams)
    {
      FOUR_C_THROW(
          "Number of provided parameters does not match number of expected parameters for function "
          "with curve number {} ({})!",
          functnr, functionname.c_str());
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double Mat::ElchSingleMat::compute_diffusion_coefficient(
    const double concentration, const double temperature) const
{
  // D(c)
  double diffusionCoefficient =
      compute_diffusion_coefficient_concentration_dependent(concentration);

  // D(c)*D(T)
  diffusionCoefficient *= compute_temperature_dependent_scale_factor(temperature,
      diffusion_coefficient_temperature_scaling_funct_num(), temp_scale_function_params_diff());

  return diffusionCoefficient;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double Mat::ElchSingleMat::compute_diffusion_coefficient_concentration_dependent(
    const double concentration) const
{
  double diffusionCoefficient(0.0);

  // evaluate pre implemented concentration dependent diffusion coefficient
  if (diffusion_coefficient_concentration_dependence_funct_num() < 0)
  {
    diffusionCoefficient =
        eval_pre_defined_funct(diffusion_coefficient_concentration_dependence_funct_num(),
            concentration, diffusion_coefficient_params());
  }
  else if (diffusion_coefficient_concentration_dependence_funct_num() == 0)
  {
    FOUR_C_THROW(
        "'DIFF_COEF_CONC_DEP_FUNCT' must not be 0! Either set it to a negative value to use one of "
        "the implemented models, or set a positive value and make use of the function framework!");
  }
  // diffusion coefficient is a function of the concentration as defined in the input file
  else
  {
    diffusionCoefficient = Global::Problem::instance()
                               ->function_by_id<Core::Utils::FunctionOfTime>(
                                   diffusion_coefficient_concentration_dependence_funct_num())
                               .evaluate(concentration);
  }

  return diffusionCoefficient;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double Mat::ElchSingleMat::compute_temperature_dependent_scale_factor(const double temperature,
    const int functionNumber, const std::vector<double>& functionParams) const
{
  double temperatureDependentScaleFactor(1.0);

  if (functionNumber == 0)
  {
    // do nothing
  }
  else if (functionNumber < 0)
  {
    temperatureDependentScaleFactor =
        eval_pre_defined_funct(functionNumber, temperature, functionParams);
  }
  else if (functionNumber > 0)
  {
    temperatureDependentScaleFactor =
        Global::Problem::instance()
            ->function_by_id<Core::Utils::FunctionOfTime>(functionNumber)
            .evaluate(temperature);
  }
  else
  {
    FOUR_C_THROW(
        "You have to set a reasonable function number for the temperature dependence.\n This can "
        "be either 0 if no temperature dependence is desired, the number of the function in which "
        "you defined the temperature dependence, or a negative number representing a predefined "
        "function");
  }

  return temperatureDependentScaleFactor;
}

/*------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------*/
double Mat::ElchSingleMat::compute_concentration_derivative_of_diffusion_coefficient(
    const double concentration, const double temperature) const
{
  // derivative of diffusion coefficient w.r.t. concentration
  double diffusion_coeff_conc_deriv(0.0);

  // evaluate derivative w.r.t. concentration of pre implemented concentration dependent diffusion
  // coefficient
  if (diffusion_coefficient_concentration_dependence_funct_num() < 0)
  {
    diffusion_coeff_conc_deriv = eval_first_deriv_pre_defined_funct(
        diffusion_coefficient_concentration_dependence_funct_num(), concentration,
        diffusion_coefficient_params());
  }
  else if (diffusion_coefficient_concentration_dependence_funct_num() == 0)
  {
    FOUR_C_THROW(
        "'DIFF_COEF_CONC_DEP_FUNCT' must not be 0! Either set it to a negative value to use one of "
        "the implemented models, or set a positive value and make use of the function framework!");
  }
  // evaluate derivative w.r.t. concentration of diffusion coefficient that is a function of the
  // concentration as defined in the input file
  else
  {
    diffusion_coeff_conc_deriv = (Global::Problem::instance()
            ->function_by_id<Core::Utils::FunctionOfTime>(
                diffusion_coefficient_concentration_dependence_funct_num())
            .evaluate_derivative(concentration));
  }

  // do the temperature dependent scaling
  diffusion_coeff_conc_deriv *= compute_temperature_dependent_scale_factor(temperature,
      diffusion_coefficient_temperature_scaling_funct_num(), temp_scale_function_params_diff());

  return diffusion_coeff_conc_deriv;
}


/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
double Mat::ElchSingleMat::compute_temperature_derivative_of_diffusion_coefficient(
    const double concentration, const double temperature) const
{
  // Computation D(c)
  double diffusion_coeff_temp_deriv =
      compute_diffusion_coefficient_concentration_dependent(concentration);

  // Computation D(c)*D'(T)
  diffusion_coeff_temp_deriv *= compute_temperature_dependent_scale_factor_deriv(temperature,
      diffusion_coefficient_temperature_scaling_funct_num(), temp_scale_function_params_diff());

  return diffusion_coeff_temp_deriv;
}

/*---------------------------------------------------------------------*

 *---------------------------------------------------------------------*/
double Mat::ElchSingleMat::compute_temperature_dependent_scale_factor_deriv(
    const double temperature, const int functionNumber,
    const std::vector<double>& functionParams) const
{
  double temperatureDependentScaleFactorDeriv(0.0);

  if (functionNumber == 0)
  {
    // do nothing
  }
  else if (functionNumber < 0)
  {
    temperatureDependentScaleFactorDeriv =
        eval_first_deriv_pre_defined_funct(functionNumber, temperature, functionParams);
  }
  else if (functionNumber > 0)
  {
    temperatureDependentScaleFactorDeriv =
        Global::Problem::instance()
            ->function_by_id<Core::Utils::FunctionOfTime>(functionNumber)
            .evaluate_derivative(temperature);
  }
  else
  {
    FOUR_C_THROW(
        "You have to set a reasonable function number for the temperature dependence.\n This can "
        "be either 0 if no temperature dependence is desired or the number of the function in "
        "which you defined the temperature dependence.");
  }

  return temperatureDependentScaleFactorDeriv;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double Mat::ElchSingleMat::compute_conductivity(
    const double concentration, const double temperature) const
{
  double conductivity = compute_conductivity_concentration_dependent(concentration);

  // do the temperature dependent scaling
  conductivity *= compute_temperature_dependent_scale_factor(
      temperature, conductivity_temperature_scaling_funct_num(), temp_scale_function_params_cond());

  return conductivity;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double Mat::ElchSingleMat::compute_conductivity_concentration_dependent(
    const double concentration) const
{
  double conductivity(0.0);

  // evaluate pre implemented concentration dependent conductivity
  if (conductivity_concentration_dependence_funct_num() < 0)
  {
    conductivity = eval_pre_defined_funct(
        conductivity_concentration_dependence_funct_num(), concentration, conductivity_params());
  }
  else if (conductivity_concentration_dependence_funct_num() == 0)
  {
    FOUR_C_THROW(
        "'COND_CONC_DEP_FUNCT' must not be 0! Either set it to a negative value to use one of the "
        "implemented models, or set a positive value and make use of the function framework!");
  }
  // conductivity is a function of the concentration as defined in the input file
  else
  {
    conductivity = Global::Problem::instance()
                       ->function_by_id<Core::Utils::FunctionOfTime>(
                           conductivity_concentration_dependence_funct_num())
                       .evaluate(concentration);
  }

  return conductivity;
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
double Mat::ElchSingleMat::compute_concentration_derivative_of_conductivity(
    const double concentration, const double temperature) const
{
  // derivative of conductivity w.r.t. concentration
  double conductivity_conc_deriv(0.0);

  // evaluate derivative w.r.t. concentration of pre implemented concentration dependent
  // conductivity
  if (conductivity_concentration_dependence_funct_num() < 0)
  {
    conductivity_conc_deriv = eval_first_deriv_pre_defined_funct(
        conductivity_concentration_dependence_funct_num(), concentration, conductivity_params());
  }
  else if (conductivity_concentration_dependence_funct_num() == 0)
  {
    FOUR_C_THROW(
        "'COND_CONC_DEP_FUNCT' must not be 0! Either set it to a negative value to use one of the "
        "implemented models, or set a positive value and make use of the function framework!");
  }
  // evaluate derivative w.r.t. concentration of conductivity that is a function of the
  // concentration as defined in the input file
  else
  {
    conductivity_conc_deriv = (Global::Problem::instance()
            ->function_by_id<Core::Utils::FunctionOfTime>(
                conductivity_concentration_dependence_funct_num())
            .evaluate_derivative(concentration));
  }

  // do the temperature dependent scaling
  conductivity_conc_deriv *= compute_temperature_dependent_scale_factor(
      temperature, conductivity_temperature_scaling_funct_num(), temp_scale_function_params_cond());

  return conductivity_conc_deriv;
}


/*-----------------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------------*/
double Mat::ElchSingleMat::compute_temperature_derivative_of_conductivity(
    const double concentration, const double temperature) const
{
  // Cond(c)
  double conductivity_temp_deriv = compute_conductivity_concentration_dependent(concentration);

  // Cond(c)*Cond'(T)
  // derivative of temp dependent scaling factor
  conductivity_temp_deriv *= compute_temperature_dependent_scale_factor_deriv(
      temperature, conductivity_temperature_scaling_funct_num(), temp_scale_function_params_cond());

  return conductivity_temp_deriv;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double Mat::ElchSingleMat::eval_pre_defined_funct(
    const int functnr, const double scalar, const std::vector<double>& functparams) const
{
  double functval(0.0);

  switch (functnr)
  {
    // a0
    case CONSTANT_FUNCTION:
      functval = functparams[0];
      break;

    // a0 + a1*c
    case LINEAR_FUNCTION:
      functval = functparams[0] + functparams[1] * scalar;
      break;

    // a0 + a1*c + a2*c^2
    case QUADRATIC_FUNCTION:
      functval = functparams[0] + functparams[1] * scalar + functparams[2] * scalar * scalar;
      break;

    // a0*c^a1
    case POWER_FUNCTION:
      functval = functparams[0] * std::pow(scalar, functparams[1]);
      break;

    // conductivity
    case CONDUCT:
    {
      const double nenner = (1.0 + functparams[2] * scalar * scalar -
                             functparams[3] * scalar * scalar * scalar * scalar);
      // functparams[0]*(functparams[1]*concentration/nenner) + 0.01 -> constant level 0.01 deleted
      // since it does not have a physical meaning (28.04.2014)
      functval = functparams[0] * (functparams[1] * scalar / nenner);
      break;
    }

    // a0*c + a1*c^1.5 + a2*c^3
    case MOD_CUBIC_FUNCTION:
      functval = functparams[0] * scalar + functparams[1] * std::pow(scalar, 1.5) +
                 functparams[2] * scalar * scalar * scalar;
      break;

    // a0 + a1*c + a2*c^2 + a3*c^3
    case CUBIC_FUNCTION:
      functval = functparams[0] + functparams[1] * scalar + functparams[2] * scalar * scalar +
                 functparams[3] * scalar * scalar * scalar;
      break;

    // thermodynamic factor Nyman 2008
    case NYMAN:
    {
      const double num =
          functparams[0] + functparams[1] * scalar + functparams[2] * scalar * scalar;
      const double denom = functparams[3] + functparams[4] * scalar +
                           functparams[5] * scalar * scalar +
                           functparams[6] * scalar * scalar * scalar;
      functval = num / denom;
      break;
    }

    // linear thermodynamic factor including Debye-Hueckel theory
    // 1 + a1*0.5*c^0.5 + a2*c
    case DEBYE_HUECKEL:
      functval = 1.0 + functparams[0] * 0.5 * std::pow(scalar, 0.5) + functparams[1] * scalar;
      break;

    // conductivity: own definition which also fulfills the Kohlrausches Square root law
    case KOHLRAUSCH_SQUAREROOT:
    {
      const double num = functparams[0] * scalar + functparams[1] * std::pow(scalar, 1.5) +
                         functparams[2] * scalar * scalar +
                         functparams[3] * scalar * scalar * scalar;
      const double denom = (1.0 + functparams[4] * scalar * scalar +
                            functparams[5] * scalar * scalar * scalar * scalar);
      // functparams[0]*(functparams[1]*concentration/nenner) + 0.01 -> constant level 0.01 deleted
      // since it does not have a physical meaning (28.04.2014)
      functval = num / denom;
      break;
    }

    // conductivity as a function of concentration according to Goldin, Colclasure, Wiedemann, Kee
    // (2012) kappa = a0*c*exp(a1*c^a2)
    case GOLDIN:
    {
      // safety check
      if (scalar < 1.e-12) FOUR_C_THROW("scalar value {} is zero or negative!", scalar);

      const double exponent = functparams[1] * std::pow(scalar, functparams[2]);

      // safety check
      if (exponent > 20.)
        FOUR_C_THROW("Overflow detected during conductivity evaluation! Exponent is too large: {}",
            exponent);

      functval = functparams[0] * scalar * std::exp(exponent);

      break;
    }
    // diffusion coefficient based on a function defined in
    // Stewart, S. G. & Newman, J. The Use of UV/vis Absorption to Measure Diffusion Coefficients in
    // LiPF6 Electrolytic Solutions Journal of The Electrochemical Society, 2008, 155, F13-F16 diff
    // = a0*exp(-a1*c^a2)
    case STEWART_NEWMAN:
    {
      functval = functparams[0] * std::exp(functparams[1] * scalar);
      break;
    }
    case TDF:
    {
      // TDF based on a function defined in
      // J. Landesfeind, A. Ehrl, M. Graf, W.A. Wall, H.A. Gasteiger: Direct electrochemical
      // determination of activity coefficients in aprotic binary electrolytes TDF = 1.0 - 0.5 a1
      // sqrt(c)/(1+a2*sqrt(c)) + a2*c
      functval = 1.0 -
                 (0.5 * functparams[0] * std::pow(scalar, 0.5)) /
                     (std::pow((1 + functparams[1] * std::pow(scalar, 0.5)), 2)) +
                 functparams[2] * scalar;
      break;
    }
    case ARRHENIUS:
    {
      // Arrhenius Ansatz for temperature dependent diffusion coefficient in solids D0 *
      // exp(-Q/(R*T)) Q: activation energy, R: universal gas constant, T:temperature, D0: max
      // diffusion coefficient D0 is provided by DIFF PARA
      // functval = exp(-a0/(R * T)
      const double R = static_cast<Mat::PAR::ElchSingleMat*>(parameter())->R_;
      functval = std::exp(-functparams[0] / (R * scalar));
      break;
    }
    case INVERSE_LINEAR:
    {
      // Temperature dependent factor for electric conductivity (sigma)
      // electric conductivity is the inverse electric resistivity (rho)
      // for "small" temperature differences the temperature dependence of the specific electric
      // resistance follows rho = rho_0 * (1 + alpha(T - T_0)) sigma(c,T) = sigma(c) * sigma(T) =
      // 1/rho_0(c) * 1/(1 + alpha*(T - T_0))
      // functval = 1/(1 + a0*(T - a1))
      functval = 1.0 / (1.0 + functparams[0] * (scalar - functparams[1]));
      break;
    }
    default:
    {
      FOUR_C_THROW("Curve number {} is not implemented!", functnr);
      break;
    }
  }

  return functval;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double Mat::ElchSingleMat::eval_first_deriv_pre_defined_funct(
    const int functnr, const double scalar, const std::vector<double>& functparams) const
{
  double firstderivfunctval(0.0);

  switch (functnr)
  {
    // d/dc: a0
    case CONSTANT_FUNCTION:
      firstderivfunctval = 0.0;
      break;

    // d/dc: a0 + a1*c
    case LINEAR_FUNCTION:
      firstderivfunctval = functparams[1];
      break;

    // d/dc: a0 + a1*c + a2*c^2
    case QUADRATIC_FUNCTION:
      firstderivfunctval = functparams[1] + 2 * functparams[2] * scalar;
      break;

    // d/dc: a0 + c^a1
    case POWER_FUNCTION:
      firstderivfunctval = functparams[0] * functparams[1] * std::pow(scalar, functparams[1] - 1.0);
      break;

    // d/dc: conductivity
    case CONDUCT:
    {
      const double nenner = (1.0 + functparams[2] * scalar * scalar -
                             functparams[3] * scalar * scalar * scalar * scalar);
      const double nennernenner = nenner * nenner;
      firstderivfunctval =
          functparams[0] *
          ((functparams[1] * nenner -
               functparams[1] * scalar *
                   (2 * functparams[2] * scalar - 4 * functparams[3] * scalar * scalar * scalar)) /
              nennernenner);
      break;
    }

    // d/dc: a0*c + a1*c^1.5 + a2*c^3
    case MOD_CUBIC_FUNCTION:
      firstderivfunctval = functparams[0] + 1.5 * functparams[1] * std::pow(scalar, 0.5) +
                           3 * functparams[2] * scalar * scalar;
      break;

    // d/dc: a0 + a1*c + a2*c^2 + a3*c^3
    case CUBIC_FUNCTION:
      firstderivfunctval =
          functparams[1] + 2 * functparams[2] * scalar + 3 * functparams[3] * scalar * scalar;
      break;

    // d/dc: thermodynamic factor Nyman 2008
    case NYMAN:
    {
      const double num =
          functparams[0] + functparams[1] * scalar + functparams[2] * scalar * scalar;
      const double denom = functparams[3] + functparams[4] * scalar +
                           functparams[5] * scalar * scalar +
                           functparams[6] * scalar * scalar * scalar;
      const double denomdenom = denom * denom;
      const double derivnum = functparams[1] + 2 * functparams[2] * scalar;
      const double derivdenom =
          functparams[4] + 2 * functparams[5] * scalar + 3 * functparams[6] * scalar * scalar;
      firstderivfunctval = (derivnum * denom - num * derivdenom) / denomdenom;
      break;
    }

    // linear thermodynamic factor including Debye-Hueckel theory
    // d/dc: 1 + a1*0.5*c^0.5 + a2*c
    case DEBYE_HUECKEL:
      firstderivfunctval = functparams[0] * 0.5 * 0.5 * std::pow(scalar, -0.5) + functparams[1];
      break;

    // d/dc: conductivity: own definition which also fulfills the Kohlrausches Square root law
    case KOHLRAUSCH_SQUAREROOT:
    {
      const double num = functparams[0] * scalar + functparams[1] * std::pow(scalar, 1.5) +
                         functparams[2] * scalar * scalar +
                         functparams[3] * scalar * scalar * scalar;
      const double denom = (1.0 + functparams[4] * scalar * scalar +
                            functparams[5] * scalar * scalar * scalar * scalar);
      const double denomdenom = denom * denom;
      const double derivnum = functparams[0] + 1.5 * functparams[1] * std::pow(scalar, 0.5) +
                              2.0 * functparams[2] * scalar +
                              3.0 * functparams[3] * scalar * scalar;
      const double derivdenom =
          2.0 * functparams[4] * scalar + 4.0 * functparams[5] * scalar * scalar * scalar;
      firstderivfunctval = ((derivnum * denom - num * derivdenom) / denomdenom);
      break;
    }

    // conductivity as a function of concentration according to Goldin, Colclasure, Wiedemann, Kee
    // (2012) d/dc: kappa = a0*c*exp(a1*c^a2)
    case GOLDIN:
    {
      // safety check
      if (scalar < 1.0e-12) FOUR_C_THROW("scalar value {} is zero or negative!", scalar);

      const double exponent = functparams[1] * std::pow(scalar, functparams[2]);

      // safety check
      if (std::abs(exponent) > 20.0)
      {
        FOUR_C_THROW(
            "Overflow detected during conductivity evaluation! Absolute value of exponent is too "
            "large: {}",
            exponent);
      }

      firstderivfunctval = functparams[0] * std::exp(exponent) *
                           (1 + functparams[1] * functparams[2] * std::pow(scalar, functparams[2]));

      break;
    }
    // diffusion coefficient based on a function defined in
    // Stewart, S. G. & Newman, J. The Use of UV/vis Absorption to Measure Diffusion Coefficients in
    // LiPF6 Electrolytic Solutions Journal of The Electrochemical Society, 2008, 155, F13-F16 diff
    // = a0*exp(-a1*c^a2) deriv (diff) = a0*a1*exp(a1*c^a2)
    case STEWART_NEWMAN:
    {
      firstderivfunctval = functparams[0] * functparams[1] * std::exp(functparams[1] * scalar);
      break;
    }
    case TDF:
    {
      // TDF based on a function defined in
      // J. Landesfeind, A. Ehrl, M. Graf, W.A. Wall, H.A. Gasteiger: Direct electrochemical
      // determination of activity coefficients in aprotic binary electrolytes TDF = 1.0 - 0.5 a1
      // sqrt(c)/(1+a2*sqrt(c)) + a2*c
      firstderivfunctval = -(0.25 * functparams[0] * std::pow(scalar, -0.5)) /
                               (std::pow((1 + functparams[1] * std::pow(scalar, 0.5)), 2)) +
                           (0.5 * functparams[0] * functparams[1]) /
                               (std::pow((1 + functparams[1] * std::pow(scalar, 0.5)), 3)) +
                           functparams[2];
      break;
    }
    case ARRHENIUS:
    {
      // Arrhenius Ansatz for temperature dependent diffusion coefficient in solids D0 *
      // exp(-Q/(R*T)) Q: activation energy, R: universal gasconstant, T:temperature, D0: max
      // diffusioncoefficient D0 is provided by DIFF PARA
      const double R = static_cast<Mat::PAR::ElchSingleMat*>(parameter())->R_;
      firstderivfunctval = std::exp(-functparams[0] / (R * scalar)) * functparams[0] / R * 1.0 /
                           std::pow(scalar, 2.0);
      break;
    }
    case INVERSE_LINEAR:
    {
      // Temperature dependent factor for electric conductivity (sigma)
      // electric conductivity is the inverse  electric resistivity (rho)
      // for "small" temperature differences the temperature dependence of the specific electric
      // resistivity follows rho = rho_0 * (1 + alpha(T - T_0)) sigma(c,T) = sigma(c) * sigma(T) =
      // 1/rho_0(c) * 1/(1 + alpha*(T - T_0))
      double base = 1.0 + functparams[0] * (scalar - functparams[1]);
      firstderivfunctval = -functparams[0] * std::pow(base, -2.0);
      break;
    }
    default:
    {
      FOUR_C_THROW("Curve number {} is not implemented!", functnr);
      break;
    }
  }

  return firstderivfunctval;
}

FOUR_C_NAMESPACE_CLOSE
