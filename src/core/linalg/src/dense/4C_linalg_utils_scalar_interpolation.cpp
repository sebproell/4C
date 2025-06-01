// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_config.hpp"

#include "4C_linalg_utils_scalar_interpolation.hpp"

#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_fixedsizematrix_solver.hpp"
#include "4C_utils_exceptions.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
template <unsigned int loc_dim, unsigned int poly_order, unsigned int num_cofeffcients>
Core::LinAlg::Matrix<num_cofeffcients, 1>
Core::LinAlg::ScalarInterpolation::polynomialshapefunction(
    const Core::LinAlg::Matrix<loc_dim, 1>& location)
{
  Core::LinAlg::Matrix<num_cofeffcients, 1> polynomial(Core::LinAlg::Initialization::zero);

  switch (loc_dim)
  {
    case 3:
    {
      double x = location(0, 0);
      double y = location(1, 0);
      double z = location(2, 0);

      switch (poly_order)
      {
        case 1:  // [1 x y z xy xz yz xyz] tri-linear
          if (num_cofeffcients == 8)
          {
            polynomial(0) = 1.;
            polynomial(1) = x;
            polynomial(2) = y;
            polynomial(3) = z;
            polynomial(4) = x * y;
            polynomial(5) = x * z;
            polynomial(6) = y * z;
            polynomial(7) = x * y * z;
          }
          break;
        case 2:  // [1 x y z x^2 y^2 z^2 xy xz yz] quadratic
          if (num_cofeffcients == 10)
          {
            polynomial(0) = 1.;
            polynomial(1) = x;
            polynomial(2) = y;
            polynomial(3) = z;
            polynomial(4) = x * x;
            polynomial(5) = y * y;
            polynomial(6) = z * z;
            polynomial(7) = x * y;
            polynomial(8) = x * z;
            polynomial(9) = y * z;
          }
          break;
        default:
          FOUR_C_THROW("Unknown 3D shape function order");
          return polynomial;
      }
    }
    break;
    case 2:
    {
      double x = location(0, 0);
      double y = location(1, 0);

      switch (poly_order)
      {
        case 1:
          if (num_cofeffcients == 4)  // [1 x y xy] linear
          {
            polynomial(0) = 1.;
            polynomial(1) = x;
            polynomial(2) = y;
            polynomial(3) = x * y;
          }
          break;
        case 2:
          if (num_cofeffcients == 8)  // [1 x y xy x^2 y^2] quadartic complete
          {
            polynomial(0) = 1.;
            polynomial(1) = x;
            polynomial(2) = y;
            polynomial(3) = x * y;
            polynomial(4) = x * x;
            polynomial(5) = y * y;
            polynomial(6) = x * x * y;
            polynomial(7) = x * y * y;
          }
          else if (num_cofeffcients == 6)  // [1 x y xy x^2 y^2] quadratic incomplete
          {
            polynomial(0) = 1.;
            polynomial(1) = x;
            polynomial(2) = y;
            polynomial(3) = x * y;
            polynomial(4) = x * x;
            polynomial(5) = y * y;
          }
          break;
        default:
          FOUR_C_THROW("Unknown 2D shape function order");
          return polynomial;
      }
    }
    break;
    case 1:
    {
      double x = location(0, 0);

      switch (poly_order)
      {
        case 1:
          if (num_cofeffcients == 2)  // [1 x] linear
          {
            polynomial(0) = 1.;
            polynomial(1) = x;
          }
          break;
        case 2:
          if (num_cofeffcients == 3)  // [1 x x*x] quadratic
          {
            polynomial(0) = 1.;
            polynomial(1) = x;
            polynomial(2) = x * x;
          }
          break;
        default:
          FOUR_C_THROW("Unknown 1D shape function order");
          return polynomial;
      }
    }
    break;
    default:
      FOUR_C_THROW("Unsupported spatial dimension for shape function");
      return polynomial;
  }
  return polynomial;
}

// 1D linear
template Core::LinAlg::Matrix<2, 1> Core::LinAlg::ScalarInterpolation::polynomialshapefunction<1, 1,
    2>(const Core::LinAlg::Matrix<1, 1>&);

// 1D quadratic
template Core::LinAlg::Matrix<3, 1> Core::LinAlg::ScalarInterpolation::polynomialshapefunction<1, 2,
    3>(const Core::LinAlg::Matrix<1, 1>&);

// 2D linear
template Core::LinAlg::Matrix<4, 1> Core::LinAlg::ScalarInterpolation::polynomialshapefunction<2, 1,
    4>(const Core::LinAlg::Matrix<2, 1>&);

// 2D quadratic (incomplete)
template Core::LinAlg::Matrix<6, 1> Core::LinAlg::ScalarInterpolation::polynomialshapefunction<2, 2,
    6>(const Core::LinAlg::Matrix<2, 1>&);

// 2D quadratic (complete)
template Core::LinAlg::Matrix<8, 1> Core::LinAlg::ScalarInterpolation::polynomialshapefunction<2, 2,
    8>(const Core::LinAlg::Matrix<2, 1>&);

// 3D linear
template Core::LinAlg::Matrix<8, 1> Core::LinAlg::ScalarInterpolation::polynomialshapefunction<3, 1,
    8>(const Core::LinAlg::Matrix<3, 1>&);

// 3D quadratic
template Core::LinAlg::Matrix<10, 1> Core::LinAlg::ScalarInterpolation::polynomialshapefunction<3,
    2, 10>(const Core::LinAlg::Matrix<3, 1>&);


/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
template <unsigned int loc_dim, unsigned int poly_order, unsigned int num_cofeffcients>
Core::LinAlg::ScalarInterpolation::ScalarInterpolator<loc_dim, poly_order,
    num_cofeffcients>::ScalarInterpolator(const ScalarInterpolationType scalar_interp_type,
    const WeightingFunction weight_func, const InterpParams& interp_params)
    : scalar_interp_type_(scalar_interp_type),
      weight_func_(weight_func),
      interp_params_(interp_params)
{
  // check for valid polynomial order
  if (scalar_interp_type != ScalarInterpolationType::LOG and poly_order < 1)
    FOUR_C_THROW("ScalarInterpolator: for MLS and LOGMLS polynomial order must be at least 1!");

  // check for valid location dimension
  if (loc_dim < 1 || loc_dim > 3)
    FOUR_C_THROW("ScalarInterpolator: location dimension must be between 1 and 3!");
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
template <unsigned int loc_dim, unsigned int poly_order, unsigned int num_cofeffcients>
std::vector<double> Core::LinAlg::ScalarInterpolation::
    ScalarInterpolator<loc_dim, poly_order, num_cofeffcients>::calculate_normalized_weights(
        const std::vector<Core::LinAlg::Matrix<loc_dim, 1>>& ref_locs,
        const Core::LinAlg::Matrix<loc_dim, 1>& interp_loc, const WeightingFunction& weight_func,
        const InterpParams& interp_params)
{
  std::vector<double> weights(ref_locs.size(), 0.0);
  std::vector<Core::LinAlg::Matrix<loc_dim, 1>> ref_locs_tmp(ref_locs);

  for (auto& loc : ref_locs_tmp)
    loc.update(-1.0, interp_loc, 1.0);  // shift locations to the interpolation point

  // calculate the weight based on the chosen weighting function
  switch (weight_func)
  {
    case WeightingFunction::inversedistance:
    {
      double exponent = loc_dim;  // default exponent for inverse distance
      if (interp_params.inverse_distance_power.has_value())
        exponent = interp_params.inverse_distance_power.value();  // use the provided power

      for (unsigned int i = 0; i < ref_locs_tmp.size(); ++i)
      {
        double norm = ref_locs_tmp[i].norm2();
        if (norm < interp_params.distance_threshold)
          weights[i] = 1.0;
        else
          weights[i] = 1.0 / std::pow(norm, exponent);
      }
    }
    break;

    case WeightingFunction::exponential:
    {
      for (unsigned int i = 0; i < ref_locs_tmp.size(); ++i)
      {
        double norm = ref_locs_tmp[i].norm2();
        if (norm < interp_params.distance_threshold)
          weights[i] = 1.0;
        else
          weights[i] = std::exp(-1.0 * interp_params.exponential_decay_c * norm * norm);
      }
    }
    break;

    case WeightingFunction::unity:
    {
      for (unsigned int i = 0; i < ref_locs_tmp.size(); ++i)
        weights[i] = 1.0 / static_cast<double>(ref_locs_tmp.size());
      return weights;
    }
    break;

    default:
      FOUR_C_THROW("Unknown weighting function.");
  }

  for (unsigned int i = 0; i < weights.size(); ++i)
  {
    if (weights[i] == 1.0)
    {
      for (auto& w : weights) w = 0.0;
      weights[i] = 1.0;  // set the weight to 1.0 for the closest point
      return weights;
    }
  }

  // normalize the weights
  double sum_of_weights = 0.0;
  for (const auto& weight : weights) sum_of_weights += weight;
  for (double& weight : weights) weight /= sum_of_weights;  // normalize the weight

  // check for sum of weights = 1.0
  sum_of_weights = 0.0;
  for (const auto& weight : weights) sum_of_weights += weight;
  if (std::abs(sum_of_weights - 1.0) > 1e-8)
    FOUR_C_THROW("Sum of weights is not equal to 1.0 after normalization.");

  return weights;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
template <unsigned int loc_dim, unsigned int poly_order, unsigned int num_cofeffcients>
std::vector<double> Core::LinAlg::ScalarInterpolation::ScalarInterpolator<loc_dim, poly_order,
    num_cofeffcients>::logarthmic(const std::vector<std::vector<double>>& scalar_data,
    const std::vector<double>& weights)
{
  const size_t field_size = scalar_data[0].size();
  std::vector<double> log_sum(field_size, 0.0);

  for (size_t i = 0; i < scalar_data.size(); ++i)
  {
    const auto& row = scalar_data[i];
    if (row.size() != field_size) FOUR_C_THROW("Inconsistent vector sizes in scalar_data.");
    if (weights[i] < 0.0)
      FOUR_C_THROW("Weights must be non-negative for logarithmic interpolation.");

    for (size_t j = 0; j < field_size; ++j)
    {
      if (row[j] <= 0.0)
        FOUR_C_THROW("Logarithmic interpolation requires positive scalar data values.");
      log_sum[j] += weights[i] * std::log(row[j]);
    }
  }

  for (auto& val : log_sum) val = std::exp(val);

  return log_sum;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
template <unsigned int loc_dim, unsigned int poly_order, unsigned int num_cofeffcients>
std::vector<double> Core::LinAlg::ScalarInterpolation::ScalarInterpolator<loc_dim, poly_order,
    num_cofeffcients>::moving_least_square(const std::vector<std::vector<double>>& scalar_data,
    const std::vector<Core::LinAlg::Matrix<loc_dim, 1>>& ref_locs,
    const Core::LinAlg::Matrix<loc_dim, 1>& interp_loc, const std::vector<double>& weights)
{
  const size_t field_size = scalar_data[0].size();
  std::vector<double> interp_scalar(field_size, 0.0);

  std::vector<Core::LinAlg::Matrix<loc_dim, 1>> ref_locs_tmp(ref_locs);

  Core::LinAlg::Matrix<num_cofeffcients, num_cofeffcients> Pinv(Core::LinAlg::Initialization::zero);
  std::vector<Core::LinAlg::Matrix<num_cofeffcients, 1>> bj(field_size);
  for (unsigned int i = 0; i < ref_locs_tmp.size(); ++i)  // number of reference locations
  {
    ref_locs_tmp[i].update(-1.0, interp_loc, 1.0);  // shift locations to the interpolation point

    Core::LinAlg::Matrix<num_cofeffcients, 1> p_i =
        Core::LinAlg::ScalarInterpolation::polynomialshapefunction<loc_dim, poly_order,
            num_cofeffcients>(ref_locs_tmp[i]);

    // bj = w * p_i * scalar_data[i][j]
    for (unsigned int j = 0; j < field_size; ++j)  // number of scalar fields
      bj[j].update(weights[i] * scalar_data[i][j], p_i, 1.0);

    // P = P + w * p_i * p_i^T
    Pinv.multiply_nt(weights[i], p_i, p_i, 1.0);
  }

  // now evaluate P^-1 with solver
  Core::LinAlg::FixedSizeSerialDenseSolver<num_cofeffcients, num_cofeffcients, 1> solve_for_inverse;
  solve_for_inverse.set_matrix(Pinv);
  int err2 = solve_for_inverse.factor();
  int err = solve_for_inverse.invert();
  if ((err != 0) || (err2 != 0))
    FOUR_C_THROW(
        "ERROR: ScalarInterpolator::moving_least_square Inversion of P failed, check the "
        "dimension, poly order, and number of coefficients.");

  // now compute the interpolated scalar values
  for (unsigned int j = 0; j < field_size; ++j)  // number of scalar fields
  {
    Core::LinAlg::Matrix<num_cofeffcients, 1> yj(Core::LinAlg::Initialization::zero);
    yj.multiply(Pinv, bj[j]);  // compute the coeffiicents

    Core::LinAlg::Matrix<loc_dim, 1> interp_loc_shifted(Core::LinAlg::Initialization::zero);

    Core::LinAlg::Matrix<num_cofeffcients, 1> p_inter =
        Core::LinAlg::ScalarInterpolation::polynomialshapefunction<loc_dim, poly_order,
            num_cofeffcients>(interp_loc_shifted);

    Core::LinAlg::Matrix<1, 1> interp_value(Core::LinAlg::Initialization::zero);
    interp_value.multiply_tn(p_inter, yj);  // compute the interpolated value

    interp_scalar[j] = interp_value(0, 0);  // store the interpolated value
  }

  return interp_scalar;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
template <unsigned int loc_dim, unsigned int poly_order, unsigned int num_cofeffcients>
std::vector<double> Core::LinAlg::ScalarInterpolation::ScalarInterpolator<loc_dim, poly_order,
    num_cofeffcients>::log_moving_least_square(const std::vector<std::vector<double>>& scalar_data,
    const std::vector<Core::LinAlg::Matrix<loc_dim, 1>>& ref_locs,
    const Core::LinAlg::Matrix<loc_dim, 1>& interp_loc, const std::vector<double>& weights)
{
  std::vector<double> interp_scalar(scalar_data[0].size(), 0.0);

  std::vector<std::vector<double>> log_scalar_data(scalar_data.size());

  for (size_t i = 0; i < scalar_data.size(); ++i)
  {
    log_scalar_data[i].resize(scalar_data[i].size());
    for (size_t j = 0; j < scalar_data[i].size(); ++j)
    {
      if (scalar_data[i][j] <= 0.0)
        FOUR_C_THROW("Logarithmic moving least squares requires positive scalar data values.");
      log_scalar_data[i][j] = std::log(scalar_data[i][j]);
    }
  }

  // Use the moving least square method on the logarithmic data
  std::vector<double> log_interp_scalar =
      moving_least_square(log_scalar_data, ref_locs, interp_loc, weights);

  // Exponentiate the result to get back to the original scale
  for (size_t j = 0; j < log_interp_scalar.size(); ++j)
    interp_scalar[j] = std::exp(log_interp_scalar[j]);

  // check the interpolated scalar values
  for (size_t j = 0; j < interp_scalar.size(); ++j)
    if (interp_scalar[j] <= 0.0)
      FOUR_C_THROW("Logarithmic moving least squares resulted in non-positive scalar values.");

  // Return the interpolated scalar values
  return interp_scalar;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
template <unsigned int loc_dim, unsigned int poly_order, unsigned int num_cofeffcients>
std::vector<double> Core::LinAlg::ScalarInterpolation::ScalarInterpolator<loc_dim, poly_order,
    num_cofeffcients>::get_interpolated_scalar(const std::vector<std::vector<double>>& scalar_data,
    const std::vector<Core::LinAlg::Matrix<loc_dim, 1>>& ref_locs,
    const Core::LinAlg::Matrix<loc_dim, 1>& interp_loc)
{
  // Check for size consistency
  if (scalar_data.empty()) FOUR_C_THROW("Scalar data vector is empty.");
  if (ref_locs.empty()) FOUR_C_THROW("Reference locations vector is empty.");
  if (scalar_data.size() != ref_locs.size())
    FOUR_C_THROW("Size mismatch: scalar_data vs ref_locs.");
  if (ref_locs.size() < 1)
    FOUR_C_THROW("At least one reference location is required for interpolation.");

  auto weights =
      this->calculate_normalized_weights(ref_locs, interp_loc, weight_func_, interp_params_);

  // If one weight is 1.0, return the corresponding scalar data directly
  for (size_t i = 0; i < scalar_data.size(); ++i)
    if (weights[i] == 1.0) return scalar_data[i];

  switch (scalar_interp_type_)  // Choose the interpolation method based on the type
  {
    case Core::LinAlg::ScalarInterpolation::ScalarInterpolationType::LOG:
      return logarthmic(scalar_data, weights);
      break;
    case Core::LinAlg::ScalarInterpolation::ScalarInterpolationType::MLS:
      return moving_least_square(scalar_data, ref_locs, interp_loc, weights);
      break;
    case Core::LinAlg::ScalarInterpolation::ScalarInterpolationType::LOGMLS:
      return log_moving_least_square(scalar_data, ref_locs, interp_loc, weights);
      break;
    default:
      FOUR_C_THROW("Unknown scalar interpolation type.");
      break;
  }
}

// Explicit template instantiations
template class Core::LinAlg::ScalarInterpolation::ScalarInterpolator<1>;
template class Core::LinAlg::ScalarInterpolation::ScalarInterpolator<2>;
template class Core::LinAlg::ScalarInterpolation::ScalarInterpolator<3>;

template class Core::LinAlg::ScalarInterpolation::ScalarInterpolator<1, 1, 2>;
template class Core::LinAlg::ScalarInterpolation::ScalarInterpolator<1, 2, 3>;

template class Core::LinAlg::ScalarInterpolation::ScalarInterpolator<2, 1, 4>;
template class Core::LinAlg::ScalarInterpolation::ScalarInterpolator<2, 2, 6>;
template class Core::LinAlg::ScalarInterpolation::ScalarInterpolator<2, 2, 8>;

template class Core::LinAlg::ScalarInterpolation::ScalarInterpolator<3, 1, 8>;
template class Core::LinAlg::ScalarInterpolation::ScalarInterpolator<3, 2, 10>;

FOUR_C_NAMESPACE_CLOSE