// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_utils_symbolic_expression.hpp"

#include <benchmark/benchmark.h>

namespace
{
  using namespace FourC::Core::Utils;

  using FirstDerivativeType = Sacado::Fad::DFad<double>;
  using SecondDerivativeType = Sacado::Fad::DFad<Sacado::Fad::DFad<double>>;

  void symbolic_expression_constant(benchmark::State& state)
  {
    // A baseline benchmark which has no variables and only evaluates a constant expression.
    SymbolicExpression<double> expr("3.14");

    for (auto _ : state)
    {
      double result = expr.value();
      benchmark::DoNotOptimize(result);
    }
  }
  BENCHMARK(symbolic_expression_constant);

  void symbolic_expression_constant_folding(benchmark::State& state)
  {
    // A benchmark that evaluates a expression which can be simplified to a constant.
    // Should run as fast as the constant benchmark above.
    SymbolicExpression<double> expr("2*1.0^2 + 2.0 + 4*3.0");
    for (auto _ : state)
    {
      double result = expr.value();
      benchmark::DoNotOptimize(result);
    }
  }
  BENCHMARK(symbolic_expression_constant_folding);

  void symbolic_expression_basic(benchmark::State& state)
  {
    // Only basic arithmetic, including expensive division.
    SymbolicExpression<double> expr("2*x + y/2 - 4*z");

    for (auto _ : state)
    {
      double result = expr.value("x", 1.0, "y", 2.0, "z", 3.0);
      benchmark::DoNotOptimize(result);
    }
  }
  BENCHMARK(symbolic_expression_basic);



  void symbolic_expression_basic_native(benchmark::State& state)
  {
    // Execute the same code at the natively compiled level to compare performance.
    volatile double x = 1.0;
    volatile double y = 2.0;
    volatile double z = 3.0;

    for (auto _ : state)
    {
      volatile double result = 2 * x + y / 2 - 4 * z;
      benchmark::DoNotOptimize(result);
    }
  }
  BENCHMARK(symbolic_expression_basic_native);


  void symbolic_expression_first_derivative(benchmark::State& state)
  {
    SymbolicExpression<double> expr("2*x^2 + y + 4*z");

    for (auto _ : state)
    {
      auto result = expr.first_derivative(
          "x", FirstDerivativeType(2, 0, 1.0), "y", FirstDerivativeType(2, 1, 2.0), "z", 3.0);
      benchmark::DoNotOptimize(result);
    }
  }
  BENCHMARK(symbolic_expression_first_derivative);


  void symbolic_expression_second_derivative(benchmark::State& state)
  {
    SymbolicExpression<double> expr("2*x^2 + y + 4*z");

    for (auto _ : state)
    {
      auto result = expr.second_derivative(
          "x", SecondDerivativeType(2, 0, 1.0), "y", SecondDerivativeType(2, 1, 2.0), "z", 3.0);
      benchmark::DoNotOptimize(result);
    }
  }
  BENCHMARK(symbolic_expression_second_derivative);


  void symbolic_expression_functions(benchmark::State& state)
  {
    SymbolicExpression<double> expr("sin(x) + cos(y) * exp(z) - log10(x + y)");

    for (auto _ : state)
    {
      double result = expr.value("x", 1.0, "y", 2.0, "z", 3.0);
      benchmark::DoNotOptimize(result);
    }
  }
  BENCHMARK(symbolic_expression_functions);

}  // namespace
