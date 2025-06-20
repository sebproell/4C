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
    std::map<std::string, double> variables;

    for (auto _ : state)
    {
      double result = expr.value(variables);
      benchmark::DoNotOptimize(result);
    }
  }
  BENCHMARK(symbolic_expression_constant);

  void symbolic_expression_basic(benchmark::State& state)
  {
    SymbolicExpression<double> expr("2*x^2 + y + 4*z");
    std::map<std::string, double> variables = {{"x", 1.0}, {"y", 2.0}, {"z", 3.0}};

    for (auto _ : state)
    {
      double result = expr.value(variables);
      benchmark::DoNotOptimize(result);
    }
  }
  BENCHMARK(symbolic_expression_basic);


  void symbolic_expression_first_derivative(benchmark::State& state)
  {
    SymbolicExpression<double> expr("2*x^2 + y + 4*z");
    std::map<std::string, FirstDerivativeType> variables = {
        {"x", FirstDerivativeType(2, 0, 1.0)},
        {"y", FirstDerivativeType(2, 1, 2.0)},
    };
    std::map<std::string, double> constants = {{"z", 3.0}};

    for (auto _ : state)
    {
      auto result = expr.first_derivative(variables, constants);
      benchmark::DoNotOptimize(result);
    }
  }
  BENCHMARK(symbolic_expression_first_derivative);


  void symbolic_expression_second_derivative(benchmark::State& state)
  {
    SymbolicExpression<double> expr("2*x^2 + y + 4*z");
    std::map<std::string, SecondDerivativeType> variables = {
        {"x", SecondDerivativeType(2, 0, 1.0)},
        {"y", SecondDerivativeType(2, 1, 2.0)},
    };
    std::map<std::string, double> constants = {{"z", 3.0}};

    for (auto _ : state)
    {
      auto result = expr.second_derivative(variables, constants);
      benchmark::DoNotOptimize(result);
    }
  }
  BENCHMARK(symbolic_expression_second_derivative);


  void symbolic_expression_functions(benchmark::State& state)
  {
    SymbolicExpression<double> expr("sin(x) + cos(y) * exp(z) - log10(x + y)");
    std::map<std::string, double> variables = {{"x", 1.0}, {"y", 2.0}, {"z", 3.0}};

    for (auto _ : state)
    {
      double result = expr.value(variables);
      benchmark::DoNotOptimize(result);
    }
  }
  BENCHMARK(symbolic_expression_functions);

}  // namespace
