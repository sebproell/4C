// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_linalg_tensor.hpp"

#include <benchmark/benchmark.h>

FOUR_C_NAMESPACE_OPEN
template <std::size_t... size>
Core::LinAlg::Tensor<double, size...> get_random_tensor()
{
  Core::LinAlg::Tensor<double, size...> tensor;
  std::generate(tensor.container().begin(), tensor.container().end(),
      []() { return static_cast<double>(std::rand()) / RAND_MAX; });
  return tensor;
}

template <std::size_t size>
static void matrix_matrix_product(benchmark::State& state)
{
  Core::LinAlg::Tensor<double, size, size> A = get_random_tensor<size, size>();
  Core::LinAlg::Tensor<double, size, size> B = get_random_tensor<size, size>();
  for (auto _ : state)
  {
    Core::LinAlg::Tensor<double, size, size> result = A * B;
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(matrix_matrix_product<2>);
BENCHMARK(matrix_matrix_product<3>);
BENCHMARK(matrix_matrix_product<4>);


template <std::size_t size>
static void matrix_matrix_t_product(benchmark::State& state)
{
  Core::LinAlg::Tensor<double, size, size> A = get_random_tensor<size, size>();
  Core::LinAlg::Tensor<double, size, size> B = get_random_tensor<size, size>();
  for (auto _ : state)
  {
    Core::LinAlg::Tensor<double, size, size> result = A * Core::LinAlg::transpose(B);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(matrix_matrix_t_product<2>);
BENCHMARK(matrix_matrix_t_product<3>);
BENCHMARK(matrix_matrix_t_product<4>);


template <std::size_t size>
static void matrix_vector_product(benchmark::State& state)
{
  Core::LinAlg::Tensor<double, size, size> A = get_random_tensor<size, size>();
  Core::LinAlg::Tensor<double, size> b = get_random_tensor<size>();
  for (auto _ : state)
  {
    Core::LinAlg::Tensor<double, size> result = A * b;
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(matrix_vector_product<2>);
BENCHMARK(matrix_vector_product<3>);
BENCHMARK(matrix_vector_product<4>);


template <std::size_t size>
static void scaled_matrix_matrix_product(benchmark::State& state)
{
  Core::LinAlg::Tensor<double, size, size> A = get_random_tensor<size, size>();
  Core::LinAlg::Tensor<double, size, size> B = get_random_tensor<size, size>();
  for (auto _ : state)
  {
    Core::LinAlg::Tensor<double, size, size> result = 3.0 * A * B;
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(scaled_matrix_matrix_product<2>);
BENCHMARK(scaled_matrix_matrix_product<3>);
BENCHMARK(scaled_matrix_matrix_product<4>);

FOUR_C_NAMESPACE_CLOSE