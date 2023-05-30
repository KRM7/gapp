/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_TEST_BENCHMARK_BINARY_HPP
#define GA_TEST_BENCHMARK_BINARY_HPP

#include "encoding/binary.hpp"
#include "algorithm/algorithm.hpp"
#include "crossover/binary.hpp"
#include "mutation/binary.hpp"
#include "stop_condition/stop_condition.hpp"
#include "stop_condition/composite.hpp"
#include "problems/single_objective.hpp"
#include "benchmark_utils.hpp"

using namespace gapp;
using namespace gapp::problems;

inline void binary_sphere()
{
    BinaryGA GA(200);

    GA.algorithm(algorithm::SingleObjective{ selection::Sigma{} });
    GA.crossover_method(crossover::binary::SinglePoint{ 0.9 });
    GA.mutation_method(mutation::binary::Flip{ 0.001 });
    GA.stop_condition(stopping::FitnessValue{ { -1E-12 } });

    benchmarkSoga(GA, 1000, Sphere{ 100 });
}

inline void binary_rastrigin()
{
    BinaryGA GA(400);

    GA.algorithm(algorithm::SingleObjective{ selection::Roulette{} });
    GA.crossover_method(crossover::binary::NPoint{ 0.75, 2 });
    GA.mutation_method(mutation::binary::Flip{ 0.015 });
    GA.stop_condition(stopping::FitnessMeanStall{ 50, 0.005 } && stopping::FitnessBestStall{ 50, 0.005 });

    benchmarkSoga(GA, 1000, Rastrigin{ 10 });
}

inline void binary_rosenbrock()
{
    BinaryGA GA(300);

    GA.algorithm(algorithm::SingleObjective{ selection::Tournament{}, replacement::KeepChildren{} });
    GA.crossover_method(crossover::binary::TwoPoint{ 0.8 });
    GA.mutation_method(mutation::binary::Flip{ 0.01 });

    benchmarkSoga(GA, 1500, Rosenbrock{ 10 });
}

inline void binary_schwefel()
{
    BinaryGA GA(200);

    GA.algorithm(algorithm::SingleObjective{ selection::Rank{}, replacement::Elitism{ 10 } });
    GA.crossover_method(crossover::binary::Uniform{ 0.7 });
    GA.mutation_method(mutation::binary::Flip{ 0.01 });
    GA.stop_condition(stopping::FitnessEvals{ 200 * 1000 });

    benchmarkSoga(GA, 1500, Schwefel{ 10 });
}

inline void binary_griewank()
{
    BinaryGA GA(200);

    GA.algorithm(algorithm::SingleObjective{ selection::Sigma{}, replacement::KeepBest{} });
    GA.crossover_method(crossover::binary::TwoPoint{ 0.8 });
    GA.mutation_method(mutation::binary::Flip{ 0.02 });
    GA.stop_condition(stopping::FitnessValue{ { -0.01 } });

    benchmarkSoga(GA, 1500, Griewank{ 10 });
}

inline void binary_ackley()
{
    BinaryGA GA(250);

    GA.algorithm(algorithm::SingleObjective{ selection::Boltzmann{}, replacement::KeepBest{} });
    GA.crossover_method(crossover::binary::SinglePoint{ 0.75 });
    GA.mutation_method(mutation::binary::Flip{ 0.01 });
    GA.stop_condition(stopping::FitnessBestStall{ 50, 0.002 });

    benchmarkSoga(GA, 2500, Ackley{ 10 });
}

inline void binary_levy()
{
    BinaryGA GA(250);

    GA.algorithm(algorithm::SingleObjective{ selection::Boltzmann{}, replacement::KeepBest{} });
    GA.crossover_method(crossover::binary::TwoPoint{ 0.8 });
    GA.mutation_method(mutation::binary::Flip{ 0.03 });

    benchmarkSoga(GA, 1500, Levy{ 10 });
}

#endif // !GA_TEST_BENCHMARK_BINARY_HPP