/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_TEST_BENCHMARK_BINARY_HPP
#define GA_TEST_BENCHMARK_BINARY_HPP

#include "encoding/binary.hpp"
#include "algorithm/algorithm.hpp"
#include "crossover/binary.hpp"
#include "mutation/binary.hpp"
#include "stop_condition/stop_condition.hpp"
#include "problems/single_objective.hpp"
#include "benchmark_utils.hpp"

using namespace genetic_algorithm;
using namespace genetic_algorithm::problems;

void binary_sphere()
{
    Sphere fitness_func(100);
    BinaryGA GA(200, fitness_func.num_bits(), fitness_func);

    GA.algorithm(algorithm::SingleObjective{ selection::Sigma{} });
    GA.crossover_method(crossover::binary::SinglePoint{ 0.9 });
    GA.mutation_method(mutation::binary::Flip{ 0.001 });
    GA.stop_condition(stopping::FitnessValue{ { -1E-12 } });

    benchmarkSoga(GA, 1000, fitness_func);
}

void binary_rastrigin()
{
    Rastrigin fitness_func(10);
    BinaryGA GA(400, fitness_func.num_bits(), fitness_func);

    GA.algorithm(algorithm::SingleObjective{ selection::Roulette{} });
    GA.crossover_method(crossover::binary::NPoint{ 0.75, 2 });
    GA.mutation_method(mutation::binary::Flip{ 0.015 });
    GA.stop_condition(stopping::AND(stopping::FitnessMeanStall{ 50, 0.005 },
                                    stopping::FitnessBestStall{ 50, 0.005 }));

    benchmarkSoga(GA, 1000, fitness_func);
}

void binary_rosenbrock()
{
    Rosenbrock fitness_func(10);
    BinaryGA GA(300, fitness_func.num_bits(), fitness_func);

    GA.algorithm(algorithm::SingleObjective{ selection::Tournament{}, update::KeepChildren{} });
    GA.crossover_method(crossover::binary::TwoPoint{ 0.8 });
    GA.mutation_method(mutation::binary::Flip{ 0.01 });

    benchmarkSoga(GA, 1500, fitness_func);
}

void binary_schwefel()
{
    Schwefel fitness_func(10);
    BinaryGA GA(200, fitness_func.num_bits(), fitness_func);

    GA.algorithm(algorithm::SingleObjective{ selection::Rank{}, update::Elitism{ 10 } });
    GA.crossover_method(crossover::binary::Uniform{ 0.7 });
    GA.mutation_method(mutation::binary::Flip{ 0.01 });
    GA.stop_condition(stopping::FitnessEvals{ 200 * 1000 });

    benchmarkSoga(GA, 1500, fitness_func);
}

void binary_griewank()
{
    Griewank fitness_func(10);
    BinaryGA GA(200, fitness_func.num_bits(), fitness_func);

    GA.algorithm(algorithm::SingleObjective{ selection::Sigma{}, update::KeepBest{} });
    GA.crossover_method(crossover::binary::TwoPoint{ 0.8 });
    GA.mutation_method(mutation::binary::Flip{ 0.2 / fitness_func.num_vars() });
    GA.stop_condition(stopping::FitnessValue{ { -0.01 } });

    benchmarkSoga(GA, 1500, fitness_func);
}

void binary_ackley()
{
    Ackley fitness_func(10);
    BinaryGA GA(250, fitness_func.num_bits(), fitness_func);

    GA.algorithm(algorithm::SingleObjective{ selection::Boltzmann{}, update::KeepBest{} });
    GA.crossover_method(crossover::binary::SinglePoint{ 0.75 });
    GA.mutation_method(mutation::binary::Flip{ 0.01 });
    GA.stop_condition(stopping::FitnessBestStall{ 50, 0.002 });

    benchmarkSoga(GA, 2500, fitness_func);
}

void binary_levy()
{
    Levy fitness_func(10);
    BinaryGA GA(250, fitness_func.num_bits(), fitness_func);

    GA.algorithm(algorithm::SingleObjective{ selection::Boltzmann{}, update::KeepBest{} });
    GA.crossover_method(crossover::binary::TwoPoint{ 0.8 });
    GA.mutation_method(mutation::binary::Flip{ 0.03 });

    benchmarkSoga(GA, 1500, fitness_func);
}

#endif // !GA_TEST_BENCHMARK_BINARY_HPP