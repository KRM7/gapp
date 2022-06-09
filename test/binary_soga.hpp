/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef BINARY_SOGA_BENCHMARK_HPP
#define BINARY_SOGA_BENCHMARK_HPP

#include "../src/encoding/binary.hpp"
#include "../src/selection/selection.hpp"
#include "../src/crossover/binary.hpp"
#include "../src/mutation/binary.hpp"
#include "../src/stop_condition/stop_condition.hpp"
#include "fitness_functions.hpp"
#include "benchmark_utils.hpp"

using namespace genetic_algorithm;

void binary_rastrigin()
{
    Rastrigin fitness_func(10);
    BinaryGA GA(400, fitness_func.num_vars * fitness_func.var_bits, fitness_func);
    
    GA.selection_method(selection::single_objective::Roulette{});
    GA.crossover_method(crossover::binary::TwoPoint{ 0.75 });
    GA.mutation_method(mutation::binary::Flip{ 0.015 });
    GA.stop_condition(stopping::AND(stopping::FitnessMeanStall{ 50, 0.005 },
                                    stopping::FitnessBestStall{ 50, 0.005 }));

    benchmarkSoga(GA, 1000, fitness_func, "Rastrigin function");
}

void binary_rosenbrock()
{
    Rosenbrock fitness_func(10);
    BinaryGA GA(300, fitness_func.num_vars * fitness_func.var_bits, fitness_func);

    GA.selection_method(selection::single_objective::Tournament{});
    GA.crossover_method(crossover::binary::TwoPoint{ 0.8 });
    GA.mutation_method(mutation::binary::Flip{ 0.01 });

    benchmarkSoga(GA, 1500, fitness_func, "Rosenbrock function");
}

void binary_schwefel()
{
    Schwefel fitness_func(10);
    BinaryGA GA(200, fitness_func.num_vars * fitness_func.var_bits, fitness_func);

    GA.selection_method(selection::single_objective::Rank{});
    GA.crossover_method(crossover::binary::Uniform{ 0.7 });
    GA.mutation_method(mutation::binary::Flip{ 0.01 });
    GA.stop_condition(stopping::FitnessEvals{ 200 * 1000 });

    benchmarkSoga(GA, 1500, fitness_func, "Schwefel function");
}

void binary_griewank()
{
    Griewank fitness_func(10);
    BinaryGA GA(250, fitness_func.num_vars * fitness_func.var_bits, fitness_func);

    GA.selection_method(selection::single_objective::Sigma{});
    GA.crossover_method(crossover::binary::TwoPoint{ 0.75 });
    GA.mutation_method(mutation::binary::Flip{ 0.04 });
    GA.stop_condition(stopping::FitnessValue{ { -0.1 } });

    benchmarkSoga(GA, 2500, fitness_func, "Griewank function");
}

void binary_ackley()
{
    Ackley fitness_func(10);
    BinaryGA GA(250, fitness_func.num_vars * fitness_func.var_bits, fitness_func);

    GA.selection_method(selection::single_objective::Boltzmann{});
    GA.crossover_method(crossover::binary::SinglePoint{ 0.75 });
    GA.mutation_method(mutation::binary::Flip{ 0.04 });
    GA.stop_condition(stopping::FitnessBestStall{ 50, 0.002 });

    benchmarkSoga(GA, 2500, fitness_func, "Ackley function");
}

#endif // !BINARY_SOGA_BENCHMARK_HPP