/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef REAL_SOGA_BENCHMARK_HPP
#define REAL_SOGA_BENCHMARK_HPP

#include "../src/encoding/real.hpp"
#include "../src/algorithm/algorithm.hpp"
#include "../src/crossover/real.hpp"
#include "../src/mutation/real.hpp"
#include "../src/stop_condition/stop_condition.hpp"
#include "fitness_functions.hpp"
#include "benchmark_utils.hpp"

using namespace genetic_algorithm;

void real_rastrigin()
{
    Rastrigin fitness_func(10);

    RCGA GA(100, fitness_func.num_vars, fitness_func, fitness_func.bounds());

    GA.algorithm(selection_::Roulette{});
    GA.crossover_method(crossover::real::SimulatedBinary{ 0.6, 2.0 });
    GA.mutation_method(mutation::real::Gauss{ 0.05 });
    GA.stop_condition(stopping::FitnessValue{ { -0.01 } });

    benchmarkSoga(GA, 1000, fitness_func, "Rastrigin function");
}

void real_rosenbrock()
{
    Rosenbrock fitness_func(10);

    RCGA GA(500, fitness_func.num_vars, fitness_func, fitness_func.bounds());

    GA.algorithm(selection_::Tournament{});
    GA.crossover_method(crossover::real::BLXa{ 0.9 });
    GA.mutation_method(mutation::real::Uniform{ 1.0 / fitness_func.num_vars });
    GA.stop_condition(stopping::FitnessEvals{ 500 * 1000 });

    benchmarkSoga(GA, 2000, fitness_func, "Rosenbrock function");
}

void real_schwefel()
{
    Schwefel fitness_func(10);

    RCGA GA(500, fitness_func.num_vars, fitness_func, fitness_func.bounds());

    GA.algorithm(selection_::Sigma{});
    GA.crossover_method(crossover::real::BLXa{ 0.7 });
    GA.mutation_method(mutation::real::NonUniform{ 1.0 / fitness_func.num_vars });
    GA.stop_condition(stopping::FitnessMeanStall{ 75, 0.01 });

    benchmarkSoga(GA, 1000, fitness_func, "Schwefel function");
}

void real_griewank()
{
    Griewank fitness_func(10);

    RCGA GA(200, fitness_func.num_vars, fitness_func, fitness_func.bounds());

    GA.algorithm(selection_::Boltzmann{});
    GA.crossover_method(crossover::real::Wright{ 0.85 });
    GA.mutation_method(mutation::real::Gauss{ 0.05 });

    benchmarkSoga(GA, 1500, fitness_func, "Griewank function");
}

void real_ackley()
{
    Ackley fitness_func(10);

    RCGA GA(200, fitness_func.num_vars, fitness_func, fitness_func.bounds());

    GA.algorithm(selection_::Boltzmann{});
    GA.crossover_method(crossover::real::Arithmetic{ 0.85 });
    GA.mutation_method(mutation::real::Polynomial{ 1.0 / fitness_func.num_vars, 60.0 });
    GA.stop_condition(stopping::FitnessBestStall{ 75, 0.002 });

    benchmarkSoga(GA, 1000, fitness_func, "Ackley function");
}

#endif // !REAL_BENCHMARK_HPP