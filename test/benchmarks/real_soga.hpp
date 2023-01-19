/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_TEST_BENCHMARK_RCGA_HPP
#define GA_TEST_BENCHMARK_RCGA_HPP

#include "encoding/real.hpp"
#include "algorithm/algorithm.hpp"
#include "crossover/real.hpp"
#include "mutation/real.hpp"
#include "stop_condition/stop_condition.hpp"
#include "problems/single_objective.hpp"
#include "benchmark_utils.hpp"

using namespace genetic_algorithm;
using namespace genetic_algorithm::problems;

void real_sphere()
{
    Sphere fitness_func(10);

    RCGA GA(200, fitness_func.num_vars(), fitness_func, fitness_func.bounds());

    GA.algorithm(algorithm::SingleObjective{ selection::Tournament{} });
    GA.crossover_method(crossover::real::Arithmetic{ 0.8 });
    GA.mutation_method(mutation::real::Boundary{ 0.05 });
    GA.stop_condition(stopping::FitnessValue{ { -1E-12 } });

    benchmarkSoga(GA, 1000, fitness_func);
}

void real_rastrigin()
{
    Rastrigin fitness_func(10);

    RCGA GA(100, fitness_func.num_vars(), fitness_func, fitness_func.bounds());

    GA.algorithm(algorithm::SingleObjective{ selection::Roulette{} });
    GA.crossover_method(crossover::real::SimulatedBinary{ 0.6, 2.0 });
    GA.mutation_method(mutation::real::Gauss{ 0.05 });
    GA.stop_condition(stopping::FitnessValue{ { -0.01 } });

    benchmarkSoga(GA, 1000, fitness_func);
}

void real_rosenbrock()
{
    Rosenbrock fitness_func(10);

    RCGA GA(500, fitness_func.num_vars(), fitness_func, fitness_func.bounds());

    GA.algorithm(algorithm::SingleObjective{ selection::Tournament{} });
    GA.crossover_method(crossover::real::BLXa{ 0.9 });
    GA.mutation_method(mutation::real::Uniform{ 1.0 / fitness_func.num_vars() });
    GA.stop_condition(stopping::FitnessEvals{ 500 * 1000 });

    benchmarkSoga(GA, 2000, fitness_func);
}

void real_schwefel()
{
    Schwefel fitness_func(10);

    RCGA GA(500, fitness_func.num_vars(), fitness_func, fitness_func.bounds());

    GA.algorithm(algorithm::SingleObjective{ selection::Sigma{} });
    GA.crossover_method(crossover::real::BLXa{ 0.7 });
    GA.mutation_method(mutation::real::NonUniform{ 1.0 / fitness_func.num_vars() });
    GA.stop_condition(stopping::FitnessMeanStall{ 75, 0.01 });

    benchmarkSoga(GA, 1000, fitness_func);
}

void real_griewank()
{
    Griewank fitness_func(10);

    RCGA GA(200, fitness_func.num_vars(), fitness_func, fitness_func.bounds());

    GA.algorithm(algorithm::SingleObjective{ selection::Boltzmann{} });
    GA.crossover_method(crossover::real::Wright{ 0.8 });
    GA.mutation_method(mutation::real::Gauss{ 0.5 / fitness_func.num_vars() });

    benchmarkSoga(GA, 1500, fitness_func);
}

void real_ackley()
{
    Ackley fitness_func(10);

    RCGA GA(200, fitness_func.num_vars(), fitness_func, fitness_func.bounds());

    GA.algorithm(algorithm::SingleObjective{ selection::Boltzmann{} });
    GA.crossover_method(crossover::real::Arithmetic{ 0.85 });
    GA.mutation_method(mutation::real::Polynomial{ 1.0 / fitness_func.num_vars(), 60.0 });
    GA.stop_condition(stopping::FitnessBestStall{ 75, 0.002 });

    benchmarkSoga(GA, 1000, fitness_func);
}

void real_levy()
{
    Levy fitness_func(10);

    RCGA GA(200, fitness_func.num_vars(), fitness_func, fitness_func.bounds());

    GA.algorithm(algorithm::SingleObjective{ selection::Boltzmann{} });
    GA.crossover_method(crossover::real::Wright{ 0.85 });
    GA.mutation_method(mutation::real::NonUniform{ 1.0 / fitness_func.num_vars() });
    GA.stop_condition(stopping::FitnessValue{ { -1E-8 } });

    benchmarkSoga(GA, 1500, fitness_func);
}

#endif // !GA_TEST_BENCHMARK_RCGA_HPP