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

using namespace gapp;
using namespace gapp::problems;

inline void real_sphere()
{
    RCGA GA(200);

    GA.algorithm(algorithm::SingleObjective{ selection::Tournament{} });
    GA.crossover_method(crossover::real::Arithmetic{ 0.8 });
    GA.mutation_method(mutation::real::Boundary{ 0.05 });
    GA.stop_condition(stopping::FitnessValue{ { -1E-12 } });

    benchmarkSoga(GA, 1000, Sphere{ 10 });
}

inline void real_rastrigin()
{
    RCGA GA(100);

    GA.algorithm(algorithm::SingleObjective{ selection::Roulette{} });
    GA.crossover_method(crossover::real::SimulatedBinary{ 0.6, 2.0 });
    GA.mutation_method(mutation::real::Gauss{ 0.05 });
    GA.stop_condition(stopping::FitnessValue{ { -0.01 } });

    benchmarkSoga(GA, 1000, Rastrigin{ 10 });
}

inline void real_rosenbrock()
{
    RCGA GA(500);

    GA.algorithm(algorithm::SingleObjective{ selection::Tournament{} });
    GA.crossover_method(crossover::real::BLXa{ 0.9 });
    GA.mutation_method(mutation::real::Uniform{ 0.1 });
    GA.stop_condition(stopping::FitnessEvals{ 500 * 1000 });

    benchmarkSoga(GA, 2000, Rosenbrock{ 10 });
}

inline void real_schwefel()
{
    RCGA GA(500);

    GA.algorithm(algorithm::SingleObjective{ selection::Sigma{} });
    GA.crossover_method(crossover::real::BLXa{ 0.7 });
    GA.mutation_method(mutation::real::NonUniform{ 0.1 });
    GA.stop_condition(stopping::FitnessMeanStall{ 75, 0.01 });

    benchmarkSoga(GA, 1000, Schwefel{ 10 });
}

inline void real_griewank()
{
    RCGA GA(200);

    GA.algorithm(algorithm::SingleObjective{ selection::Boltzmann{} });
    GA.crossover_method(crossover::real::Wright{ 0.8 });
    GA.mutation_method(mutation::real::Gauss{ 0.05 });

    benchmarkSoga(GA, 1500, Griewank{ 10 });
}

inline void real_ackley()
{
    RCGA GA(200);

    GA.algorithm(algorithm::SingleObjective{ selection::Boltzmann{} });
    GA.crossover_method(crossover::real::Arithmetic{ 0.85 });
    GA.mutation_method(mutation::real::Polynomial{ 0.1, 60.0 });
    GA.stop_condition(stopping::FitnessBestStall{ 75, 0.002 });

    benchmarkSoga(GA, 1000, Ackley{ 10 });
}

inline void real_levy()
{
    RCGA GA(200);

    GA.algorithm(algorithm::SingleObjective{ selection::Boltzmann{} });
    GA.crossover_method(crossover::real::Wright{ 0.85 });
    GA.mutation_method(mutation::real::NonUniform{ 0.1 });
    GA.stop_condition(stopping::FitnessValue{ { -1E-8 } });

    benchmarkSoga(GA, 1500, Levy{ 10 });
}

#endif // !GA_TEST_BENCHMARK_RCGA_HPP