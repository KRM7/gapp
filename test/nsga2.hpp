/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef NSGA2_BENCHMARK_HPP
#define NSGA2_BENCHMARK_HPP

#include "../src/encoding/real.hpp"
#include "../src/algorithm/algorithm.hpp"
#include "../src/crossover/real.hpp"
#include "../src/mutation/real.hpp"
#include "fitness_functions.hpp"
#include "benchmark_utils.hpp"

using namespace genetic_algorithm;

void nsga2_kur()
{
    KUR fitness_func(3);

    RCGA GA(fitness_func.num_vars, fitness_func, fitness_func.bounds());

    GA.population_size(100);
    GA.algorithm(algorithm::NSGA2{});
    GA.crossover_method(crossover::real::SimulatedBinary{ 0.8 });
    GA.mutation_method(mutation::real::Gauss{ 1.0 / fitness_func.bounds().size() });

    benchmarkMoga(GA, 250, "NSGA2", "KUR");
}

void nsga2_zdt2()
{
    ZDT2 fitness_func(30);
    
    RCGA GA(fitness_func.num_vars, fitness_func, fitness_func.bounds());

    GA.population_size(100);
    GA.algorithm(algorithm::NSGA2{});
    GA.crossover_method(crossover::real::SimulatedBinary{ 0.8 });
    GA.mutation_method(mutation::real::Gauss{ 1.0 / fitness_func.bounds().size() });

    benchmarkMoga(GA, 250, "NSGA2", "ZDT2");
}

void nsga2_zdt3()
{
    ZDT3 fitness_func(30);

    RCGA GA(fitness_func.num_vars, fitness_func, fitness_func.bounds());

    GA.population_size(100);
    GA.algorithm(algorithm::NSGA2{});
    GA.crossover_method(crossover::real::SimulatedBinary{ 0.8 });
    GA.mutation_method(mutation::real::Gauss{ 1.0 / fitness_func.bounds().size() });

    benchmarkMoga(GA, 250, "NSGA2", "ZDT3");
}

void nsga2_zdt6()
{
    ZDT6 fitness_func(10);
    
    RCGA GA(fitness_func.num_vars, fitness_func, fitness_func.bounds());

    GA.population_size(100);
    GA.algorithm(algorithm::NSGA2{});
    GA.crossover_method(crossover::real::SimulatedBinary{ 0.8 });
    GA.mutation_method(mutation::real::Gauss{ 1.0 / fitness_func.bounds().size() });

    benchmarkMoga(GA, 250, "NSGA2", "ZDT6");
}

void nsga2_dtlz1()
{
    DTLZ1 fitness_func(7, 3);

    RCGA GA(fitness_func.num_vars, fitness_func, fitness_func.bounds());

    GA.population_size(100);
    GA.algorithm(algorithm::NSGA2{});
    GA.crossover_method(crossover::real::SimulatedBinary{ 0.9, 15.0 });
    GA.mutation_method(mutation::real::Uniform{ 1.0 / fitness_func.bounds().size() });

    benchmarkMoga(GA, 1500, "NSGA2", "DTLZ1");
}

void nsga2_dtlz2()
{
    DTLZ2 fitness_func(12, 3);

    RCGA GA(fitness_func.num_vars, fitness_func, fitness_func.bounds());

    GA.population_size(100);
    GA.algorithm(algorithm::NSGA2{});
    GA.crossover_method(crossover::real::SimulatedBinary{ 0.9, 15.0 });
    GA.mutation_method(mutation::real::Uniform{ 1.0 / fitness_func.bounds().size() });

    benchmarkMoga(GA, 1500, "NSGA2", "DTLZ2");
}


#endif // !NSGA2_BENCHMARK_HPP