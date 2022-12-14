/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_TEST_BENCHMARK_NSGA3_HPP
#define GA_TEST_BENCHMARK_NSGA3_HPP

#include "encoding/real.hpp"
#include "encoding/binary.hpp"
#include "algorithm/algorithm.hpp"
#include "crossover/real.hpp"
#include "crossover/binary.hpp"
#include "mutation/real.hpp"
#include "mutation/binary.hpp"
#include "stop_condition/stop_condition.hpp"
#include "problems/many_objective.hpp"
#include "benchmark_utils.hpp"

using namespace genetic_algorithm;
using namespace genetic_algorithm::problems;

template<typename Problem>
void benchmark_real_nsga3(const Problem& problem, size_t generations, size_t population_size = 100)
{
    RCGA GA(population_size, problem.num_vars(), problem, problem.bounds());

    GA.algorithm(algorithm::NSGA3{});
    GA.crossover_method(crossover::real::SimulatedBinary{ 0.9 });
    GA.mutation_method(mutation::real::Uniform{ 1.0 / problem.num_vars() });

    benchmarkMoga(GA, generations, "NSGA3", problem);
}

template<typename Problem>
void benchmark_binary_nsga3(const Problem& problem, size_t generations, size_t population_size = 100)
{
    BinaryGA GA(population_size, problem.num_vars(), problem);

    GA.algorithm(algorithm::NSGA3{});
    GA.crossover_method(crossover::binary::TwoPoint{ 0.8 });
    GA.mutation_method(mutation::binary::Flip{ 1.0 / problem.num_vars() });

    benchmarkMoga(GA, generations, "NSGA3", problem);
}

inline void benchmark_nsga3_zdt(size_t generations = 250, size_t population_size = 100)
{
    benchmark_real_nsga3(Kursawe{}, generations, population_size);
    benchmark_real_nsga3(ZDT1{}, generations, population_size);
    benchmark_real_nsga3(ZDT2{}, generations, population_size);
    benchmark_real_nsga3(ZDT3{}, generations, population_size);
    benchmark_real_nsga3(ZDT4{}, generations, population_size);
    benchmark_binary_nsga3(ZDT5{}, generations, population_size);
    benchmark_real_nsga3(ZDT6{}, generations, population_size);
}

inline void benchmark_nsga3_dtlz(size_t generations = 1000, size_t population_size = 100, size_t dim = 3)
{
    benchmark_real_nsga3(DTLZ1{ dim }, generations, population_size);
    benchmark_real_nsga3(DTLZ2{ dim }, generations, population_size);
    benchmark_real_nsga3(DTLZ3{ dim }, generations, population_size);
    benchmark_real_nsga3(DTLZ4{ dim }, generations, population_size);
    benchmark_real_nsga3(DTLZ5{ dim }, generations, population_size);
    benchmark_real_nsga3(DTLZ6{ dim }, generations, population_size);
    benchmark_real_nsga3(DTLZ7{ dim }, generations, population_size);
}

#endif // !GA_TEST_BENCHMARK_NSGA3_HPP