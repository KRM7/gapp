/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_TEST_BENCHMARK_NSGA2_HPP
#define GAPP_TEST_BENCHMARK_NSGA2_HPP

#include "encoding/real.hpp"
#include "encoding/binary.hpp"
#include "algorithm/algorithm.hpp"
#include "crossover/real.hpp"
#include "crossover/binary.hpp"
#include "mutation/real.hpp"
#include "mutation/binary.hpp"
#include "stop_condition/stop_condition.hpp"
#include "problems/multi_objective.hpp"
#include "benchmark_utils.hpp"

using namespace gapp;
using namespace gapp::problems;

template<typename Problem>
void benchmark_real_nsga2(Problem problem, size_t generations, size_t population_size = 100)
{
    RCGA GA(population_size);

    GA.algorithm(algorithm::NSGA2{});
    GA.crossover_method(crossover::real::SimulatedBinary{ 0.9 });
    GA.mutation_method(mutation::real::Uniform{ 1.0 / problem.num_vars() });

    benchmarkMoga(GA, generations, "NSGA2", std::move(problem));
}

template<typename Problem>
void benchmark_binary_nsga2(Problem problem, size_t generations, size_t population_size = 100)
{
    BinaryGA GA(population_size);

    GA.algorithm(algorithm::NSGA2{});
    GA.crossover_method(crossover::binary::TwoPoint{ 0.8 });
    GA.mutation_method(mutation::binary::Flip{ 1.0 / problem.num_vars() });

    benchmarkMoga(GA, generations, "NSGA2", std::move(problem));
}

inline void benchmark_nsga2_zdt(size_t generations = 250, size_t population_size = 100)
{
    benchmark_real_nsga2(Kursawe{}, generations, population_size);
    benchmark_real_nsga2(ZDT1{}, generations, population_size);
    benchmark_real_nsga2(ZDT2{}, generations, population_size);
    benchmark_real_nsga2(ZDT3{}, generations, population_size);
    benchmark_real_nsga2(ZDT4{}, generations, population_size);
    benchmark_binary_nsga2(ZDT5{}, generations, population_size);
    benchmark_real_nsga2(ZDT6{}, generations, population_size);
}

inline void benchmark_nsga2_dtlz(size_t generations = 1000, size_t population_size = 100, size_t dim = 3)
{
    benchmark_real_nsga2(DTLZ1{ dim }, generations, population_size);
    benchmark_real_nsga2(DTLZ2{ dim }, generations, population_size);
    benchmark_real_nsga2(DTLZ3{ dim }, generations, population_size);
    benchmark_real_nsga2(DTLZ4{ dim }, generations, population_size);
    benchmark_real_nsga2(DTLZ5{ dim }, generations, population_size);
    benchmark_real_nsga2(DTLZ6{ dim }, generations, population_size);
    benchmark_real_nsga2(DTLZ7{ dim }, generations, population_size);
}

#endif // !GAPP_TEST_BENCHMARK_NSGA2_HPP
