/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_TEST_BENCHMARK_INTEGER_HPP
#define GA_TEST_BENCHMARK_INTEGER_HPP

#include "encoding/integer.hpp"
#include "algorithm/algorithm.hpp"
#include "crossover/integer.hpp"
#include "mutation/integer.hpp"
#include "stop_condition/stop_condition.hpp"
#include "problems/integer.hpp"
#include "benchmark_utils.hpp"

using namespace genetic_algorithm;
using namespace genetic_algorithm::problems;

void integer_hello()
{
    StringFinder fitness_func("HELLO WORLD!");

    IntegerGA GA(100, fitness_func.num_vars(), fitness_func, 96);

    GA.algorithm(algorithm::SingleObjective{ selection::Tournament{} });
    GA.crossover_method(crossover::integer::TwoPoint{ 0.85 });
    GA.mutation_method(mutation::integer::Uniform{ 0.01 });
    //+32
    benchmarkInt(GA, 500, fitness_func);
}

void integer_sentence()
{
    StringFinder fitness_func("Lorem ipsum dolor sit amet, consectetur adipiscing elit. Pellentesque gravida ut ipsum at tincidunt.");

    IntegerGA GA(250, fitness_func.num_vars(), fitness_func, 96);

    GA.algorithm(algorithm::SingleObjective{ selection::Boltzmann{} });
    GA.crossover_method(crossover::integer::Uniform{ 0.8 });
    GA.mutation_method(mutation::integer::Uniform{ 5.0 / 250 });
    // +32
    benchmarkInt(GA, 1000, fitness_func);
}

#endif // !GA_TEST_BENCHMARK_INTEGER_HPP