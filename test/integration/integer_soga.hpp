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

using namespace gapp;
using namespace gapp::problems;

inline void integer_hello()
{
    IntegerGA GA(100);

    GA.algorithm(algorithm::SingleObjective{ selection::Tournament{} });
    GA.crossover_method(crossover::integer::TwoPoint{ 0.85 });
    GA.mutation_method(mutation::integer::Uniform{ 0.01 });
    
    benchmarkInt(GA, 500, StringFinder{ "Hello World!" });
}

inline void integer_sentence()
{
    IntegerGA GA(250);

    GA.algorithm(algorithm::SingleObjective{ selection::Boltzmann{} });
    GA.crossover_method(crossover::integer::Uniform{ 0.8 });
    GA.mutation_method(mutation::integer::Uniform{ 5.0 / 250 });
    
    benchmarkInt(GA, 1000, StringFinder{ "Lorem ipsum dolor sit amet, consectetur adipiscing elit. Pellentesque gravida ut ipsum at tincidunt." });
}

#endif // !GA_TEST_BENCHMARK_INTEGER_HPP