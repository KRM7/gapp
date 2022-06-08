#ifndef INTEGER_BENCHMARK_HPP
#define INTEGER_BENCHMARK_HPP

#include <vector>
#include <chrono>
#include <iostream>

#include "../src/encoding/integer.hpp"
#include "../src/selection/selection.hpp"
#include "../src/crossover/integer.hpp"
#include "../src/mutation/integer.hpp"
#include "fitness_functions.hpp"
#include "benchmark_utils.hpp"

using namespace genetic_algorithm;

void integer_hello()
{
    MatchString fitness_func("HELLO WORLD!");

    IntegerGA GA(100, fitness_func.num_vars(), fitness_func, 96);

    GA.selection_method(selection::single_objective::Tournament{});
    GA.crossover_method(crossover::integer::TwoPoint{ 0.85 });
    GA.mutation_method(mutation::integer::Uniform{ 0.01 });
    //+32
    benchmarkSoga(GA, 500, fitness_func, "Hello");
}

void integer_sentence()
{
    MatchString fitness_func("Lorem ipsum dolor sit amet, consectetur adipiscing elit. Pellentesque gravida ut ipsum at tincidunt.");

    IntegerGA GA(250, fitness_func.num_vars(), fitness_func, 96);

    GA.selection_method(selection::single_objective::Boltzmann{});
    GA.crossover_method(crossover::integer::Uniform{ 0.8 });
    GA.mutation_method(mutation::integer::Uniform{ 5.0 / 250 });
    // +32
    benchmarkSoga(GA, 1000, fitness_func, "Lorem");
}

#endif // !INTEGER_BENCHMARK_HPP