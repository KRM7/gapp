/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_TEST_BENCHMARK_PERM_HPP
#define GA_TEST_BENCHMARK_PERM_HPP

#include "encoding/permutation.hpp"
#include "algorithm/algorithm.hpp"
#include "crossover/permutation.hpp"
#include "mutation/permutation.hpp"
#include "stop_condition/stop_condition.hpp"
#include "problems/travelling_salesman.hpp"
#include "benchmark_utils.hpp"

using namespace genetic_algorithm;
using namespace genetic_algorithm::problems;

void perm_tsp52()
{
    TSP52 tsp52;

    PermutationGA GA(tsp52, 500);

    GA.algorithm(algorithm::SingleObjective{ selection::Sigma{} });
    GA.crossover_method(crossover::perm::Edge{ 0.9 });
    GA.mutation_method(mutation::perm::Inversion{ 0.05 });

    benchmarkTSP(GA, 1250, tsp52);
}

void perm_tsp76()
{
    TSP76 tsp76;

    PermutationGA GA(tsp76, 400);

    GA.algorithm(algorithm::SingleObjective{ selection::Tournament{} });
    GA.crossover_method(crossover::perm::Order1{ 0.9 });
    GA.mutation_method(mutation::perm::Shift{ 0.1 });

    benchmarkTSP(GA, 1000, tsp76);
}

void perm_tsp124()
{
    TSP124 tsp124;

    PermutationGA GA(tsp124, 500);

    GA.algorithm(algorithm::SingleObjective{ selection::Boltzmann{} });
    GA.crossover_method(crossover::perm::Position{ 0.9 });
    GA.mutation_method(mutation::perm::Swap3{ 0.4 });

    benchmarkTSP(GA, 1500, tsp124);
}

void perm_tsp152()
{
    TSP152 tsp152;

    PermutationGA GA(tsp152, 500);

    GA.algorithm(algorithm::SingleObjective{ selection::Tournament{} });
    GA.crossover_method(crossover::perm::PMX{ 0.9 });
    GA.mutation_method(mutation::perm::Shift{ 0.6 });

    benchmarkTSP(GA, 1250, tsp152);
}

void perm_tsp226()
{
    TSP226 tsp226;

    PermutationGA GA(tsp226, 500);

    GA.algorithm(algorithm::SingleObjective{ selection::Roulette{} });
    GA.crossover_method(crossover::perm::Cycle{ 0.9 });
    GA.mutation_method(mutation::perm::Shuffle{ 0.2 });

    benchmarkTSP(GA, 1250, tsp226);
}

void perm_tsp299()
{
    TSP299 tsp299;

    PermutationGA GA(tsp299, 500);

    GA.algorithm(algorithm::SingleObjective{ selection::Boltzmann{} });
    GA.crossover_method(crossover::perm::Order2{ 0.9 });
    GA.mutation_method(mutation::perm::Inversion{ 0.3 });

    benchmarkTSP(GA, 1000, tsp299);
}

void perm_tsp439()
{
    TSP439 tsp439;

    PermutationGA GA(tsp439, 500);

    GA.algorithm(algorithm::SingleObjective{ selection::Boltzmann{} });
    GA.crossover_method(crossover::perm::Order2{ 0.9 });
    GA.mutation_method(mutation::perm::Inversion{ 0.3 });

    benchmarkTSP(GA, 1000, tsp439);
}

#endif // !GA_TEST_BENCHMARK_PERM_HPP