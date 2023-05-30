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

using namespace gapp;
using namespace gapp::problems;

inline void perm_tsp52()
{
    PermutationGA GA(500);

    GA.algorithm(algorithm::SingleObjective{ selection::Sigma{} });
    GA.crossover_method(crossover::perm::Edge{ 0.9 });
    GA.mutation_method(mutation::perm::Inversion{ 0.5 });

    benchmarkTSP(GA, 1250, TSP52{});
}

inline void perm_tsp76()
{
    PermutationGA GA(400);

    GA.algorithm(algorithm::SingleObjective{ selection::Tournament{} });
    GA.crossover_method(crossover::perm::Order1{ 0.9 });
    GA.mutation_method(mutation::perm::Shift{ 0.5 });

    benchmarkTSP(GA, 1000, TSP76{});
}

inline void perm_tsp124()
{
    PermutationGA GA(500);

    GA.algorithm(algorithm::SingleObjective{ selection::Boltzmann{} });
    GA.crossover_method(crossover::perm::Position{ 0.9 });
    GA.mutation_method(mutation::perm::Swap3{ 0.5 });

    benchmarkTSP(GA, 1500, TSP124{});
}

inline void perm_tsp152()
{
    PermutationGA GA(500);

    GA.algorithm(algorithm::SingleObjective{ selection::Tournament{} });
    GA.crossover_method(crossover::perm::PMX{ 0.9 });
    GA.mutation_method(mutation::perm::Shift{ 0.5 });

    benchmarkTSP(GA, 1250, TSP152{});
}

inline void perm_tsp226()
{
    PermutationGA GA(500);

    GA.algorithm(algorithm::SingleObjective{ selection::Roulette{} });
    GA.crossover_method(crossover::perm::Cycle{ 0.9 });
    GA.mutation_method(mutation::perm::Shuffle{ 0.7 });

    benchmarkTSP(GA, 1250, TSP226{});
}

inline void perm_tsp299()
{
    PermutationGA GA(500);

    GA.algorithm(algorithm::SingleObjective{ selection::Boltzmann{} });
    GA.crossover_method(crossover::perm::Order2{ 0.9 });
    GA.mutation_method(mutation::perm::Inversion{ 0.3 });

    benchmarkTSP(GA, 1000, TSP299{});
}

inline void perm_tsp439()
{
    PermutationGA GA(500);

    GA.algorithm(algorithm::SingleObjective{ selection::Boltzmann{} });
    GA.crossover_method(crossover::perm::Order2{ 0.9 });
    GA.mutation_method(mutation::perm::Inversion{ 0.3 });

    benchmarkTSP(GA, 1000, TSP439{});
}

#endif // !GA_TEST_BENCHMARK_PERM_HPP