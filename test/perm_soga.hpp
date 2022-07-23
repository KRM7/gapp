/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef PERM_SOGA_BENCHMARK_HPP
#define PERM_SOGA_BENCHMARK_HPP

#include "../src/encoding/permutation.hpp"
#include "../src/algorithm/algorithm.hpp"
#include "../src/crossover/permutation.hpp"
#include "../src/mutation/permutation.hpp"
#include "fitness_functions.hpp"
#include "benchmark_utils.hpp"

using namespace std;
using namespace genetic_algorithm;

void perm_tsp52()
{
    TSP tsp52("test/tsp_data/tsp52.txt");

    PermutationGA GA(tsp52.num_vars(), tsp52);

    GA.population_size(500);
    GA.algorithm(selection::Sigma{});
    GA.crossover_method(crossover::perm::Edge{ 0.9 });
    GA.mutation_method(mutation::perm::Inversion{ 0.05 });

    benchmarkSoga(GA, 1250, tsp52, "TSP52");
}

void perm_tsp76()
{
    TSP tsp76("test/tsp_data/tsp76.txt");

    PermutationGA GA(tsp76.num_vars(), tsp76);

    GA.population_size(400);
    GA.algorithm(selection::Tournament{});
    GA.crossover_method(crossover::perm::Order1{ 0.9 });
    GA.mutation_method(mutation::perm::Shift{ 0.1 });

    benchmarkSoga(GA, 1000, tsp76, "TSP76");
}

void perm_tsp124()
{
    TSP tsp124("test/tsp_data/tsp124.txt");

    PermutationGA GA(tsp124.num_vars(), tsp124);

    GA.population_size(500);
    GA.algorithm(selection::Boltzmann{});
    GA.crossover_method(crossover::perm::Position{ 0.9 });
    GA.mutation_method(mutation::perm::Swap3{ 0.4 });

    benchmarkSoga(GA, 1500, tsp124, "TSP124");
}

void perm_tsp152()
{
    TSP tsp152("test/tsp_data/tsp152.txt");

    PermutationGA GA(tsp152.num_vars(), tsp152);

    GA.population_size(500);
    GA.algorithm(selection::Tournament{});
    GA.crossover_method(crossover::perm::PMX{ 0.9 });
    GA.mutation_method(mutation::perm::Shift{ 0.6 });

    benchmarkSoga(GA, 1250, tsp152, "TSP152");
}

void perm_tsp226()
{
    TSP tsp226("test/tsp_data/tsp226.txt");

    PermutationGA GA(tsp226.num_vars(), tsp226);

    GA.population_size(500);
    GA.algorithm(selection::Roulette{});
    GA.crossover_method(crossover::perm::Cycle{ 0.9 });
    GA.mutation_method(mutation::perm::Shuffle{ 0.2 });

    benchmarkSoga(GA, 1250, tsp226, "TSP226");
}

void perm_tsp299()
{
    TSP tsp299("test/tsp_data/tsp299.txt");

    PermutationGA GA(tsp299.num_vars(), tsp299);

    GA.population_size(500);
    GA.algorithm(selection::Boltzmann{});
    GA.crossover_method(crossover::perm::Order2{ 0.9 });
    GA.mutation_method(mutation::perm::Inversion{ 0.3 });

    benchmarkSoga(GA, 1000, tsp299, "TSP299");
}

void perm_tsp439()
{
    TSP tsp439("test/tsp_data/tsp439.txt");

    PermutationGA GA(tsp439.num_vars(), tsp439);

    GA.population_size(500);
    GA.algorithm(selection::Boltzmann{});
    GA.crossover_method(crossover::perm::Order2{ 0.9 });
    GA.mutation_method(mutation::perm::Inversion{ 0.3 });

    benchmarkSoga(GA, 1000, tsp439, "TSP439");
}

#endif // !PERM_SOGA_BENCHMARK_HPP