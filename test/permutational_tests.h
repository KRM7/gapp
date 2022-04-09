/* Benchmark/test functions for the permutationGA. */

#ifndef PERMUTATIONAL_TESTS_H
#define PERMUTATIONAL_TESTS_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <chrono>

#include "../src/algorithms/permutation_ga.hpp"
#include "../src/selection/selection.hpp"
#include "../src/crossover/permutation.hpp"
#include "../src/mutation/permutation.hpp"
#include "fitness_functions.h"
#include "utils.h"

using namespace std;
using namespace genetic_algorithm;

void perm52Test()
{
    /* Init GA. */
    TSP tsp52("test/tsp_data/tsp52.txt");

    PermutationGA GA(tsp52.num_vars(), tsp52);

    /* Set some optional parameters. */
    GA.population_size(500);
    GA.selection_method(selection::single_objective::Sigma{});
    GA.crossover_method(crossover::perm::Edge{ 0.9 });
    GA.mutation_method(mutation::perm::Inversion{ 0.05 });

    /* Run the GA with a timer. */
    auto tbegin = chrono::high_resolution_clock::now();
    auto sols = GA.run(1250);
    auto tend = chrono::high_resolution_clock::now();

    auto duration = chrono::duration_cast<chrono::microseconds>(tend - tbegin).count();
    double time_spent = double(duration) / 1E+6;

    /* Print the results. */
    cout << "\n\nThe number of optimal sols found for the TSP52: " << sols.size() << "\n";
    cout << "The length of the shortest route found: " << -sols[0].fitness[0] << " (best is " << -tsp52.optimal_value() << ").\n";
    cout << "Number of fitness evals: " << GA.num_fitness_evals() << "\n";
    cout << "Time taken: " << time_spent << " s\n\n";
}

void perm76Test()
{
    /* Init GA. */
    TSP tsp76("test/tsp_data/tsp76.txt");

    PermutationGA GA(tsp76.num_vars(), tsp76);

    /* Set some optional parameters. */
    GA.population_size(400);
    GA.selection_method(selection::single_objective::Tournament{});
    GA.crossover_method(crossover::perm::Order1{ 0.9 });
    GA.mutation_method(mutation::perm::Shift{ 0.1 });

    /* Run the GA with a timer. */
    auto tbegin = chrono::high_resolution_clock::now();
    auto sols = GA.run(1000);
    auto tend = chrono::high_resolution_clock::now();

    auto duration = chrono::duration_cast<chrono::microseconds>(tend - tbegin).count();
    double time_spent = double(duration) / 1E+6;

    /* Print the results. */
    cout << "\n\nThe number of optimal sols found for the TSP76: " << sols.size() << "\n";
    cout << "The length of the shortest route found: " << -sols[0].fitness[0] << " (best is " << -tsp76.optimal_value() << ").\n";
    cout << "Number of fitness evals: " << GA.num_fitness_evals() << "\n";
    cout << "Time taken: " << time_spent << " s\n\n";
}

void perm124Test()
{
    /* Init GA. */
    TSP tsp124("test/tsp_data/tsp124.txt");

    PermutationGA GA(tsp124.num_vars(), tsp124);

    /* Set some optional parameters. */
    GA.population_size(500);
    GA.selection_method(selection::single_objective::Boltzmann{});
    GA.crossover_method(crossover::perm::Position{ 0.9 });
    GA.mutation_method(mutation::perm::Swap3{ 0.4 });

    /* Run the GA with a timer. */
    auto tbegin = chrono::high_resolution_clock::now();
    auto sols = GA.run(1500);
    auto tend = chrono::high_resolution_clock::now();

    auto duration = chrono::duration_cast<chrono::microseconds>(tend - tbegin).count();
    double time_spent = double(duration) / 1E+6;

    /* Print the results. */
    cout << "\n\nThe number of optimal sols found for the TSP124: " << sols.size() << "\n";
    cout << "The length of the shortest route found: " << -sols[0].fitness[0] << " (best is " << -tsp124.optimal_value() << ").\n";
    cout << "Number of fitness evals: " << GA.num_fitness_evals() << "\n";
    cout << "Time taken: " << time_spent << " s\n\n";
}

void perm152Test()
{
    /* Init GA. */
    TSP tsp152("test/tsp_data/tsp152.txt");

    PermutationGA GA(tsp152.num_vars(), tsp152);

    /* Set some optional parameters. */
    GA.population_size(500);
    GA.selection_method(selection::single_objective::Tournament{});
    GA.crossover_method(crossover::perm::PMX{ 0.9 });
    GA.mutation_method(mutation::perm::Shift{ 0.6 });

    /* Run the GA with a timer. */
    auto tbegin = chrono::high_resolution_clock::now();
    auto sols = GA.run(1250);
    auto tend = chrono::high_resolution_clock::now();

    auto duration = chrono::duration_cast<chrono::microseconds>(tend - tbegin).count();
    double time_spent = double(duration) / 1E+6;

    /* Print the results. */
    cout << "\n\nThe number of optimal sols found for the TSP152: " << sols.size() << "\n";
    cout << "The length of the shortest route found: " << -sols[0].fitness[0] << " (best is " << -tsp152.optimal_value() << ").\n";
    cout << "Number of fitness evals: " << GA.num_fitness_evals() << "\n";
    cout << "Time taken: " << time_spent << " s\n\n";
}

void perm226Test()
{
    /* Init GA. */
    TSP tsp226("test/tsp_data/tsp226.txt");

    PermutationGA GA(tsp226.num_vars(), tsp226);

    /* Set some optional parameters. */
    GA.population_size(500);
    GA.selection_method(selection::single_objective::Roulette{});
    GA.crossover_method(crossover::perm::Cycle{ 0.9 });
    GA.mutation_method(mutation::perm::Shuffle{ 0.2 });

    /* Run the GA with a timer. */
    auto tbegin = chrono::high_resolution_clock::now();
    auto sols = GA.run(1250);
    auto tend = chrono::high_resolution_clock::now();

    auto duration = chrono::duration_cast<chrono::microseconds>(tend - tbegin).count();
    double time_spent = double(duration) / 1E+6;

    /* Print the results. */
    cout << "\n\nThe number of optimal sols found for the TSP226: " << sols.size() << "\n";
    cout << "The length of the shortest route found: " << -sols[0].fitness[0] << " (best is " << -tsp226.optimal_value() << ").\n";
    cout << "Number of fitness evals: " << GA.num_fitness_evals() << "\n";
    cout << "Time taken: " << time_spent << " s\n\n";
}

void perm299Test()
{
    /* Init GA. */
    TSP tsp299("test/tsp_data/tsp299.txt");

    PermutationGA GA(tsp299.num_vars(), tsp299);

    /* Set some optional parameters. */
    GA.population_size(500);
    GA.selection_method(selection::single_objective::Boltzmann{});
    GA.crossover_method(crossover::perm::Order2{ 0.9 });
    GA.mutation_method(mutation::perm::Inversion{ 0.3 });

    /* Run the GA with a timer. */
    auto tbegin = chrono::high_resolution_clock::now();
    auto sols = GA.run(1000);
    auto tend = chrono::high_resolution_clock::now();

    auto duration = chrono::duration_cast<chrono::microseconds>(tend - tbegin).count();
    double time_spent = double(duration) / 1E+6;

    /* Print the results. */
    cout << "\n\nThe number of optimal sols found for the TSP299: " << sols.size() << "\n";
    cout << "The length of the shortest route found: " << -sols[0].fitness[0] << " (best is " << -tsp299.optimal_value() << ").\n";
    cout << "Number of fitness evals: " << GA.num_fitness_evals() << "\n";
    cout << "Time taken: " << time_spent << " s\n\n";
}

void perm439Test()
{
    /* Init GA. */
    TSP tsp439("test/tsp_data/tsp439.txt");

    PermutationGA GA(tsp439.num_vars(), tsp439);

    /* Set some optional parameters. */
    GA.population_size(500);
    GA.selection_method(selection::single_objective::Boltzmann{});
    GA.crossover_method(crossover::perm::Order2{ 0.9 });
    GA.mutation_method(mutation::perm::Inversion{ 0.3 });

    /* Run the GA with a timer. */
    auto tbegin = chrono::high_resolution_clock::now();
    auto sols = GA.run(1000);
    auto tend = chrono::high_resolution_clock::now();

    auto duration = chrono::duration_cast<chrono::microseconds>(tend - tbegin).count();
    double time_spent = double(duration) / 1E+6;

    /* Print the results. */
    cout << "\n\nThe number of optimal sols found for the TSP439: " << sols.size() << "\n";
    cout << "The length of the shortest route found: " << -sols[0].fitness[0] << " (best is " << -tsp439.optimal_value() << ").\n";
    cout << "Number of fitness evals: " << GA.num_fitness_evals() << "\n";
    cout << "Time taken: " << time_spent << " s\n\n";
}

#endif // !PERMUTATIONAL_TESTS_H