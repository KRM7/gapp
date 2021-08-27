/* Benchmark/test functions for the permutationGA. */

#ifndef PERMUTATIONAL_TESTS_H
#define PERMUTATIONAL_TESTS_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <chrono>

#include "../src/permutation_ga.h"
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
    GA.crossover_rate(0.9);
    GA.mutation_rate(0.05);
    GA.selection_method(PermutationGA::SogaSelection::sigma);
    GA.crossover_method(PermutationGA::CrossoverMethod::edge);
    GA.mutation_method(PermutationGA::MutationMethod::inversion);
    GA.max_gen(1250);

    /* Run the GA with a timer. */
    auto tbegin = chrono::high_resolution_clock::now();
    auto sols = GA.run();
    auto tend = chrono::high_resolution_clock::now();

    auto duration = chrono::duration_cast<chrono::microseconds>(tend - tbegin).count();
    double time_spent = double(duration) / 1E+6;

    /* Print the results. */
    cout << "\n\nThe number of optimal sols found for the TSP52: " << sols.size() << "\n";
    cout << "The length of the shortest route found: " << -sols[0].fitness[0] << " (best is " << -tsp52.optimal_value() << ").\n";
    cout << "Number of fitness evals: " << GA.num_fitness_evals() << "\n";
    cout << "Time taken: " << time_spent << " s\n\n";

    /* Print the results. */
    //displayStats(GA.soga_history());
}

void perm124Test()
{
    /* Init GA. */
    TSP tsp124("test/tsp_data/tsp124.txt");

    PermutationGA GA(tsp124.num_vars(), tsp124);

    /* Set some optional parameters. */
    GA.population_size(500);
    GA.crossover_rate(0.95);
    GA.mutation_rate(0.5);
    GA.selection_method(PermutationGA::SogaSelection::boltzmann);
    GA.crossover_method(PermutationGA::CrossoverMethod::order);
    GA.mutation_method(PermutationGA::MutationMethod::swap);
    GA.max_gen(1500);

    /* Run the GA with a timer. */
    auto tbegin = chrono::high_resolution_clock::now();
    auto sols = GA.run();
    auto tend = chrono::high_resolution_clock::now();

    auto duration = chrono::duration_cast<chrono::microseconds>(tend - tbegin).count();
    double time_spent = double(duration) / 1E+6;

    /* Print the results. */
    cout << "\n\nThe number of optimal sols found for the TSP124: " << sols.size() << "\n";
    cout << "The length of the shortest route found: " << -sols[0].fitness[0] << " (best is " << -tsp124.optimal_value() << ").\n";
    cout << "Number of fitness evals: " << GA.num_fitness_evals() << "\n";
    cout << "Time taken: " << time_spent << " s\n\n";

    /* Print the results. */
    //displayStats(GA.soga_history());
}

void perm226Test()
{
    /* Init GA. */
    TSP tsp226("test/tsp_data/tsp226.txt");

    PermutationGA GA(tsp226.num_vars(), tsp226);

    /* Set some optional parameters. */
    GA.population_size(500);
    GA.crossover_rate(0.9);
    GA.mutation_rate(0.6);
    GA.selection_method(PermutationGA::SogaSelection::roulette);
    GA.crossover_method(PermutationGA::CrossoverMethod::cycle);
    GA.mutation_method(PermutationGA::MutationMethod::scramble);
    GA.max_gen(1000);

    /* Run the GA with a timer. */
    auto tbegin = chrono::high_resolution_clock::now();
    auto sols = GA.run();
    auto tend = chrono::high_resolution_clock::now();

    auto duration = chrono::duration_cast<chrono::microseconds>(tend - tbegin).count();
    double time_spent = double(duration) / 1E+6;

    /* Print the results. */
    cout << "\n\nThe number of optimal sols found for the TSP226: " << sols.size() << "\n";
    cout << "The length of the shortest route found: " << -sols[0].fitness[0] << " (best is " << -tsp226.optimal_value() << ").\n";
    cout << "Number of fitness evals: " << GA.num_fitness_evals() << "\n";
    cout << "Time taken: " << time_spent << " s\n\n";

    /* Print the results. */
    //displayStats(GA.soga_history());
}

void perm439Test()
{
    /* Init GA. */
    TSP tsp439("test/tsp_data/tsp439.txt");

    PermutationGA GA(tsp439.num_vars(), tsp439);

    /* Set some optional parameters. */
    GA.population_size(500);
    GA.crossover_rate(0.9);
    GA.mutation_rate(0.3);
    GA.selection_method(PermutationGA::SogaSelection::boltzmann);
    GA.crossover_method(PermutationGA::CrossoverMethod::order);
    GA.mutation_method(PermutationGA::MutationMethod::inversion);
    GA.max_gen(1000);

    /* Run the GA with a timer. */
    auto tbegin = chrono::high_resolution_clock::now();
    auto sols = GA.run();
    auto tend = chrono::high_resolution_clock::now();

    auto duration = chrono::duration_cast<chrono::microseconds>(tend - tbegin).count();
    double time_spent = double(duration) / 1E+6;

    /* Print the results. */
    cout << "\n\nThe number of optimal sols found for the TSP439: " << sols.size() << "\n";
    cout << "The length of the shortest route found: " << -sols[0].fitness[0] << " (best is " << -tsp439.optimal_value() << ").\n";
    cout << "Number of fitness evals: " << GA.num_fitness_evals() << "\n";
    cout << "Time taken: " << time_spent << " s\n\n";

    /* Print the results. */
    //displayStats(GA.soga_history());
}

#endif // !PERMUTATIONAL_TESTS_H