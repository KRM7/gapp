/* Simple example showing the usage of the permutationGA. */

#include "../src/genetic_algorithm.h"
#include "../test/fitness_functions.h" /* For the fitness function that will be used. */

#include <cstdio>
#include <cstddef>
#include <iostream>
#include <iomanip>
#include <chrono>

using namespace std;
using namespace chrono;
using namespace genetic_algorithm;

int main()
{
    /*
    * Define the fitness function that will be used in the GA. (Assuming fitness maximization.)
    * The fitness function must take a vector of geneType (size_t for the permutationGA),
    * and return the fitness vector (a vector with 1 element for single-objective problems).
    * In this example, it will be a TSP with 439 nodes, already used in the benchmarks.
    */
    TSP tsp439("../test/tsp_data/tsp439.txt");

    /* Create the GA using the appropriate chromosome length (number of nodes) and the fitness function. */
    size_t num_nodes = tsp439.num_vars();
    PermutationGA GA(num_nodes, tsp439);

    /* Set the other parameters of the GA. */
    GA.population_size(500);
    GA.crossover_rate(0.9);
    GA.mutation_rate(0.3);  /* This is a per-chromosome mutation rate for the permutationGAs. */
    GA.max_gen(1000);

    GA.selection_method(PermutationGA::SogaSelection::boltzmann);
    GA.crossover_method(PermutationGA::CrossoverMethod::order);
    GA.mutation_method(PermutationGA::MutationMethod::inversion);

    GA.repairFunction = nullptr;  /* Could be something like 2-opt for example (nullptr by default). */

    /* Run the GA. */
    auto tbegin = high_resolution_clock::now();
    auto sols = GA.run();   /* This will take some time even with optimizations. */
    auto tend = high_resolution_clock::now();
    double time_spent = double(duration_cast<microseconds>(tend - tbegin).count()) / 1E+6;

    /* Print the results. */
    cout << "The number of optimal sols found for the TSP439: " << sols.size() << "\n";
    cout << "The length of the shortest route found: " << -sols[0].fitness[0] << " (theoretical best is " << -tsp439.optimal_value() << ").\n";
    cout << "The number of fitness function evals performed: " << GA.num_fitness_evals() << "\n";
    cout << "Time taken: " << time_spent << " s\n\n";

    getchar();
    return 0;
}