/* Example showing the usage of the NSGA-III algorithm for multi-objective optimization. */

#include "../src/real_ga.h"             /* For the real-encoded genetic algorithm class. */
#include "../test/fitness_functions.h"  /* For the fitness function that will be used in the example. */
#include "../test/utils.h"              /* For writing the results to a file. */

#include <cstdio>
#include <cstddef>
#include <iostream>
#include <iomanip>
#include <utility>
#include <vector>

using namespace std;
using namespace genetic_algorithm;

int main()
{
    /*
    * Define the fitness function that will be used in the GA (assuming fitness maximization).
    * The fitness function must be thread-safe.
    * The fitness function must take a vector of geneType (based on the encoding, double in this case),
    * and return the fitness vector.
    * In this example, it will be the DTLZ1 function with 7 variables and 3 objectives, already used in the benchmarks.
    */
    size_t num_vars = 7;
    size_t num_obj = 3;
    DTLZ1 dtlz1Function(num_vars, num_obj);

    /* Set the lower and upper limit of each gene for the RCGA. Same for every gene (variable) in this example. */
    pair<double, double> limit = { dtlz1Function.lbound(), dtlz1Function.ubound() };
    vector<pair<double, double>> limits(dtlz1Function.num_vars, limit);

    /*
    * Create the GA using the appropriate chromosome length and the fitness function.
    * This example uses real-encoding, but the NSGA-II algorithm is not dependant on any particular encoding.
    */
    size_t chrom_len = dtlz1Function.num_vars;
    RCGA GA(chrom_len, dtlz1Function, limits);

    GA.mode(RCGA::Mode::multi_objective_decomp);    /* Use the NSGA-III algorithm. (The encoding type doesn't matter.) */

    /* Set the other parameters of the GA. */
    GA.population_size(100);
    GA.crossover_rate(0.9);
    /* Use the default mutation rate (1/chrom_len) */
    GA.max_gen(750);

    GA.crossover_method(RCGA::CrossoverMethod::simulated_binary);
    GA.sim_binary_crossover_param(15.0);

    GA.mutation_method(RCGA::MutationMethod::random);

    /* Keep every pareto-optimal solution. */
    GA.archive_optimal_solutions = true;


    /* Run the algorithm. */
    GA.run();


    /* Print the results. */
    cout << "The number of pareto-optimal solutions found for the DTLZ1 problem with the NSGA-III: " << GA.solutions().size() << "\n";
    cout << "The number of fitness function evals: " << GA.num_fitness_evals() << "\n";

    /* Write the fitness values of the solutions to a file for plotting. */
    writeResultsToFile(GA.population(), "mo_results/nsga3_dtlz1_last.txt");   /* Last generation. */
    writeResultsToFile(GA.solutions(), "mo_results/nsga3_dtlz1_sols.txt");    /* Every pareto-optimal solution. */


    getchar();
    return 0;
}