/* Simple example showing the usage of the NSGA-II. */

#include "../src/genetic_algorithm.h"
#include "../test/fitness_functions.h"  /* For the fitness function that will be used. */
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
    * Define the fitness function that will be used in the GA. (Assuming fitness maximization.)
    * The fitness function must take a vector of geneType (based on the encoding, in this case double),
    * and return the fitness vector.
    * In this example, it will be the Kursawe function with 3 variables and 2 objectives, already used in the benchmarks.
    */
    size_t num_vars = 3;
    KUR kursaweFunction(num_vars);

    /* Set the lower and upper limit of each gene for the RCGA. */
    pair<double, double> limit = { kursaweFunction.lbound(), kursaweFunction.ubound() };
    vector<pair<double, double>> limits(kursaweFunction.num_vars, limit);

    /* 
    * Create the GA using the appropriate chromosome length and the fitness function.
    * This example uses real-encoding, but the NSGA-II algorithm is not dependant on any particular encoding.
    */
    size_t chrom_len = kursaweFunction.num_vars;
    RCGA GA(chrom_len, kursaweFunction, limits);

    GA.mode(RCGA::Mode::multi_objective_sorting); /* Use the NSGA-II algorithm. (The encoding type doesn't matter.) */

    /* Set the other parameters of the GA. */
    GA.population_size(100);
    GA.crossover_rate(0.8);
    //use default mutation rate (1/chrom_len)
    GA.max_gen(250);

    GA.crossover_method(RCGA::CrossoverMethod::simulated_binary);
    GA.mutation_method(RCGA::MutationMethod::gauss);
    GA.gauss_mutation_param(2.0);

    /* Set an early-stop condition. */
    GA.stop_condition(RCGA::StopCondition::fitness_evals);
    GA.max_fitness_evals(20000);

    GA.archive_optimal_solutions = true;    /* Keep every pareto-optimal solution. */

    /* Run the GA. */
    auto sols = GA.run();

    /* Print the results. */
    cout << "The number of pareto-optimal solutions found for the KUR problem with the NSGA-II: " << sols.size() << "\n";
    cout << "The number of fitness function evals: " << GA.num_fitness_evals() << "\n";

    /* Write the fitness values of the solutions to a file for plotting. */
    writeResultsToFile(GA.population(), "mo_results/nsga2_kur_last.txt");   /* Last generation. */
    writeResultsToFile(sols, "mo_results/nsga2_kur_sols.txt");              /* Every pareto-optimal solution. */
    

    getchar();
    return 0;
}