/* Simple example showing the usage of the RCGA (real-coded genetic algorithm). */

#include "../src/real_ga.h"
#include "../test/fitness_functions.h" /* For the fitness function that will be used. */

#include <cstdio>
#include <cstddef>
#include <iostream>
#include <iomanip>
#include <utility>

using namespace std;
using namespace genetic_algorithm;

int main()
{
    /*
    * Define the fitness function that will be used in the GA. (Assuming fitness maximization.)
    * The fitness function must take a vector of geneType (double for the RCGA),
    * and return the fitness vector (a vector with 1 element for single-objective problems).
    * In this example, it will be the Rastrigin function in 10 dimensions, already used in the benchmarks.
    */
    size_t num_vars = 10;
    Rastrigin rastriginFunction(num_vars);

    /* In the case of the real-coded GA, the lower and upper limit of each gene must be set. */
    pair<double, double> limit = { rastriginFunction.lbound(), rastriginFunction.ubound() };
    vector<pair<double, double>> limits(rastriginFunction.num_vars, limit);

    /* Create the GA using the appropriate chromosome length and the fitness function. */
    RCGA GA(num_vars, rastriginFunction, limits);

    /* Set the other parameters of the GA. */
    GA.population_size(100);
    GA.crossover_rate(0.6);
    GA.mutation_rate(0.05); /* The mutation rate is set to 1/chrom_len by default. */
    GA.max_gen(1000);

    GA.selection_method(RCGA::SogaSelection::boltzmann);
    //GA.boltzmann_temps(0.1, 4.0);

    GA.crossover_method(RCGA::CrossoverMethod::simulated_binary);
    GA.sim_binary_crossover_param(4.0);

    GA.mutation_method(RCGA::MutationMethod::gauss);
    GA.gauss_mutation_param(2.0);

    /* Set an early stop-condition. */
    GA.stop_condition(RCGA::StopCondition::fitness_value);
    GA.fitness_threshold({ -0.005 });

    /* Run the GA. */
    auto sols = GA.run();

    /* Print the results. */
    auto num_evals = GA.num_fitness_evals();

    cout << "The results of the algorithm:\n";
    cout << " The number of fitness function evals: " << num_evals << "\n";
    cout << " The best solutions found:\n" << scientific << setprecision(2);
    for (const auto& sol : sols)
    {
        cout << "  f(x) = " << -sol.fitness[0] << " at x = (";
        for (const auto& gene : sol.chromosome)
        {
            cout << gene << ", ";
        }
        cout << ")\n";
    }

    getchar();
    return 0;
}