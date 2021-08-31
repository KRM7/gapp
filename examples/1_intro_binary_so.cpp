/* Example showing the basics through a single-objective genetic algorithm using binary encoding. */

//#include "../src/genetic_algorithm.h  /* This would include everything from the library. */
#include "../src/binary_ga.h"           /* For the BinaryGA class. */
#include "../test/fitness_functions.h"  /* For the fitness function used, and for converting the binary chromosomes to reals. */

#include <cstdio>
#include <cstddef>
#include <algorithm>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace genetic_algorithm;

int main()
{
    /* 
    * Define the fitness function that will be used in the GA. (Assuming fitness maximization.)
    * The fitness function can be anything that std::function can store (eg. functions, lambdas, etc.), and it must be thread-safe.
    * The fitness function must take a vector of geneType (char for the binaryGA), and return the fitness vector (a vector with 1 element
    * for single-objective problems).
    * In this example, it will be the Rastrigin function in 10 dimensions, already used in the benchmarks.
    */
    size_t num_vars = 10;
    Rastrigin rastriginFunction(num_vars);

    /* Create the GA using the appropriate chromosome length and the fitness function. */
    size_t chrom_len = rastriginFunction.num_vars * rastriginFunction.var_bits;
    BinaryGA GA(chrom_len, rastriginFunction);

    /* Set the algorithm used (based on the number of objectives). */
    GA.mode(BinaryGA::Mode::single_objective);  /* Single-objective is the default, so this isn't neccesary in this case. */

    /* Run the GA using the default parameters. */
    auto solutions = GA.run();                  /* Returns a vector of candidates (the optimal solutions). */

    /* Display results. */
    cout << "The results of the first run:\n";
    cout << " Best fitness: " << scientific << setprecision(2) << -solutions[0].fitness[0] << "\n\n";


    /* Set the parameters of the GA for better results. */

    /* Basic parameters. */
    GA.population_size(400);
    GA.crossover_rate(0.75);
    GA.mutation_rate(0.015);    /* The mutation rate is set to 1/chrom_len by default. */
    GA.max_gen(1500);           /* The GA will always stop when reaching the maximum number of generations. */


    /* 
    * The selection method can be changed for the single-objective algorithm.
    * The multi-objective algorithms always use their predefined selection methods.
    * Some selection methods have additional parameters, which can also be set (in this case the size of the tournaments),
    * but the default values are usually fine.
    */
    GA.selection_method(BinaryGA::SogaSelection::tournament);
    //GA.tournament_size(2);    /* 2 by default. */

    /*
    * Select the crossover method method used.
    * Note that in this case instead of using the n-point crossover with 2 crossover points, we could also set
    * the crossover method to two-point.
    */
    GA.crossover_method(BinaryGA::CrossoverMethod::n_point);
    GA.num_crossover_points(2);
    //GA.crossover_method(BinaryGA::CrossoverMethod::two_point);

    /*
    * It is possible to set an early-stop condition (optional).
    * Some of the stop conditions only work for the single-objective algorithm.
    */
    GA.stop_condition(BinaryGA::StopCondition::fitness_mean_stall); /* The stall stop conditions only work for SOGAs. */
    GA.stall_gen_count(50);
    GA.stall_threshold(0.005);

    /* A different stop condition example. */
    //GA.stop_condition(BinaryGA::StopCondition::fitness_evals);
    //GA.max_fitness_evals(10000);

    
    /* Some other settings. */

    /*
    * If archive_optimal_solutions is set to true, the GA will save every pareto optimal solutions, not just the ones in the last generation.
    * False is the default value.
    * This doesn't do much for the single-objective algorithms.
    */
    GA.archive_optimal_solutions = false;
    
    /*
    * If the fitness function doesn't change over time (always returns the same fitness vector for a given chromosome),
    * some unnecesary fitness function evaluations can be saved.
    * This is turned on by default (changing_fitness_func = false).
    */
    GA.changing_fitness_func = false;

    /*
    * It's possible to use a preset initial population instead of randomly generating one (which is the default behaviour).
    * The number of Candidates in the preset vector does not change the population size used in the algorithm from the population size that was set:
    * If there are less candidates in the preset initial population than the population size set, then the remaining
    * candidates will be generated randomly, while candidates will be discarded if the preset population is larger
    * than the population size set.
    */
    //GA.presetInitialPopulation(/* insert initial candidate vector */);


    /* Run the GA again with these new settings. The return value can also be ignored. */
    GA.run();


    /* Print the results of the second run. */
    auto sols = GA.solutions();
    size_t num_evals = GA.num_fitness_evals();
    //auto hist = GA.soga_history();    /* Contains stats of fitness values from each generation (min, max, mean, SD) of the SOGA. */

    cout << "The results of the second run of the algorithm:\n";
    cout << " The number of fitness function evals: " << num_evals << "\n";
    cout << " The best solutions found:\n";
    for (const auto& sol : sols)
    {
        /* Convert the binary chromosomes of the solutions. */
        vector<double> chrom = convertToReals(sol.chromosome, rastriginFunction.var_bits, rastriginFunction.intval());
        for_each(chrom.begin(), chrom.end(), [&rastriginFunction](double& var) { var += rastriginFunction.lbound(); });

        cout << "  f(x) = " << -sol.fitness[0] << " at x = (";
        for (const auto& gene : chrom)
        {
            cout << gene << ", ";
        }
        cout << ")\n";
    }

    getchar();
    return 0;
}