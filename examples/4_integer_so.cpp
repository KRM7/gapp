/* Simple example showing the usage of the integerGA. */

#include "../src/integer_ga.h"      /* For the integer encoded genetic algorithm class. */

#include <cstdio>
#include <cstddef>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

using namespace std;
using namespace genetic_algorithm;

int main()
{
    /*
    * Define the fitness function that will be used in the GA (assuming fitness maximization).
    * The fitness function must be thread-safe.
    * The fitness function must take a vector of geneType (size_t for the integerGA),
    * and return the fitness vector (a vector with 1 element for single-objective problems).
    * In this example, the goal is to generate a predefined string using the integerGA.
    */
    auto fitnessFunction = [](const vector<size_t>& chrom) -> vector<double>
    {
        static const string target = "Lorem ipsum dolor sit amet, consectetur adipiscing elit. Pellentesque gravida ut ipsum at tincidunt.";

        double fitness = 0.0;
        for (size_t i = 0; i < chrom.size(); i++)
        {
            fitness += double(char(chrom[i] + 32) == target[i]);
        }

        return vector<double>(1, fitness);
    };

    /* Create the GA using the appropriate chromosome length (length of the string) and the fitness function. */
    size_t str_len = 100;
    size_t base = 96;                               /* Number of values a gene can take (number of valid characters in this case). */
    IntegerGA GA(str_len, fitnessFunction, base);
    //GA.mode(IntegerGA::Mode::single_objective);   /* Default. */


    /* Set other parameters of the GA. */
    GA.population_size(250);
    GA.crossover_rate(0.8);
    /* The default mutation rate (1/chrom_len) is used. */
    GA.swap_rate(0.3);          /* Per-candidate probability of a single-swap mutation. Unique to the IntegerGA. */
    GA.inversion_rate(0.2);     /* Per-candidate probability of an inversion during mutation. Unique to the IntegerGA. */
    GA.max_gen(1000);

    GA.selection_method(IntegerGA::SogaSelection::tournament);
    GA.crossover_method(IntegerGA::CrossoverMethod::uniform);
    /* The default mutation method is used. */

    /* Set an early-stop condition (optional). */
    GA.stop_condition(IntegerGA::StopCondition::fitness_value); /* Stop once a solution dominates the reference fitness vector (fitness maximization). */
    GA.fitness_threshold({ 99.9 });                             /* We know the best (highest) possible fitness value is 100. */


    /* Run the GA. */
    auto sols = GA.run();


    /* Print the results. */
    cout << "The best solutions found are:\n";
    for (const auto& sol : sols)
    {
        cout << " ";
        for (const auto& gene : sol.chromosome)
        {
            cout << char(gene + 32);
        }
        cout << " Fitness value: " << sol.fitness[0] << "\n";
    }


    getchar();
    return 0;
}