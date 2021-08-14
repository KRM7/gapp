/* Example showing the usage of user-defined genetic operators in the algorithms. */

#include "../src/genetic_algorithm.h"
#include "../test/fitness_functions.h"  /* For the fitness function that will be used in the example. */

#include <cstdio>
#include <cstddef>
#include <iostream>
#include <iomanip>
#include <algorithm>

using namespace std;
using namespace genetic_algorithm;

/* 
* The example will show user-defined genetic operators used with a single-objective binary GA.
* (The genetic operators implemented here for the GA are already implemented and usable in the BinaryGA.)
* The genetic operators that can be defined by the user are:
*   selection (only for SOGAs)
*   crossover
*   mutation
*   repair (not shown here)
*/

BinaryGA::Candidate mySelection(const BinaryGA::Population&);
BinaryGA::CandidatePair myCrossover(const BinaryGA::Candidate&, const BinaryGA::Candidate&, double);
void myMutation(BinaryGA::Candidate&, double);

int main()
{
    /*
    * Define the fitness function that will be used in the GA. (Assuming fitness maximization.)
    * The fitness function must take a vector of geneType (char for the BinaryGA),
    * and return the fitness vector (a vector with 1 element for single-objective problems).
    * In this example, it will be the Rastrigin function in 10 dimensions, already used in the benchmarks and earlier examples.
    */
    size_t num_vars = 10;
    Rastrigin rastriginFunction(num_vars);

    /* Create the GA using the appropriate chromosome length and the fitness function. */
    size_t chrom_len = rastriginFunction.num_vars * rastriginFunction.var_bits;
    BinaryGA GA(chrom_len, rastriginFunction);

    /* Set the parameters of the GA. */
    GA.population_size(400);
    GA.crossover_rate(0.75);
    GA.mutation_rate(0.015);    /* The mutation rate is set to 1/chrom_len by default. */
    GA.max_gen(500);

    /* Set the genetic operators to the custom ones defined below. */
    GA.selection_method(mySelection);
    GA.crossover_method(myCrossover);
    GA.mutation_method(myMutation);

    /* Run the GA. */
    auto sols = GA.run();

    /* Print results. */
    cout << " The number of fitness function evals performed: " << GA.num_fitness_evals() << "\n";
    cout << " The best solutions found:\n";
    for (const auto& sol : sols)
    {
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

/*
* Simple binary tournament selection.
* The selection operator can only by user-defined for the single-objective GAs, the multiobjective algorithms (NSGA-II and NSGA-III)
* always use their own selection operators.
*/
BinaryGA::Candidate mySelection(const BinaryGA::Population& pop)
{
    size_t idx1 = rng::generateRandomIdx(pop.size());
    size_t idx2 = rng::generateRandomIdx(pop.size());

    return (pop[idx1].fitness > pop[idx2].fitness) ? pop[idx1] : pop[idx2]; /* The algorithm assumes fitness maximization. */
}

/* Simple single-point crossover for binary chromosomes. */
BinaryGA::CandidatePair myCrossover(const BinaryGA::Candidate& parent1, const BinaryGA::Candidate& parent2, double crossover_rate)
{
    BinaryGA::Candidate child1(parent1), child2(parent2);   /* Children are same as parents (is_evaluated = true). */

    /* Perform crossover with crossover_rate probability. */
    if (rng::generateRandomDouble() < crossover_rate)   /* The rng namespace has functions for random number generation if needed. */
    {
        size_t cx_point = rng::generateRandomInt(size_t{ 1 }, parent1.chromosome.size() - 1);
        for (size_t i = 0; i < cx_point; i++)
        {
            child1.chromosome[i] = parent2.chromosome[i];
            child2.chromosome[i] = parent1.chromosome[i];
        }
        child1.is_evaluated = false;    /* The children were changed, they will need to be evaluated. */
        child2.is_evaluated = false;
    }

    return make_pair(child1, child2);
}

/* Simple mutation for binary chromosomes. (The genes are encoded as char types for the binary GA.) */
void myMutation(BinaryGA::Candidate& child, double mutation_rate)
{
    for (auto& gene : child.chromosome)
    {
        /* Flip each gene with mutation_rate probability. */
        if (rng::generateRandomDouble() < mutation_rate)
        {
            gene = bool(gene) ? char(0) : char(1);
            child.is_evaluated = false;
        }
    }
}