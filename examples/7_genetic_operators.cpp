/* Example showing the usage of the genetic operators in the GAs. */

#include "gapp.hpp"
#include <algorithm>
#include <iostream>
#include <cassert>

using namespace gapp;

class MyCrossover : public crossover::Crossover<PermutationGene>
{
public:
    using Crossover::Crossover;

    CandidatePair<GeneType> crossover(const GA<GeneType>&, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const override
    {
        auto child1 = parent1;
        auto child2 = parent2;

        // perform the crossover ...

        return { std::move(child1), std::move(child2) };
    }
};

class MyMutation : public mutation::Mutation<PermutationGene>
{
public:
    using Mutation::Mutation;

    void mutate(const GA<GeneType>&, const Candidate<GeneType>&, Chromosome<GeneType>& chromosome) const override
    {
        if (rng::randomReal() < mutation_rate())
        {
            std::reverse(chromosome.begin(), chromosome.end());
        }
    }
};

int main()
{
    PermutationGA ga;
    ga.solve(problems::TSP52{}); // run using the default crossover and mutation methods

    // using other crossover, mutation operators

    ga.crossover_method(crossover::perm::Edge{}); // using the default crossover probability
    ga.mutation_method(mutation::perm::Inversion{ /* mutation_rate = */ 0.3 });

    std::cout << "The default crossover probability is "<< ga.crossover_rate() << ".\n";

    // changing the crossover and mutation probabilities

    ga.crossover_method(crossover::perm::Edge{ /* crossover_rate = */ 0.92 });
    std::cout << "The crossover probability is " << ga.crossover_rate() << ".\n";

    ga.crossover_rate(0.71);
    ga.mutation_rate(0.1);

    std::cout << "The crossover probability is " << ga.crossover_rate() << ".\n";
    std::cout << "The mutation probability is " << ga.mutation_rate() << ".\n";

    // user defined crossover and mutation methods

    ga.crossover_method(MyCrossover{ /* crossover_rate = */ 0.123 });
    ga.mutation_method(MyMutation{ /* mutation_rate = */ 0.456 });

    // using a repair function

    ga.repair_function([](const GA<PermutationGene>&, const Chromosome<PermutationGene>& chrom)
    {
        auto new_chrom = chrom;
        std::swap(new_chrom.front(), new_chrom.back());
        return new_chrom;
    });
}
