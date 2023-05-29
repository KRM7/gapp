/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CROSSOVER_PERMUTATION_HPP
#define GA_CROSSOVER_PERMUTATION_HPP

#include "crossover_base.decl.hpp"
#include "../population/candidate.hpp"
#include "../encoding/gene_types.hpp"

/** Predefined crossover operators for the permutation encoded %GA. */
namespace genetic_algorithm::crossover::perm
{
    /**
    * Order (OX1) crossover operator for the permutation encoded %GA.
    * 
    * In order to create a child, a range of genes is randomly selected from parent1 and copied directly
    * to the child into the same position, while the remaining genes are filled in from parent2 in the
    * order they appear in, starting at the end of the randomly selected range.
    * 
    * The second child is created by repeating this process with the roles of the two parents swapped.
    * The same range of genes is used for the directly copied genes.
    */
    class Order1 final : public Crossover<PermutationGene>
    {
    public:
        using Crossover::Crossover;
    private:
        CandidatePair<GeneType> crossover(const GA<GeneType>& ga, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const override;
    };

    /**
    * Order based (OX2) crossover operator for the permutation encoded %GA.
    * This crossover operator is a slightly modified version of the Order1 operator.
    * 
    * In order to create a child, a range of genes is randomly selected from parent1 and copied directly
    * to the child into the same position, while the remaining genes are filled in from parent 2 in the
    * order they appear in, starting at the beginning of the chromosome.
    * 
    * The second child is created by repeating this process with the roles of the two parents swapped.
    * The same range of genes is used for the directly copied genes.
    */
    class Order2 final : public Crossover<PermutationGene>
    {
    public:
        using Crossover::Crossover;
    private:
        CandidatePair<GeneType> crossover(const GA<GeneType>& ga, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const override;
    };

    /**
    * %Position/position-based (POS) crossover operator for the permutation encoded %GA.
    * This crossover operator is a modification of the Order1 operator.
    *
    * In order to create a child, a random number of positions are selected randomly from parent1
    * (instead of the continuous range selected in the Order1 crossover operator), then these genes are copied
    * directly from parent1 to the child into the same positons. The remaining genes which are still missing from
    * the child are copied from parent2 in the order they appear in, starting at the beginning of the chromosome.
    * 
    * The second child is created by repeating this process with the roles of the two parents swapped,
    * but using the same positions for direct copying that were used to create the first child.
    */
    class Position final : public Crossover<PermutationGene>
    {
    public:
        using Crossover::Crossover;
    private:
        CandidatePair<GeneType> crossover(const GA<GeneType>& ga, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const override;
    };

    /**
    * %Cycle (CX) crossover operator for the permutation encoded %GA.
    * 
    * This operator works by identifying cycles of genes between the parent chromosomes,
    * and building the 2 child solutions from these cycles.
    * Each of the genes in the children appears in the same position in one of the parents.
    */
    class Cycle final : public Crossover<PermutationGene>
    {
    public:
        using Crossover::Crossover;
    private:
        CandidatePair<GeneType> crossover(const GA<GeneType>& ga, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const override;
    };

    /**
    * %Edge assembly/recombination (EAX) crossover operator for the permutation encoded %GA.
    * 
    * The children are created from the parents by trying to keep as many edges
    * present in the parents as possible, and not introducing new edges into the children.
    * This crossover operator is significantly slower than the other implemented operators,
    * but produces good results.
    */
    class Edge final : public Crossover<PermutationGene>
    {
    public:
        using Crossover::Crossover;
    private:
        CandidatePair<GeneType> crossover(const GA<GeneType>& ga, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const override;
    };

    /**
    * Partially mapped (%PMX) crossover operator for the permutation encoded %GA.
    * 
    * Similar to the Order1 crossover, a random range of genes is selected from parent1 and copied
    * directly intto the same positions of the child chromosome, and the remaining genes not yet in
    * the child are filled in from parent2 using a different method from the one used in the Order1
    * crossover.
    * 
    * The second child is created by performing the same process with the roles of the 2 parents swapped.
    */
    class PMX final : public Crossover<PermutationGene>
    {
    public:
        using Crossover::Crossover;
    private:
        CandidatePair<GeneType> crossover(const GA<GeneType>& ga, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const override;
    };

} // namespace genetic_algorithm::crossover::perm

#endif // !GA_CROSSOVER_PERMUTATION_HPP