/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_MUTATION_PERMUTATION_HPP
#define GA_MUTATION_PERMUTATION_HPP

#include "mutation_base.hpp"

/** Predefined mutation operators for the permutation encoded genetic algorithm. */
namespace genetic_algorithm::mutation::perm
{
    /**
    * Inversion mutation operator for the permutation encoded genetic algorithm.
    * Each individual is mutated with the specified mutation probability.
    * In the mutated individuals, a random range of genes is chosen and reversed.
    */
    class Inversion : public Mutation<size_t>
    {
    public:
        using Mutation::Mutation;
    private:
        void mutate(const GA<size_t>& ga, Candidate<size_t>& candidate) const override;
    };

    /**
    * Single swap/swap2 mutation operator for the permutation encoded genetic algorithm.
    * Each individual is mutated with the set mutation probability.
    * Two randomly selected genes are swapped in the mutated individuals.
    */
    class Swap2 final : public Mutation<size_t>
    {
    public:
        using Mutation::Mutation;
    private:
        void mutate(const GA<size_t>& ga, Candidate<size_t>& candidate) const override;
    };

    /**
    * Swap-3 mutation operator for the permutation encoded genetic algorithm.
    * Each individual is mutated with the set mutation probability.
    * Three randomly selected genes are swapped in the mutated individuals (0-1-2)->(2-0-1).
    */
    class Swap3 final : public Mutation<size_t>
    {
    public:
        using Mutation::Mutation;
    private:
        void mutate(const GA<size_t>& ga, Candidate<size_t>& candidate) const override;
    };

    /**
    * Shuffle/scramble mutation operator for the permutation encoded genetic algorithm.
    * Each individual is mutated with the set mutation probability.
    * In the mutated candidates, a random range of genes is selected and then shuffled.
    */
    class Shuffle final : public Mutation<size_t>
    {
    public:
        using Mutation::Mutation;
    private:
        void mutate(const GA<size_t>& ga, Candidate<size_t>& candidate) const override;
    };

    /**
    * Shift/slide mutation operator for the permutation encoded genetic algorithm.
    * Each individual is mutated with the set mutation proability.
    * In the mutated candidates, a random range of genes is selected and then moved to a
    * different position in the chromosome.
    */
    class Shift final : public Mutation<size_t>
    {
    public:
        using Mutation::Mutation;
    private:
        void mutate(const GA<size_t>& ga, Candidate<size_t>& candidate) const override;
    };

} // namespace genetic_algorithm::mutation::perm

#endif // !GA_MUTATION_PERMUTATION_HPP