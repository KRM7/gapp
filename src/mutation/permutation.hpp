/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_MUTATION_PERMUTATION_HPP
#define GA_MUTATION_PERMUTATION_HPP

#include "mutation_base.hpp"
#include "../encoding/gene_types.hpp"
#include "../utility/bounded_value.hpp"

/** Predefined mutation operators for the permutation encoded genetic algorithm. */
namespace genetic_algorithm::mutation::perm
{
    /**
    * %Inversion mutation operator for the permutation encoded %GA.
    * 
    * Each individual is mutated with the specified mutation probability.
    * In the mutated individuals, a randomly selected range of genes are reversed.
    * 
    * The operator has a single parameter (@p range_max) that specifies the maximum
    * length of the reversed ranges relative to the chromosome length.
    */
    class Inversion final : public Mutation<PermutationGene>
    {
    public:
        /**
        * Create an inversion mutation operator.
        * 
        * @param pm The mutation probability. Must be in the closed interval [0.0, 1.0].
        * @param range_max The maximum length of the reversed ranges. Must be in the closed interval [0.0, 1.0].
        */
        constexpr explicit Inversion(Probability pm, Normalized<double> range_max = 0.75) noexcept :
            Mutation(pm), range_max_(range_max)
        {}

        /**
        * Set the maximum length of the ranges that can be selected to be reversed
        * by the operator. The parameter specifies the maximum range length relative
        * to the overall chromosome length of a candidate.
        * 
        * @param rm The maximum length of the reversed ranges. Must be in the closed interval [0.0, 1.0].
        */
        constexpr void range_max(Normalized<double> rm) noexcept { range_max_ = rm; }

        /** @returns The maximum length of the reversed ranges. */
        [[nodiscard]]
        constexpr double range_max() const noexcept { return range_max_; }

    private:
        void mutate(const GA<GeneType>& ga, Candidate<GeneType>& candidate) const override;

        Normalized<double> range_max_;
    };

    /**
    * Single swap/swap2 mutation operator for the permutation encoded %GA.
    * 
    * Each candidate solution is mutated with the set mutation probability. In the mutated
    * candidates, two distinct genes are randomly selected and then swapped.
    */
    class Swap2 final : public Mutation<PermutationGene>
    {
    public:
        using Mutation::Mutation;
    private:
        void mutate(const GA<GeneType>& ga, Candidate<GeneType>& candidate) const override;
    };

    /**
    * Swap-3 mutation operator for the permutation encoded %GA.
    * 
    * Each candidate solution is mutated with the set mutation probability.
    * In the mutated candidates, 3 distinct genes are randomly selected and then reordered
    * as: (a-b-c) -> (c-a-b).
    */
    class Swap3 final : public Mutation<PermutationGene>
    {
    public:
        using Mutation::Mutation;
    private:
        void mutate(const GA<GeneType>& ga, Candidate<GeneType>& candidate) const override;
    };

    /**
    * %Shuffle/scramble mutation operator for the permutation encoded %GA.
    * 
    * Each candidate solution is mutated with the set mutation probability.
    * In the mutated candidates, a random range of genes is selected and then randomly shuffled.
    * 
    * The operator has a single parameter (@p range_max) that specifies the maximum
    * length of the shuffled ranges relative to the chromosome length.
    * 
    * @note There is a possbility that the shuffled chromosome will be the same as the original
    *   chromosome, so the probability of a chromosome being changed won't be exactly equal to the
    *   set mutation probability (it will be slightly lower).
    */
    class Shuffle final : public Mutation<PermutationGene>
    {
    public:
        /**
        * Create a shuffle mutation operator.
        *
        * @param pm The mutation probability. Must be in the closed interval [0.0, 1.0].
        * @param range_max The maximum length of the shuffled ranges. Must be in the closed interval [0.0, 1.0].
        */
        constexpr explicit Shuffle(Probability pm, Normalized<double> range_max = 0.5) noexcept :
            Mutation(pm), range_max_(range_max)
        {}

        /**
        * Set the maximum length of the ranges that can be selected to be shuffled
        * by the operator. The parameter specifies the maximum range length relative
        * to the overall chromosome length of a candidate.
        *
        * @param rm The maximum length of the shuffled ranges. Must be in the closed interval [0.0, 1.0].
        */
        constexpr void range_max(Normalized<double> rm) noexcept { range_max_ = rm; }

        /** @returns The maximum length of the shuffled ranges. */
        [[nodiscard]]
        constexpr double range_max() const noexcept { return range_max_; }

    private:
        void mutate(const GA<GeneType>& ga, Candidate<GeneType>& candidate) const override;

        Normalized<double> range_max_;
    };

    /**
    * %Shift/slide mutation operator for the permutation encoded %GA.
    * 
    * Each candidate solution is mutated with the set mutation proability.
    * In the mutated candidates, a random range of genes is selected and then moved to a
    * different position in the chromosome.
    * 
    * The operator has a single parameter (@p range_max) that specifies the maximum
    * length of the moved ranges relative to the chromosome length.
    */
    class Shift final : public Mutation<PermutationGene>
    {
    public:
        /**
        * Create a shift mutation operator.
        *
        * @param pm The mutation probability. Must be in the closed interval [0.0, 1.0].
        * @param range_max The maximum length of the moved ranges. Must be in the closed interval [0.0, 1.0].
        */
        constexpr explicit Shift(Probability pm, Normalized<double> range_max = 0.75) noexcept :
            Mutation(pm), range_max_(range_max)
        {}

        /**
        * Set the maximum length of the ranges that can be selected to be moved
        * by the operator. The parameter specifies the maximum range length relative
        * to the overall chromosome length of a candidate.
        *
        * @param rm The maximum length of the moved ranges. Must be in the closed interval [0.0, 1.0].
        */
        constexpr void range_max(Normalized<double> rm) noexcept { range_max_ = rm; }

        /** @returns The maximum length of the moved ranges. */
        [[nodiscard]]
        constexpr double range_max() const noexcept { return range_max_; }

    private:
        void mutate(const GA<GeneType>& ga, Candidate<GeneType>& candidate) const override;

        Normalized<double> range_max_;
    };

} // namespace genetic_algorithm::mutation::perm

#endif // !GA_MUTATION_PERMUTATION_HPP