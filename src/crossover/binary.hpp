/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CROSSOVER_BINARY_HPP
#define GA_CROSSOVER_BINARY_HPP

#include "crossover_base.decl.hpp"
#include "../population/candidate.hpp"
#include "../encoding/gene_types.hpp"
#include "../utility/bounded_value.hpp"
#include <cstddef>

/** Predefined crossover operators for the binary encoded genetic algorithm. */
namespace genetic_algorithm::crossover::binary
{
    /**
    * Standard single-point crossover operator for the binary encoded %GA.
    * 
    * A random position is selected in the chromosomes as the crossover point,
    * and the genes before this crossover point are swapped between the parents
    * in order to create the child solutions.
    */
    class SinglePoint final : public Crossover<BinaryGene>
    {
    public:
        using Crossover::Crossover;
    private:
        CandidatePair<GeneType> crossover(const GA<GeneType>& ga, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const override;
    };

    /**
    * Two-point crossover operator for the binary encoded %GA.
    * 
    * 2 random points are selected in the chromosomes as the crossover points,
    * and the genes between these 2 crossover points are swapped between the parents in
    * order to create the child solutions.
    * This operation is effectively the same as performing 2 consecutive single-point
    * crossovers on the parents.
    */
    class TwoPoint final : public Crossover<BinaryGene>
    {
    public:
        using Crossover::Crossover;
    private:
        CandidatePair<GeneType> crossover(const GA<GeneType>& ga, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const override;
    };

    /**
    * General N-point crossover operator for the binary encoded %GA.
    * 
    * N random points are selected in the chromosomes as the crossover points for
    * performing the crossover.
    * This operation is effectively the same as performing N consecutive single-point
    * crossovers on the parents to generate the child solutions.
    */
    class NPoint final : public Crossover<BinaryGene>
    {
    public:
        /**
        * Create an N-point crossover operator.
        * 
        * @param n The number of crossover points. Must be at least 1.
        */
        constexpr explicit NPoint(Positive<size_t> n) noexcept :
            n_(n)
        {}

        /**
        * Create an N-point crossover operator.
        *
        * @param pc The crossover probability. Must be in the closed interval [0.0, 1.0].
        * @param n The number of crossover points. Must be at least 1.
        */
        constexpr NPoint(Probability pc, Positive<size_t> n) noexcept :
            Crossover(pc), n_(n)
        {}

        /**
        * Set the number of crossover points used in for the crossovers.
        * The number of crossover points can't be 0, and all values greater than the chromosome
        * length will be treated the same, as if they are equal to the chromosome length.
        * 
        * @param n The number of crossover points to use.
        */
        constexpr void num_crossover_points(Positive<size_t> n) noexcept { n_ = n; }

        /** @returns The number of crossover points used. */
        [[nodiscard]]
        constexpr size_t num_crossover_points() const noexcept { return n_; };

    private:
        CandidatePair<GeneType> crossover(const GA<GeneType>& ga, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const override;

        Positive<size_t> n_;
    };

    /**
    * Uniform crossover operator for the binary encoded %GA.
    * 
    * Each pair of genes of the chromosomes are swapped with a set probability
    * between the parents to create the child solutions.
    */
    class Uniform final : public Crossover<BinaryGene>
    {
    public:
        /** Create a uniform crossover operator using the default crossover and swap rates. */
        constexpr Uniform() noexcept = default;

        /**
        * Create a uniform crossover operator.
        *
        * @param pc The crossover probability. Must be in the closed interval [0.0, 1.0].
        * @param swap_prob The probability of swapping each pair of genes between the 2 parents.
        *   Must be in the closed interval [0.0, 1.0].
        */
        constexpr explicit Uniform(Probability pc, Probability swap_prob = 0.5) noexcept :
            Crossover(pc), ps_(swap_prob)
        {}

        /**
        * Set the swap probability used in the crossovers.
        * The swap probability is the probability of swapping a given pair of genes
        * between the parents.
        *
        * @param swap_prob The probability of swapping each pair of genes between the 2 parents.
        *   Must be in the closed interval [0.0, 1.0].
        */
        constexpr void swap_probability(Probability ps) noexcept { ps_ = ps; }

        /** @returns The swap probability used for the crossovers. */
        [[nodiscard]]
        constexpr Probability swap_probability() const noexcept { return ps_; }

    private:
        CandidatePair<GeneType> crossover(const GA<GeneType>& ga, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const override;

        Probability ps_ = 0.5;
    };

} // namespace genetic_algorithm::crossover::binary

#endif // !GA_CROSSOVER_BINARY_HPP