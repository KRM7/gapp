/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CROSSOVER_INTEGER_HPP
#define GA_CROSSOVER_INTEGER_HPP

#include "crossover_base.decl.hpp"
#include "../population/candidate.hpp"
#include "../encoding/gene_types.hpp"
#include "../utility/bounded_value.hpp"
#include <cstddef>

/** Predefined crossover operators for the integer encoded genetic algorithms (IntegerGA). */
namespace genetic_algorithm::crossover::integer
{
    /**
    * Standard single-point crossover operator for the integer encoded algorithms.
    * A random crossover point (locus) is selected and the genes before the locus are swapped
    * between the parents to create the children.
    */
    class SinglePoint final : public Crossover<IntegerGene>
    {
    public:
        using Crossover::Crossover;
    private:
        CandidatePair<GeneType> crossover(const GA<GeneType>& ga, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const override;
    };

    /**
    * Two-point crossover operator for the integer encoded algorithms.
    * Two random crossover points are selected, and the genes between the 2 point are
    * swapped between the parents in order to create the children.
    * (Works as if 2 consecutive single-point crossovers were performed on the parents.)
    */
    class TwoPoint final : public Crossover<IntegerGene>
    {
    public:
        using Crossover::Crossover;
    private:
        CandidatePair<GeneType> crossover(const GA<GeneType>& ga, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const override;
    };

    /**
    * General N-point crossover operator for the integer encoded algorithms. \n
    * The result of the crossover is equivalent to performing N consecutive single-point
    * crossovers on the parents using different crossover points.
    */
    class NPoint final : public Crossover<IntegerGene>
    {
    public:
        /**
        * Create an N-point crossover operator.
        *
        * @param n The number of crossover points. Must be greater than 0.
        */
        constexpr explicit NPoint(Positive<size_t> n) noexcept :
            n_(n)
        {}

        /**
        * Create an N-point crossover operator.
        *
        * @param pc The crossover probability. Must be in the closed interval [0.0, 1.0].
        * @param n The number of crossover points. Must be greater than 0.
        */
        constexpr NPoint(Probability pc, Positive<size_t> n) noexcept :
            Crossover(pc), n_(n)
        {}

        /**
        * Set the number of crossover points used for the operator to @p n. \n
        * The number of crossover points can't be 0, and all values greater than the chromosome
        * length are treated the same (as if n == chrom_len).
        *
        * @param n The number of crossover points.
        */
        constexpr void num_crossover_points(Positive<size_t> n) noexcept { n_ = n; }

        /** @returns The number of crossover points set. */
        [[nodiscard]]
        constexpr size_t num_crossover_points() const noexcept { return n_; };

    private:
        CandidatePair<GeneType> crossover(const GA<GeneType>& ga, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const override;

        Positive<size_t> n_;
    };

    /**
    * Uniform crossover operator for the integer encoded algorithms. \n
    * Each pair of genes of the chromosomes are swapped with a set probability between the parents to create the children.
    */
    class Uniform final : public Crossover<IntegerGene>
    {
    public:
        /** Create a uniform crossover operator. */
        constexpr Uniform() noexcept = default;

        /**
        * Create a uniform crossover operator.
        *
        * @param pc The crossover probability. Must be in the closed interval [0.0, 1.0].
        * @param ps The probability of swapping each pair of genes between the 2 parents. Must be in the closed interval [0.0, 1.0].
        */
        constexpr explicit Uniform(Probability pc, Probability ps = 0.5) noexcept :
            Crossover(pc), ps_(ps)
        {}

        /**
        * Set the probability of swapping each pair of genes between the parents during
        * the crossovers to @p ps. This value must be in the closed interval [0.0, 1.0].
        *
        * @param swap_prob The probability of swapping each pair of genes between the 2 parents.
        */
        constexpr void swap_probability(Probability ps) noexcept { ps_ = ps; }

        /** @returns The swap probability set. */
        [[nodiscard]]
        constexpr Probability swap_probability() const noexcept { return ps_; }

    private:
        CandidatePair<GeneType> crossover(const GA<GeneType>& ga, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const override;

        Probability ps_ = 0.5;
    };

} // namespace genetic_algorithm::crossover::integer

#endif // !GA_CROSSOVER_INTEGER_HPP