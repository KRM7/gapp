/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CROSSOVER_BINARY_HPP
#define GA_CROSSOVER_BINARY_HPP

#include "crossover_base.decl.hpp"
#include "../population/candidate.hpp"
#include "../encoding/gene_types.hpp"

/** Predefined crossover operators for the binary encoded genetic algorithms (BinaryGA). */
namespace genetic_algorithm::crossover::binary
{
    /**
    * Standard single-point crossover operator for the binary encoded algorithms. \n
    * A random crossover point (locus) is selected and the genes before the locus are swapped
    * between the parents to create the children.
    */
    class SinglePoint final : public Crossover<BinaryGene>
    {
    public:
        using Crossover::Crossover;
    private:
        CandidatePair<GeneType> crossover(const GA<GeneType>& ga, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const override;
    };

    /**
    * Two-point crossover operator for the binary encoded algorithms. \n
    * Two random crossover points are selected, and the genes between the 2 point are
    * swapped between the parents in order to create the children. \n
    * (Works as if 2 consecutive single-point crossovers were performed on the parents.)
    */
    class TwoPoint final : public Crossover<BinaryGene>
    {
    public:
        using Crossover::Crossover;
    private:
        CandidatePair<GeneType> crossover(const GA<GeneType>& ga, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const override;
    };

    /**
    * General N-point crossover operator for the binary encoded algorithms. \n
    * The result of the crossover is equivalent to performing N consecutive single-point
    * crossovers on the parents using different crossover points.
    */
    class NPoint final : public Crossover<BinaryGene>
    {
    public:
        /**
        * Create an N-point crossover operator.
        * 
        * @param n The number of crossover points. Must be greater than 0.
        */
        explicit NPoint(size_t n);

        /**
        * Create an N-point crossover operator.
        *
        * @param pc The crossover probability. Must be in the closed interval [0.0, 1.0].
        * @param n The number of crossover points. Must be greater than 0.
        */
        NPoint(Probability pc, size_t n);

        /**
        * Set the number of crossover points used for the operator to @p n. \n
        * The number of crossover points can't be 0, and all values greater than the chromosome
        * length are treated the same (as if n == chrom_len).
        * 
        * @param n The number of crossover points.
        */
        void num_crossover_points(size_t n);

        /** @returns The number of crossover points set. */
        [[nodiscard]]
        size_t num_crossover_points() const noexcept { return n_; };

    private:
        CandidatePair<GeneType> crossover(const GA<GeneType>& ga, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const override;

        size_t n_;
    };

    /**
    * Uniform crossover operator for the binary encoded algorithms. \n
    * Each pair of genes of the chromosomes are swapped with a set probability between the parents to create the children.
    */
    class Uniform final : public Crossover<BinaryGene>
    {
    public:
        /**
        * Create a uniform crossover operator.
        *
        * @param pc The crossover probability. Must be in the closed interval [0.0, 1.0].
        * @param swap_prob The probability of swapping each pair of genes between the 2 parents. Must be in the closed interval [0.0, 1.0].
        */
        explicit Uniform(Probability pc, Probability swap_prob = 0.5) noexcept;

        /**
        * Set the probability of swapping each pair of genes between the parents during
        * the crossovers to @p ps. This value must be in the closed interval [0.0, 1.0].
        * 
        * @param swap_prob The probability of swapping each pair of genes between the 2 parents.
        */
        void swap_probability(Probability ps);

        /** @returns The swap probability set. */
        [[nodiscard]]
        Probability swap_probability() const noexcept { return ps_; }

    private:
        CandidatePair<GeneType> crossover(const GA<GeneType>& ga, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const override;

        Probability ps_;
    };

} // namespace genetic_algorithm::crossover::binary

#endif // !GA_CROSSOVER_BINARY_HPP