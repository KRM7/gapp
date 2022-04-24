/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CROSSOVER_INTEGER_HPP
#define GA_CROSSOVER_INTEGER_HPP

#include "crossover_base.decl.hpp"
#include "../population/candidate.hpp"
#include <cstddef>

/** Predefined crossover operators for the integer encoded genetic algorithms (IntegerGA). */
namespace genetic_algorithm::crossover::integer
{
    /**
    * Standard single-point crossover operator for the integer encoded algorithms.
    * A random crossover point (locus) is selected and the genes before the locus are swapped
    * between the parents to create the children.
    */
    class SinglePoint final : public Crossover<size_t>
    {
    public:
        using Crossover::Crossover;
    private:
        CandidatePair<GeneType> crossover(const GaInfo& ga, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const override;
    };

    /**
    * Two-point crossover operator for the integer encoded algorithms.
    * Two random crossover points are selected, and the genes between the 2 point are
    * swapped between the parents in order to create the children.
    * (Works as if 2 consecutive single-point crossovers were performed on the parents.)
    */
    class TwoPoint final : public Crossover<size_t>
    {
    public:
        using Crossover::Crossover;
    private:
        CandidatePair<GeneType> crossover(const GaInfo& ga, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const override;
    };

    /**
    * Uniform crossover operator for the integer encoded algorithms.
    * Each pair of genes of the chromosomes are swapped with 0.5 probability between the parents to create the children.
    */
    class Uniform final : public Crossover<size_t>
    {
    public:
        using Crossover::Crossover;
    private:
        CandidatePair<GeneType> crossover(const GaInfo& ga, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const override;
    };

} // namespace genetic_algorithm::crossover::integer

#endif // !GA_CROSSOVER_INTEGER_HPP