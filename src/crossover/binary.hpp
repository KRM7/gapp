/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CROSSOVER_BINARY_HPP
#define GA_CROSSOVER_BINARY_HPP

#include "crossover_base.decl.hpp"
#include "../population/candidate.hpp"

/** Predefined crossover operators for the binary encoded genetic algorithms (BinaryGA). */
namespace genetic_algorithm::crossover::binary
{
    /**
    * Standard single-point crossover operator for the binary encoded algorithms. \n
    * A random crossover point (locus) is selected and the genes before the locus are swapped
    * between the parents to create the children.
    */
    class SinglePoint final : public Crossover<char>
    {
    public:
        using Crossover::Crossover;
    private:
        CandidatePair<char> crossover(const GaInfo& ga, const Candidate<char>& parent1, const Candidate<char>& parent2) const override;
    };

    /**
    * Two-point crossover operator for the binary encoded algorithms. \n
    * Two random crossover points are selected, and the genes between the 2 point are
    * swapped between the parents in order to create the children. \n
    * (Works as if 2 consecutive single-point crossovers were performed on the parents.)
    */
    class TwoPoint final : public Crossover<char>
    {
    public:
        using Crossover::Crossover;
    private:
        CandidatePair<char> crossover(const GaInfo& ga, const Candidate<char>& parent1, const Candidate<char>& parent2) const override;
    };

    /**
    * Uniform crossover operator for the binary encoded algorithms. \n
    * Every gene of the chromosomes is swapped with 0.5 probability between the parents to create the children.
    */
    class Uniform final : public Crossover<char>
    {
    public:
        using Crossover::Crossover;
    private:
        CandidatePair<char> crossover(const GaInfo& ga, const Candidate<char>& parent1, const Candidate<char>& parent2) const override;
    };

} // namespace genetic_algorithm::crossover::binary

#endif // !GA_CROSSOVER_BINARY_HPP