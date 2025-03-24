/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_ENCODING_INTEGER_HPP
#define GAPP_ENCODING_INTEGER_HPP

#include "gene_types.hpp"
#include "../core/ga_base.hpp"
#include "../core/candidate.hpp"
#include "../crossover/integer.hpp"
#include "../mutation/integer.hpp"
#include <cstddef>

namespace gapp
{
    template<>
    struct GaTraits<IntegerGene>
    {
        using DefaultCrossover = crossover::integer::TwoPoint;
        using DefaultMutation = mutation::integer::Uniform;

        static constexpr Probability defaultMutationRate(size_t chrom_len) noexcept { return 1.0 / chrom_len; }

        static Chromosome<IntegerGene> randomChromosome(size_t chrom_len, BoundsView<IntegerGene> bounds);
    };

    /**
    * Integer-encoded genetic algorithm class. This is the main solver
    * that should be used for integer-encoded objective functions.
    * 
    * Similar to the binary-encoded %GA, but the values of a genes can be any integer
    * in a closed interval, not just 0 or 1. The concrete interval is specified by the gene bounds.
    * 
    * @see BinaryGA
    */
    class IntegerGA final : public GA<IntegerGene>
    {
    public:
        using GA::GA;
    };

} // namespace gapp

#endif // !GAPP_ENCODING_INTEGER_HPP
