/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_ENCODING_REAL_HPP
#define GAPP_ENCODING_REAL_HPP

#include "gene_types.hpp"
#include "../core/ga_base.hpp"
#include "../core/ga_traits.hpp"
#include "../core/candidate.hpp"
#include "../crossover/real.hpp"
#include "../mutation/real.hpp"
#include <cstddef>

namespace gapp
{
    template<>
    struct GaTraits<RealGene>
    {
        using DefaultCrossover = crossover::real::Wright;
        using DefaultMutation = mutation::real::Gauss;

        static constexpr Probability defaultMutationRate(size_t chrom_len) noexcept { return 1.0 / chrom_len; }

        static Chromosome<RealGene> randomChromosome(size_t chrom_len, BoundsView<RealGene> bounds);
    };

    /**
    * Real-encoded genetic algorithm class. This is the main solver
    * that should be used for real-encoded objective functions.
    */
    class RCGA final : public GA<RealGene>
    {
    public:
        using GA::GA;
    };

} // namespace gapp

#endif // !GAPP_ENCODING_REAL_HPP
