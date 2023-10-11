/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ENCODING_REAL_HPP
#define GA_ENCODING_REAL_HPP

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
    };

    /**
    * Real-encoded genetic algorithm class. This is the main solver
    * that should be used for real-encoded objective functions.
    */
    class RCGA final : public GA<RealGene>
    {
    public:
        using GA::GA;
    private:
        Candidate<GeneType> generateCandidate() const override;
    };

} // namespace gapp

#endif // !GA_ENCODING_REAL_HPP