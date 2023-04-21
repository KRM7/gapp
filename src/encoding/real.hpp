/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ENCODING_REAL_HPP
#define GA_ENCODING_REAL_HPP

#include "gene_types.hpp"
#include "../core/ga_base.hpp"
#include "../core/ga_traits.hpp"
#include "../population/candidate.hpp"
#include "../crossover/real.hpp"
#include "../mutation/real.hpp"
#include <cstddef>

namespace genetic_algorithm
{
    template<>
    struct GaTraits<RealGene>
    {
        using DefaultCrossover = crossover::real::Wright;
        using DefaultMutation = mutation::real::Gauss;

        static constexpr Probability defaultMutationRate(size_t chrom_len) noexcept { return 1.0 / chrom_len; }
    };

    /**
    * Real coded genetic algorithm. \n
    * Each gene of the chromosomes is a floating point value.
    */
    class RCGA final : public GA<RealGene>
    {
    public:
        using GA::GA;
    private:
        Candidate<GeneType> generateCandidate() const override;
    };

} // namespace genetic_algorithm

#endif // !GA_ENCODING_REAL_HPP