/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "permutation.hpp"
#include "crossover_base.hpp"
#include "crossover_dtl.hpp"
#include "../population/candidate.hpp"
#include "../utility/rng.hpp"
#include "../utility/utility.hpp"
#include <vector>
#include <utility>
#include <cstddef>

namespace gapp::crossover::perm
{
    auto Order1::crossover(const GA<GeneType>&, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const -> CandidatePair<GeneType>
    {
        GA_ASSERT(parent1.chromosome.size() == parent2.chromosome.size(), "Mismatching parent chromosome lengths.");
        /* The genes of the parent chromosomes must be unique. */

        const size_t chrom_len = parent1.chromosome.size();

        if (chrom_len < 2) return { parent1, parent2 };

        const size_t length = rng::randomInt(1_sz, chrom_len - 1_sz);
        const size_t first = rng::randomInt(0_sz, chrom_len - length);
        const size_t last = first + length;

        auto child1 = dtl::order1CrossoverImpl(parent1, parent2, first, last);
        auto child2 = dtl::order1CrossoverImpl(parent2, parent1, first, last);

        return { std::move(child1), std::move(child2) };
    }

    auto Order2::crossover(const GA<GeneType>&, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const -> CandidatePair<GeneType>
    {
        GA_ASSERT(parent1.chromosome.size() == parent2.chromosome.size(), "Mismatching parent chromosome lengths.");
        /* The genes of the parent chromosomes must be unique. */

        const size_t chrom_len = parent1.chromosome.size();

        if (chrom_len < 2) return { parent1, parent2 };

        const size_t length = rng::randomInt(1_sz, chrom_len - 1_sz);
        const size_t first = rng::randomInt(0_sz, chrom_len - length);
        const size_t last = first + length;

        auto child1 = dtl::order2CrossoverImpl(parent1, parent2, first, last);
        auto child2 = dtl::order2CrossoverImpl(parent2, parent1, first, last);

        return { std::move(child1), std::move(child2) };
    }

    auto Position::crossover(const GA<GeneType>&, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const -> CandidatePair<GeneType>
    {
        GA_ASSERT(parent1.chromosome.size() == parent2.chromosome.size(), "Mismatching parent chromosome lengths.");
        /* The genes of the parent chromosomes must be unique. */

        const size_t chrom_len = parent1.chromosome.size();

        if (chrom_len < 2) return { parent1, parent2 };

        const size_t ns = rng::randomInt(1_sz, chrom_len - 1_sz);
        const auto idxs = rng::sampleUnique(0_sz, chrom_len, ns);

        auto child1 = dtl::positionCrossoverImpl(parent1, parent2, idxs);
        auto child2 = dtl::positionCrossoverImpl(parent2, parent1, idxs);

        return { std::move(child1), std::move(child2) };
    }

    auto Cycle::crossover(const GA<GeneType>&, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const -> CandidatePair<GeneType>
    {
        GA_ASSERT(parent1.chromosome.size() == parent2.chromosome.size(), "Mismatching parent chromosome lengths.");
        /* The genes of the parent chromosomes must be unique. */

        const size_t chrom_len = parent1.chromosome.size();

        if (chrom_len < 2) return { parent1, parent2 };

        return dtl::cycleCrossoverImpl(parent1, parent2);
    }

    auto Edge::crossover(const GA<GeneType>&, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const -> CandidatePair<GeneType>
    {
        GA_ASSERT(parent1.chromosome.size() == parent2.chromosome.size(), "Mismatching parent chromosome lengths.");
        /* The genes of the parent chromosomes must be unique. */

        const size_t chrom_len = parent1.chromosome.size();

        if (chrom_len < 2) return { parent1, parent2 };

        auto child1 = dtl::edgeCrossoverImpl(parent1, parent2);
        auto child2 = dtl::edgeCrossoverImpl(parent2, parent1);

        return { std::move(child1), std::move(child2) };
    }

    auto PMX::crossover(const GA<GeneType>&, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const -> CandidatePair<GeneType>
    {
        GA_ASSERT(parent1.chromosome.size() == parent2.chromosome.size(), "Mismatching parent chromosome lengths.");
        /* The genes of the parent chromosomes must be unique. */

        const size_t chrom_len = parent1.chromosome.size();
        
        if (chrom_len < 2) return { parent1, parent2 };

        const size_t range_len = rng::randomInt(1_sz, chrom_len - 1);

        const size_t first = rng::randomInt(0_sz, chrom_len - range_len);
        const size_t last = first + range_len;

        auto child1 = dtl::pmxCrossoverImpl(parent1, parent2, first, last);
        auto child2 = dtl::pmxCrossoverImpl(parent2, parent1, first, last);

        return { std::move(child1), std::move(child2) };
    }

} // namespace gapp::crossover::perm