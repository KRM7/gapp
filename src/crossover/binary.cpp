/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "binary.hpp"
#include "crossover_dtl.hpp"
#include "../core/candidate.hpp"
#include "../utility/rng.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <vector>
#include <utility>
#include <cstddef>

namespace gapp::crossover::binary
{
    auto SinglePoint::crossover(const GA<GeneType>&, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const -> CandidatePair<GeneType>
    {
        GAPP_ASSERT(parent1.chromosome.size() == parent2.chromosome.size(), "Mismatching parent chromosome lengths.");

        const size_t chrom_len = parent1.chromosome.size();
        const size_t crossover_point = rng::randomInt(0_sz, chrom_len);

        return dtl::singlePointCrossoverImpl(parent1, parent2, crossover_point);
    }

    auto TwoPoint::crossover(const GA<GeneType>&, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const -> CandidatePair<GeneType>
    {
        GAPP_ASSERT(parent1.chromosome.size() == parent2.chromosome.size(), "Mismatching parent chromosome lengths.");

        const size_t chrom_len = parent1.chromosome.size();
        const size_t crossover_point1 = rng::randomInt(0_sz, chrom_len);
        const size_t crossover_point2 = rng::randomInt(0_sz, chrom_len);

        return dtl::twoPointCrossoverImpl(parent1, parent2, { crossover_point1, crossover_point2 });
    }

    auto NPoint::crossover(const GA<GeneType>&, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const -> CandidatePair<GeneType>
    {
        GAPP_ASSERT(parent1.chromosome.size() == parent2.chromosome.size(), "Mismatching parent chromosome lengths.");

        const size_t chrom_len = parent1.chromosome.size();
        const size_t num_crossover_points = std::min(size_t(n_), chrom_len);

        auto crossover_points = rng::sampleUnique(0_sz, chrom_len, num_crossover_points);

        return dtl::nPointCrossoverImpl(parent1, parent2, std::move(crossover_points));
    }

    auto Uniform::crossover(const GA<GeneType>&, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const -> CandidatePair<GeneType>
    {
        GAPP_ASSERT(parent1.chromosome.size() == parent2.chromosome.size(), "Mismatching parent chromosome lengths.");

        const size_t chrom_len = parent1.chromosome.size();
        const size_t num_swapped = rng::randomBinomial(chrom_len, ps_);
        const auto swapped_indices = rng::sampleUnique(0_sz, chrom_len, num_swapped);

        Candidate child1{ parent1 }, child2{ parent2 };

        for (const auto& idx : swapped_indices)
        {
            using std::swap;
            swap(child1.chromosome[idx], child2.chromosome[idx]);
        }

        return { std::move(child1), std::move(child2) };
    }

} // namespace gapp::crossover::binary