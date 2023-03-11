/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "integer.hpp"
#include "crossover_dtl.hpp"
#include "../core/ga_base.hpp"
#include "../population/candidate.hpp"
#include "../utility/rng.hpp"
#include "../utility/bounded_value.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <vector>
#include <utility>
#include <stdexcept>
#include <cstddef>

namespace genetic_algorithm::crossover::integer
{
    auto SinglePoint::crossover(const GA<GeneType>&, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const -> CandidatePair<GeneType>
    {
        const size_t chrom_len = parent1.chromosome.size();
        const size_t crossover_point = rng::randomInt(0_sz, chrom_len);

        return dtl::singlePointCrossoverImpl(parent1, parent2, crossover_point);
    }

    auto TwoPoint::crossover(const GA<GeneType>&, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const -> CandidatePair<GeneType>
    {
        const size_t chrom_len = parent1.chromosome.size();

        return dtl::twoPointCrossoverImpl(parent1, parent2, { rng::randomInt(0_sz, chrom_len), rng::randomInt(0_sz, chrom_len) });
    }

    auto NPoint::crossover(const GA<GeneType>&, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const -> CandidatePair<GeneType>
    {
        if (parent1.chromosome.size() != parent2.chromosome.size())
        {
            GA_THROW(std::invalid_argument, "The parent chromosomes must be the same length for the n-point crossover.");
        }

        const size_t chrom_len = parent1.chromosome.size();
        const size_t num_cx_points = std::min(size_t(n_), chrom_len);

        auto cx_points = rng::sampleUnique(0_sz, chrom_len, num_cx_points);

        return dtl::nPointCrossoverImpl(parent1, parent2, std::move(cx_points));
    }

    auto Uniform::crossover(const GA<GeneType>&, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const -> CandidatePair<GeneType>
    {
        if (parent1.chromosome.size() != parent2.chromosome.size())
        {
            GA_THROW(std::invalid_argument, "The parent chromosomes must be the same length for the uniform crossover.");
        }
        
        const size_t chrom_len = parent1.chromosome.size();
        const size_t num_swapped = rng::randomBinomial(chrom_len, ps_);
        const auto swapped_indices = rng::sampleUnique(0_sz, chrom_len, num_swapped);

        Candidate child1 = parent1;
        Candidate child2 = parent2;

        for (const auto& idx : swapped_indices)
        {
            using std::swap;
            swap(child1.chromosome[idx], child2.chromosome[idx]);
        }

        return { std::move(child1), std::move(child2) };
    }

} // namespace genetic_algorithm::crossover::integer