/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "binary.hpp"
#include "mutation_dtl.hpp"
#include "../utility/rng.hpp"
#include "../utility/utility.hpp"
#include <cstddef>

namespace genetic_algorithm::mutation::binary
{
    void Flip::mutate(const GaInfo&, Candidate<GeneType>& candidate) const
    {
        size_t flip_count = dtl::approxMutateCnt(candidate.chromosome.size(), pm_);
        auto flipped_indices = rng::sampleUnique(0_sz, candidate.chromosome.size(), flip_count);

        for (const auto& idx : flipped_indices)
        {
            candidate.chromosome[idx] = !bool(candidate.chromosome[idx]);
        }
    }

} // namespace genetic_algorithm::mutation::binary