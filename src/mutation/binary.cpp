/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "binary.hpp"
#include "mutation_dtl.hpp"
#include "../utility/rng.hpp"
#include <cstddef>

namespace genetic_algorithm::mutation::binary
{
    void Flip::mutate(const GaInfo&, Candidate<GeneType>& candidate) const
    {
        size_t flip_cnt = dtl::approxMutateCnt(candidate.chromosome.size(), pm_);
        auto flipped_idxs = rng::sampleUnique(candidate.chromosome.size(), flip_cnt);

        for (const auto& idx : flipped_idxs)
        {
            candidate.chromosome[idx] = !bool(candidate.chromosome[idx]);
        }
    }

} // namespace genetic_algorithm::mutation::binary