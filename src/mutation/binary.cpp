/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "binary.hpp"
#include "../utility/rng.hpp"
#include "../utility/utility.hpp"
#include <cstddef>

namespace genetic_algorithm::mutation::binary
{
    void Flip::mutate(const GA<GeneType>&, Candidate<GeneType>& candidate) const
    {
        const size_t flip_count = rng::randomBinomial(candidate.chromosome.size(), mutation_rate());
        const auto flipped_indices = rng::sampleUnique(0_sz, candidate.chromosome.size(), flip_count);

        for (const auto& idx : flipped_indices)
        {
            candidate.chromosome[idx] = GeneType(!bool(candidate.chromosome[idx]));
        }
    }

} // namespace genetic_algorithm::mutation::binary