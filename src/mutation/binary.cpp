/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "binary.hpp"
#include "../core/ga_base.hpp"
#include "../core/candidate.hpp"
#include "../utility/rng.hpp"
#include "../utility/utility.hpp"
#include <vector>
#include <cstddef>

namespace gapp::mutation::binary
{
    void Flip::mutate(const GA<GeneType>&, const Candidate<GeneType>&, Chromosome<GeneType>& chromosome) const
    {
        const size_t flip_count = rng::randomBinomial(chromosome.size(), mutation_rate());
        const auto flipped_indices = rng::sampleUnique(0_sz, chromosome.size(), flip_count);

        for (const auto& idx : flipped_indices)
        {
            chromosome[idx] ^= 1;
        }
    }

} // namespace gapp::mutation::binary