/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "binary.hpp"
#include "../core/candidate.hpp"
#include "../core/ga_info.hpp"
#include "../utility/rng.hpp"
#include "../utility/utility.hpp"
#include <vector>
#include <cstddef>

namespace gapp::mutation::binary
{
    void Flip::initialize(const GaInfo& ga)
    {
        random_binomial_.init(ga.chrom_len<GeneType>(), mutation_rate());
    }

    void Flip::mutate(const GaInfo&, const Candidate<GeneType>&, Chromosome<GeneType>& chromosome) const
    {
        const size_t flip_count = random_binomial_(chromosome.size(), mutation_rate());
        const auto flipped_indices = rng::sampleUnique(0_sz, chromosome.size(), flip_count);

        for (const auto& idx : flipped_indices)
        {
            chromosome[idx] ^= 1;
        }
    }

} // namespace gapp::mutation::binary
