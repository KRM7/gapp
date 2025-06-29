/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "binary.hpp"
#include "../core/candidate.hpp"
#include "../utility/rng.hpp"
#include "../utility/distribution.hpp"
#include <algorithm>
#include <functional>

namespace gapp
{
    Chromosome<BinaryGene> GaTraits<BinaryGene>::randomChromosome(size_t chrom_len)
    {
        rng::uniform_bool_distribution bool_dist;

        Chromosome<BinaryGene> chrom(chrom_len);
        std::generate(chrom.begin(), chrom.end(), std::bind_front(bool_dist, rng::prng));

        return chrom;
    }

} // namespace gapp
