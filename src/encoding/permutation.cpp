/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "permutation.hpp"
#include "../core/candidate.hpp"
#include "../utility/rng.hpp"
#include <algorithm>
#include <numeric>

namespace gapp
{
    Chromosome<PermutationGene> GaTraits<PermutationGene>::randomChromosome(size_t chrom_len)
    {
        Chromosome<PermutationGene> chrom(chrom_len);
        std::iota(chrom.begin(), chrom.end(), PermutationGene{ 0 });
        std::shuffle(chrom.begin(), chrom.end(), rng::prng);

        return chrom;
    }

} // namespace gapp
