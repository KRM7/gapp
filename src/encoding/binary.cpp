/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "binary.hpp"
#include "../core/candidate.hpp"
#include "../utility/rng.hpp"
#include <algorithm>

namespace gapp
{
    Chromosome<BinaryGene> GaTraits<BinaryGene>::randomChromosome(size_t chrom_len)
    {
        Chromosome<BinaryGene> chrom(chrom_len);
        std::generate(chrom.begin(), chrom.end(), rng::randomBool);

        return chrom;
    }

} // namespace gapp
