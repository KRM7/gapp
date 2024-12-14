/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "integer.hpp"
#include "../core/candidate.hpp"
#include "../utility/rng.hpp"

namespace gapp
{
    Chromosome<IntegerGene> GaTraits<IntegerGene>::randomChromosome(size_t chrom_len, const BoundsVector<IntegerGene>& bounds)
    {
        GAPP_ASSERT(chrom_len == bounds.size(), "The size of the bounds vector must match the chromosome length.");

        Chromosome<IntegerGene> chrom(chrom_len);
        for (size_t idx = 0; idx < chrom_len; idx++)
        {
            chrom[idx] = rng::randomInt(bounds[idx].lower(), bounds[idx].upper());
        }

        return chrom;
    }

} // namespace gapp
