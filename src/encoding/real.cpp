/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "real.hpp"
#include "../core/candidate.hpp"
#include "../utility/rng.hpp"
#include "../utility/utility.hpp"

namespace gapp
{
    Chromosome<RealGene> GaTraits<RealGene>::randomChromosome(size_t chrom_len, const BoundsVector<RealGene>& bounds)
    {
        GAPP_ASSERT(chrom_len == bounds.size(), "The size of the bounds vector must match the chromosome length.");

        Chromosome<RealGene> chrom(chrom_len);
        for (size_t i = 0; i < chrom.size(); i++)
        {
            chrom[i] = rng::randomReal(bounds[i].lower(), bounds[i].upper());
        }

        return chrom;
    }

} // namespace gapp
