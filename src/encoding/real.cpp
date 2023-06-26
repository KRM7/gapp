/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "real.hpp"
#include "../core/ga_base.hpp"
#include "../population/candidate.hpp"
#include "../utility/rng.hpp"
#include "../utility/utility.hpp"
#include <vector>

namespace gapp
{
    auto RCGA::generateCandidate() const -> Candidate<GeneType>
    {
        GAPP_ASSERT(chrom_len() == gene_bounds().size(), "The size of the bounds vector must match the chromosome length.");

        Candidate<GeneType> solution(chrom_len());
        for (size_t i = 0; i < solution.chromosome.size(); i++)
        {
            solution.chromosome[i] = rng::randomReal(gene_bounds()[i].lower(), gene_bounds()[i].upper());
        }

        return solution;
    }

} // namespace gapp