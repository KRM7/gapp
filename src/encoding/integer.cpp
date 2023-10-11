/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "integer.hpp"
#include "../core/candidate.hpp"
#include "../utility/rng.hpp"
#include <vector>

namespace gapp
{
    auto IntegerGA::generateCandidate() const -> Candidate<GeneType>
    {
        Candidate<GeneType> solution(chrom_len());

        for (size_t idx = 0; idx < chrom_len(); idx++)
        {
            solution.chromosome[idx] = rng::randomInt(gene_bounds()[idx].lower(), gene_bounds()[idx].upper());
        }

        return solution;
    }

} // namespace gapp