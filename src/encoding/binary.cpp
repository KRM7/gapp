/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "binary.hpp"
#include "../core/ga_base.hpp"
#include "../population/candidate.hpp"
#include "../utility/rng.hpp"
#include <algorithm>
#include <vector>

namespace genetic_algorithm
{
    auto BinaryGA::generateCandidate() const -> Candidate<GeneType>
    {
        Candidate<GeneType> solution(chrom_len());
        std::generate(solution.chromosome.begin(), solution.chromosome.end(), rng::randomBool);

        return solution;
    }

} // namespace genetic_algorithm