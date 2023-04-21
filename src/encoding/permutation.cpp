/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "permutation.hpp"
#include "../core/ga_base.hpp"
#include "../population/candidate.hpp"
#include "../utility/rng.hpp"
#include <algorithm>
#include <numeric>
#include <vector>

namespace genetic_algorithm
{
    auto PermutationGA::generateCandidate() const -> Candidate<GeneType>
    {
        Candidate<GeneType> solution(chrom_len());
        std::iota(solution.chromosome.begin(), solution.chromosome.end(), GeneType{ 0 });
        std::shuffle(solution.chromosome.begin(), solution.chromosome.end(), rng::prng);

        return solution;
    }

} // namespace genetic_algorithm