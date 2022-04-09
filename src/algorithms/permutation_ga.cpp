/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "permutation_ga.hpp"
#include "../utility/rng.hpp"
#include <algorithm>
#include <numeric>
#include <cassert>

namespace genetic_algorithm
{
    PermutationGA::PermutationGA(size_t chrom_len, FitnessFunction fitnessFunction)
        : GA(chrom_len, std::move(fitnessFunction))
    {
    }

    PermutationGA::Candidate PermutationGA::generateCandidate() const
    {
        assert(chrom_len_ > 0);

        Candidate sol;

        Chromosome chrom(chrom_len_);
        std::iota(chrom.begin(), chrom.end(), GeneType{ 0 });
        std::shuffle(chrom.begin(), chrom.end(), rng::prng);

        sol.chromosome = chrom;

        return sol;
    }

} // namespace genetic_algorithm