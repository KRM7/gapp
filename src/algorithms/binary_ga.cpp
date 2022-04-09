/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "binary_ga.hpp"
#include "../utility/rng.hpp"
#include <utility>
#include <cassert>
#include <cstdlib>

namespace genetic_algorithm
{
    BinaryGA::BinaryGA(size_t chrom_len, FitnessFunction fitness_function)
        : GA(chrom_len, std::move(fitness_function))
    {
    }

    BinaryGA::Candidate BinaryGA::generateCandidate() const
    {
        assert(chrom_len_ > 0);

        Candidate sol;
        sol.chromosome.reserve(chrom_len_);
        for (size_t i = 0; i < chrom_len_; i++)
        {
            sol.chromosome.push_back(GeneType{ rng::randomBool() });
        }

        return sol;
    }

} // namespace genetic_algorithm