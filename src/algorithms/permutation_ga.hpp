/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_PERMUTATION_GA_HPP
#define GA_PERMUTATION_GA_HPP

#include "ga_base.hpp"
#include <cstddef>

namespace genetic_algorithm
{
    /**
    * Genetic algorithm that uses permutational encoding. \n
    * The genes of the chromosomes are all unique unsigned integers on [0, chrom_len-1].
    */
    class PermutationGA : public GA<size_t>
    {
    public:
        /**
        * Basic contructor for the PermutationGA.
        *
        * @param chrom_len The number of genes in the chromosomes.
        * @param fitness_function The fitness function used in the algorithm to find the maximum of.
        */
        PermutationGA(size_t chrom_len, FitnessFunction fitnessFunction);

    private:
        Candidate generateCandidate() const override;
    };

} // namespace genetic_algorithm


/* IMPLEMENTATION */

#include "../rng.hpp"
#include <algorithm>
#include <numeric>
#include <vector>
#include <stdexcept>
#include <cassert>

namespace genetic_algorithm
{
    inline PermutationGA::PermutationGA(size_t chrom_len, FitnessFunction fitnessFunction)
        : GA(chrom_len, fitnessFunction)
    {
    }

    inline PermutationGA::Candidate PermutationGA::generateCandidate() const
    {
        assert(chrom_len_ > 0);

        Candidate sol;

        std::vector<size_t> chrom(chrom_len_);
        std::iota(chrom.begin(), chrom.end(), 0U);
        std::shuffle(chrom.begin(), chrom.end(), rng::prng);

        sol.chromosome = chrom;

        return sol;
    }

} // namespace genetic_algorithm


#endif // !GA_PERMUTATION_GA_H