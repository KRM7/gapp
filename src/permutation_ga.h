/*
*  MIT License
*
*  Copyright (c) 2021 Krisztián Rugási
*
*  Permission is hereby granted, free of charge, to any person obtaining a copy
*  of this softwareand associated documentation files (the "Software"), to deal
*  in the Software without restriction, including without limitation the rights
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*  copies of the Software, and to permit persons to whom the Software is
*  furnished to do so, subject to the following conditions:
*
*  The above copyright noticeand this permission notice shall be included in all
*  copies or substantial portions of the Software.
*
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
*  SOFTWARE.
*/

/**
* This file contains the permutations GA class.
*
* @file permutation_ga.h
*/

#ifndef GA_PERMUTATION_GA_H
#define GA_PERMUTATION_GA_H

#include <cstddef>

#include "base_ga.h"

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
        PermutationGA(size_t chrom_len, fitnessFunction_t fitnessFunction);

    private:
        Candidate generateCandidate() const override;
    };

} // namespace genetic_algorithm


/* IMPLEMENTATION */

#include <algorithm>
#include <numeric>
#include <vector>
#include <stdexcept>
#include <cassert>

#include "rng.hpp"

namespace genetic_algorithm
{
    inline PermutationGA::PermutationGA(size_t chrom_len, fitnessFunction_t fitnessFunction)
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