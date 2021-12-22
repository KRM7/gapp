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
* This file contains the real encoded genetic algorithm class.
*
* @file real_ga.h
*/

#ifndef GA_RCGA_H
#define GA_RCGA_H

#include <vector>
#include <utility>

#include "base_ga.h"

namespace genetic_algorithm
{
    /**
    * Standard genetic algorithm that uses real encoding. \n
    * Each gene of the chromosomes is a real value.
    */
    class RCGA : public GA<double>
    {
    public:
        /**
        * For the gene boundaries. \n
        * For example: { {gene1_min, gene1_max}, {gene2_min, gene2_max}, ... }
        */
        using limits_t = std::vector<std::pair<double, double>>;

        /**
        * Basic constructor for the RCGA.
        *
        * @param chrom_len The number of real genes in the chromosomes of the candidates.
        * @param fitness_function The fitness function to find the maximum of with the algorithm.
        * @param bounds The boundaries of the real coded genes (their min and max values).
        */
        RCGA(size_t chrom_len, fitnessFunction_t fitnessFunction, limits_t bounds);

        /**
        * Sets the boundaries of the real-coded genes. \n
        * Each element must contain the lower and upper bounds of the corresponding real gene. (The min and max values of the gene.) \n
        * The number of elements must be the same as the length of the chromosomes, and the lower bounds must not be higher than the upper bounds. \n
        * Eg. in the case of chromosomes of 2 values where both genes must be between -1 and 1: \n
        * limits = { {-1.0, 1.0}, {-1.0, 1.0} }
        *
        * @param limits The boundaries of the real genes.
        */
        void limits(limits_t limits);
        [[nodiscard]] limits_t limits() const;

    private:
        limits_t limits_;

        Candidate generateCandidate() const override;
    };

} // namespace genetic_algorithm


/* IMPLEMENTATION */

#include <algorithm>
#include <utility>
#include <stdexcept>
#include <cstddef>
#include <cassert>

#include "rng.hpp"

namespace genetic_algorithm
{
    inline RCGA::RCGA(size_t chrom_len, fitnessFunction_t fitnessFunction, limits_t bounds)
        : GA(chrom_len, fitnessFunction), limits_(bounds)
    {
        if (bounds.size() != chrom_len)
        {
            throw std::invalid_argument("The size of the bounds must be the same as the number of genes.");
        }
        if (std::any_of(bounds.begin(), bounds.end(), [](std::pair<double, double> b) { return b.first > b.second; }))
        {
            throw std::invalid_argument("The lower bound must be lower than the upper bound for each gene.");
        }

    }

    inline void RCGA::limits(limits_t limits)
    {
        if (limits.size() != chrom_len_)
        {
            throw std::invalid_argument("The number of limits must be equal to the chromosome length.");
        }
        if (std::any_of(limits.begin(), limits.end(), [](std::pair<double, double> b) {return b.first > b.second; }))
        {
            throw std::invalid_argument("The lower bound must be lower than the upper bound for each gene.");
        }

        limits_ = limits;
    }

    inline RCGA::limits_t RCGA::limits() const
    {
        return limits_;
    }

    inline RCGA::Candidate RCGA::generateCandidate() const
    {
        assert(chrom_len_ > 0);
        assert(chrom_len_ == limits_.size());

        Candidate sol;
        sol.chromosome.reserve(chrom_len_);
        for (size_t i = 0; i < chrom_len_; i++)
        {
            sol.chromosome.push_back(rng::randomReal(limits_[i].first, limits_[i].second));
        }

        return sol;
    }

} // namespace genetic_algorithm

#endif // !GA_RCGA_H