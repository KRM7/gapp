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
* This file contains the integer coded genetic algorithm class.
*
* @file integer_ga.h
*/

#ifndef GA_INTEGER_GA_H
#define GA_INTEGER_GA_H

#include <cstddef>

#include "base_ga.h"

namespace genetic_algorithm
{
    /**
    * Integer coded GA. \n
    * Same as @ref BinaryGA, but the genes of the chromosomes can be any integer on [0, base], not just 0 or 1. \n
    * It also uses a slightly different mutation function with swaps and inversions.
    */
    class IntegerGA : public GA<size_t>
    {
    public:
        /**
        * Basic contructor for the IntegerGA.
        *
        * @param chrom_len The number of genes in each chromosome.
        * @param fitness_function The fitness function used in the algorithm.
        * @param base The number of values a gene can take. Must be > 1. If 2, same as the @ref BinaryGA.
        */
        IntegerGA(size_t chrom_len, FitnessFunction fitnessFunction, size_t base);

        /**
        * Sets the number of values a gene can take to @p base. \n
        * The value of the base must be at least 2, and the GA is essentially the same as the
        * BinaryGA if the base is set to 2.
        *
        * @param base The number of values a gene can be.
        */
        void base(size_t base);
        [[nodiscard]] size_t base() const;

    private:
        size_t base_ = 4;

        Candidate generateCandidate() const override;
    };

} // namespace genetic_algorithm


/* IMPLEMENTATION */

#include <stdexcept>
#include <cassert>

#include "rng.hpp"

namespace genetic_algorithm
{
    inline IntegerGA::IntegerGA(size_t chrom_len, FitnessFunction fitnessFunction, size_t base)
        : GA(chrom_len, fitnessFunction), base_(base)
    {
        if (base < 2) throw std::invalid_argument("The base must be at least 2.");
    }

    inline void IntegerGA::base(size_t base)
    {
        if (base < 2) throw std::invalid_argument("The base must be at least 2.");

        base_ = base;
    }

    inline size_t IntegerGA::base() const
    {
        return base_;
    }

    inline IntegerGA::Candidate IntegerGA::generateCandidate() const
    {
        assert(chrom_len_ > 0);
        assert(base_ > 1);

        Candidate sol;
        sol.chromosome.reserve(chrom_len_);
        for (size_t i = 0; i < chrom_len_; i++)
        {
            sol.chromosome.push_back(rng::randomInt(size_t{ 0 }, base_ - 1));
        }

        return sol;
    }

} // namespace genetic_algorithm

#endif // !GA_INTEGER_GA_H