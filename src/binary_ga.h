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
* This file contains the binary genetic algorithm class.
*
* @file binary_ga.h
*/

#ifndef GA_BINARY_GA_H
#define GA_BINARY_GA_H

#include <cstddef>
#include "base_ga.h"
#include "candidate.h"

namespace genetic_algorithm
{
    /**
    * Standard genetic algorithm with binary encoding. \n
    * (The binary genes are encoded as char types.)
    */
    class BinaryGA : public GA<char>
    {
    public:
        /**
        * Basic contructor for the binary GA.
        *
        * @param chrom_len The length of the binary chromosomes.
        * @param fitness_function The fitness function to find the maximum of in the algorithm.
        */
        BinaryGA(size_t chrom_len, fitnessFunction_t fitness_function);

    private:
        Candidate generateCandidate() const override;
    };

} // namespace genetic_algorithm


/* IMPLEMENTATION */

#include <cassert>
#include <cstdlib>

#include "rng.hpp"

namespace genetic_algorithm
{
    inline BinaryGA::BinaryGA(size_t chrom_len, fitnessFunction_t fitness_function)
        : GA(chrom_len, fitness_function) 
    {
    }

    inline BinaryGA::Candidate BinaryGA::generateCandidate() const
    {
        assert(chrom_len_ > 0);

        Candidate sol;
        sol.chromosome.reserve(chrom_len_);
        for (size_t i = 0; i < chrom_len_; i++)
        {
            sol.chromosome.push_back(char(rng::randomBool()));
        }

        return sol;
    }

} // namespace genetic_algorithm

#endif // !GA_BINARY_GA_H