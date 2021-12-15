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
        * Possible mutation operators that can be used in the BinaryGA. \n
        * Set the mutation method used in the algorithm with @ref mutation_method. \n
        * The function used for the mutations with the custom method can be set with @ref setMutationFunction.
        */
        enum class MutationMethod
        {
            standard,        /**< Standard mutation operator used in binary coded genetic algorithms. */
            custom           /**< Custom mutation operator defined by the user. @see setMutationFunction */
        };

        /**
        * Basic contructor for the binary GA.
        *
        * @param chrom_len The length of the binary chromosomes.
        * @param fitness_function The fitness function to find the maximum of in the algorithm.
        */
        BinaryGA(size_t chrom_len, fitnessFunction_t fitness_function);

        /**
        * Sets the mutation function used in the algorithm to @f.
        * @see MutationMethod
        *
        * @param method The mutation function to use.
        */
        void mutation_method(mutationFunction_t f);

        /**
        * Sets the mutation method used in the algorithm to @p method.
        * @see MutationMethod
        *
        * @param method The mutation method to use.
        */
        void mutation_method(MutationMethod method);
        [[nodiscard]] MutationMethod mutation_method() const;

    private:

        MutationMethod mutation_method_ = MutationMethod::standard;

        Candidate generateCandidate() const override;
        void mutate(Candidate& child) const override;

        static void standardMutate(Candidate& child, double pm);

    };

} // namespace genetic_algorithm


/* IMPLEMENTATION */

#include <algorithm>
#include <unordered_set>
#include <utility>
#include <stdexcept>
#include <cmath>
#include <cassert>
#include <cstdlib>

#include "rng.h"

namespace genetic_algorithm
{
    inline BinaryGA::BinaryGA(size_t chrom_len, fitnessFunction_t fitness_function)
        : GA(chrom_len, fitness_function) 
    {
    }

    inline void BinaryGA::mutation_method(mutationFunction_t f)
    {
        if (f == nullptr) throw std::invalid_argument("The function used for the crossovers can't be a nullptr.");

        mutation_method_ = MutationMethod::custom;
        customMutate = f;
    }

    inline void BinaryGA::mutation_method(MutationMethod method)
    {
        if (static_cast<size_t>(method) > 1) throw std::invalid_argument("Invalid mutation method selected.");

        mutation_method_ = method;
    }

    inline BinaryGA::MutationMethod BinaryGA::mutation_method() const
    {
        return mutation_method_;
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

    inline void BinaryGA::mutate(Candidate& child) const
    {
        switch (mutation_method_)
        {
            case MutationMethod::standard:
                standardMutate(child, mutation_rate_);
                break;
            case MutationMethod::custom:
                customMutate(child, mutation_rate_);
                break;
            default:
                assert(false);    /* Invalid mutation method. Shouldnt get here. */
                std::abort();
        }
    }

    inline void BinaryGA::standardMutate(Candidate& child, double pm)
    {
        assert(0.0 <= pm && pm <= 1.0);

        /* Calc number of mutated genes. */
        double mean = child.chromosome.size() * pm;
        double SD = child.chromosome.size() * pm * (1.0 - pm);

        size_t mutation_count = size_t(std::round(rng::randomNormal(mean, SD)));
        mutation_count = std::clamp(mutation_count, size_t{ 0 }, child.chromosome.size());

        /* The child will (very likely) be changed, and will need to be evaluated. */
        if (mutation_count > 0) child.is_evaluated = false;

        /* Flip mutation_count number of random genes. One gene may be flipped multiple times, but uncommon for long chromosomes. */
        for (size_t i = 0; i < mutation_count; i++)
        {
            size_t idx = rng::randomIdx(child.chromosome.size());
            child.chromosome[idx] = !child.chromosome[idx];
        }
    }

} // namespace genetic_algorithm

#endif // !GA_BINARY_GA_H