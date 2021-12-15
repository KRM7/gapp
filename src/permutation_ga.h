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
        * Possible mutation methods that can be used in the PermutationGA. \n
        * Includes commonly used mutation operators in permutations GAs, but a custom mutation function can
        * also be used to perform the mutations with the custom option. \n
        * Set the mutation method used in the algorithm with @ref mutation_method.
        */
        enum class MutationMethod
        {
            swap,         /**< Single-swap mutation operator. Uses no parameters. */
            scramble,     /**< Scramble mutation operator. Uses no parameters. */
            inversion,    /**< Inversion mutation operator. Uses no parameters. */
            custom        /**< Custom mutation function defined by the user. @see setMutationFunction */
        };

        /**
        * Basic contructor for the PermutationGA.
        *
        * @param chrom_len The number of genes in the chromosomes.
        * @param fitness_function The fitness function used in the algorithm to find the maximum of.
        */
        PermutationGA(size_t chrom_len, fitnessFunction_t fitnessFunction);

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

        MutationMethod mutation_method_ = MutationMethod::inversion;

        Candidate generateCandidate() const override;
        void mutate(Candidate& child) const override;

        static void swapMutate(Candidate& child, double pm);
        static void scrambleMutate(Candidate& child, double pm);
        static void inversionMutate(Candidate& child, double pm);
    };

} // namespace genetic_algorithm


/* IMPLEMENTATION */

#include <algorithm>
#include <random>
#include <vector>
#include <unordered_set>
#include <utility>
#include <tuple>
#include <stdexcept>
#include <cassert>
#include <cstdlib>

#include "rng.h"

namespace genetic_algorithm
{
    inline PermutationGA::PermutationGA(size_t chrom_len, fitnessFunction_t fitnessFunction)
        : GA(chrom_len, fitnessFunction)
    {
    }

    inline void PermutationGA::mutation_method(mutationFunction_t f)
    {
        if (f == nullptr) throw std::invalid_argument("The function used for the crossovers can't be a nullptr.");

        mutation_method_ = MutationMethod::custom;
        customMutate = f;
    }


    inline void PermutationGA::mutation_method(MutationMethod method)
    {
        if (static_cast<size_t>(method) > 3) throw std::invalid_argument("Invalid mutation method selected.");

        mutation_method_ = method;
    }

    inline PermutationGA::MutationMethod PermutationGA::mutation_method() const
    {
        return mutation_method_;
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

    inline void PermutationGA::mutate(Candidate& child) const
    {
        switch (mutation_method_)
        {
            case MutationMethod::swap:
                swapMutate(child, mutation_rate_);
                break;
            case MutationMethod::scramble:
                scrambleMutate(child, mutation_rate_);
                break;
            case MutationMethod::inversion:
                inversionMutate(child, mutation_rate_);
                break;
            case MutationMethod::custom:
                customMutate(child, mutation_rate_);
                break;
            default:
                assert(false);    /* Invalid mutation method. Shouldn't get here. */
                std::abort();
        }
    }

    inline void PermutationGA::swapMutate(Candidate& child, double pm)
    {
        assert(0.0 <= pm && pm <= 1.0);

        /* Perform mutation with pm probability. */
        if (rng::randomReal() <= pm)
        {
            /* r1 and r2 might be the same index, but its rare for long chromosomes. */
            size_t r1 = rng::randomIdx(child.chromosome.size());
            size_t r2 = rng::randomIdx(child.chromosome.size());

            std::swap(child.chromosome[r1], child.chromosome[r2]);

            /* If the indices are different, the child was changed and will need evaluation. */
            if (r1 != r2) child.is_evaluated = false;
        }
    }

    inline void PermutationGA::scrambleMutate(Candidate& child, double pm)
    {
        assert(0.0 <= pm && pm <= 1.0);

        /* Perform mutation with pm probability. */
        if (rng::randomReal() <= pm)
        {
            /* Pick a random range of genes. The bounds may be the same, but its rare for long chromosomes. */
            size_t r1 = rng::randomIdx(child.chromosome.size());
            size_t r2 = rng::randomIdx(child.chromosome.size());
            auto [idx1, idx2] = std::minmax(r1, r2);

            std::shuffle(child.chromosome.begin() + idx1, child.chromosome.begin() + idx2 + 1, rng::prng);

            /* If the indices are different, the child was very likely changed and will need evaluation. */
            if (r1 != r2) child.is_evaluated = false;
        }
    }

    inline void PermutationGA::inversionMutate(Candidate& child, double pm)
    {
        assert(0.0 <= pm && pm <= 1.0);

        /* Perform mutation with pm probability. */
        if (rng::randomReal() <= pm)
        {
            /* Pick a random range of genes. The bounds of the range may be the same, but its rare for long chromosomes. */
            size_t r1 = rng::randomIdx(child.chromosome.size());
            size_t r2 = rng::randomIdx(child.chromosome.size());
            auto [idx1, idx2] = std::minmax(r1, r2);

            std::reverse(child.chromosome.begin() + idx1, child.chromosome.begin() + idx2 + 1);

            /* If the indices are different, the child was changed and will need evaluation. */
            if (r1 != r2) child.is_evaluated = false;
        }
    }

} // namespace genetic_algorithm


#endif // !GA_PERMUTATION_GA_H