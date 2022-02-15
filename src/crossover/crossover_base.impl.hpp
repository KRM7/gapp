/*
*  MIT License
*
*  Copyright (c) 2021 Krisztián Rugási
*
*  Permission is hereby granted, free of charge, to any person obtaining a copy
*  of this software and associated documentation files (the "Software"), to deal
*  in the Software without restriction, including without limitation the rights
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*  copies of the Software, and to permit persons to whom the Software is
*  furnished to do so, subject to the following conditions:
*
*  The above copyright notice and this permission notice shall be included in all
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

#ifndef GA_CROSSOVER_BASE_IMPL_HPP
#define GA_CROSSOVER_BASE_IMPL_HPP

#include "crossover_base.decl.hpp"
#include "../base_ga.h"
#include "../rng.hpp"

#include <algorithm>
#include <stdexcept>
#include <cassert>

namespace genetic_algorithm::crossover
{
    template<gene T>
    inline Crossover<T>::Crossover(double pc)
    {
        crossover_rate(pc);
    }

    template<gene T>
    inline void Crossover<T>::crossover_rate(double pc)
    {
        if (!(0.0 <= pc && pc <= 1.0))
        {
            throw std::invalid_argument("The crossover probability must be in the closed range [0.0, 1.0]");
        }

        pc_ = pc;
    }

    template<gene T>
    inline CandidatePair<T> Crossover<T>::operator()(const GA<T>& ga, const Candidate<T>& parent1, const Candidate<T>& parent2) const
    {
        assert(0.0 <= pc_ && pc_ <= 1.0);

        /* Only need to perform the crossover with the set pc probability. Return early with (1 - pc) probability. */
        if (rng::randomReal() >= pc_)
        {
            return { parent1, parent2 };
        }

        /*
        * If the parents are the same, the crossover doesn't need to be performed.
        * This assumes that with 2 parents that have the same chromosomes, the children would be the same
        * as the parents, which is true for every crossover operator implemented, but
        * could be an issue for user defined crossovers.
        */
        if (parent1 == parent2)
        {
            return { parent1, parent2 };
        }

        /* Perform the actual crossover. */
        auto [child1, child2] = crossover(ga, parent1, parent2);

        child1.is_evaluated = false;
        child2.is_evaluated = false;

        /*
        * Check if either of the children are the same as one of the parents.
        * (This can happen in edge cases even if the parents are different.)
        * If one of the children are the same as one of the parents, then the fitness function
        * evaluation for that child can be skipped (if the fitness function is the same) by assigning
        * it the same fitness as the parent.
        */
        if (child1 == parent1)
        {
            child1.fitness = parent1.fitness;
            child1.is_evaluated = true;
        }
        else if (child1 == parent2)
        {
            child1.fitness = parent2.fitness;
            child1.is_evaluated = true;
        }
        if (child2 == parent1)
        {
            child2.fitness = parent1.fitness;
            child2.is_evaluated = true;
        }
        else if (child2 == parent2)
        {
            child2.fitness = parent2.fitness;
            child2.is_evaluated = true;
        }

        return { child1, child2 };
    }


    template<gene T>
    inline BoundedCrossover<T>::BoundedCrossover(const std::vector<std::pair<T, T>>& bounds, double pc)
        : Crossover<T>(pc)
    {
        this->bounds(bounds);
    }

    template<gene T>
    inline void BoundedCrossover<T>::bounds(const std::vector<std::pair<T, T>>& bounds)
    {
        if (std::any_of(bounds.begin(), bounds.end(), [](std::pair<double, double> bound) {return bound.first > bound.second; }))
        {
            throw std::invalid_argument("The lower bound must be lower than the upper bound for each gene.");
        }

        bounds_ = bounds;
    }

} // namespace genetic_algorithm::crossover

#endif // !GA_CROSSOVER_BASE_IMPL_HPP