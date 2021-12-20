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

#ifndef GA_MUTATION_BASE_IMPL_HPP
#define GA_MUTATION_BASE_IMPL_HPP

#include "mutation_base.decl.hpp"
#include "../base_ga.h"

#include <algorithm>
#include <stdexcept>

namespace genetic_algorithm::mutation
{
    template<regular_hashable GeneType>
    inline Mutation<GeneType>::Mutation(double pm)
    {
        mutation_rate(pm);
    }

    template<regular_hashable GeneType>
    inline void Mutation<GeneType>::mutation_rate(double pm)
    {
        if (!(0.0 <= pm && pm <= 1.0))
        {
            throw std::invalid_argument("The mutation probability must be in the closed range [0.0, 1.0]");
        }

        pm_ = pm;
    }

    template<regular_hashable GeneType>
    inline void Mutation<GeneType>::operator()(const GA<GeneType>& ga, Candidate<GeneType>& candidate) const
    {
        auto old_chrom = candidate.chromosome;

        mutate(ga, candidate);

        if (candidate.chromosome != old_chrom)
        {
            candidate.is_evaluated = false;
        }
    }

    template<regular_hashable GeneType>
    inline BoundedMutation<GeneType>::BoundedMutation(const std::vector<std::pair<GeneType, GeneType>>& bounds, double pm) :
        Mutation<GeneType>(pm)
    {
        this->bounds(bounds);
    }

    template<regular_hashable GeneType>
    inline void BoundedMutation<GeneType>::bounds(const std::vector<std::pair<GeneType, GeneType>>& bounds)
    {
        if (std::any_of(bounds.begin(), bounds.end(), [](std::pair<double, double> bound) {return bound.first > bound.second; }))
        {
            throw std::invalid_argument("The lower bound must be lower than the upper bound for each gene.");
        }

        bounds_ = bounds;
    }

} // namespace genetic_algorithm::mutation

#endif // !GA_MUTATION_BASE_IMPL_HPP