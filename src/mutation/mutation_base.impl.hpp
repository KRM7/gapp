/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_MUTATION_BASE_IMPL_HPP
#define GA_MUTATION_BASE_IMPL_HPP

#include "mutation_base.decl.hpp"
#include "../algorithms/ga_base.hpp"

#include <algorithm>
#include <stdexcept>

namespace genetic_algorithm::mutation
{
    template<Gene T>
    inline Mutation<T>::Mutation(double pm)
    {
        mutation_rate(pm);
    }

    template<Gene T>
    inline void Mutation<T>::mutation_rate(double pm)
    {
        if (!(0.0 <= pm && pm <= 1.0))
        {
            throw std::invalid_argument("The mutation probability must be in the closed range [0.0, 1.0]");
        }

        pm_ = pm;
    }

    template<Gene T>
    inline void Mutation<T>::operator()(const GA<T>& ga, Candidate<T>& candidate) const
    {
        auto old_chrom = candidate.chromosome;

        mutate(ga, candidate);

        if (candidate.chromosome != old_chrom)
        {
            candidate.is_evaluated = false;
        }
    }

    template<Gene T>
    inline BoundedMutation<T>::BoundedMutation(const std::vector<std::pair<T, T>>& bounds, double pm) :
        Mutation<T>(pm)
    {
        this->bounds(bounds);
    }

    template<Gene T>
    inline void BoundedMutation<T>::bounds(const std::vector<std::pair<T, T>>& bounds)
    {
        if (std::any_of(bounds.begin(), bounds.end(), [](std::pair<double, double> bound) {return bound.first > bound.second; }))
        {
            throw std::invalid_argument("The lower bound must be lower than the upper bound for each gene.");
        }

        bounds_ = bounds;
    }

} // namespace genetic_algorithm::mutation

#endif // !GA_MUTATION_BASE_IMPL_HPP