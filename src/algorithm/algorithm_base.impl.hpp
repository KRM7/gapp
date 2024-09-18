﻿/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_ALGORITHM_ALGORITHM_BASE_IMPL_HPP
#define GAPP_ALGORITHM_ALGORITHM_BASE_IMPL_HPP

#include "algorithm_base.decl.hpp"
#include "../core/ga_info.hpp"
#include "../core/ga_base.decl.hpp"
#include "../core/population.hpp"
#include "../utility/algorithm.hpp"
#include "../utility/functional.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <iterator>
#include <cstddef>

namespace gapp::algorithm
{
    template<typename T>
    const Candidate<T>& Algorithm::select(const GA<T>& ga, const Population<T>& pop, const FitnessMatrix& fmat) const
    {
        GAPP_ASSERT(ga.population_size() == pop.size());
        GAPP_ASSERT(pop.size() == fmat.size());

        const size_t selected_idx = selectImpl(ga, fmat);

        GAPP_ASSERT(selected_idx < pop.size(), "An invalid index was returned by selectImpl().");

        return pop[selected_idx];
    }

    template<typename T>
    Population<T> Algorithm::nextPopulation(const GA<T>& ga, Population<T> parents, Population<T> children)
    {
        GAPP_ASSERT(ga.population_size() == parents.size());
        GAPP_ASSERT(ga.population_size() <= children.size());

        parents.reserve(parents.size() + children.size());
        std::move(children.begin(), children.end(), std::back_inserter(parents));
        const FitnessMatrix fmat = detail::toFitnessMatrix(parents);

        const auto next_indices = nextPopulationImpl(ga, fmat);

        GAPP_ASSERT(std::all_of(next_indices.begin(), next_indices.end(), detail::less_than(parents.size())),
                  "An invalid index was returned by nextPopulationImpl().");

        return detail::select(std::move(parents), next_indices);
    }

    template<typename T>
    Candidates<T> Algorithm::optimalSolutions(const GA<T>& ga, const Population<T>& pop) const
    {
        GAPP_ASSERT(ga.population_size() == pop.size());

        const auto optimal_indices = optimalSolutionsImpl(ga);

        GAPP_ASSERT(std::all_of(optimal_indices.begin(), optimal_indices.end(), detail::less_than(pop.size())),
                  "An invalid index was returned by optimalSolutionsImpl().");

        auto optimal_sols = optimal_indices.empty() ?
            detail::findParetoFront(pop) :
            detail::select(pop, optimal_indices);

        return optimal_sols;
    }

} // namespace gapp::algorithm

#endif //!GAPP_ALGORITHM_ALGORITHM_BASE_IMPL_HPP
