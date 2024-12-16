/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

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
    void Algorithm::initialize(const GA<T>& ga)
    {
        GAPP_ASSERT(ga.population_size() == ga.population().size());

        initializeImpl(ga);
    }

    template<typename T>
    void Algorithm::prepareSelections(const GA<T>& ga, const Population<T>& pop)
    {
        GAPP_ASSERT(ga.population_size() == pop.size());

        prepareSelectionsImpl(ga, pop);
    }

    template<typename T>
    const Candidate<T>& Algorithm::select(const GA<T>& ga, const Population<T>& pop) const
    {
        GAPP_ASSERT(ga.population_size() == pop.size());

        const size_t selected_idx = selectImpl(ga, pop);

        GAPP_ASSERT(selected_idx < pop.size(), "An invalid index was returned by selectImpl().");

        return pop[selected_idx];
    }

    template<typename T>
    Population<T> Algorithm::nextPopulation(const GA<T>& ga, Population<T> parents, Population<T> children)
    {
        GAPP_ASSERT(ga.population_size() == parents.size());
        GAPP_ASSERT(ga.population_size() <= children.size());

        parents.insert(parents.end(), std::move_iterator(children.begin()), std::move_iterator(children.end()));

        const auto next_indices = nextPopulationImpl(ga, parents);

        GAPP_ASSERT(std::all_of(next_indices.begin(), next_indices.end(), detail::less_than(parents.size())),
                    "An invalid index was returned by nextPopulationImpl().");

        GAPP_ASSERT(next_indices.size() == ga.population_size(),
                    "The number of indices returned by nextPopulationImpl() is incorrect.");

        return detail::select(std::move(parents), next_indices);
    }

    template<typename T>
    Candidates<T> Algorithm::optimalSolutions(const GA<T>& ga, const Population<T>& pop) const
    {
        GAPP_ASSERT(ga.population_size() == pop.size());

        const auto optimal_indices = optimalSolutionsImpl(ga, pop);

        GAPP_ASSERT(!optimal_indices.empty(),
                    "No optimal solutions were returned by optimalSolutionsImpl().");

        GAPP_ASSERT(std::all_of(optimal_indices.begin(), optimal_indices.end(), detail::less_than(pop.size())),
                    "An invalid index was returned by optimalSolutionsImpl().");

        return detail::select(pop, optimal_indices);
    }

    inline small_vector<size_t> Algorithm::optimalSolutionsImpl(const GaInfo& ga, const PopulationView&) const
    {
        return detail::findParetoFront(ga.fitness_matrix());
    }

} // namespace gapp::algorithm

#endif //!GAPP_ALGORITHM_ALGORITHM_BASE_IMPL_HPP
