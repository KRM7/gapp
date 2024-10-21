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
    const Candidate<T>& Algorithm::select(const GA<T>& ga, const PopulationView& pop) const
    {
        GAPP_ASSERT(ga.population_size() == pop.size());

        const CandidateInfo& selected = selectImpl(ga, pop);

        GAPP_ASSERT(std::any_of(pop.begin(), pop.end(), detail::reference_to(selected)),
                    "An invalid candidate was returned by selectImpl().");

        return static_cast<const Candidate<T>&>(selected);
    }

    template<typename T>
    Population<T> Algorithm::nextPopulation(const GA<T>& ga, Population<T> parents, Population<T> children)
    {
        GAPP_ASSERT(ga.population_size() == parents.size());
        GAPP_ASSERT(ga.population_size() <= children.size());

        parents.insert(parents.end(), std::move_iterator(children.begin()), std::move_iterator(children.end()));

        const CandidatePtrVec next_pop = nextPopulationImpl(ga, parents);

        GAPP_ASSERT(std::all_of(next_pop.begin(), next_pop.end(), detail::points_into(parents)),
                    "An invalid candidate was returned by nextPopulationImpl().");

        return detail::map(next_pop, [](const CandidateInfo* candidate_info)
        {
            const auto* candidate = static_cast<const Candidate<T>*>(candidate_info);
            return std::move(*const_cast<Candidate<T>*>(candidate)); // NOLINT(*const-cast)
        });
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
