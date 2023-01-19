/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_ALGORITHM_BASE_IMPL_HPP
#define GA_ALGORITHM_ALGORITHM_BASE_IMPL_HPP

#include "algorithm_base.decl.hpp"
#include "../population/population.hpp"
#include "../utility/algorithm.hpp"
#include "../utility/functional.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <cstddef>

namespace genetic_algorithm::algorithm
{
    template<Gene T>
    auto Algorithm::select(const GaInfo& ga, const Population<T>& pop, const FitnessMatrix& fmat) const -> const Candidate<T>&
    {
        const size_t selected_idx = selectImpl(ga, fmat);

        if (selected_idx >= pop.size())
        {
            GA_THROW(std::logic_error, "An invalid candidate was selected by the algorithm.");
        }

        return pop[selected_idx];
    }

    template<Gene T>
    auto Algorithm::nextPopulation(const GaInfo& ga, Population<T>&& parents, Population<T>&& children) -> Population<T>
    {
        const size_t popsize = parents.size();

        parents.reserve(parents.size() + children.size());
        std::move(children.begin(), children.end(), std::back_inserter(parents));
        const FitnessMatrix fmat = detail::toFitnessMatrix(parents);

        const auto next_indices = nextPopulationImpl(ga, fmat.begin(), fmat.begin() + popsize, fmat.end());

        if (std::any_of(next_indices.begin(), next_indices.end(), detail::greater_eq_than(parents.size())))
        {
            GA_THROW(std::logic_error, "An invalid candidate was selected for the next population by the algorithm.");
        }

        return detail::select(std::move(parents), next_indices);
    }

    template<Gene T>
    auto Algorithm::optimalSolutions(const GaInfo& ga, const Population<T>& pop) const -> Candidates<T>
    {
        const auto optimal_indices = optimalSolutionsImpl(ga);

        if (std::any_of(optimal_indices.begin(), optimal_indices.end(), detail::greater_eq_than(pop.size())))
        {
            GA_THROW(std::logic_error, "An invalid optimal solution index was returned by optimalSolutionsImpl.");
        }

        auto optimal_sols = optimal_indices.empty() ?
            detail::findParetoFront(pop) :
            detail::select(pop, optimal_indices);

        return optimal_sols;
    }

} // namespace genetic_algorithm::algorithm

#endif //!GA_ALGORITHM_ALGORITHM_BASE_IMPL_HPP
