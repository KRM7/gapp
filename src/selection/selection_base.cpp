/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "selection_base.hpp"
#include "../algorithms/ga_info.hpp"
#include "../utility/algorithm.hpp"
#include <algorithm>
#include <cassert>

namespace genetic_algorithm::selection
{
    void Selection::init(const GaInfo&) 
    {
    }

    std::vector<size_t> Selection::nextPopulation(const GaInfo& ga, FitnessMatrix& fitness_matrix)
    {
        assert(std::all_of(fitness_matrix.begin(), fitness_matrix.end(), [](const FitnessVector& sol) { return !sol.empty(); }));

        auto selected_indices =
        detail::partial_argsort(fitness_matrix.begin(), fitness_matrix.begin() + ga.population_size(), fitness_matrix.end(),
        [](const FitnessVector& lhs, const FitnessVector& rhs)
        {
            return detail::paretoCompareLess(rhs, lhs);
        });
        selected_indices.resize(ga.population_size());

        return selected_indices;
    }

} // namespace genetic_algorithm::selection