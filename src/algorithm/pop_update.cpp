/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "pop_update.hpp"
#include "../core/ga_info.hpp"
#include "../population/population.hpp"
#include "../utility/algorithm.hpp"
#include "../utility/math.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <numeric>
#include <cassert>

namespace genetic_algorithm::update
{
    std::vector<size_t> KeepChildren::operator()(const GaInfo& ga,
                                                 FitnessMatrix::const_iterator,
                                                 FitnessMatrix::const_iterator,
                                                 FitnessMatrix::const_iterator)
    {
        return detail::index_vector(ga.population_size(), ga.population_size());
    }

    std::vector<size_t> Elitism::operator()(const GaInfo& ga,
                                            FitnessMatrix::const_iterator first,
                                            FitnessMatrix::const_iterator children_first,
                                            [[maybe_unused]] FitnessMatrix::const_iterator last)
    {
        assert(ga.population_size() >= n_);

        const auto sorted_parent_indices = detail::partial_argsort(first, first + n_, children_first,
        [](const FitnessVector& lhs, const FitnessVector& rhs) noexcept
        {
            return math::paretoCompareLess(rhs, lhs); // descending
        });

        std::vector<size_t> indices(ga.population_size());
        std::copy(sorted_parent_indices.begin(), sorted_parent_indices.begin() + n_, indices.begin());
        std::iota(indices.begin() + n_, indices.end(), ga.population_size());

        return indices;
    }

    std::vector<size_t> KeepBest::operator()(const GaInfo& ga,
                                             FitnessMatrix::const_iterator first,
                                             FitnessMatrix::const_iterator children_first,
                                             FitnessMatrix::const_iterator last)
    {
        assert(size_t(children_first - first) == ga.population_size());

        auto sorted_indices = detail::partial_argsort(first, children_first, last,
        [](const FitnessVector& lhs, const FitnessVector& rhs) noexcept
        {
            return math::paretoCompareLess(rhs, lhs); // descending
        });
        sorted_indices.resize(ga.population_size());

        return sorted_indices;
    }

} // namespace genetic_algorithm::update