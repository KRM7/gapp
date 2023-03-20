/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "soga_replacement.hpp"
#include "../core/ga_info.hpp"
#include "../population/population.hpp"
#include "../utility/algorithm.hpp"
#include "../utility/math.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <numeric>
#include <utility>

namespace genetic_algorithm::replacement
{
    std::vector<size_t> KeepChildren::nextPopulationImpl(const GaInfo& ga, FitnessMatrix::const_iterator, FitnessMatrix::const_iterator, FitnessMatrix::const_iterator)
    {
        return detail::index_vector(ga.population_size(), ga.population_size());
    }

    std::vector<size_t> Elitism::nextPopulationImpl(const GaInfo& ga, FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator children_first, FitnessMatrix::const_iterator)
    {
        const size_t n = std::min(n_, ga.population_size());

        const auto sorted_parent_indices = detail::partial_argsort(first, first + n, children_first,
        [](const auto& lhs, const auto& rhs) noexcept
        {
            return math::paretoCompareLess(rhs, lhs); // descending
        });

        std::vector<size_t> indices(ga.population_size());
        std::copy(sorted_parent_indices.begin(), sorted_parent_indices.begin() + n, indices.begin());
        std::iota(indices.begin() + n, indices.end(), ga.population_size());

        return indices;
    }

    std::vector<size_t> KeepBest::nextPopulationImpl(const GaInfo& ga, FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator children_first, FitnessMatrix::const_iterator last)
    {
        GA_ASSERT(size_t(children_first - first) == ga.population_size());

        auto sorted_indices = detail::partial_argsort(first, children_first, last,
        [](const auto& lhs, const auto& rhs) noexcept
        {
            return math::paretoCompareLess(rhs, lhs); // descending
        });
        sorted_indices.resize(ga.population_size());

        return sorted_indices;
    }


    Lambda::Lambda(ReplacementCallable f) noexcept
    {
        GA_ASSERT(f, "The population replacement method can't be a nullptr.");

        replacement_ = std::move(f);
    }

    std::vector<size_t> Lambda::nextPopulationImpl(const GaInfo& ga, FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator children_first, FitnessMatrix::const_iterator last)
    {
        GA_ASSERT(replacement_);

        return replacement_(ga, first, children_first, last);
    }

} // namespace genetic_algorithm::replacement