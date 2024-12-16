/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "soga_replacement.hpp"
#include "../core/ga_info.hpp"
#include "../core/population.hpp"
#include "../utility/algorithm.hpp"
#include "../utility/math.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <numeric>
#include <utility>

namespace gapp::replacement
{
    small_vector<size_t> KeepChildren::nextPopulationImpl(const GaInfo& ga, const PopulationView&)
    {
        return detail::index_vector(ga.population_size(), ga.population_size());
    }

    small_vector<size_t> Elitism::nextPopulationImpl(const GaInfo& ga, const PopulationView& pop)
    {
        GAPP_ASSERT(pop.size() >= 2 * ga.population_size());

        const size_t elite_count = std::min(n_, ga.population_size());

        const auto parents_first = pop.begin();
        const auto sorted_last   = pop.begin() + elite_count;
        const auto parents_last  = pop.begin() + ga.population_size();

        auto indices = detail::partial_argsort(parents_first, sorted_last, parents_last,
        [](const auto& lhs, const auto& rhs) noexcept
        {
            return math::paretoCompareLess(rhs.fitness, lhs.fitness); // descending
        });

        indices.resize(ga.population_size());
        std::iota(indices.begin() + elite_count, indices.end(), ga.population_size());

        return indices;
    }

    small_vector<size_t> KeepBest::nextPopulationImpl(const GaInfo& ga, const PopulationView& pop)
    {
        GAPP_ASSERT(pop.size() >= ga.population_size());

        auto indices = detail::partial_argsort(pop.begin(), pop.begin() + ga.population_size(), pop.end(),
        [](const auto& lhs, const auto& rhs) noexcept
        {
            return math::paretoCompareLess(rhs.fitness, lhs.fitness); // descending
        });

        indices.resize(ga.population_size());

        return indices;
    }


    Lambda::Lambda(ReplacementCallable f) noexcept :
        replacement_(std::move(f))
    {
        GAPP_ASSERT(replacement_, "The population replacement method can't be a nullptr.");
    }

    small_vector<size_t> Lambda::nextPopulationImpl(const GaInfo& ga, const PopulationView& pop)
    {
        GAPP_ASSERT(replacement_);
        return replacement_(ga, pop);
    }

} // namespace gapp::replacement
