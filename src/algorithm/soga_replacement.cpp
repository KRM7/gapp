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
    CandidatePtrVec KeepChildren::nextPopulationImpl(const GaInfo& ga, const PopulationView& pop)
    {
        GAPP_ASSERT(pop.size() >= 2 * ga.population_size());

        CandidatePtrVec next_pop(ga.population_size());
        const auto children_first = pop.begin() + ga.population_size();
        const auto children_last  = pop.begin() + 2 * ga.population_size();

        std::transform(children_first, children_last, next_pop.begin(), [](const auto& sol) { return &sol; });

        return next_pop;
    }

    CandidatePtrVec Elitism::nextPopulationImpl(const GaInfo& ga, const PopulationView& pop)
    {
        GAPP_ASSERT(pop.size() >= 2 * ga.population_size());

        const size_t elite_count = std::min(n_, ga.population_size());

        const auto parents_first = pop.begin();
        const auto sorted_last   = pop.begin() + elite_count;
        const auto parents_last  = pop.begin() + ga.population_size();

        const auto sorted_parent_indices = detail::partial_argsort(parents_first, sorted_last, parents_last,
        [](const auto& lhs, const auto& rhs) noexcept
        {
            return math::paretoCompareLess(rhs.fitness, lhs.fitness); // descending
        });

        CandidatePtrVec next_pop(ga.population_size());

        for (size_t i = 0; i < elite_count; i++)
        {
            next_pop[i] = &pop[sorted_parent_indices[i]];
        }

        for (size_t i = elite_count; i < ga.population_size(); i++)
        {
            next_pop[i] = &pop[ga.population_size() + i - elite_count];
        }

        return next_pop;
    }

    CandidatePtrVec KeepBest::nextPopulationImpl(const GaInfo& ga, const PopulationView& pop)
    {
        GAPP_ASSERT(pop.size() >= ga.population_size());

        const auto indices = detail::partial_argsort(pop.begin(), pop.begin() + ga.population_size(), pop.end(),
        [](const auto& lhs, const auto& rhs) noexcept
        {
            return math::paretoCompareLess(rhs.fitness, lhs.fitness); // descending
        });

        CandidatePtrVec next_pop(ga.population_size());

        for (size_t i = 0; i < ga.population_size(); i++)
        {
            next_pop[i] = &pop[indices[i]];
        }

        return next_pop;
    }


    Lambda::Lambda(ReplacementCallable f) noexcept :
        replacement_(std::move(f))
    {
        GAPP_ASSERT(replacement_, "The population replacement method can't be a nullptr.");
    }

    CandidatePtrVec Lambda::nextPopulationImpl(const GaInfo& ga, const PopulationView& pop)
    {
        GAPP_ASSERT(replacement_);
        return replacement_(ga, pop);
    }

} // namespace gapp::replacement
