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
    small_vector<size_t> KeepChildren::nextPopulationImpl(const GaInfo& ga, const FitnessMatrix&)
    {
        return detail::index_vector(ga.population_size(), ga.population_size());
    }

    small_vector<size_t> Elitism::nextPopulationImpl(const GaInfo& ga, const FitnessMatrix& fmat)
    {
        GAPP_ASSERT(fmat.size() >= 2 * ga.population_size());

        const size_t elite_count = std::min(n_, ga.population_size());

        const auto sorted_parent_indices = detail::partial_argsort(fmat.begin(), fmat.begin() + elite_count, fmat.begin() + ga.population_size(),
        [](const auto& lhs, const auto& rhs) noexcept
        {
            return math::paretoCompareLess(rhs, lhs); // descending
        });

        small_vector<size_t> indices(ga.population_size());
        std::copy(sorted_parent_indices.begin(), sorted_parent_indices.begin() + elite_count, indices.begin());
        std::iota(indices.begin() + elite_count, indices.end(), ga.population_size());

        return indices;
    }

    small_vector<size_t> KeepBest::nextPopulationImpl(const GaInfo& ga, const FitnessMatrix& fmat)
    {
        GAPP_ASSERT(fmat.size() >= ga.population_size());

        auto sorted_indices = detail::partial_argsort(fmat.begin(), fmat.begin() + ga.population_size(), fmat.end(),
        [](const auto& lhs, const auto& rhs) noexcept
        {
            return math::paretoCompareLess(rhs, lhs); // descending
        });
        sorted_indices.resize(ga.population_size());

        return sorted_indices;
    }


    Lambda::Lambda(ReplacementCallable f) noexcept
    {
        GAPP_ASSERT(f, "The population replacement method can't be a nullptr.");

        replacement_ = std::move(f);
    }

    small_vector<size_t> Lambda::nextPopulationImpl(const GaInfo& ga, const FitnessMatrix& fmat)
    {
        GAPP_ASSERT(replacement_);

        return replacement_(ga, fmat);
    }

} // namespace gapp::replacement