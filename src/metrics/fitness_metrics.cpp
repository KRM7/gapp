/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include "fitness_metrics.hpp"
#include "pop_stats.hpp"
#include "../core/ga_info.hpp"

namespace gapp::metrics
{
    void FitnessMin::update(const GaInfo& ga)
    {
        const auto& fitness_matrix = ga.fitness_matrix();
        data_.append_row(detail::minFitness(fitness_matrix.begin(), fitness_matrix.end()));
    }

    void FitnessMax::update(const GaInfo& ga)
    {
        const auto& fitness_matrix = ga.fitness_matrix();
        data_.append_row(detail::maxFitness(fitness_matrix.begin(), fitness_matrix.end()));
    }

    void FitnessMean::update(const GaInfo& ga)
    {
        const auto& fitness_matrix = ga.fitness_matrix();
        data_.append_row(detail::fitnessMean(fitness_matrix.begin(), fitness_matrix.end()));
    }

    void FitnessVariance::update(const GaInfo& ga)
    {
        const auto& fitness_matrix = ga.fitness_matrix();
        data_.append_row(detail::fitnessVariance(fitness_matrix.begin(), fitness_matrix.end()));
    }

    void FitnessStdDev::update(const GaInfo& ga)
    {
        const auto& fitness_matrix = ga.fitness_matrix();
        data_.append_row(detail::fitnessStdDev(fitness_matrix.begin(), fitness_matrix.end()));
    }

} // namespace gapp::metrics
