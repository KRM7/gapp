/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "stop_condition.hpp"
#include "../core/ga_info.hpp"
#include "../population/population.hpp"
#include "../metrics/pop_stats.hpp"
#include "../utility/math.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <vector>
#include <stdexcept>

namespace genetic_algorithm::stopping
{
    static bool metricImproved(FitnessVector& best_so_far, const FitnessVector& new_val, double delta) noexcept
    {
        GA_ASSERT(best_so_far.size() == new_val.size());

        bool improved = false;
        for (size_t i = 0; i < new_val.size(); i++)
        {
            if (new_val[i] >= (best_so_far[i] + delta))
            {
                best_so_far[i] = new_val[i];
                improved = true;
                /* No break because the entire fitness vector needs to be updated. */
            }
        }

        return improved;
    }


    bool FitnessEvals::stop_condition(const GaInfo& ga)
    {
        return (ga.num_fitness_evals() >= max_fitness_evals_);
    }


    bool FitnessValue::stop_condition(const GaInfo& ga)
    {
        GA_ASSERT(ga.num_objectives() == fitness_threshold_.size(),
                  "The size of the fitness threshold vector must match the number of objectives.");

        const auto& fitness_matrix = ga.fitness_matrix();

        return std::any_of(fitness_matrix.begin(), fitness_matrix.end(),
        [&](const auto& fitness_vector) noexcept
        {
            return !math::paretoCompareLess(fitness_vector, fitness_threshold_);
        });
    }


    void FitnessMeanStall::initialize(const GaInfo& ga)
    {
        reset();
        best_fitness_mean_ = FitnessVector(ga.num_objectives(), -math::inf<double>);
    }

    bool FitnessMeanStall::stop_condition(const GaInfo& ga)
    {
        const FitnessMatrix& fitness_matrix = ga.fitness_matrix();
        const FitnessVector current_mean = detail::fitnessMean(fitness_matrix.begin(), fitness_matrix.end());

        const bool improved = metricImproved(best_fitness_mean_, current_mean, delta_);
        improved ? reset() : (void)--cntr_;

        return cntr_ == 0;
    }


    void FitnessBestStall::initialize(const GaInfo& ga)
    {
        reset();
        fitness_max_ = FitnessVector(ga.num_objectives(), -math::inf<double>);
    }

    bool FitnessBestStall::stop_condition(const GaInfo& ga)
    {
        const FitnessMatrix& fitness_matrix = ga.fitness_matrix();
        const FitnessVector current_max = detail::maxFitness(fitness_matrix.begin(), fitness_matrix.end());

        const bool improved = metricImproved(fitness_max_, current_max, delta_);
        improved ? reset() : (void)--cntr_;

        return cntr_ == 0;
    }

} // namespace genetic_algorithm::stopping