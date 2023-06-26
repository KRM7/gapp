/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include "distribution_metrics.hpp"
#include "pop_stats.hpp"
#include "../core/ga_info.hpp"
#include "../population/population.hpp"
#include "../utility/math.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <functional>

namespace gapp::metrics
{
    using math::Point;

    std::span<const double> NadirPoint::value_at(size_t generation) const noexcept
    {
        GAPP_ASSERT(generation < data_.size());

        return data_[generation];
    }

    void NadirPoint::initialize(const GaInfo& ga)
    {
        data_.clear();
        data_.reserve(ga.max_gen(), ga.num_objectives());
    }

    void NadirPoint::update(const GaInfo& ga)
    {
        GAPP_ASSERT(ga.population_size() > 0);

        const auto& fitness_matrix = ga.fitness_matrix();
        data_.append_row(detail::findNadirPoint(fitness_matrix));
    }


    Hypervolume::Hypervolume(FitnessVector ref_point) noexcept :
        ref_point_(std::move(ref_point))
    {}

    double Hypervolume::value_at(size_t generation) const noexcept
    {
        GAPP_ASSERT(generation < data_.size());

        return data_[generation];
    }

    void Hypervolume::initialize(const GaInfo& ga)
    {
        GAPP_ASSERT(ref_point_.size() == ga.num_objectives());

        data_.clear();
        data_.reserve(ga.max_gen());
    }

    void Hypervolume::update(const GaInfo& ga)
    {
        GAPP_ASSERT(ref_point_.size() == ga.num_objectives());

        const auto& fitness_matrix = ga.fitness_matrix();
        data_.push_back(detail::hypervolume(fitness_matrix, ref_point_));
    }


    double AutoHypervolume::value_at(size_t generation) const noexcept
    {
        GAPP_ASSERT(generation < data_.size());

        return data_[generation];
    }

    void AutoHypervolume::initialize(const GaInfo& ga)
    {
        data_.clear();
        data_.reserve(ga.max_gen());

        ideal_points_.clear();
        ideal_points_.reserve(ga.max_gen(), ga.num_objectives());

        worst_point_ = FitnessVector(ga.num_objectives(), math::inf<double>);
    }

    void AutoHypervolume::update(const GaInfo& ga)
    {
        const auto& fitness_matrix = ga.fitness_matrix();

        FitnessVector worst_point = detail::minFitness(fitness_matrix.begin(), fitness_matrix.end());
        FitnessVector ideal_point = detail::maxFitness(fitness_matrix.begin(), fitness_matrix.end());

        FitnessVector prev_worst_point = worst_point_;
        worst_point_ = detail::elementwise_min(std::move(worst_point_), worst_point);

        if (worst_point_ != prev_worst_point)
        {
            for (size_t i = 0; i < data_.size(); i++)
            {
                const double old_vol = math::volumeBetween(ideal_points_[i], prev_worst_point);
                const double new_vol = math::volumeBetween(ideal_points_[i], worst_point_);

                data_[i] += (new_vol - old_vol);
            }
        }

        data_.push_back(detail::hypervolume(fitness_matrix, worst_point_));
        ideal_points_.append_row(std::move(ideal_point));
    }

} // namespace gapp::metrics
