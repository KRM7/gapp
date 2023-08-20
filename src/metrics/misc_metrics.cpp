/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include "misc_metrics.hpp"
#include "../core/ga_info.hpp"
#include "../utility/utility.hpp"
#include <utility>
#include <cstddef>

namespace gapp::metrics
{
    void FitnessEvaluations::initialize(const GaInfo& ga)
    {
        data_.clear();
        data_.reserve(ga.max_gen());
        sum_ = 0;
    }

    void FitnessEvaluations::update(const GaInfo& ga)
    {
        const size_t old_sum = std::exchange(sum_, ga.num_fitness_evals());

        data_.push_back(sum_ - old_sum);
    }

} // namespace gapp::metrics
