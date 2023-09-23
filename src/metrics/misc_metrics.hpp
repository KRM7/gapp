/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_METRICS_MISC_METRICS_HPP
#define GA_METRICS_MISC_METRICS_HPP

#include "monitor.hpp"
#include <vector>
#include <cstddef>

namespace gapp::metrics
{
    /** Record the number of fitness function evaluations performed in each generation. */
    class FitnessEvaluations final : public Monitor<FitnessEvaluations, std::vector<size_t>>
    {
        void initialize(const GaInfo& ga) override;
        void update(const GaInfo& ga) override;

        size_t sum_ = 0;
    };

} // namespace gapp::metrics

#endif // !GA_METRICS_MISC_METRICS_HPP
