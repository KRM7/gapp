/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_METRICS_FITNESS_METRICS_HPP
#define GAPP_METRICS_FITNESS_METRICS_HPP

#include "monitor.hpp"
#include "../core/ga_info.hpp"
#include "../core/candidate.hpp"
#include "../utility/utility.hpp"
#include <span>
#include <cstddef>

namespace gapp::metrics
{
    /* The base class used for the objective-wise fitness metrics. */
    template<typename Derived>
    class FitnessMonitor : public Monitor<Derived, FitnessMatrix>
    {
    private:
        void initialize(const GaInfo& ga) override
        {
            this->data_.clear();
            this->data_.reserve(ga.max_gen(), ga.num_objectives());
        }
    };

    /** Record the objective-wise minimum of the fitness values in the population. */
    class FitnessMin final : public FitnessMonitor<FitnessMin>
    {
    private:
        void update(const GaInfo& ga) override;
    };

    /** Record the objective-wise maximum of the fitness values in the population. */
    class FitnessMax final : public FitnessMonitor<FitnessMax>
    {
    private:
        void update(const GaInfo& ga) override;
    };

    /** Record the objective-wise mean of the fitness values in the population. */
    class FitnessMean final : public FitnessMonitor<FitnessMean>
    {
    private:
        void update(const GaInfo& ga) override;
    };

    /** Record the objective-wise variance of the fitness values in the population. */
    class FitnessVariance final : public FitnessMonitor<FitnessVariance>
    {
    private:
        void update(const GaInfo& ga) override;
    };

    /** Record the objective-wise standard deviation of the fitness values in the population. */
    class FitnessStdDev final : public FitnessMonitor<FitnessStdDev>
    {
    private:
        void update(const GaInfo& ga) override;
    };

} // namespace gapp::metrics

#endif // !GAPP_METRICS_FITNESS_METRICS_HPP
