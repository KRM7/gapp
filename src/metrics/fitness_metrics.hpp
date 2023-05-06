/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_METRICS_FITNESS_METRICS_HPP
#define GA_METRICS_FITNESS_METRICS_HPP

#include "monitor.hpp"
#include "../core/ga_info.hpp"
#include "../population/candidate.hpp"
#include "../utility/utility.hpp"
#include <span>
#include <cstddef>

namespace genetic_algorithm::metrics
{
    /* The base class used for the objective-wise fitness metrics. */
    template<typename Derived>
    class FitnessMonitor : public Monitor<Derived, FitnessMatrix>
    {
    public:
        std::span<const double> value_at(size_t generation) const noexcept
        {
            GA_ASSERT(generation < this->data_.size());

            return this->data_[generation];
        }
    private:
        void initialize(const GaInfo& ga) override
        {
            this->data_.clear();
            this->data_.reserve(ga.max_gen(), ga.num_objectives());
        }
    };

    /** Record the minimum of the fitness values in the population for each objective. */
    class FitnessMin final : public FitnessMonitor<FitnessMin>
    {
    private:
        void update(const GaInfo& ga) override;
    };

    /** Record the maximum of the fitness values in the population for each objective. */
    class FitnessMax final : public FitnessMonitor<FitnessMax>
    {
    private:
        void update(const GaInfo& ga) override;
    };

    /** Record the mean of the fitness values in the population for each objective. */
    class FitnessMean final : public FitnessMonitor<FitnessMean>
    {
    private:
        void update(const GaInfo& ga) override;
    };

    /** Record the variance of the fitness values in the population for each objective. */
    class FitnessVariance final : public FitnessMonitor<FitnessVariance>
    {
    private:
        void update(const GaInfo& ga) override;
    };

    /** Record the standard deviation of the fitness values in the population for each objective. */
    class FitnessStdDev final : public FitnessMonitor<FitnessStdDev>
    {
    private:
        void update(const GaInfo& ga) override;
    };

} // namespace genetic_algorithm::metrics

#endif // !GA_METRICS_FITNESS_METRICS_HPP
