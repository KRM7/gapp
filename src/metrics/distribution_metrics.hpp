/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_METRICS_DISTRIBUTION_METRICS_HPP
#define GA_METRICS_DISTRIBUTION_METRICS_HPP

#include "monitor.hpp"
#include "../population/candidate.hpp"
#include <span>

namespace gapp::metrics
{
    /**
    * Record the nadir point of the population's fitness values in each generation.
    * This metric is intended for multi-objective problems, but it also works for
    * single-objective problems, where it will be equivalent to the FitnessMax metric.
    */
    class NadirPoint final : public Monitor<NadirPoint, FitnessMatrix>
    {
        void initialize(const GaInfo& ga) override;
        void update(const GaInfo& ga) override;
    };

    /**
    * Record the hypervolume of the population's fitness values in each generation
    * relative to some reference point. The coordinates of this reference point should be
    * less than any fitness value it will be compared to
    * (ie. the worst point of the objective space).
    * 
    * This metric is intended for multi-objective problems, but it also works for
    * single-objective ones.
    * 
    * @note This metric can be computationally expensive for large populations and dimensions.
    */
    class Hypervolume final : public Monitor<Hypervolume, std::vector<double>>
    {
    public:
        /**
        * Create a hypervolume metric.
        * 
        * @param ref_point The reference point that will be used to calculate the hypervolume.
        *   The size of this point should match the number of objectives of the fitness functions,
        *   and it should be dominated by every point in the objective space that it will be
        *   compared to.
        */
        explicit Hypervolume(FitnessVector ref_point) noexcept;

        /** @returns The reference point used for computing the hypervolumes. */
        const FitnessVector& ref_point() const noexcept { return ref_point_; }

    private:
        void initialize(const GaInfo& ga) override;
        void update(const GaInfo& ga) override;

        FitnessVector ref_point_;
    };

    /**
    * Record the hypervolume of the population's fitness values in each generation.
    * The reference point used for the calculation of the hypervolumes is determined automatically
    * as the objective-wise worst point encountered by the %GA throughout the run.
    * 
    * While the reference point will be updated throughout a run, the old hypervolume
    * values will also be updated along with it, so every generation's hypervolume is
    * computed relative to the same reference point.
    *
    * This metric is intended for multi-objective problems, but it also works for
    * single-objective ones.
    * 
    * @note This metric can be computationally expensive for large populations and dimensions.
    */
    class AutoHypervolume final : public Monitor<AutoHypervolume, std::vector<double>>
    {
    public:
        /** @returns The reference point used for computing the hypervolumes. */
        const FitnessVector& ref_point() const noexcept { return worst_point_; }

    private:
        void initialize(const GaInfo& ga) override;
        void update(const GaInfo& ga) override;

        FitnessVector worst_point_;
        FitnessMatrix ideal_points_;
    };

} // namespace gapp::metrics

#endif // !GA_METRICS_DISTRIBUTION_METRICS_HPP
