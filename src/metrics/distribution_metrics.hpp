/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_METRICS_DISTRIBUTION_METRICS_HPP
#define GA_METRICS_DISTRIBUTION_METRICS_HPP

#include "monitor.hpp"
#include "../population/candidate.hpp"
#include <span>

namespace genetic_algorithm::metrics
{
    /**
    * Record the nadir point of the population's fitness values in each generation. \n
    * Intended for multi-objective problems, but works and is equivalent to FitnessMax for single-objective problems.
    */
    class NadirPoint final : public Monitor<NadirPoint, FitnessMatrix>
    {
    public:
        std::span<const double> value_at(size_t generation) const noexcept;
    private:
        void initialize(const GaInfo& ga) override;
        void update(const GaInfo& ga) override;
    };

    /**
    * Record the hypervolume of the population's fitness values in each generation
    * relative to some reference point. \n
    * The coordinates of the reference point should be less than any fitness value
    * it will be compared to (ie. the worst point of the objective space). \n
    * 
    * Intended for multi-objective problems, but also works for single-objective ones. \n
    * This metric can be computationally expensive for large populations and dimensions.
    */
    class Hypervolume final : public Monitor<Hypervolume, std::vector<double>>
    {
    public:
        /**
        * Create a hypervolume metric.
        * 
        * @param ref_point The reference point that will be used to calculate the hypervolume.
        *   This point should be dominated by every point in the objective space.
        */
        explicit Hypervolume(FitnessVector ref_point) noexcept;

        /** @returns The reference point used for computing the hypervolumes. */
        const FitnessVector ref_point() const noexcept { return ref_point_; }

        double value_at(size_t generation) const noexcept;

    private:
        void initialize(const GaInfo& ga) override;
        void update(const GaInfo& ga) override;

        FitnessVector ref_point_;
    };

    /**
    * Record the hypervolume of the population's fitness values in each generation.
    * The reference point used for the calculation of the hypervolumes is determined automatically
    * as the objective-wise worst point encountered by the algorithm during the run. \n
    * Every generation's hypervolume is computed relative to the same reference point. \n
    *
    * Intended for multi-objective problems, but also works for single-objective ones. \n
    * This metric can be computationally expensive for large populations and dimensions.
    */
    class AutoHypervolume final : public Monitor<AutoHypervolume, std::vector<double>>
    {
    public:
        /** @returns The reference point used for computing the hypervolumes. */
        const FitnessVector ref_point() const noexcept { return worst_point_; }

        double value_at(size_t generation) const noexcept;
    private:
        void initialize(const GaInfo& ga) override;
        void update(const GaInfo& ga) override;

        FitnessVector worst_point_;
        FitnessMatrix ideal_points_;
    };

} // namespace genetic_algorithm::metrics

#endif // !GA_METRICS_DISTRIBUTION_METRICS_HPP
