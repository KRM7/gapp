/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_STOP_CONDITION_DECL_HPP
#define GA_STOP_CONDITION_DECL_HPP

#include "stop_condition_base.hpp"
#include "../core/population.hpp"
#include "../utility/bounded_value.hpp"
#include <cstddef>

namespace gapp::stopping
{
    /**
    * Early stop condition based on the number of fitness function evaluations performed.
    * 
    * The %GA will stop early if a set maximum number of objective function evaluations
    * have been performed.
    * 
    * @note The stop condition is only checked once at the end of each generation,
    *   so the actual number objective function evaluations might be higher than
    *   the limit that was set for this stop condition.
    */
    class FitnessEvals final : public StopCondition
    {
    public:
        /**
        * Create an early-stop condition based on the number of objective
        * function evaluations performed.
        *
        * @param max_fitness_evals The maximum number of fitness function evaluations
        *   that should be performed.
        */
        constexpr explicit FitnessEvals(size_t max_fitness_evals) noexcept :
            max_fitness_evals_(max_fitness_evals)
        {}

        /**
        * Set the maximum number of fitness function evaluations allowed in the %GA.
        * The actual number of evaluations might be somewhat higher, as the stop
        * condition is only checked once at the end of each generation.
        *
        * @param max_fitness_evals The maximum number of objective function evaluations
        *   that should be performed.
        */
        constexpr void max_fitness_evals(size_t max_fitness_evals) noexcept { max_fitness_evals_ = max_fitness_evals; }

        /** @returns The maximum number of fitness function evaluations allowed. */
        [[nodiscard]]
        constexpr size_t max_fitness_evals() const noexcept { return max_fitness_evals_; }

    private:
        size_t max_fitness_evals_;

        bool stop_condition(const GaInfo& ga) override;
    };

    /**
    * Early stop condition based on the fitness of the best solution discovered so far.
    * 
    * The %GA will stop early if a solution is found that is equal to or better than a
    * fitness threshold vector (assuming maximization).
    */
    class FitnessValue final : public StopCondition
    {
    public:
        /**
        * Create an early-stop condition that is based on reaching a fitness threshold.
        * Assumes fitness maximization.
        *
        * @param fitness_threshold The fitness threshold vector used. The size of this vector
        *   should be the same as the size of the fitness vectors of the population 
        *   (ie. the number of objectives should match).
        */
        explicit FitnessValue(FitnessVector fitness_threshold) noexcept :
            fitness_threshold_(std::move(fitness_threshold))
        {}

        /**
        * Set the fitness threshold vector.
        *
        * @param fitness_threshold The fitness threshold vector used. The size of this vector
        *   should be the same as the size of the fitness vectors of the population
        *   (ie. the number of objectives should match).
        */
        void fitness_threshold(FitnessVector threshold) noexcept { fitness_threshold_ = std::move(threshold); }

        /** @returns The fitness threshold vector used. */
        [[nodiscard]]
        const FitnessVector& fitness_threshold() const& noexcept { return fitness_threshold_; }

    private:
        FitnessVector fitness_threshold_;

        bool stop_condition(const GaInfo& ga) override;
    };

    /**
    * Early stop condition based on the mean fitness vector of the population.
    * 
    * The mean fitness values of the population are calculated separately for each objective,
    * and the %GA is stopped if the value hasn't improved in any of the objectives
    * for a certain number of generations.
    * 
    * For multi-objective problems, the mean fitness vector is considered improved if it's better
    * in at least one objective than the old one (assuming fitness maximization).
    */
    class FitnessMeanStall final : public StopCondition
    {
    public:
        /** Create an early-stop condition based on the mean fitness values of the population. */
        FitnessMeanStall() noexcept :
            patience_(0), delta_(1E-6), cntr_(0)
        {}

        /**
        * Create an early-stop condition based on the mean fitness values of the population.
        *
        * @param patience The number of generations to wait without stopping even if there is no improvement.
        * @param delta The minimum fitness difference considered an improvement. Negative values may also be used.
        */
        explicit FitnessMeanStall(size_t patience, double delta = 1E-6) noexcept :
            patience_(patience), delta_(delta), cntr_(0)
        {}

        /**
        * Set the patience value used for the stop condition.
        * 
        * @note The patiance counter is reset every time there is an improvement in the
        *   mean fitness vector.
        *
        * @param patience The number of generations to wait without stopping if there is no improvement.
        */
        constexpr void patience(size_t patience) noexcept { patience_ = patience; reset(); }

        /** @returns The current patience value of the stop condition. */
        [[nodiscard]]
        constexpr size_t patience() const noexcept { return patience_; }

        /**
        * Set the delta parameter of the stop condition.
        * The same delta value is used for every objective in the case of multi-objective problems.
        *
        * @param delta The minimum fitness difference considered an improvement.
        *   Negative values may also be used.
        */
        constexpr void delta(double delta) noexcept { delta_ = delta; }

        /** @returns The value of the delta parameter. */
        [[nodiscard]]
        constexpr double delta() const noexcept { return delta_; }

        void initialize(const GaInfo& ga) override;

    private:
        FitnessVector best_fitness_mean_;
        size_t patience_;
        double delta_;
        size_t cntr_;

        constexpr void reset() noexcept { cntr_ = patience_ + 1; }
        bool stop_condition(const GaInfo& ga) override;
    };

    /**
    * Early stop condition based on the best fitness vector of the population.
    * 
    * The best fitness values of the population are calculated separately for each objective,
    * and the %GA is stopped if the best fitness value hasn't improved in any of the
    * objectives for a certain number of generations.
    * 
    * For multi-objective problems, the best fitness vector is considered improved if it's better
    * in at least one coordinate than the previous best (assuming fitness maximization).
    */
    class FitnessBestStall final : public StopCondition
    {
    public:
        /** Create an early-stop condition based on the best fitness values of the population. */
        FitnessBestStall() noexcept :
            patience_(0), delta_(1E-6), cntr_(0)
        {}

        /**
        * Create an early-stop condition based on the best fitness values of the population.
        *
        * @param patience The number of generations to wait before stopping even if there is no improvement.
        * @param delta The minimum fitness difference considered an improvement. Negative values may also be used.
        */
        explicit FitnessBestStall(size_t patience, double delta = 1E-6) noexcept :
            patience_(patience), delta_(delta), cntr_(0)
        {}

        /**
        * Set the patience value used for the stop condition.
        * 
        * @note The patiance counter is reset every time there is an improvement in the
        *   mean fitness vector.
        *
        * @param patience The number of generations to wait before stopping if there is no improvement.
        */
        constexpr void patience(size_t patience) noexcept { patience_ = patience; reset(); }

       /** @returns The current patience value of the stop condition. */
        [[nodiscard]]
        constexpr size_t patience() const noexcept { return patience_; }

        /**
        * Set the delta parameter of the stop condition.
        * The same delta is used for every objective in the case of multi-objective problems.
        *
        * @param delta The minimum fitness difference considered an improvement.
        *   Negative values may also be used.
        */
        constexpr void delta(double delta) noexcept { delta_ = delta; }

        /** @returns The value of the delta parameter. */
        [[nodiscard]]
        constexpr double delta() const noexcept { return delta_; }

        void initialize(const GaInfo& ga) override;

    private:
        FitnessVector fitness_max_;
        size_t patience_;
        double delta_;
        size_t cntr_;

        constexpr void reset() noexcept { cntr_ = patience_ + 1; }
        bool stop_condition(const GaInfo& ga) override;
    };

    /**
    * This early-stop condition will always evaluate to false, so no early stopping
    * will be used by the %GA.
    * This can be used to ensure that the %GA will only stop upon reaching the
    * maximum number of generations set for it.
    */
    class NoEarlyStop final : public StopCondition
    {
    public:
        using StopCondition::StopCondition;
    private:
        bool stop_condition(const GaInfo&) noexcept override { return false; };
    };

} // namespace gapp::stopping

#endif // !GA_STOP_CONDITION_DECL_HPP