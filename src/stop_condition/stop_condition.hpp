/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_STOP_CONDITION_DECL_HPP
#define GA_STOP_CONDITION_DECL_HPP

#include "stop_condition_base.hpp"
#include "../population/population.hpp"
#include "../utility/bounded_value.hpp"
#include <cstddef>

/* Early stop conditions for the genetic algorithms. */
namespace genetic_algorithm::stopping
{
    using FitnessVector = detail::FitnessVector;
    using FitnessMatrix = detail::FitnessMatrix;

    /**
    * Early stop condition based on the number of fitness function evaluations performed. \n
    * The algorithm will stop early if the set maximum number of objective function evaluations
    * have been performed before reaching the maximum number of generations.
    * The stop-condition is only checked at the end of each generation, so the number of
    * actual fitness function evaluations might be somewhat higher than the limit that was set.
    */
    class FitnessEvals final : public StopCondition
    {
    public:
        /**
        * Create a stop condition based on the number of objective function evaluations performed.
        *
        * @param max_fitness_evals The maximum number of fitness function evaluations to perform.
        */
        constexpr explicit FitnessEvals(size_t max_fitness_evals) noexcept :
            max_fitness_evals_(max_fitness_evals)
        {}

        /**
        * Set the maximum number of fitness function evaluations allowed in the algorithm. \n
        * The actual number of fitness function evaluations might be somewhat higher as the stop
        * condition is only checked once per generation.
        *
        * @param max_fitness_evals The maximum number of fitness function evaluations to perform.
        */
        constexpr void max_fitness_evals(size_t max_fitness_evals) noexcept { max_fitness_evals_ = max_fitness_evals; }

        /** @returns The number of fitness function evaluations allowed. */
        [[nodiscard]]
        constexpr size_t max_fitness_evals() const noexcept { return max_fitness_evals_; }

    private:
        size_t max_fitness_evals_;

        bool stop_condition(const GaInfo& ga) override;
    };

    /**
    * Early stop condition based on the fitness of the best solution discovered so far. \n
    * The algorithm will stop if the best solution's fitness vector is equal to or dominates the set
    * fitness threshold vector (assuming fitness maximization).
    */
    class FitnessValue final : public StopCondition
    {
    public:
        /**
        * Create a stop condition based on reaching a set fitness threshold.
        *
        * @param fitness_threshold The fitness threshold vector used for checking the stop condition (assuming fitness maximization).
        *   The size of this vector should be the same as the size of the fitness vectors (ie. the number of objectives).
        */
        explicit FitnessValue(FitnessVector fitness_threshold) noexcept :
            fitness_threshold_(std::move(fitness_threshold))
        {}

        /**
        * Set the fitness threshold vector used when evaluating the stop condition.
        *
        * @param threshold The fitness threshold at which the algorithm will be stopped (assuming fitness maximization).
        *   The size of this vector should be the same as the size of the fitness vectors (ie. the number of objectives).
        */
        void fitness_threshold(FitnessVector threshold) noexcept { fitness_threshold_ = std::move(threshold); }

        /** @returns The fitness threshold vector of the operator. */
        [[nodiscard]]
        const FitnessVector& fitness_threshold() const& noexcept { return fitness_threshold_; }

    private:
        FitnessVector fitness_threshold_;

        bool stop_condition(const GaInfo& ga) override;
    };

    /**
    * Early stop condition based on the mean fitness vector of the population. \n The mean fitness
    * values of the population are calculated along each fitness dimension, and the algorithm
    * is stopped if the mean fitness value hasn't improved for a set number of generations. \n
    * For multi-objective problems, the mean fitness vector is considered better if it's better
    * in at least one coordinate than the old one (assuming fitness maximization).
    */
    class FitnessMeanStall final : public StopCondition
    {
    public:
        /** Create a stop condition based on the mean fitness values of the population. */
        constexpr FitnessMeanStall() noexcept :
            patience_(0), delta_(1E-6), cntr_(0)
        {}

        /**
        * Create a stop condition based on the mean fitness values of the population.
        *
        * @param patience The number of generations to wait without stopping even if there is no improvement.
        * @param delta The minimum fitness difference considered an improvement.
        */
        constexpr explicit FitnessMeanStall(size_t patience, double delta = 1E-6) noexcept :
            patience_(patience), delta_(delta), cntr_(0)
        {}

        /**
        * Set the patience value used for the stop condition.
        *
        * @param patience The number of generations to wait without stopping if there is no improvement.
        */
        constexpr void patience(size_t patience) noexcept { patience_ = patience; reset(); }

        /** @returns The currently set patience value for the stop condition. */
        [[nodiscard]]
        constexpr size_t patience() const noexcept { return patience_; }

        /**
        * Sets the delta value used for the stop condition. \n
        * The same delta value is used for every fitness coordinate in multi-objective problems. \n
        * Negative delta values may also be used.
        *
        * @param delta The minimum fitness difference considered an improvement.
        */
        constexpr void delta(double delta) noexcept { delta_ = delta; }

        /** @returns The minimum fitness improvement value. */
        [[nodiscard]]
        constexpr double delta() const noexcept { return delta_; }

    private:
        size_t patience_;
        double delta_;
        size_t cntr_;
        FitnessVector best_fitness_mean_;

        constexpr void reset() noexcept { cntr_ = patience_ + 1; }
        bool stop_condition(const GaInfo& ga) override;
    };

    /**
    * Early stop condition based on the best fitness vector of the population. \n The best fitness
    * values of the population are calculated for each fitness dimension, and the algorithm
    * is stopped if the best fitness value hasn't improved for a set number of generations. \n
    * For multi-objective problems, the best fitness vector is considered better if it's better
    * in at least one coordinate than the previous best (assuming fitness maximization).
    */
    class FitnessBestStall final : public StopCondition
    {
    public:
        /** Create a stop condition based on the best fitness values of the population. */
        constexpr FitnessBestStall() noexcept :
            patience_(0), delta_(1E-6), cntr_(0)
        {}

        /**
        * Create a stop condition based on the best fitness values of the population.
        *
        * @param patience The number of generations to wait without stopping even if there is no improvement.
        * @param delta The minimum fitness difference considered an improvement.
        */
        constexpr explicit FitnessBestStall(size_t patience, double delta = 1E-6) noexcept :
            patience_(patience), delta_(delta), cntr_(0)
        {}

        /**
        * Set the patience value used for the stop condition.
        *
        * @param patience The number of generations to wait without stopping if there is no improvement.
        */
        constexpr void patience(size_t patience) noexcept { patience_ = patience; reset(); }

        /** @returns The currently set patience value for the stop condition. */
        [[nodiscard]]
        constexpr size_t patience() const noexcept { return patience_; }

        /**
        * Set the delta value used for the stop condition. \n
        * The same delta is used for every fitness coordinate in multi-objective problems. \n
        * Negative delta values may also be used.
        *
        * @param delta The minimum fitness difference considered an improvement.
        */
        constexpr void delta(double delta) noexcept { delta_ = delta; }

        /** @returns The currently set minimum fitness improvement value. */
        [[nodiscard]]
        constexpr double delta() const noexcept { return delta_; }

    private:
        size_t patience_;
        double delta_;
        size_t cntr_;
        FitnessVector best_fitness_max_;

        constexpr void reset() noexcept { cntr_ = patience_ + 1; }
        bool stop_condition(const GaInfo& ga) override;
    };

    /**
    * This stop condition always evaluates to false, so no early stopping
    * will be used by the algorithm. \n
    * The algorithm will only stop when reaching the maximum number of generations.
    */
    class NoEarlyStop final : public StopCondition
    {
    public:
        using StopCondition::StopCondition;
    private:
        bool stop_condition(const GaInfo&) noexcept override { return false; };
    };

} // namespace genetic_algorithm::stopping

#endif // !GA_STOP_CONDITION_DECL_HPP