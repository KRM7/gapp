/*
*  MIT License
*
*  Copyright (c) 2021 Krisztián Rugási
*
*  Permission is hereby granted, free of charge, to any person obtaining a copy
*  of this software and associated documentation files (the "Software"), to deal
*  in the Software without restriction, including without limitation the rights
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*  copies of the Software, and to permit persons to whom the Software is
*  furnished to do so, subject to the following conditions:
*
*  The above copyright notice and this permission notice shall be included in all
*  copies or substantial portions of the Software.
*
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
*  SOFTWARE.
*/

/** Implementations of early-stop conditions for the genetic algorithms. */

#ifndef GA_STOP_CONDITION_DECL_HPP
#define GA_STOP_CONDITION_DECL_HPP

#include "stop_condition_base.hpp"
#include "../candidate.hpp"
#include "../concepts.hpp"

#include <vector>
#include <cstddef>

/* Early stop conditions for the genetic algorithms. */
namespace genetic_algorithm::stopping
{
    /**
    * Early stop condition based on the number of fitness function evaluations performed.
    * The algorithm will stop if the set maximum amount of objective function evaluations
    * have been performed before reaching the maximum number of generations.
    * The stop-condition is only checked at the end of each generation, so the number of
    * actual fitness function evaluations might be higher than the limit that was set.
    */
    template<gene T>
    class FitnessEvals final : public StopCondition<T>
    {
    public:
        /**
        * Create a stop condition based on the number of fitness evaluations performed.
        *
        * @param max_fitness_evals The maximum number of fitness function evaluations to perform.
        */
        explicit FitnessEvals(size_t max_fitness_evals);

        /**
        * Create a stop condition with the specified parameters.
        *
        * @param ga The context.
        * @param max_fitness_evals The maximum number of fitness function evaluations to perform.
        */
        FitnessEvals(const GA<T>& ga, size_t max_fitness_evals);

        /** Evaluate the stop condition and return true if the genetic algorithm should stop. */
        bool operator()(const GA<T>& ga) override;

        /**
        * Set the maximum number of fitness function evaluations allowed in the algorithm.
        * The actual number of fitness function evaluations might be somewhat higher as the stop
        * condition is only checked once per generation.
        *
        * @param max_fitness_evals The maximum number of fitness function evaluations to perform.
        */
        void max_fitness_evals(size_t max_fitness_evals);

        /** @returns The number of fitness function evaluations allowed. */
        [[nodiscard]] size_t max_fitness_evals() const noexcept { return max_fitness_evals_; }

    private:
        size_t max_fitness_evals_;
    };

    /**
    * Early stop condition based on the fitness of the best solution discovered so far.
    * The algorithm will stop if the best solutions fitness vector dominates the set
    * fitness threshold vector.
    * (Assuming fitness maximization.)
    */
    template<gene T>
    class FitnessValue final : public StopCondition<T>
    {
    public:
        /**
        * Create a stop condition based on reaching a fitness threshold.
        *
        * @param fitness_threshold The fitness threshold vector used for checking the stop condition.
        */
        explicit FitnessValue(const std::vector<double>& fitness_threshold);

        /**
        * Create a stop condition with the specified parameters.
        *
        * @param ga The context.
        * @param fitness_threshold The fitness threshold vector used for checking the stop condition.
        */
        FitnessValue(const GA<T>& ga, const std::vector<double>& fitness_threshold);

        /** Evaluate the stop condition and return true if the genetic algorithm should stop. */
        bool operator()(const GA<T>& ga) override;

        /**
        * Sets the fitness threshold vector used when evaluating the stop condition.
        *
        * @param threshold The fitness threshold at which the algorithm will be stopped (assuming fitness maximization).
        */
        void fitness_threshold(const std::vector<double>& threshold);

        /** @returns The currently set fitness threshold vector. */
        const std::vector<double>& fitness_threshold() const noexcept { return fitness_threshold_; }

    private:
        std::vector<double> fitness_threshold_;
    };

    /**
    * Early stop condition based on the mean fitness vector of the population. The mean fitness
    * values of the population are calculated along each fitness dimension, and the algorithm
    * is stopped if the mean fitness value hasn't improved for a set number of generations.
    * For multi-objective problems, the mean fitness vector is considered better if it's better
    * in at least one coordinate. (Assuming fitness maximization.)
    */
    template<gene T>
    class FitnessMeanStall final : public StopCondition<T>
    {
    public:
        /**
        * Create a stop condition based on the mean fitness values of the population.
        *
        * @param patience The number of generations to wait without stopping even if there is no improvement.
        * @param delta The minimum fitness difference considered an improvement.
        */
        explicit FitnessMeanStall(size_t patience = 0, double delta = 1E-6);

        /**
        * Create a stop condition with the specified parameters.
        *
        * @param ga The context.
        * @param patience The number of generations to wait without stopping even if there is no improvement.
        * @param delta The minimum fitness difference considered an improvement.
        */
        explicit FitnessMeanStall(const GA<T>& ga, size_t patience = 0, double delta = 1E-6);

        /** Evaluate the stop condition and return true if the genetic algorithm should stop. */
        bool operator()(const GA<T>& ga) override;

        /**
        * Sets the patience value used for the stop condition.
        *
        * @param patience The number of generations to wait without stopping if there is no improvement.
        */
        void patience(size_t patience);

        /** @returns The currently set patience value for the stop condition. */
        [[nodiscard]] size_t patience() const noexcept { return patience_; }

        /**
        * Sets the delta value used for the stop condition.
        * The same delta is used for each fitness coordinate in multi-objective problems.
        *
        * @param delta The minimum fitness difference considered an improvement.
        */
        void delta(double delta);

        /** @returns The currently set minimum fitness improvement value. */
        [[nodiscard]] double delta() const noexcept { return delta_; }

    private:
        size_t patience_;
        double delta_;
        size_t cntr_;
        std::vector<double> best_fitness_mean_;

        void resetCntr();
    };

    /**
    * Early stop condition based on the best fitness vector of the population. The best fitness
    * values of the population are calculated along each fitness dimension, and the algorithm
    * is stopped if the best fitness value hasn't improved for a set number of generations.
    * For multi-objective problems, the best fitness vector is considered better if it's better
    * in at least one coordinate than the previous best. (Assuming fitness maximization.)
    */
    template<gene T>
    class FitnessBestStall final : public StopCondition<T>
    {
    public:
        /**
        * Create a stop condition based on the best fitness values of the population.
        *
        * @param patience The number of generations to wait without stopping even if there is no improvement.
        * @param delta The minimum fitness difference considered an improvement.
        */
        explicit FitnessBestStall(size_t patience = 0, double delta = 1E-6);

        /**
        * Create a stop condition with the specified parameters.
        *
        * @param ga The context.
        * @param patience The number of generations to wait without stopping even if there is no improvement.
        * @param delta The minimum fitness difference considered an improvement.
        */
        explicit FitnessBestStall(const GA<T>& ga, size_t patience = 0, double delta = 1E-6);

        /** Evaluate the stop condition and return true if the genetic algorithm should stop. */
        bool operator()(const GA<T>& ga) override;

        /**
        * Sets the patience value used for the stop condition.
        *
        * @param patience The number of generations to wait without stopping if there is no improvement.
        */
        void patience(size_t patience);

        /** @returns The currently set patience value for the stop condition. */
        [[nodiscard]] size_t patience() const noexcept { return patience_; }

        /**
        * Sets the delta value used for the stop condition.
        * The same delta is used for each fitness coordinate in multi-objective problems.
        *
        * @param delta The minimum fitness difference considered an improvement.
        */
        void delta(double delta);

        /** @returns The currently set minimum fitness improvement value. */
        [[nodiscard]] double delta() const noexcept { return delta_; }

    private:
        size_t patience_;
        double delta_;
        size_t cntr_;
        std::vector<double> best_fitness_max_;

        void resetCntr();
    };

}

#endif // !GA_STOP_CONDITION_DECL_HPP