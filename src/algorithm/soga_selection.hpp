/* Copyright (c) 2022 Kriszti�n Rug�si. Subject to the MIT License. */

#ifndef GA_ALGORITHM_SOGA_SELECTION_HPP
#define GA_ALGORITHM_SOGA_SELECTION_HPP

#include "../core/ga_info.hpp"
#include "../population/population.hpp"
#include <vector>
#include <utility>
#include <functional>
#include <concepts>
#include <type_traits>
#include <cstddef>

namespace genetic_algorithm
{
    class GaInfo;
}

    /**
    * Roulette selection operator for single-objective optimization, assuming fitness maximization.
    * The probability of selecting an individual from the population is proportional to it's fitness value.
    *
    * The selection algorithm is slightly modified so that it also works with negative fitness values.
    * Has no parameters.
    */
    class Roulette
    {
    public:
        void initialize(const GaInfo&) {}
        void prepareSelections(const GaInfo& ga, const FitnessMatrix& fmat);
        size_t select(const GaInfo& ga, const FitnessMatrix& fmat);

    private:
        std::vector<double> cdf_;
    };

    /**
    * Tournament selection operator for single-objective optimization.
    * When performing a selection, a number of individuals are randomly chosen (from a uniform distribution)
    * from the population, and the best one is selected from these (assuming fitness maximization).
    *
    * The number of candidates in the tournaments is determined by the size parameter of the operator.
    */
    class Tournament
    {
    public:
        /**
        * Create a tournament selection operator for a single-objective ga.
        *
        * @param size The size of the tournaments.
        */
        explicit Tournament(size_t size = 2);

        /**
        * Sets the number of individuals that participate in a tournament.
        * Must be at least 2.
        *
        * @param size The size of the tournaments during tournament selection.
        */
        void size(size_t size);

        /** @returns The tournament size used. */
        [[nodiscard]]
        size_t size() const noexcept { return tourney_size_; }

        void initialize(const GaInfo&) {}
        void prepareSelections(const GaInfo& ga, const FitnessMatrix& fmat);
        size_t select(const GaInfo& ga, const FitnessMatrix& fmat);

    private:
        size_t tourney_size_;
        std::vector<double> fvec_;
    };

    /**
    * Rank selection operator for single-objective optimization.
    * The individuals of the population are assigned selection weights between a
    * minimum and maximum value based on their rank in the population relative to
    * other individuals (assuming fitness maximization).
    */
    class Rank
    {
    public:
        /**
        * Create a rank selection operator with the specified weight limits.
        * The @p min_weight must be on the closed interval [0.0, @p max_weight].
        * The @p max_weight must be greater than @p min_weight.
        *
        * @param min_weight The selection weight assigned to the worst individual of the population.
        * @param max_weight The selection weight assigned to the best individual of the population.
        */
        explicit Rank(double min_weight = 0.1, double max_weight = 1.1);

        /**
        * Set the minimum selection weight.
        * Must be on the closed interval [0.0, max_weight].
        *
        * @param min_weight The selection weight assigned to the worst individual of the population.
        */
        void min_weight(double min_weight);

        /** @returns The minimum selection weight. */
        [[nodiscard]]
        double min_weight() const noexcept { return min_weight_; }

        /**
        * Set the maximum selection weight.
        * Must be greater than the minimum weight.
        *
        * @param max_weight The selection weight assigned to the best individual of the population.
        */
        void max_weight(double max_weight);

        /** @returns The maximum selection weight. */
        [[nodiscard]]
        double max_weight() const noexcept { return max_weight_; }

        /**
        * Sets the minimum and maximum selection weights used.
        * The @p min_weight must be on the closed interval [0.0, @p max_weight].
        * The @p max_weight must be greater than @p min_weight.
        *
        * @param min_weight The selection weight assigned to the worst individual of the population.
        * @param max_weight The selection weight assigned to the best individual of the population.
        */
        void weights(double min_weight, double max_weight);

        /** @returns The minimum and maximum selection weights. */
        [[nodiscard]]
        std::pair<double, double> weights() const noexcept { return { min_weight_, max_weight_ }; }

        void initialize(const GaInfo&) {}
        void prepareSelections(const GaInfo& ga, const FitnessMatrix& fmat);
        size_t select(const GaInfo& ga, const FitnessMatrix& fmat);

    private:
        double min_weight_;
        double max_weight_;
        std::vector<double> cdf_;
    };

    /**
    * Sigma scaling selection operator for single-objective optimization.
    * The probability of selecting an individual from the population is proportional to
    * its scaled fitness value, which is: f' = (f - f_mean) / (scale * f_sd)
    */
    class Sigma
    {
    public:
        /**
        * Create a sigma scaling selection operator with the specified scaling.
        * The @p scale must be on the closed interval [1.0, DBL_MAX].
        *
        * @param scale The scaling parameter used.
        */
        explicit Sigma(double scale = 3.0);

        /**
        * Sets the scaling parameter used.
        * Must be on the closed interval [1.0, DBL_MAX].
        *
        * @param scale The scaling parameter used.
        */
        void scale(double scale);

        /** @returns The scale parameter. */
        [[nodiscard]]
        double scale() const noexcept { return scale_; }

        void initialize(const GaInfo&) {}
        void prepareSelections(const GaInfo& ga, const FitnessMatrix& fmat);
        size_t select(const GaInfo& ga, const FitnessMatrix& fmat);

    private:
        double scale_;
        std::vector<double> cdf_;
    };

    /**
    * Boltzmann selection operator for single-objective optimization.
    * The selection pressure throughout the algorithm is controlled by a temperature function.
    * In the early generations, the temperature is high, resulting in low selection pressure,
    * in the later generations the temperature decreases, resulting in higher selection
    * pressure.
    */
    class Boltzmann
    {
    public:
        using TemperatureFunction = std::function<double(size_t, size_t)>;  /**< The type of the temperature function. */
        
        /**
        * Create a Boltzmann selection operator with the specified temperature function.
        *
        * @param f The temperature function used by the operator.
        */
        explicit Boltzmann(TemperatureFunction f = boltzmannDefaultTemp);

        /**
        * Sets the temperature function used.
        * The temperature functions signature should be:
        *   double f(size_t current_generation, size_t max_generation)
        * and return the current temperature.
        *
        * @param f The temperature function.
        */
        void temperature_function(TemperatureFunction f); 
        
        void initialize(const GaInfo&) {}
        void prepareSelections(const GaInfo& ga, const FitnessMatrix& fmat);
        size_t select(const GaInfo& ga, const FitnessMatrix& fmat);
    private:
        TemperatureFunction temperature_;
        std::vector<double> cdf_;

        static double boltzmannDefaultTemp(size_t gen, size_t max_gen) noexcept;
    };

    static_assert(Selection<Roulette>);
    static_assert(Selection<Tournament>);
    static_assert(Selection<Rank>);
    static_assert(Selection<Sigma>);
    static_assert(Selection<Boltzmann>);

} // namespace genetic_algorithm::selection

#endif // !GA_ALGORITHM_SOGA_SELECTION_HPP