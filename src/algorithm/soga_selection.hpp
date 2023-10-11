/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_SOGA_SELECTION_HPP
#define GA_ALGORITHM_SOGA_SELECTION_HPP

#include "selection_base.hpp"
#include "../core/population.hpp"
#include "../utility/bounded_value.hpp"
#include <vector>
#include <utility>
#include <functional>
#include <cstddef>

namespace gapp
{
    class GaInfo;

} // namespace gapp

namespace gapp::selection
{
    /**
    * %Roulette selection operator for the single-objective algorithm.
    *
    * The probability of selecting a candidate from the population is proportional
    * to it's fitness value. The operator assumes maximization, so candidates with
    * higher fitness values will have a higher probability of being selected.
    * 
    *   \f[ p_i = \frac{f_i}{\sum_i f_i} \f]
    *
    * This operator is a modified version of the standard roulette-selection method,
    * so it also works if there are negative fitness values in the population, but
    * all of the fitness must be finite for the operator to work correctly.
    */
    class Roulette final : public Selection
    {
    private:
        void prepareSelectionsImpl(const GaInfo& ga, const FitnessMatrix& fmat) override;
        size_t selectImpl(const GaInfo& ga, const FitnessMatrix& fmat) const override;

        std::vector<double> cdf_;
    };

    /**
    * %Tournament selection operator for the single-objective algorithm.
    * 
    * When performing a selection, the operator selects a set number of candidate solutions
    * randomly from the population using a uniform distribution, and then the best one is
    * selected from these. The number of candidates initially picked is controlled by the
    * @p size parameter of the operator.
    * 
    * The operator assumes fitness maximization, and works with arbitrary fitness values,
    * including infinite fitness values.
    */
    class Tournament final : public Selection
    {
    public:
        /**
        * Create a tournament selection operator.
        *
        * @param size The tournament size to use. Must be at least 1.
        */
        constexpr explicit Tournament(Positive<size_t> size = 2) noexcept :
            tourney_size_(size)
        {}

        /**
        * Set the number of candidates that will be picked for a tournament.
        * 
        * If the tourney size is 1, the selection operator is equivalent to
        * randomly selecting a candidate from a uniform distribution.
        *
        * @param size The size of the tournaments during tournament selection. Must be at least 1.
        */
        constexpr void size(Positive<size_t> size) noexcept { tourney_size_ = size; }

        /** @returns The tournament size used. */
        [[nodiscard]]
        constexpr size_t size() const noexcept { return tourney_size_; }

    private:
        size_t selectImpl(const GaInfo& ga, const FitnessMatrix& fmat) const override;

        Positive<size_t> tourney_size_;
    };

    /**
    * %Rank selection operator for the single-objective algorithm.
    * 
    * The individuals of the population are assigned selection weights between a
    * minimum and maximum value based on their rank in the population relative to
    * other individuals, assuming fitness maximization.
    * The selection probabilities are then determined based on these weights.
    * 
    * The operator works with arbitrary fitness values, infinite values are also allowed
    * to be present in the fitness matrix of the population.
    */
    class Rank final : public Selection
    {
    public:
        /**
        * Create a rank selection operator using the specified weight limits.
        *
        * @param min_weight The selection weight assigned to the worst individual of the population.
        *   Must be in the closed interval [0.0, @p max_weight].
        * @param max_weight The selection weight assigned to the best individual of the population.
        *   Must be greater than @p min_weight.
        */
        explicit Rank(NonNegative<double> min_weight = 0.1, NonNegative<double> max_weight = 1.1) noexcept;

        /**
        * Sets the minimum and maximum selection weights used.
        *
        * @param min_weight The selection weight assigned to the worst individual of the population.
        *   Must be in the closed interval [0.0, @p max_weight].
        * @param max_weight The selection weight assigned to the best individual of the population.
        *   Must be greater than @p min_weight.
        */
        void weights(NonNegative<double> min_weight, NonNegative<double> max_weight) noexcept;

        /** @returns The minimum and maximum selection weights used. */
        [[nodiscard]]
        std::pair<double, double> weights() const noexcept { return { min_weight_, max_weight_ }; }

        /** @returns The minimum selection weight used. */
        [[nodiscard]]
        double min_weight() const noexcept { return min_weight_; }

        /** @returns The maximum selection weight used. */
        [[nodiscard]]
        double max_weight() const noexcept { return max_weight_; }

    private:
        void prepareSelectionsImpl(const GaInfo& ga, const FitnessMatrix& fmat) override;
        size_t selectImpl(const GaInfo& ga, const FitnessMatrix& fmat) const override;

        std::vector<double> cdf_;
        NonNegative<double> min_weight_;
        NonNegative<double> max_weight_;
    };

    /**
    * %Sigma scaling selection operator for the single-objective algorithm.
    * 
    * The fitness values of the population are scaled based on the mean and the
    * standard deviation of the fitness values in the population, and the probability
    * of selecting a candidate will be proportional to its scaled fitness value.
    * The operator has a parameter (scale, or S) that controls how the values are
    * scaled. Smaller values of the parameter will emphasize the differences between
    * the fitness values of the candidates.
    * 
    *   \f[ f_i' = \frac{ f_i - f_{\textrm{mean}} }{ S\ f_{\textrm{sd}} } \f]
    * 
    * The operator assumes fitness maximization, and all of the fitness values of
    * the population must be finite.
    */
    class Sigma final : public Selection
    {
    public:
        /**
        * Create a sigma scaling selection operator.
        *
        * @param scale The scaling parameter to use. Must be greater than 0.
        */
        explicit Sigma(Positive<double> scale = 3.0) noexcept :
            scale_(scale)
        {}

        /**
        * Set the scaling parameter used.
        * 
        * The value of this parameter determines how the fitness values are
        * scaled. Smaller values of the parameter will emphasize the differences
        * between the candidates, meaning that even candidates with small differences
        * in their fitnesses can have large differences in their selection probabilities.
        * Larger values will lead to the candidates having a more equal probability of being
        * selected regardless of the differences in fitnesses.
        *
        * @param scale The scaling parameter to use. Must be greater than 0.
        */
        void scale(Positive<double> scale) noexcept { scale_ = scale; }

        /** @returns The scaling parameter used. */
        [[nodiscard]]
        double scale() const noexcept { return scale_; }

    private:
        void prepareSelectionsImpl(const GaInfo& ga, const FitnessMatrix& fmat) override;
        size_t selectImpl(const GaInfo& ga, const FitnessMatrix& fmat) const override;

        std::vector<double> cdf_;
        Positive<double> scale_;
    };

    /**
    * %Boltzmann selection operator forthe single-objective algorithm.
    * 
    * The fitness values of the candidates are scaled based on the overall fitness
    * values of the population, and the probability of selecting a candidate will be
    * proportional to its scaled fitness value.
    * How the fitness values are scaled changes over time in a run (from generation to
    * generation) based on a temperature function. In the early generations this temperature
    * value will be high, leading to the candidates having close to equal probabilities of
    * being selected. The temperature value will decrease over the generations, and in the
    * later generations even small differences in the fitness values of the candidates will
    * lead to large differences in their selection probabilities.
    * 
    *   \f[ p_i = \frac{ e^{\frac{f_i}{T(t)}} }{ \sum_i e^{\frac{f_i}{T(t)}} } \f]
    * 
    * The operator assumes fitness maximization, and all of the fitness values of
    * the population must be finite.
    */
    class Boltzmann final : public Selection
    {
    public:
        /**
        * The general callable type that can be used as a temperature function.
        * The function should return the temperature in the given generation, with its signature being: \n
        * \t    double f(size_t current_generation, size_t max_generation)
        */
        using TemperatureFunction = std::function<double(size_t, size_t)>;
        
        /**
        * Create a %Boltzmann selection operator.
        *
        * @param f The temperature function to use. Can't be a nullptr.
        */
        explicit Boltzmann(TemperatureFunction f = boltzmannDefaultTemp) noexcept;
        
    private:
        void prepareSelectionsImpl(const GaInfo& ga, const FitnessMatrix& fmat) override;
        size_t selectImpl(const GaInfo& ga, const FitnessMatrix& fmat) const override;

        static double boltzmannDefaultTemp(size_t gen, size_t max_gen) noexcept;

        std::vector<double> cdf_;
        TemperatureFunction temperature_;
    };

    /*
    * Wraps a callable with the right signature so that it can be used as a selection
    * method in the single-objective algorithms.
    */
    class Lambda final : public Selection
    {
    public:
        using SelectionCallable = std::function<size_t(const GaInfo&, const FitnessMatrix&)>;

        explicit Lambda(SelectionCallable f) noexcept;

    private:
        size_t selectImpl(const GaInfo& ga, const FitnessMatrix& fmat) const override;

        SelectionCallable selection_;
    };

} // namespace gapp::selection

#endif // !GA_ALGORITHM_SOGA_SELECTION_HPP