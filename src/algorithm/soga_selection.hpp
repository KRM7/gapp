/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_SOGA_SELECTION_HPP
#define GA_ALGORITHM_SOGA_SELECTION_HPP

#include "selection_base.hpp"
#include "../population/population.hpp"
#include "../utility/bounded_value.hpp"
#include <vector>
#include <utility>
#include <functional>
#include <cstddef>

namespace genetic_algorithm
{
    class GaInfo;

} // namespace genetic_algorithm

namespace genetic_algorithm::selection
{
    using detail::FitnessVector;
    using detail::FitnessMatrix;

    /**
    * Roulette selection operator for single-objective optimization, assuming fitness maximization.
    * The probability of selecting an individual from the population is proportional to it's fitness value. \n
    *
    * All of the fitness values must be finite for the roulette selection operator to work
    * (the algorithm is modified so that it also works with negative fitness values). \n
    * 
    * Has no parameters.
    */
    class Roulette final : public Selection
    {
    private:
        void prepareSelectionsImpl(const GaInfo& ga, const FitnessMatrix& fmat) override;
        size_t selectImpl(const GaInfo& ga, const FitnessMatrix& fmat) const override;

        std::vector<double> cdf_;
    };

    /**
    * Tournament selection operator for single-objective optimization. \n
    * When performing a selection, a number of individuals are randomly chosen (from a uniform distribution)
    * from the population, and the best one is selected from these (assuming fitness maximization). \n
    * 
    * The operator works with arbitrary fitness values, infinite values are also allowed
    * in the fitness matrix.
    *
    * The number of candidates in the tournaments is determined by the size parameter of the operator.
    */
    class Tournament final : public Selection
    {
    public:
        /**
        * Create a tournament selection operator for the single-objective GA.
        *
        * @param size The size of the tournaments. Must be greater than 0.
        */
        constexpr explicit Tournament(Positive<size_t> size = 2) noexcept :
            tourney_size_(size)
        {}

        /**
        * Set the number of individuals that participate in a tournament. \n
        * If the tourney size is 1, the selection operator is equivalent to
        * randomly selecting a candidate from a uniform distribution.
        *
        * @param size The size of the tournaments during tournament selection. Must be greater than 0.
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
    * Rank selection operator for single-objective optimization. \n
    * The individuals of the population are assigned selection weights between a
    * minimum and maximum value based on their rank in the population relative to
    * other individuals (assuming fitness maximization).
    * 
    * The operator works with arbitrary fitness values, infinite values are also allowed
    * in the fitness matrix.
    */
    class Rank final : public Selection
    {
    public:
        /**
        * Create a rank selection operator with the specified weight limits. \n
        * The @p min_weight must be on the closed interval [0.0, @p max_weight].
        * The @p max_weight must be greater than @p min_weight.
        *
        * @param min_weight The selection weight assigned to the worst individual of the population.
        * @param max_weight The selection weight assigned to the best individual of the population.
        */
        explicit Rank(NonNegative<double> min_weight = 0.1, NonNegative<double> max_weight = 1.1) noexcept;

        /**
        * Sets the minimum and maximum selection weights used. \n
        * The @p min_weight must be on the closed interval [0.0, @p max_weight].
        * The @p max_weight must be greater than @p min_weight.
        *
        * @param min_weight The selection weight assigned to the worst individual of the population.
        * @param max_weight The selection weight assigned to the best individual of the population.
        */
        void weights(NonNegative<double> min_weight, NonNegative<double> max_weight) noexcept;

        /** @returns The minimum and maximum selection weights. */
        [[nodiscard]]
        std::pair<double, double> weights() const noexcept { return { min_weight_, max_weight_ }; }

        /** @returns The minimum selection weight. */
        [[nodiscard]]
        double min_weight() const noexcept { return min_weight_; }

        /** @returns The maximum selection weight. */
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
    * Sigma scaling selection operator for single-objective optimization. \n
    * The probability of selecting an individual from the population is proportional to
    * its scaled fitness value, which is: f' = (f - f_mean) / (scale * f_sd)
    */
    class Sigma final : public Selection
    {
    public:
        /**
        * Create a sigma scaling selection operator with the specified scaling.
        *
        * @param scale The scaling parameter used. Must be greater than 0.0.
        */
        explicit Sigma(Positive<double> scale = 3.0) noexcept :
            scale_(scale)
        {}

        /**
        * Set the scaling parameter used.
        *
        * @param scale The scaling parameter of the selection.
        */
        void scale(Positive<double> scale) noexcept { scale_ = scale; }

        /** @returns The scale parameter. */
        [[nodiscard]]
        double scale() const noexcept { return scale_; }

    private:
        void prepareSelectionsImpl(const GaInfo& ga, const FitnessMatrix& fmat) override;
        size_t selectImpl(const GaInfo& ga, const FitnessMatrix& fmat) const override;

        std::vector<double> cdf_;
        Positive<double> scale_;
    };

    /**
    * Boltzmann selection operator for single-objective optimization. \n
    * The selection pressure throughout the algorithm is controlled by a temperature function.
    * In the early generations, the temperature is high, resulting in low selection pressure,
    * in the later generations the temperature decreases, resulting in higher selection
    * pressure.
    */
    class Boltzmann final : public Selection
    {
    public:
        /**
        * The type of the temperature function. \n
        * The function should return the temperature in the given generation, with its signature being: \n
        *   double f(size_t current_generation, size_t max_generation) \n
        */
        using TemperatureFunction = std::function<double(size_t, size_t)>;
        
        /**
        * Create a Boltzmann selection operator with the specified temperature function.
        *
        * @param f The temperature function used by the operator.
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

} // namespace genetic_algorithm::selection

#endif // !GA_ALGORITHM_SOGA_SELECTION_HPP