/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_SELECTION_DECL_HPP
#define GA_SELECTION_DECL_HPP

#include "selection_base.hpp"
#include "selection_dtl.hpp"
#include "../population/candidate.hpp"

#include <vector>
#include <utility>
#include <functional>

namespace genetic_algorithm::selection
{
    /**
    * Roulette selection operator for single-objective optimization, assuming fitness maximization.
    * The probability of selecting an individual from the population is proportional to it's fitness value.
    * 
    * The selection algorithm is slightly modified so that it also works with negative fitness values.
    * Has no parameters.
    */
    template<Gene T>
    class Roulette final : public Selection<T>
    {
    public:
        using Selection<T>::Selection;

        void prepare(const GA<T>& ga, const Population<T>& pop) override;
        Candidate<T> select(const GA<T>& ga, const Population<T>& pop) override;

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
    template<Gene T>
    class Tournament final : public Selection<T>
    {
    public:
        /**
        * Create a tournament selection operator for a single-objective ga.
        * 
        * @param size The size of the tournaments.
        */
        explicit Tournament(size_t size = 2);

        /**
        * Create a tournament selection operator for a single-objective ga.
        *
        * @param ga The algorithm that uses the selection operator.
        * @param size The size of the tournaments.
        */
        explicit Tournament(const GA<T>& ga, size_t size = 2);

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

        void prepare(const GA<T>& ga, const Population<T>& pop) override;
        Candidate<T> select(const GA<T>& ga, const Population<T>& pop) override;

    private:
        size_t tourney_size_;
    };

    /**
    * Rank selection operator for single-objective optimization.
    * The individuals of the population are assigned selection weights between a
    * minimum and maximum value based on their rank in the population relative to
    * other individuals (assuming fitness maximization).
    */
    template<Gene T>
    class Rank final : public Selection<T>
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
        * Create a rank selection operator with the specified weight limits.
        * The @p min_weight must be on the closed interval [0.0, @p max_weight].
        * The @p max_weight must be greater than @p min_weight.
        * 
        * @param ga The algorithm that uses the selection method.
        * @param min_weight The selection weight assigned to the worst individual of the population.
        * @param max_weight The selection weight assigned to the best individual of the population.
        */
        explicit Rank(const GA<T>& ga, double min_weight = 0.1, double max_weight = 1.1);

        /**
        * Set the minimum selection weight.
        * Must be on the closed interval [0.0, max_weight].
        * 
        * @param min_weight The selection weight assigned to the worst individual of the population.
        */
        void min_weight(double min_weight);

        /** @returns The minimum selection weight. */
        [[nodiscard]]
        double min_weigth() const noexcept { return min_weight_; }

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

        void prepare(const GA<T>& ga, const Population<T>& pop) override;
        Candidate<T> select(const GA<T>& ga, const Population<T>& pop) override;

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
    template<Gene T>
    class Sigma final : public Selection<T>
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
        * Create a sigma scaling selection operator with the specified scaling.
        * The @p scale must be on the closed interval [1.0, DBL_MAX].
        *
        * @param ga The algorithm that uses the selection operator.
        * @param scale The scaling parameter used.
        */
        explicit Sigma(const GA<T>& ga, double scale = 3.0);

        /**
        * Sets the scaling parameter used.
        * Must be on the closed interval [1.0, DBL_MAX].
        *
        * @param scale The scaling parameter used.
        */
        void scale(double scale);
        [[nodiscard]]

        /** @returns The scale parameter. */
        double scale() const noexcept { return scale_; }

        void prepare(const GA<T>& ga, const Population<T>& pop) override;
        Candidate<T> select(const GA<T>& ga, const Population<T>& pop) override;

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
    template<Gene T>
    class Boltzmann final : public Selection<T>
    {
    public:
        using TemperatureFunction = std::function<double(size_t, size_t)>;  /** The type of the temperature function. */

        /** 
        * Create a Boltzmann selection operator with the specified temperature function.
        * 
        * @param f The temperature function used by the operator.
        */
        explicit Boltzmann(TemperatureFunction f = dtl::boltzmannDefaultTemp);

        /**
        * Create a Boltzmann selection operator with the specified temperature function.
        *
        * @param ga The algorithm that uses the selection operator.
        * @param f The temperature function used by the operator.
        */
        explicit Boltzmann(const GA<T>& ga, TemperatureFunction f = dtl::boltzmannDefaultTemp);

        /**
        * Sets the temperature function used.
        * The temperature functions signature should be:
        *   double f(size_t current_generation, size_t max_generation)
        * and return the current temperature.
        *
        * @param f The temperature function.
        */
        void temperature_function(TemperatureFunction f);

        void prepare(const GA<T>& ga, const Population<T>& pop) override;
        Candidate<T> select(const GA<T>& ga, const Population<T>& pop) override;

    private:
        TemperatureFunction temperature_;
        std::vector<double> cdf_;
    };

    /**
    * Non-dominated sorting genetic algorithm (NSGA-II) for multi-objective optimization.
    */
    template<Gene T>
    class NSGA2 final : public Selection<T>
    {
    public:
        using Selection<T>::Selection;

        void init(const GA<T>& ga) override;
        void prepare(const GA<T>& ga, const Population<T>& pop) override;
        Candidate<T> select(const GA<T>& ga, const Population<T>& pop) override;
        Population<T> nextPopulation(const GA<T>& ga, Population<T>& old_pop, CandidateVec<T>& children) override;

    private:
        std::vector<size_t> ranks_;
        std::vector<double> dists_;

        /* Returns true if Pop[lidx] is better than Pop[ridx]. */
        bool crowdedCompare(size_t lidx, size_t ridx) const;
    };

    /**
    * NSGA-III algorithm for many-objective optimization.
    */
    template<Gene T>
    class NSGA3 final : public Selection<T>
    {
    public:
        using Selection<T>::Selection;

        void init(const GA<T>& ga) override;
        void prepare(const GA<T>& ga, const Population<T>& pop) override;
        Candidate<T> select(const GA<T>& ga, const Population<T>& pop) override;
        Population<T> nextPopulation(const GA<T>& ga, Population<T>& old_pop, CandidateVec<T>& children) override;

    private:
        using Point = std::vector<double>;

        struct CandidateInfo
        {
            size_t rank = 0;
            size_t ref_idx = 0;
            double ref_dist = 0.0;
            size_t niche_count = 0;
        };

        std::vector<CandidateInfo> sol_props_;
        std::vector<Point> ref_points_;
        std::vector<size_t> ref_niche_counts_;

        Point ideal_point_;
        Point nadir_point_;
        std::vector<Point> extreme_points_;

        static void updateIdealPoint(Point& ideal_point, const Population<T>& pop);
        static std::vector<Point> initExtremePoints(const Population<T>& pop, const Point& ideal_point);
        static void updateExtremePoints(std::vector<Point>& extreme_points, const Population<T>& pop, const Point& ideal_point);
        static Point nadirPoint(const std::vector<Point>& extreme_points);

        /* Find the closest reference point to each candidate after normalization, and their distances. */
        void associatePopWithRefs(const Population<T>& pop);

        /* Return the niche counts of the ref points and assign niche counts to the candidates. */
        static std::vector<size_t> calcNicheCounts(const GA<T>& ga, std::vector<CandidateInfo>& props);

        /* Returns true if Pop[lidx] is better than Pop[ridx]. */
        bool nichedCompare(size_t lidx, size_t ridx) const;
    };


    template<typename T>
    Roulette(T) -> Roulette<typename T::GeneType>;

    template<typename T>
    NSGA2(T) -> NSGA2<typename T::GeneType>;

    template<typename T>
    NSGA3(T) -> NSGA3<typename T::GeneType>;

} // namespace genetic_algorithm::selection

#endif // !GA_SELECTION_DECL_HPP