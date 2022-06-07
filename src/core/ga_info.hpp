/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_GA_INFO_HPP
#define GA_GA_INFO_HPP

#include "../population/population.hpp"
#include <vector>
#include <atomic>
#include <memory>
#include <concepts>
#include <cstddef>

namespace genetic_algorithm
{
    namespace selection
    {
        class Selection;

        template<typename T>
        concept SelectionMethod = requires
        {
            requires std::derived_from<T, Selection>;
            requires std::copy_constructible<T>;
        };
    }
    namespace stopping
    {
        class StopCondition;

        template<typename T>
        concept StopMethod = requires
        {
            requires std::derived_from<T, StopCondition>;
            requires std::copy_constructible<T>;
        };
    }

    class GaInfo
    {
    public:
        GaInfo(size_t chrom_len);
        GaInfo(size_t population_size, size_t chrom_len);

        GaInfo(GaInfo&&) noexcept;
        GaInfo& operator=(GaInfo&&) noexcept;
        virtual ~GaInfo();

        GaInfo(const GaInfo&) = delete;
        GaInfo& operator=(const GaInfo&) = delete;

        using FitnessMatrix = detail::FitnessMatrix;
        using StopConditionFunction = std::function<bool(const GaInfo&)>;

        /**
        * Should be set to false if the fitness function does not change over time. \n
        * (The fitness function will always return the same fitness values for a given chromosome.) \n
        * Used to eliminate unnecesary objective function evaluations.
        */
        bool dynamic_fitness = false;

        /**
        * All pareto optimal optimal solutions found by the algorithm will be stored in the solutions,
        * not just the ones in the current population if set to true. \n
        */
        bool archive_optimal_solutions = false;

        /**
        * Sets the length of the chromosomes (number of genes) of the Candidate solutions used in the algorithm to @p len. \n
        * The chromosome length must be at least 1.
        *
        * @param len The length of the chromosomes.
        */
        void chrom_len(size_t len);

        /** @returns The chromosome length used for the candidates of the population. */
        [[nodiscard]]
        size_t chrom_len() const noexcept;

        /**
        * Sets the number of Candidates used in the population to @p size. \n
        * The population size must be at least 1.
        *
        * @param size The size of the populations.
        */
        void population_size(size_t size);

        /** @returns The population size of the algorithm. */
        [[nodiscard]]
        size_t population_size() const noexcept;

        /** @returns The maximum number of generations set for the algorithm. */
        [[nodiscard]]
        size_t max_gen() const noexcept;

        /** @returns The number of objectives. Determined by the algorithm and returns 0 before the start of the algorithm. */
        [[nodiscard]]
        size_t num_objectives() const noexcept;

        /** @returns The fitness matrix of the population. */
        [[nodiscard]]
        const FitnessMatrix& fitness_matrix() const&;

        /** @returns The number of fitness evaluations performed so far by the algorithm. */
        [[nodiscard]]
        size_t num_fitness_evals() const noexcept;

        /** @returns The current value of the generation counter. */
        [[nodiscard]]
        size_t generation_cntr() const noexcept;

        /**
        * Sets the crossover rate of the crossover operator used by the algorithm to @p pc.
        *
        * @param pc The crossover probability. Must be in the closed interval [0.0, 1.0].
        */
        virtual void crossover_rate(double pc) = 0;

        /** @returns The current crossover rate set for the crossover operator. */
        [[nodiscard]]
        virtual double crossover_rate() const noexcept = 0;

        /**
        * Sets the mutation rate of the mutation operator used by the algorithm to @p pm.
        *
        * @param pm The mutation probability. Must be in the closed interval [0.0, 1.0].
        */
        virtual void mutation_rate(double pm) = 0;

        /** @returns The current mutation rate set for the mutation operator. */
        [[nodiscard]]
        virtual double mutation_rate() const noexcept = 0;

        /**
        * Set the selection method used in the algorithm to @p f. \n
        * The selection method also determines the type of the algorithm (single- or multi-objective).
        *
        * @param method The selection method used in the single-objective algorithm.
        */
        template<selection::SelectionMethod F>
        void selection_method(const F& f);

        /**
        * Set the selection method used in the algorithm to @p f. \n
        * The selection method also determines the type of the algorithm (single- or multi-objective).
        *
        * @param method The selection method used in the single-objective algorithm.
        */
        template<selection::SelectionMethod F>
        void selection_method(std::unique_ptr<F>&& f);

        /** @returns The current selection method used by the algorithm. */
        template<selection::SelectionMethod F = selection::Selection>
        [[nodiscard]]
        F& selection_method() &;

        /**
        * Set an early-stop condition for the genetic algorithm. \n
        * The algorithm will always stop when reaching the maximum generations set, regardless of the early-stop
        * condition set here.
        *
        * @param f The early-stop method the algorithm should use.
        */
        template<stopping::StopMethod F>
        void stop_condition(const F& f);

        /**
        * Set an early-stop condition for the genetic algorithm. \n
        * The algorithm will always stop when reaching the maximum generations set, regardless of the early-stop
        * condition set here.
        *
        * @param f The early-stop method the algorithm should use.
        */
        template<stopping::StopMethod F>
        void stop_condition(std::unique_ptr<F>&& f);

        /**
        * Set the early-stop condition used by the algorithm to @p f. \n
        * The algorithm will always stop when reaching thhe maximum generations set, regardless
        * of this early-stop condition.
        *
        * @param f The function used to check for the early-stop condition.
        */
        void stop_condition(StopConditionFunction f);

        /** @returns The current stop condition used by the algorithm. */
        template<stopping::StopMethod F = stopping::StopCondition>
        [[nodiscard]]
        F& stop_condition() &;

    protected:

        FitnessMatrix fitness_matrix_;

        std::unique_ptr<selection::Selection> selection_;
        std::unique_ptr<stopping::StopCondition> stop_condition_;

        std::atomic<size_t> num_fitness_evals_ = 0;
        size_t generation_cntr_ = 0;
        size_t num_objectives_ = 0;

        size_t chrom_len_;
        size_t population_size_ = 100;
        size_t max_gen_ = 500;

        bool can_continue_ = false;

        void max_gen(size_t max_gen);
        void num_objectives(size_t n);
    };

} // namespace genetic_algorithm


/* IMPLEMENTATION */

#include "../selection/selection_base.hpp"
#include "../stop_condition/stop_condition_base.hpp"
#include <utility>
#include <stdexcept>

namespace genetic_algorithm
{
    template<selection::SelectionMethod F>
    void GaInfo::selection_method(const F& f)
    {
        selection_ = std::make_unique<F>(f);
        can_continue_ = false;
    }

    template<selection::SelectionMethod F>
    void GaInfo::selection_method(std::unique_ptr<F>&& f)
    {
        if (!f)
        {
            throw std::invalid_argument("The selection method can't be a nullptr.");
        }

        selection_ = std::move(f);
        can_continue_ = false;
    }

    template<selection::SelectionMethod F>
    F& GaInfo::selection_method() &
    {
        return dynamic_cast<F&>(*selection_);
    }

    template<stopping::StopMethod F>
    void GaInfo::stop_condition(const F& f)
    {
        stop_condition_ = std::make_unique<F>(f);
    }

    template<stopping::StopMethod F>
    void GaInfo::stop_condition(std::unique_ptr<F>&& f)
    {
        if (!f)
        {
            throw std::invalid_argument("The stop condition can't be a nullptr.");
        }

        stop_condition_ = std::move(f);
    }

    template<stopping::StopMethod F>
    F& GaInfo::stop_condition() &
    {
        return dynamic_cast<F&>(*stop_condition_);
    }

} // namespace genetic_algorithm

#endif // !GA_GA_INFO_HPP