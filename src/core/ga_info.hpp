/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CORE_GA_INFO_HPP
#define GA_CORE_GA_INFO_HPP

#include "../population/population.hpp"
#include "../selection/selection_base.fwd.hpp"
#include "../stop_condition/stop_condition_base.fwd.hpp"
#include <vector>
#include <atomic>
#include <memory>
#include <cstddef>

namespace genetic_algorithm
{
    /**
    * Base class that all GAs are derived from. \n
    * Contains all of the general information about a GA that does not depend
    * on the encoding (gene) type. \n
    * Move-only.
    */
    class GaInfo
    {
    public:
        /**
        * Construct a genetic algorithm.
        * 
        * @param chrom_len The size of the Chromosomes in the Population (each solution has a single Chromosome).
        */
        GaInfo(size_t chrom_len);

        /**
        * Construct a genetic algorithm.
        * 
        * @param population_size The number of solutions in the Population.
        * @param chrom_len The size of the Chromosomes in the Population (each solution has a single Chromosome).
        */
        GaInfo(size_t population_size, size_t chrom_len);

        /**
        * Type returned by fitness_matrix(), used to represent the fitness matrix of the population. \n
        * Each element of the matrix is the fitness vector of the corresponding solution of the population. \n
        * Eg. fmat[0] is the fitness vector of the first member of the population.
        */
        using FitnessMatrix = detail::FitnessMatrix;

        /**
        * The type of the elements of a FitnessMatrix. \n
        * Contains fitness values along each objective axis.
        */
        using FitnessVector = detail::FitnessVector;

        /**
        * The type of the stop condition function if a function or lambda is used instead of a functor
        * derived from StopCondition. \n
        * The function should return true if the algorithm should stop.
        */
        using StopConditionFunction = std::function<bool(const GaInfo&)>;

        /**
        * Should be set to false if the fitness function does not change while running the algorithm. \n
        * (The fitness function will always return the same fitness values for a given Chromosome.) \n
        * Used to eliminate unnecessary objective function evaluations.
        */
        bool dynamic_fitness = false;

        /**
        * All pareto optimal optimal Candidates found by the algorithm will be stored and kept in the solutions,
        * regardless of when they were founds,
        * not just the optimal solutions of the last Population if this is set to true. \n
        */
        bool keep_all_optimal_solutions = false;

        /**
        * Set the size of the chromosomes (number of genes) of the Candidate solutions used in the algorithm. \n
        * The chromosome length must be at least 1.
        *
        * @param len The length of the Chromosomes.
        */
        void chrom_len(size_t len);

        /** @returns The Chromosome length used for the Candidates of the Population. */
        [[nodiscard]]
        size_t chrom_len() const noexcept;

        /**
        * Set the number of Candidates used in the Population. \n
        * Must be at least 1.
        *
        * @param size The number of Candidates in a Population.
        */
        void population_size(size_t size);

        /** @returns The number of Candidates in the Population. */
        [[nodiscard]]
        size_t population_size() const noexcept;

        /** @returns The maximum number of generations set for the algorithm. */
        [[nodiscard]]
        size_t max_gen() const noexcept;

        /** @returns The number of objectives as determined by the algorithm based on the fitness function. */
        [[nodiscard]]
        size_t num_objectives() const noexcept;

        /**
        * @returns The fitness matrix of the population.
        * Each element of the matrix is the fitness vector of the corresponding solution in the population. \n
        * Eg. fmat[0] is the fitness vector of the first member of the population.
        */
        [[nodiscard]]
        const FitnessMatrix& fitness_matrix() const& noexcept;

        /** @returns The number of fitness evaluations performed by the algorithm. This value is updated after every objective function evaluation. */
        [[nodiscard]]
        size_t num_fitness_evals() const noexcept;

        /** @returns The current generation's number. */
        [[nodiscard]]
        size_t generation_cntr() const noexcept;
        
        /**
        * Set the crossover rate of the crossover operator used by the algorithm.
        *
        * @param pc The crossover probability. Must be in the closed interval [0.0, 1.0].
        */
        virtual void crossover_rate(double pc) = 0;

        /** @returns The crossover rate set for the crossover operator. */
        [[nodiscard]]
        virtual double crossover_rate() const noexcept = 0;

        /**
        * Set the mutation rate of the mutation operator used by the algorithm.
        *
        * @param pm The mutation probability. Must be in the closed interval [0.0, 1.0].
        */
        virtual void mutation_rate(double pm) = 0;

        /** @returns The mutation rate set for the mutation operator. */
        [[nodiscard]]
        virtual double mutation_rate() const noexcept = 0;

        /**
        * Set the selection method used in the algorithm. \n
        * The selection method will also determine the type of the algorithm (single- or multi-objective),
        * so it must be consistent with the size of the fitness vectors returned by the fitness function.
        * @see Selection
        *
        * @param method The selection method used in the algorithm.
        */
        template<typename F>
        requires selection::SelectionMethod<F> && std::is_final_v<F>
        void selection_method(F&& f);

        /**
        * Set the selection method used in the algorithm. \n
        * The selection method will also determine the type of the algorithm (single- or multi-objective),
        * so it must be consistent with the size of the fitness vectors returned by the fitness function.
        * @see Selection
        *
        * @param method The selection method used in the algorithm.
        */
        template<selection::SelectionMethod F>
        void selection_method(std::unique_ptr<F>&& f);

        /** @returns The selection method used by the algorithm, cast to type @p F. */
        template<selection::SelectionMethod F = selection::Selection>
        [[nodiscard]]
        F& selection_method() &;

        /**
        * Set an early-stop condition for the genetic algorithm. \n
        * The algorithm will always stop when reaching the maximum generations set, regardless of the stop
        * condition set here.
        * @see StopCondition
        *
        * @param f The StopCondition the algorithm should use.
        */
        template<typename F>
        requires stopping::StopMethod<F> && std::is_final_v<F>
        void stop_condition(F&& f);

        /**
        * Set an early-stop condition for the genetic algorithm. \n
        * The algorithm will always stop when reaching the maximum generations set, regardless of the stop
        * condition set here.
        * @see StopCondition
        *
        * @param f The StopCondition the algorithm should use.
        */
        template<stopping::StopMethod F>
        void stop_condition(std::unique_ptr<F>&& f);

        /**
        * Set an early-stop condition for the genetic algorithm. \n
        * The algorithm will always stop when reaching the maximum generations set, regardless of the stop
        * condition set here.
        * @see StopCondition
        * @see StopConditionFunction
        *
        * @param f The function used to check for the early-stop condition.
        */
        void stop_condition(StopConditionFunction f);

        /** @returns The stop condition used by the algorithm, cast to type @p F. */
        template<stopping::StopMethod F = stopping::StopCondition>
        [[nodiscard]]
        F& stop_condition() &;


        /* Move-only. */
        GaInfo(const GaInfo&)            = delete;
        GaInfo& operator=(const GaInfo&) = delete;

        GaInfo(GaInfo&&) noexcept;
        GaInfo& operator=(GaInfo&&) noexcept;

        virtual ~GaInfo();

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

        inline static constexpr size_t DEFAULT_POPSIZE = 100;
    };

} // namespace genetic_algorithm


/* IMPLEMENTATION */

#include "../selection/selection_base.hpp"
#include "../stop_condition/stop_condition_base.hpp"
#include <utility>
#include <stdexcept>

namespace genetic_algorithm
{
    template<typename F>
    requires selection::SelectionMethod<F> && std::is_final_v<F>
    void GaInfo::selection_method(F&& f)
    {
        selection_method(std::make_unique<F>(std::forward<F>(f)));
    }

    template<selection::SelectionMethod F>
    void GaInfo::selection_method(std::unique_ptr<F>&& f)
    {
        if (!f) throw std::invalid_argument("The selection method can't be a nullptr.");

        selection_ = std::move(f);
        can_continue_ = false;
    }

    template<selection::SelectionMethod F>
    F& GaInfo::selection_method() &
    {
        return dynamic_cast<F&>(*selection_);
    }

    template<typename F>
    requires stopping::StopMethod<F>&& std::is_final_v<F>
    void GaInfo::stop_condition(F&& f)
    {
        stop_condition(std::make_unique<F>(std::forward<F>(f)));
    }

    template<stopping::StopMethod F>
    void GaInfo::stop_condition(std::unique_ptr<F>&& f)
    {
        if (!f) throw std::invalid_argument("The stop condition can't be a nullptr.");

        stop_condition_ = std::move(f);
    }

    template<stopping::StopMethod F>
    F& GaInfo::stop_condition() &
    {
        return dynamic_cast<F&>(*stop_condition_);
    }

} // namespace genetic_algorithm

#endif // !GA_CORE_GA_INFO_HPP