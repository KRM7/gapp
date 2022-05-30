/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_GA_BASE_DECL_HPP
#define GA_GA_BASE_DECL_HPP

#include "ga_info.hpp"
#include "../population/candidate.hpp"
#include "../utility/concepts.hpp"
#include <algorithm>
#include <vector>
#include <unordered_set>
#include <utility>
#include <functional>
#include <type_traits>
#include <atomic>
#include <cstddef>
#include <memory>
#include <concepts>

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
    namespace crossover
    {
        template<Gene T>
        class Crossover;

        template<typename T, typename G>
        concept CrossoverMethod = requires
        {
            requires Gene<G>;
            requires std::derived_from<T, Crossover<G>>;
            requires std::copy_constructible<T>;
        };
    }
    namespace mutation
    {
        template<Gene T>
        class Mutation;

        template<typename T, typename G>
        concept MutationMethod = requires
        {
            requires Gene<G>;
            requires std::derived_from<T, Mutation<G>>;
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

    /**
    * Base GA class.
    *
    * @tparam geneType The type of the genes in the candidates' chromosomes.
    */
    template<Gene T, typename Derived>
    class GA : public GaInfo
    {
    public:
        using GeneType = T;                             /**< The gene type used in the algorithm. */
        using Candidate = Candidate<GeneType>;          /**< The candidate type used in the algorithm. */
        using Chromosome = std::vector<GeneType>;       /**< The type of the chromosomes used in the algorithm. */
        using CandidatePair = CandidatePair<GeneType>;  /**< A pair of Candidates. */
        using Candidates = std::vector<Candidate>;      /**< A vector of Candidates. */
        using Population = std::vector<Candidate>;      /**< The population type used in the algorithm. */

        using FitnessFunction = std::function<std::vector<double>(const Chromosome&)>;      /**< The type of the fitness function. */
        using CrossoverFunction = std::function<CandidatePair(const GaInfo&, const Candidate&, const Candidate&)>;
        using MutationFunction = std::function<void(const GaInfo&, Candidate&)>;
        using StopConditionFunction = std::function<bool(const GaInfo&)>;
        using RepairFunction = std::function<Chromosome(const Chromosome&)>;                /**< The type of the repair function. */
        using CallbackFunction = std::function<void(const GA<T, Derived>&)>;

        /**
        * Create a genetic algorithm.
        *
        * @param chrom_len The length of the chromosomes (number of genes).
        * @param fitness_function The fitness function to find the maximum of.
        */
        GA(size_t chrom_len, FitnessFunction fitness_function);

        /**
        * Create a genetic algorithm.
        *
        * @param population_size The number of candidates in the population.
        * @param chrom_len The length of the chromosomes (number of genes).
        * @param fitness_function The fitness function to find the maximum of.
        */
        GA(size_t population_size, size_t chrom_len, FitnessFunction fitness_function);

        /**
        * Run the genetic algorithm for the set number of generations. \n
        * The algorithm will always stop when reaching the maximum number of
        * generations have been reached.
        *
        * @param num_generations The maximum number of generations.
        * @returns The optimal solutions.
        */
        [[maybe_unused]] const Candidates& run(size_t num_generations = 500);

        /**
        * Continue running the genetic algorithm for the set number of generations. \n
        * Equivalent to calling run if the algorithm hasn't been run before, or if the selection
        * method was changed.
        *
        * @param num_generations The number of generations to run the algorithm for.
        * @returns The optimal solutions.
        */
        [[maybe_unused]] const Candidates& continueFor(size_t num_generations);

        /** @returns The pareto optimal solutions found by the algorithm. */
        [[nodiscard]] const Candidates& solutions() const noexcept;

        /** @returns The current population of the algorithm. Not the same as the solutions. */
        [[nodiscard]] const Population& population() const noexcept;

        /** @returns The fitness matrix of the current population of the algorithm. */
        [[nodiscard]] FitnessMatrix fitness_matrix() const override final;

        /**
        * Set the initial population to be used in the algorithm to @p pop instead of randomly generating it. \n
        * If @p pop is empty, the initial population will be randomly generated. \n
        * If the preset population's size is not equal to the population size set, either additional randomly generated
        * Candidates will be added to fill out the initial population, or some Candidates will be discarded from the preset.
        *
        * @param pop The initial population to use in the algorithm.
        */
        void initial_population(const Population& pop);

        /**
        * Set the fitness function used by the algorithm to @p f. \n
        * The fitness function should return a vector with a size equal to the number of objectives.
        *
        * @param f The fitness function to find the maximum of.
        */
        void fitness_function(FitnessFunction f);

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
        [[nodiscard]] F& selection_method();

        /**
        * Set the crossover method the algorithm will use to @p f.
        * 
        * @param f The crossover method for the algorithm.
        */
        template<crossover::CrossoverMethod<T> F>
        void crossover_method(const F& f);

        /**
        * Set the crossover method the algorithm will use to @p f.
        *
        * @param f The crossover method for the algorithm.
        */
        template<crossover::CrossoverMethod<T> F>
        void crossover_method(std::unique_ptr<F>&& f);

        /**
        * Set the crossover method the algorithm will use to @p f.
        *
        * @param f The crossover method for the algorithm.
        */
        void crossover_method(CrossoverFunction f);

        /** @returns The current crossover method used by the algorithm. */
        template<crossover::CrossoverMethod<T> F = crossover::Crossover<GeneType>>
        F& crossover_method();

        /**
        * Sets the crossover rate used in the algorithm to @p pc.
        *
        * @param pc The crossover probability. Must be in the closed interval [0.0, 1.0].
        */
        void crossover_rate(double pc);

        /** @returns The current crossover rate set. */
        [[nodiscard]] double crossover_rate() const noexcept;

        /**
        * Set the mutation method the algorithm will use to @p f.
        *
        * @param f The crossover method for the algorithm.
        */
        template<mutation::MutationMethod<T> F>
        void mutation_method(const F& f);

        /**
        * Set the mutation method the algorithm will use to @p f.
        *
        * @param f The crossover method for the algorithm.
        */
        template<mutation::MutationMethod<T> F>
        void mutation_method(std::unique_ptr<F>&& f);

        /**
        * Set the mutation method the algorithm will use to @p f.
        *
        * @param f The crossover method for the algorithm.
        */
        void mutation_method(MutationFunction f);

        /** @returns The current mutation method used by the algorithm. */
        template<mutation::MutationMethod<T> F = mutation::Mutation<GeneType>>
        F& mutation_method();

        /**
        * Sets the mutation rate used in the algorithm to @p pm.
        *
        * @param pm The mutation probability. Must be in the closed interval [0.0, 1.0].
        */
        void mutation_rate(double pm);

        /** @returns The current crossover rate set. */
        [[nodiscard]] double mutation_rate() const noexcept;

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
        [[nodiscard]] F& stop_condition();

        /**
        * Set a repair function for the genetic algorithm. \n
        * This function will be applied to each chromosome after the mutations. \n
        * This can be used for example to perform a local search. \n
        * 
        * No repair function is used in the algorithm by default, set to nullptr if
        * no repair function should be used.
        * 
        * @param f The repair function the algorithm will use.
        */
        void repair_function(const RepairFunction& f);

        CallbackFunction endOfGenerationCallback = nullptr;

    protected:

        Population population_;
        Candidates solutions_;
        Population initial_population_; /* Only if preset pop is used. */

        FitnessFunction fitness_function_;
        std::unique_ptr<selection::Selection> selection_;
        std::unique_ptr<crossover::Crossover<GeneType>> crossover_;
        std::unique_ptr<mutation::Mutation<GeneType>> mutation_;
        std::unique_ptr<stopping::StopCondition> stop_condition_;
        RepairFunction repair_ = nullptr;

        bool can_continue_ = false;

        void initialize();
        size_t getNumObjectives(const FitnessFunction& f) const;
        Population generateInitialPopulation() const;
        void evaluate(Population& pop);
        Population nextPopulation(Population& pop, Population& children) const;
        void updateOptimalSolutions(Candidates& optimal_sols, const Population& pop) const;
        void repair(Population& pop) const;
        bool stopCondition() const;
        void advance();

    private:

        constexpr Derived& derived() noexcept;
        constexpr const Derived& derived() const noexcept;
        Candidate generateCandidate() const;
    };

    /** Genetic algorithm types. */
    template<typename T>
    concept GeneticAlgorithm = detail::DerivedFromSpecializationOf<T, GA>;

} // namespace genetic_algorithm

#endif // !GA_GA_BASE_DECL_HPP