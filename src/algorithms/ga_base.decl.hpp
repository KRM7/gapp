/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_GA_BASE_DECL_HPP
#define GA_GA_BASE_DECL_HPP

#include <algorithm>
#include <vector>
#include <unordered_set>
#include <utility>
#include <functional>
#include <type_traits>
#include <atomic>
#include <cstddef>
#include <memory>
#include "../candidate.hpp"
#include "../concepts.hpp"


namespace genetic_algorithm
{
    namespace selection
    {
        template<gene T>
        class Selection;
    }
    namespace crossover
    {
        template<gene T>
        class Crossover;
    }
    namespace mutation
    {
        template<gene T>
        class Mutation;
    }
    namespace stopping
    {
        template<gene T>
        class StopCondition;
    }

    /**
    * Base GA class.
    *
    * @tparam geneType The type of the genes in the candidates' chromosomes.
    */
    template<typename geneType>
    class GA
    {
    public:

        /** Structure containing stats of the single-objective algorithm. */
        struct History
        {
            std::vector<double> fitness_mean;   /**< The mean fitness values of each generation. */
            std::vector<double> fitness_sd;     /**< The standard deviation of the fitness values of each generation. */
            std::vector<double> fitness_min;    /**< The lowest fitness value in each generation. */
            std::vector<double> fitness_max;    /**< The highest fitness value in each generation. */

            void clear() noexcept;
            void reserve(size_t new_capacity);
            void add(double mean, double sd, double min, double max);
        };

        using GeneType = geneType;
        using Candidate = Candidate<GeneType>;  /**< The candidates used in the algorithm, each representing a solution to the problem. */

        using Chromosome = std::vector<GeneType>;                                           /**< . */
        using CandidatePair = CandidatePair<GeneType>;                                      /**< . */
        using CandidateVec = std::vector<Candidate>;                                        /**< . */
        using CandidateSet = std::unordered_set<Candidate, CandidateHasher<GeneType>>;      /**< . */
        using Population = std::vector<Candidate>;                                          /**< . */

        using FitnessFunction = std::function<std::vector<double>(const Chromosome&)>;      /**< The type of the fitness function. */
        using SelectionFunction = std::function<Candidate(const Population&)>;              /**< The type of the selection function. */
        using MutationFunction = std::function<void(Candidate&, double)>;                   /**< The type of the mutation function. */
        using RepairFunction = std::function<Chromosome(const Chromosome&)>;                /**< The type of the repair function. */
        using CallbackFunction = std::function<void(const GA&)>;

        /**
        * Should be set to false if the fitness function does not change over time. \n
        * (The fitness function will always return the same value for a given chromosome.) \n
        * Used to eliminate unnecesary fitness evaluations.
        */
        bool changing_fitness_func = false;

        /**
        * All pareto optimal optimal solutions found in the algorithm will be stored in the solutions,
        * not just the ones in the current population if this is set to true. \n
        * Setting it to false can speed up the algorithm.
        */
        bool archive_optimal_solutions = false;

        /**
        * The repair function applied to each Candidate of the population after the mutations if it isn't a nullptr. \n
        * This can be used to perform local search after the mutations, implementing a memetic algorithm.
        */
        RepairFunction repairFunction = nullptr;

        CallbackFunction endOfGenerationCallback = nullptr;

        /**
        * Standard constructor for the GA.
        *
        * @param chrom_len The length of the chromosomes (number of genes).
        * @param fitness_function The fitness function to find the maximum of.
        */
        GA(size_t chrom_len, FitnessFunction fitness_function);

        virtual ~GA() = default;

        /**
        * Runs the genetic algorithm with the selected settings.
        *
        * @returns The optimal solutions.
        */
        [[maybe_unused]] const CandidateVec& run();


        /** @returns A vector of the pareto optimal solutions found while running the algorithm. */
        [[nodiscard]] const CandidateVec& solutions() const;

        /** @returns The number of fitness evaluations performed while running the algorithm. */
        [[nodiscard]] size_t num_fitness_evals() const;

        /** @returns The current value of the generation counter. */
        [[nodiscard]] size_t generation_cntr() const;

        /** @returns The population of the final generation in the algorithm. */
        [[nodiscard]] const Population& population() const;

        /** @returns A History object containing stats from each generation of the single objective genetic algorithm. */
        [[nodiscard]] const History& soga_history() const;

        /**
        * Sets the length of the chromosomes (number of genes) of the Candidate solutions used in the algorithm to @p len. \n
        * The chromosome length must be at least 1.
        *
        * @param len The length of the chromosomes.
        */
        void chrom_len(size_t len);
        [[nodiscard]] size_t chrom_len() const;

        /**
        * Sets the number of Candidates used in the population to @p size. \n
        * The population size must be at least 1.
        *
        * @param size The size of the populations.
        */
        void population_size(size_t size);
        [[nodiscard]] size_t population_size() const;

        /**
        * Sets the maximum number of generations the algorithm runs for to @p max_gen. The
        * algorithm will always stop when this generation has been reached regardless of what stop
        * condition was set, but it can stop earlier when another stop condition is selected.
        * @see stop_condition @see StopCondition \n
        * The value of @p max_gen must be at least 1.
        *
        * @param max_gen The maximum number of generations.
        */
        void max_gen(size_t max_gen);
        [[nodiscard]] size_t max_gen() const;

        /**
        * Sets the initial population to be used in the algorithm to @p pop instead of randomly generating it. \n
        * If @p pop is empty, the initial population will be randomly generated. \n
        * If the preset population's size is not equal to the population size set, either additional randomly generated
        * Candidates will be added to fill out the initial population, or some Candidates will be discarded from the preset.
        *
        * @param pop The initial population to use in the algorithm.
        */
        void presetInitialPopulation(const Population& pop);

        /**
        * Sets the fitness function used by the algorithm to @p f. \n
        * The fitness function should return a vector whose size is equal to the number of objectives, and
        * each element of the vector should be finite.
        *
        * @param f The fitness function to find the maximum of.
        */
        void setFitnessFunction(FitnessFunction f);

        [[nodiscard]] size_t num_objectives() const noexcept;


        /* CROSSOVER */

        template<typename CrossoverType>
        //requires std::derived_from<CrossoverType, crossover::Crossover<GeneType>> && std::copy_constructible<CrossoverType>
        void crossover_method(const CrossoverType& f);

        void crossover_method(std::unique_ptr<crossover::Crossover<GeneType>>&& f);

        template<typename CrossoverType = crossover::Crossover<GeneType>>
        //requires std::derived_from<CrossoverType, crossover::Crossover<GeneType>>
        CrossoverType& crossover_method() const;


        /* MUTATION */

        template<typename MutationType>
        //requires std::derived_from<MutationType, mutation::Mutation<GeneType>>&& std::copy_constructible<MutationType>
        void mutation_method(const MutationType& f);

        void mutation_method(std::unique_ptr<mutation::Mutation<GeneType>>&& f);

        template<typename MutationType = mutation::Mutation<GeneType>>
        //requires std::derived_from<MutationType, mutation::Crossover<GeneType>>
        MutationType& mutation_method() const;


        /* STOP CONDITION */

        /* The algorithm always stops when the set max_gen generation has been reached regardless of the stop condition
        set. */
        template<typename StopType>
        //requires std::derived_from<StopType, stopping::StopCondition<GeneType>> && std::copy_constructible<StopType>
        void stop_condition(const StopType& f);

        void stop_condition(std::unique_ptr<stopping::StopCondition<GeneType>>&& f);
        
        template<typename StopType = stopping::StopCondition<GeneType>>
        //requires std::derived_from<StopType, stopping::StopCondition<GeneType>>
        StopType& stop_condition() const;


        /* SELECTION METHOD */

        /**
        * Sets the selection method used in the single-objective algorithm to @p method. \n
        * The selection method set is ignored in the other algorithm types. @see Mode
        *
        * @param method The selection method used in the single-objective algorithm.
        */
        template<typename SelectionType>
        //requires std::derived_from<SelectionType, selection::Selection<GeneType>> && std::copy_constructible<SelectionType>
        void selection_method(const SelectionType& f);

        void selection_method(std::unique_ptr<selection::Selection<GeneType>>&& f);

        template<typename SelectionType = selection::Selection<GeneType>>
        SelectionType& selection_method() const;

    protected:

        Population population_;
        size_t generation_cntr_ = 0;
        size_t num_objectives_ = 0;        /* Determined from the fitness function. */

        /* Results. */
        CandidateVec solutions_;
        std::atomic<size_t> num_fitness_evals_ = 0;
        History soga_history_;

        /* Basic parameters of the GA. */
        size_t chrom_len_;
        size_t population_size_ = 100;
        size_t max_gen_ = 500;

        /* Initial population settings. */
        Population initial_population_preset_;

        /* User supplied functions used in the GA. All of these are optional except for the fitness function. */
        FitnessFunction fitnessFunction;
        //selectionFunction_t customSelection = nullptr;

        std::unique_ptr<selection::Selection<GeneType>> selection_;
        std::unique_ptr<crossover::Crossover<GeneType>> crossover_;
        std::unique_ptr<mutation::Mutation<GeneType>> mutation_;
        std::unique_ptr<stopping::StopCondition<GeneType>> stop_condition_;


        /* General functions for the genetic algorithms. */

        void init();
        virtual Candidate generateCandidate() const = 0;
        Population generateInitialPopulation() const;
        bool stopCondition();
        void evaluate(Population& pop);
        void updateOptimalSolutions(CandidateVec& optimal_sols, const Population& pop) const;
        void repair(Population& pop) const;
        void updateStats(const Population& pop);

    };

} // namespace genetic_algorithm

#endif // !GA_GA_BASE_DECL_HPP