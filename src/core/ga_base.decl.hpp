/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CORE_GA_BASE_DECL_HPP
#define GA_CORE_GA_BASE_DECL_HPP

#include "ga_info.hpp"
#include "fitness_function.hpp"
#include "../population/candidate.hpp"
#include "../mutation/mutation_base.fwd.hpp"
#include "../crossover/crossover_base.fwd.hpp"
#include <algorithm>
#include <vector>
#include <utility>
#include <functional>
#include <type_traits>
#include <concepts>
#include <memory>
#include <cstddef>

namespace genetic_algorithm
{
    /**
    * Base GA class. \n
    * Contains everything of the genetic algorithm, except for the initialize() and
    * generateCandidate() functions, which should be implemented in the derived classes.
    *
    * @tparam Gene The type of the genes in the candidates' chromosomes.
    */
    template<Gene T>
    class GA : public GaInfo
    {
    public:
        struct GeneBounds { T lower; T upper; };        /**< The type used to represent the lower and upper bounds of a gene. */
        
        using GeneType = T;                             /**< The gene type used in the chromosomes. */
        using BoundsVector = std::vector<GeneBounds>;   /**< A vector of gene bounds. */

        /**
        * The general callable type that can be used as a crossover method in the algorithms. \n
        * Takes two candidate solutions (parents), and returns two candidates generated from these
        * parent solutions (children). \n
        * 
        * The returned children's chromosomes don't have to be the same length as the parents, but
        * all of the other operators used should be able to handle different chromosome lengths in this case
        * (fitness function, mutation, repair).
        */
        using CrossoverCallable = std::function<CandidatePair<T>(const GA&, const Candidate<T>&, const Candidate<T>&)>;

        /**
        * The general callable type that can be used as a mutation operator in the algorithm. \n
        * Takes a candidate solution, and mutates this solution. \n
        * 
        * The function is allowed to change the candidate's chromosome's length, but
        * all of the other operators used should be able to handle different chromosome lengths in this
        * case (fitness function, crossover, repair).
        */
        using MutationCallable = std::function<void(const GA&, Candidate<T>&)>;

        /**
        * The general callable type that can be used as a repair function in the algorithm. \n
        * Takes a candidate solution, and performs some operation on it. \n
        * 
        * This function is allowed the change the candidate's chromosome's length, but
        * all of the other operators used should be able to handle different chromosome lengths in this
        * case (fitness function, crossover, mutation).
        */
        using RepairCallable = std::function<Chromosome<T>(const GA&, const Chromosome<T>&)>;


        /**
        * Create a genetic algorithm.
        *
        * @param fitness_function The fitness function used. @see FitnessFunction
        * @param population_size The number of candidates in the population.
        */
        GA(std::unique_ptr<FitnessFunction<T>> fitness_function, size_t population_size);

        /**
        * Set the fitness function used by the algorithm.
        * @see FitnessFunction
        *
        * @param f The fitness function to use in the algorithm.
        */
        template<typename F>
        requires FitnessFunctionType<F, T> && std::is_final_v<F>
        void fitness_function(F f);

        /**
        * Set the fitness function used by the algorithm.
        * @see FitnessFunction
        *
        * @param f The fitness function to use in the algorithm.
        */
        void fitness_function(std::unique_ptr<FitnessFunction<T>> f);

        /** @returns The fitness function used in the algorithm. */
        [[nodiscard]]
        FitnessFunction<T>& fitness_function() & noexcept;

        /** @returns The fitness function used in the algorithm. */
        [[nodiscard]]
        const FitnessFunction<T>& fitness_function() const& noexcept;

        /** @returns The chromosome length used for the candidates of the population. */
        [[nodiscard]]
        size_t chrom_len() const noexcept final;

        /** @returns True if variable chromosome lengths are allowed and used. */
        [[nodiscard]]
        bool variable_chrom_len() const noexcept final;

        /** @returns The number of objectives of the fitness function. */
        [[nodiscard]]
        size_t num_objectives() const noexcept final;

        /** @returns True if a dynamic fitness function is used. */
        [[nodiscard]]
        bool dynamic_fitness() const noexcept final;


        /** @returns The lower and upper bounds of each of the chromosomes' genes (the ranges are inclusive). */
        [[nodiscard]]
        const BoundsVector& gene_bounds() const noexcept;

        /**
        * Set the crossover method the algorithm will use. \n
        * The crossover function should be thread-safe if parallel execution is enabled (enabled by default).
        * 
        * @see Crossover
        * 
        * @param f The crossover method used by the algorithm.
        */
        template<typename F>
        requires crossover::CrossoverType<F, T> && std::is_final_v<F>
        void crossover_method(F f);

        /**
        * Set the crossover method the algorithm will use. \n
        * The crossover function should be thread-safe if parallel execution is enabled (enabled by default).
        * 
        * @see Crossover
        *
        * @param f The crossover method used by the algorithm. Can't be a nullptr.
        */
        void crossover_method(std::unique_ptr<crossover::Crossover<T>> f);

        /**
        * Set the crossover method the algorithm will use. \n
        * The crossover function should be thread-safe if parallel execution is enabled (enabled by default).
        * 
        * @see CrossoverCallable
        * @see Crossover
        *
        * @param f The crossover method used by the algorithm.
        */
        void crossover_method(CrossoverCallable f);

        /** @returns The crossover operator used by the algorithm. */
        [[nodiscard]]
        crossover::Crossover<T>& crossover_method() & noexcept;

        /** @returns The crossover operator used by the algorithm. */
        [[nodiscard]]
        const crossover::Crossover<T>& crossover_method() const& noexcept;

        /**
        * Set the crossover rate of the crossover operator used by the algorithm.
        *
        * @param pc The crossover probability. Must be in the closed interval [0.0, 1.0].
        */
        void crossover_rate(Probability pc) noexcept final;

        /** @returns The crossover rate set for the crossover operator. */
        [[nodiscard]]
        Probability crossover_rate() const noexcept final;

        /**
        * Set the mutation method the algorithm will use. \n
        * The mutation function should be thread-safe if parallel execution is enabled (enabled by default).
        * 
        * @see Mutation
        *
        * @param f The crossover function that will be used by the algorithm.
        */
        template<typename F>
        requires mutation::MutationType<F, T> && std::is_final_v<F>
        void mutation_method(F f);

        /**
        * Set the mutation method the algorithm will use. \n
        * The mutation function should be thread-safe if parallel execution is enabled (enabled by default).
        * 
        * @see Mutation
        *
        * @param f The crossover function that will be used by the algorithm. Can't be a nullptr.
        */
        void mutation_method(std::unique_ptr<mutation::Mutation<T>> f);

        /**
        * Set the mutation method the algorithm will use. \n
        * The mutation function should be thread-safe if parallel execution is enabled (enabled by default).
        * 
        * @see MutationCallable
        * @see Mutation
        *
        * @param f The crossover function that will be used by the algorithm.
        */
        void mutation_method(MutationCallable f);

        /** @returns The mutation operator used by the algorithm. */
        [[nodiscard]]
        mutation::Mutation<T>& mutation_method() & noexcept;

        /** @returns The mutation operator used by the algorithm. */
        [[nodiscard]]
        const mutation::Mutation<T>& mutation_method() const& noexcept;

        /**
        * Set the mutation rate of the mutation operator.
        *
        * @param pm The mutation probability. Must be in the closed interval [0.0, 1.0].
        */
        void mutation_rate(Probability pm) noexcept final;

        /** @returns The mutation rate of the mutation operator. */
        [[nodiscard]]
        Probability mutation_rate() const noexcept final;

        /**
        * Set a repair function for the algorithm. \n
        * This function will be applied to each chromosome after the mutations have been performed. \n
        * This can be used to perform a local search for example. \n
        * 
        * No repair function is used in the algorithm by default, set to nullptr if no repair function should be used.
        * 
        * The repair function should be thread-safe if parallel execution is enabled (enabled by default).
        * 
        * @param f The repair function the algorithm will use.
        */
        void repair_function(RepairCallable f);

        /** 
        * @returns The pareto optimal solutions found by the algorithm.
        * These are just the optimal solutions of the last generation's population if
        * keep_all_optimal_sols is not set, otherwise it contains every optimal solution
        * found during the run, updated in every generation.
        */
        [[nodiscard]]
        const Candidates<T>& solutions() const noexcept { return solutions_; }

        /**
        * @returns The current population of the algorithm.
        * This is not the same as the solutions, this will always be the same size as the population_size set,
        * and may also contain non pareto-optimal solutions.
        */
        [[nodiscard]]
        const Population<T>& population() const noexcept { return population_; }

        /**
        * Run the genetic algorithm. \n
        *
        * @param initial_population The population to use as the first generation of the algorithm.
        * @returns The pareto-optimal solutions found. @see solutions
        */
        Candidates<T> run(const Population<T>& initial_population = {});

        /**
        * Run the genetic algorithm for a set number of generations. \n
        * The algorithm will always stop when reaching this number of generations.
        *
        * @param num_generations The maximum number of generations.
        * @returns The pareto-optimal solutions found. @see solutions
        */
        Candidates<T> run(size_t num_generations);

        /**
        * Run the genetic algorithm for a set number of generations, using the specified initial population. \n
        * The algorithm will always stop when reaching this set number of generations. \n
        * Te initial population may have any number of solutions, but only the first popsize elements
        * will be used if it's greater than the population size set for the algorithm. If it has less
        * elements than the population size, the rest of the candidates will be generated using
        * the generateCandidate method.
        *
        * @param initial_population The population to use as the first generation of the algorithm.
        * @param num_generations The maximum number of generations.
        * @returns The pareto-optimal solutions found. @see solutions
        */
        Candidates<T> run(const Population<T>& initial_population, size_t num_generations);

        /**
        * Continue running the genetic algorithm for the set number of generations. \n
        * The algorithm can only continue if run has already been called at least once,
        * and the selection_method hasn't been changed since then.
        * (Otherwise this function will behave the same as calling run().)
        *
        * @param num_generations The number of generations to run the algorithm for.
        * @returns The optimal solutions.
        */
        Candidates<T> continueFor(size_t num_generations);

    protected:

        BoundsVector bounds_;

    private:

        Population<T> population_;
        Candidates<T> solutions_;

        std::unique_ptr<FitnessFunction<T>> fitness_function_;
        std::unique_ptr<crossover::Crossover<T>> crossover_;
        std::unique_ptr<mutation::Mutation<T>> mutation_;
        RepairCallable repair_ = nullptr;

        virtual void initialize() = 0;
        virtual Candidate<T> generateCandidate() const = 0;

        void initializeAlgorithm(const Population<T>& initial_population);
        Population<T> generatePopulation(size_t pop_size, const Population<T>& initial_population) const;
        void prepareSelections() const;
        const Candidate<T>& select() const;
        CandidatePair<T> crossover(const Candidate<T>& parent1, const Candidate<T>& parent2) const;
        void mutate(Candidate<T>& sol) const;
        void repair(Candidate<T>& sol) const;
        void updatePopulation(Population<T>&& children);
        bool stopCondition() const;

        void evaluate(Candidate<T>& sol);
        void updateOptimalSolutions(Candidates<T>& optimal_sols, const Population<T>& pop) const;

        void advance();


        bool hasValidFitness(const Candidate<T>& sol) const noexcept;
        bool hasValidChromosome(const Candidate<T>& sol) const noexcept;
        bool isValidEvaluatedPopulation(const Population<T>& pop) const;
        bool isValidUnevaluatedPopulation(const Population<T>& pop) const;
        bool fitnessMatrixIsSynced() const;


        /* Make the protected members of GaInfo private. */
        using GaInfo::fitness_matrix_;
        using GaInfo::algorithm_;
        using GaInfo::stop_condition_;
        using GaInfo::population_size_;
        using GaInfo::max_gen_;
        using GaInfo::generation_cntr_;
        using GaInfo::num_fitness_evals_;
        using GaInfo::keep_all_optimal_sols_;
        using GaInfo::is_initialized_;
    };

    /** Genetic algorithm types. */
    template<typename T>
    concept GeneticAlgorithmType = detail::DerivedFromSpecializationOf<T, GA>;

} // namespace genetic_algorithm

#endif // !GA_CORE_GA_BASE_DECL_HPP