/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CORE_GA_BASE_DECL_HPP
#define GA_CORE_GA_BASE_DECL_HPP

#include "ga_info.hpp"
#include "ga_traits.hpp"
#include "fitness_function.hpp"
#include "../encoding/gene_types.hpp"
#include "../population/candidate.hpp"
#include "../stop_condition/stop_condition.hpp"
#include "../utility/bounded_value.hpp"
#include "../utility/type_traits.hpp"
#include <algorithm>
#include <vector>
#include <utility>
#include <functional>
#include <type_traits>
#include <concepts>
#include <memory>
#include <cstddef>

namespace genetic_algorithm::crossover
{
    template<typename T>
    class Crossover;

} // namespace genetic_algorithm::crossover

namespace genetic_algorithm::mutation
{
    template<typename T>
    class Mutation;

} // namespace genetic_algorithm::mutation

namespace genetic_algorithm
{
    /**
    * Base GA class. \n
    * Contains everything of the genetic algorithm, except for the initialize() and
    * generateCandidate() functions, which should be implemented in the derived classes.
    *
    * @tparam Gene The type of the genes in the candidates' chromosomes.
    */
    template<typename T>
    class GA : public GaInfo
    {
    public:
        /** The gene type used for the chromosomes. */
        using GeneType = T;

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
        * Create a genetic algorithm using the default genetic operators. \n
        * The algorithm used will be deduced from the number of objectives of the fitness function
        * (single- or multi-objective) along with the mutation probability used for the mutation operator.
        *
        * @param population_size The number of candidates in the population. Must be at least 1.
        */
        GA(Positive<size_t> population_size = DEFAULT_POPSIZE);

        /**
        * Create a genetic algorithm using the default genetic operators. \n
        * The mutation probability used for the mutation operator will be deduced from the fitness function.
        *
        * @param population_size The number of candidates in the population. Must be at least 1.
        * @param algorithm The algorithm to use. Can't be a nullptr.
        */
        GA(Positive<size_t> population_size, std::unique_ptr<algorithm::Algorithm> algorithm);

        /**
        * Create a genetic algorithm. \n
        * The algorithm used will be deduced from the number of objectives of the fitness function
        * (single- or multi-objective).
        *
        * @param population_size The number of candidates in the population. Must be at least 1.
        * @param crossover The crossover operator to use in the algorithm. Can't be a nullptr.
        * @param mutation The mutation operator to use in the algorithm. Can't be a nullptr.
        * @param stop_condition The early-stop condition to use. Can't be a nullptr.
        */
        GA(Positive<size_t> population_size,
           std::unique_ptr<crossover::Crossover<T>> crossover,
           std::unique_ptr<mutation::Mutation<T>> mutation,
           std::unique_ptr<stopping::StopCondition> stop_condition = std::make_unique<stopping::NoEarlyStop>());

        /**
        * Create a genetic algorithm.
        *
        * @param population_size The number of candidates in the population. Must be at least 1.
        * @param algorithm The algorithm to use. Can't be a nullptr.
        * @param crossover The crossover operator to use in the algorithm. Can't be a nullptr.
        * @param mutation The mutation operator to use in the algorithm. Can't be a nullptr.
        * @param stop_condition The early-stop condition to use. Can't be a nullptr.
        */
        GA(Positive<size_t> population_size,
           std::unique_ptr<algorithm::Algorithm> algorithm,
           std::unique_ptr<crossover::Crossover<T>> crossover,
           std::unique_ptr<mutation::Mutation<T>> mutation,
           std::unique_ptr<stopping::StopCondition> stop_condition = std::make_unique<stopping::NoEarlyStop>());

        /**
        * Create a genetic algorithm using the default genetic operators. \n
        * The mutation probability used for the mutation operator will be deduced from the fitness function.
        *
        * @param population_size The number of candidates in the population. Must be at least 1.
        * @param algorithm The algorithm to use. Can't be a nullptr.
        */
        template<typename AlgorithmType>
        requires std::derived_from<AlgorithmType, algorithm::Algorithm> && std::is_final_v<AlgorithmType>
        GA(Positive<size_t> population_size, AlgorithmType algorithm);

        /**
        * Create a genetic algorithm. \n
        * The algorithm used will be deduced from the number of objectives of the fitness function
        * (single- or multi-objective).
        *
        * @param population_size The number of candidates in the population. Must be at least 1.
        * @param crossover The crossover operator to use in the algorithm.
        * @param mutation The mutation operator to use in the algorithm.
        * @param stop_condition The early-stop condition to use.
        */
        template<typename CrossoverType, typename MutationType, typename StoppingType = stopping::NoEarlyStop>
        requires std::derived_from<CrossoverType, crossover::Crossover<T>> && std::is_final_v<CrossoverType> &&
                 std::derived_from<MutationType, mutation::Mutation<T>> && std::is_final_v<MutationType> &&
                 std::derived_from<StoppingType, stopping::StopCondition> && std::is_final_v<StoppingType>
        GA(Positive<size_t> population_size, CrossoverType crossover, MutationType mutation, StoppingType stop_condition = {});

        /**
        * Create a genetic algorithm.
        *
        * @param population_size The number of candidates in the population. Must be at least 1.
        * @param algorithm The algorithm to use.
        * @param crossover The crossover operator to use in the algorithm.
        * @param mutation The mutation operator to use in the algorithm.
        * @param stop_condition The early-stop condition to use.
        */
        template<typename AlgorithmType, typename CrossoverType, typename MutationType, typename StoppingType = stopping::NoEarlyStop>
        requires std::derived_from<AlgorithmType, algorithm::Algorithm> && std::is_final_v<AlgorithmType> &&
                 std::derived_from<CrossoverType, crossover::Crossover<T>> && std::is_final_v<CrossoverType> &&
                 std::derived_from<MutationType, mutation::Mutation<T>> && std::is_final_v<MutationType> &&
                 std::derived_from<StoppingType, stopping::StopCondition> && std::is_final_v<StoppingType>
        GA(Positive<size_t> population_size, AlgorithmType algorithm, CrossoverType crossover, MutationType mutation, StoppingType stop_condition = {});


        /** @returns The chromosome length used for the candidates of the population. */
        [[nodiscard]]
        size_t chrom_len() const noexcept final;

        /** @returns True if variable chromosome lengths are allowed and used. */
        [[nodiscard]]
        bool variable_chrom_len() const noexcept final;

        /** @returns True if a dynamic fitness function is used. */
        [[nodiscard]]
        bool dynamic_fitness() const noexcept final;


        /** @returns The lower and upper bounds of each of the chromosomes' genes (the ranges are inclusive). */
        [[nodiscard]]
        const BoundsVector<T>& gene_bounds() const noexcept requires (is_bounded<T>);

        /**
        * Set the crossover method the algorithm will use. \n
        * The crossover function should be thread-safe if parallel execution is enabled (enabled by default).
        * 
        * @see Crossover
        * 
        * @param f The crossover method used by the algorithm.
        */
        template<typename F>
        requires std::derived_from<F, crossover::Crossover<T>> && std::is_final_v<F>
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
        requires std::derived_from<F, mutation::Mutation<T>> && std::is_final_v<F>
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
        * No repair function is used in the algorithm by default, set it to nullptr if no repair function should be used.
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
        const Candidates<T>& solutions() const& noexcept { return solutions_; }

        /**
        * @returns The current population of the algorithm.
        * This is not the same as the solutions, this will always be the same size as the population_size set,
        * and may also contain non pareto-optimal solutions.
        */
        [[nodiscard]]
        const Population<T>& population() const& noexcept { return population_; }


        /********************* NON-BOUNDED ALGORITHM SOLVE *********************/

        /**
        * Run the genetic algorithm for a set number of generations for a fitness function, using the specified initial population. \n
        * The initial population may have any number of solutions, but only the first popsize elements
        * will be used if it's greater than the population size set for the algorithm. If it has less
        * elements than the population size, the rest of the candidates will be generated using
        * the generateCandidate method.
        *
        * @param fitness_function The fitness function the algorithm should find the maximum of. Can't be a nullptr.
        * @param generations The maximum number of generations. The algorithm will always stop when reaching this, even if another stop condition was set.
        * @param initial_population The population to use as the first generation of the algorithm.
        * @returns The pareto-optimal solutions found (this is not the final population).
        */
        Candidates<T> solve(std::unique_ptr<FitnessFunctionBase<T>> fitness_function, size_t generations, Population<T> initial_population = {}) requires (!is_bounded<T>);

        /**
        * Run the genetic algorithm for a fitness function, using the specified initial population. \n
        * The initial population may have any number of solutions, but only the first popsize elements
        * will be used if it's greater than the population size set for the algorithm. If it has less
        * elements than the population size, the rest of the candidates will be generated using
        * the generateCandidate method.
        *
        * @param fitness_function The fitness function the algorithm should find the maximum of. Can't be a nullptr.
        * @param initial_population The population to use as the first generation of the algorithm.
        * @returns The pareto-optimal solutions found (this is not the final population).
        */
        Candidates<T> solve(std::unique_ptr<FitnessFunctionBase<T>> fitness_function, Population<T> initial_population = {}) requires (!is_bounded<T>);

        /**
        * Run the genetic algorithm for a fitness function, using the specified initial population. \n
        * The initial population may have any number of solutions, but only the first popsize elements
        * will be used if it's greater than the population size set for the algorithm. If it has less
        * elements than the population size, the rest of the candidates will be generated using
        * the generateCandidate method.
        *
        * @param fitness_function The fitness function the algorithm should find the maximum of.
        * @param initial_population The population to use as the first generation of the algorithm.
        * @returns The pareto-optimal solutions found (this is not the final population).
        */
        template<typename F>
        requires (!is_bounded<T> && std::derived_from<F, FitnessFunctionBase<T>> && std::is_final_v<F>)
        Candidates<T> solve(F fitness_function, Population<T> initial_population = {});

        /**
        * Run the genetic algorithm for a set number of generations for a fitness function, using the specified initial population. \n
        * The initial population may have any number of solutions, but only the first popsize elements
        * will be used if it's greater than the population size set for the algorithm. If it has less
        * elements than the population size, the rest of the candidates will be generated using
        * the generateCandidate method.
        *
        * @param fitness_function The fitness function the algorithm should find the maximum of. Can't be a nullptr.
        * @param generations The maximum number of generations. The algorithm will always stop when reaching this, even if another stop condition was set.
        * @param initial_population The population to use as the first generation of the algorithm.
        * @returns The pareto-optimal solutions found (this is not the final population).
        */
        template<typename F>
        requires (!is_bounded<T> && std::derived_from<F, FitnessFunctionBase<T>> && std::is_final_v<F>)
        Candidates<T> solve(F fitness_function, size_t generations, Population<T> initial_population = {});


        /********************* BOUNDED ALGORITHM SOLVE *********************/

        /**
        * Run the genetic algorithm for a set number of generations for a fitness function, using the specified initial population. \n
        * The initial population may have any number of solutions, but only the first popsize elements
        * will be used if it's greater than the population size set for the algorithm. If it has less
        * elements than the population size, the rest of the candidates will be generated using
        * the generateCandidate method.
        *
        * @param fitness_function The fitness function the algorithm should find the maximum of. Can't be a nullptr.
        * @param bounds The lower and upper bounds of each of the chromosomes' genes.
        * @param generations The maximum number of generations. The algorithm will always stop when reaching this, even if another stop condition was set.
        * @param initial_population The population to use as the first generation of the algorithm.
        * @returns The pareto-optimal solutions found (this is not the final population).
        */
        Candidates<T> solve(std::unique_ptr<FitnessFunctionBase<T>> fitness_function, BoundsVector<T> bounds, size_t generations, Population<T> initial_population = {}) requires (is_bounded<T>);

        /**
        * Run the genetic algorithm for a set number of generations for a fitness function, using the specified initial population. \n
        * The initial population may have any number of solutions, but only the first popsize elements
        * will be used if it's greater than the population size set for the algorithm. If it has less
        * elements than the population size, the rest of the candidates will be generated using
        * the generateCandidate method.
        *
        * @param fitness_function The fitness function the algorithm should find the maximum of. Can't be a nullptr.
        * @param bounds The lower and upper bounds of a gene. These bounds will be used for every gene of the chromosomes.
        * @param generations The maximum number of generations. The algorithm will always stop when reaching this, even if another stop condition was set.
        * @param initial_population The population to use as the first generation of the algorithm.
        * @returns The pareto-optimal solutions found (this is not the final population).
        */
        Candidates<T> solve(std::unique_ptr<FitnessFunctionBase<T>> fitness_function, Bounds<T> bounds, size_t generations, Population<T> initial_population = {}) requires (is_bounded<T>);

        /**
        * Run the genetic algorithm for a set number of generations for a fitness function, using the specified initial population. \n
        * The initial population may have any number of solutions, but only the first popsize elements
        * will be used if it's greater than the population size set for the algorithm. If it has less
        * elements than the population size, the rest of the candidates will be generated using
        * the generateCandidate method.
        *
        * @param fitness_function The fitness function the algorithm should find the maximum of. Can't be a nullptr.
        * @param bounds The lower and upper bounds of each of the chromosomes' genes.
        * @param generations The maximum number of generations. The algorithm will always stop when reaching this, even if another stop condition was set.
        * @param initial_population The population to use as the first generation of the algorithm.
        * @returns The pareto-optimal solutions found (this is not the final population).
        */
        Candidates<T> solve(std::unique_ptr<FitnessFunctionBase<T>> fitness_function, BoundsVector<T> bounds, Population<T> initial_population = {}) requires (is_bounded<T>);

        /**
        * Run the genetic algorithm for a set number of generations for a fitness function, using the specified initial population. \n
        * The initial population may have any number of solutions, but only the first popsize elements
        * will be used if it's greater than the population size set for the algorithm. If it has less
        * elements than the population size, the rest of the candidates will be generated using
        * the generateCandidate method.
        *
        * @param fitness_function The fitness function the algorithm should find the maximum of. Can't be a nullptr.
        * @param bounds The lower and upper bounds of a gene. These bounds will be used for every gene of the chromosomes.
        * @param generations The maximum number of generations. The algorithm will always stop when reaching this, even if another stop condition was set.
        * @param initial_population The population to use as the first generation of the algorithm.
        * @returns The pareto-optimal solutions found (this is not the final population).
        */
        Candidates<T> solve(std::unique_ptr<FitnessFunctionBase<T>> fitness_function, Bounds<T> bounds, Population<T> initial_population = {}) requires (is_bounded<T>);

        /**
        * Run the genetic algorithm for a set number of generations for a fitness function, using the specified initial population. \n
        * The initial population may have any number of solutions, but only the first popsize elements
        * will be used if it's greater than the population size set for the algorithm. If it has less
        * elements than the population size, the rest of the candidates will be generated using
        * the generateCandidate method.
        *
        * @param fitness_function The fitness function the algorithm should find the maximum of.
        * @param bounds The lower and upper bounds of each of the chromosomes' genes.
        * @param generations The maximum number of generations. The algorithm will always stop when reaching this, even if another stop condition was set.
        * @param initial_population The population to use as the first generation of the algorithm.
        * @returns The pareto-optimal solutions found (this is not the final population).
        */
        template<typename F>
        requires (is_bounded<T> && std::derived_from<F, FitnessFunctionBase<T>> && std::is_final_v<F>)
        Candidates<T> solve(F fitness_function, BoundsVector<T> bounds, Population<T> initial_population = {});

        /**
        * Run the genetic algorithm for a set number of generations for a fitness function, using the specified initial population. \n
        * The initial population may have any number of solutions, but only the first popsize elements
        * will be used if it's greater than the population size set for the algorithm. If it has less
        * elements than the population size, the rest of the candidates will be generated using
        * the generateCandidate method.
        *
        * @param fitness_function The fitness function the algorithm should find the maximum of.
        * @param bounds The lower and upper bounds of a gene. These bounds will be used for every gene of the chromosomes.
        * @param generations The maximum number of generations. The algorithm will always stop when reaching this, even if another stop condition was set.
        * @param initial_population The population to use as the first generation of the algorithm.
        * @returns The pareto-optimal solutions found (this is not the final population).
        */
        template<typename F>
        requires (is_bounded<T> && std::derived_from<F, FitnessFunctionBase<T>> && std::is_final_v<F>)
        Candidates<T> solve(F fitness_function, Bounds<T> bounds, Population<T> initial_population = {});

        /**
        * Run the genetic algorithm for a set number of generations for a fitness function, using the specified initial population. \n
        * The initial population may have any number of solutions, but only the first popsize elements
        * will be used if it's greater than the population size set for the algorithm. If it has less
        * elements than the population size, the rest of the candidates will be generated using
        * the generateCandidate method.
        *
        * @param fitness_function The fitness function the algorithm should find the maximum of.
        * @param bounds The lower and upper bounds of each of the chromosomes' genes.
        * @param generations The maximum number of generations. The algorithm will always stop when reaching this, even if another stop condition was set.
        * @param initial_population The population to use as the first generation of the algorithm.
        * @returns The pareto-optimal solutions found (this is not the final population).
        */
        template<typename F>
        requires (is_bounded<T> && std::derived_from<F, FitnessFunctionBase<T>> && std::is_final_v<F>)
        Candidates<T> solve(F fitness_function, BoundsVector<T> bounds, size_t generations, Population<T> initial_population = {});

        /**
        * Run the genetic algorithm for a set number of generations for a fitness function, using the specified initial population. \n
        * The initial population may have any number of solutions, but only the first popsize elements
        * will be used if it's greater than the population size set for the algorithm. If it has less
        * elements than the population size, the rest of the candidates will be generated using
        * the generateCandidate method.
        *
        * @param fitness_function The fitness function the algorithm should find the maximum of. Can't be a nullptr.
        * @param bounds The lower and upper bounds of a gene. These bounds will be used for every gene of the chromosomes.
        * @param generations The maximum number of generations. The algorithm will always stop when reaching this, even if another stop condition was set.
        * @param initial_population The population to use as the first generation of the algorithm.
        * @returns The pareto-optimal solutions found (this is not the final population).
        */
        template<typename F>
        requires (is_bounded<T> && std::derived_from<F, FitnessFunctionBase<T>> && std::is_final_v<F>)
        Candidates<T> solve(F fitness_function, Bounds<T> bounds, size_t generations, Population<T> initial_population = {});

    private:

        using MaybeBoundsVector = std::conditional_t<is_bounded<T>, BoundsVector<T>, detail::empty_t>;

        Population<T> population_;
        Candidates<T> solutions_;

        [[no_unique_address]]
        MaybeBoundsVector bounds_;

        std::unique_ptr<FitnessFunctionBase<T>> fitness_function_;
        std::unique_ptr<crossover::Crossover<T>> crossover_;
        std::unique_ptr<mutation::Mutation<T>> mutation_;
        RepairCallable repair_ = nullptr;

        bool use_default_mutation_rate_ = false;


        virtual void initialize() {};
        virtual Candidate<T> generateCandidate() const = 0;

        void initializeAlgorithm(MaybeBoundsVector bounds, Population<T> initial_population);
        Population<T> generatePopulation(Positive<size_t> pop_size, Population<T> initial_population) const;
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


        size_t findNumberOfObjectives() const;
        std::unique_ptr<algorithm::Algorithm> defaultAlgorithm() const;
        Probability defaultMutationRate() const;

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
        using GaInfo::num_objectives_;
        using GaInfo::generation_cntr_;
        using GaInfo::num_fitness_evals_;
        using GaInfo::keep_all_optimal_sols_;
        using GaInfo::use_default_algorithm_;
    };

} // namespace genetic_algorithm

#endif // !GA_CORE_GA_BASE_DECL_HPP