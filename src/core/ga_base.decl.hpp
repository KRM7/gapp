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
    * The base class used for all of the genetic algorithms. Contains all of the
    * encoding/gene specific information of the %GA in addition to the non-encoding
    * dependent parts contained by GaInfo.
    * 
    * This class should be used as the base class of new genetic algorithms. It has
    * 2 virtual functions that should be implemented in the derived classes:
    *   - initialize (optional) : Perform some additional initialization specific to the derived class.
    *   - generateCandidate     : Returns a (randomly) generated candidate solution. Used to create the initial population.
    * 
    * @note Before declaring a class that is derived from this class, the GaTraits struct should be specialized for the
    *   particular gene type that is used. If the gene type is bounded, the is_bounded variable template should
    *   also be specialized.
    * 
    * @see GaTraits, is_bounded
    *
    * @tparam T The gene type used for the candidate's chromosomes.
    */
    template<typename T>
    class GA : public GaInfo
    {
    public:
        /** The gene type of the candidates. */
        using GeneType = T;

        /**
        * The general callable type that can be used as a crossover method in the %GA (when
        * not using a crossover method that is derived from crossover::Crossover).
        * The function takes two candidate solutions (parents), and returns two candidates generated
        * from these (children).
        * 
        * @see crossover_method
        */
        using CrossoverCallable = std::function<CandidatePair<T>(const GA&, const Candidate<T>&, const Candidate<T>&)>;

        /**
        * The general callable type that can be used as a mutation method in the %GA (when
        * not using a mutation method derived from mutation::Mutation).
        * The function takes a candidate solution, and changes this solution's chromosome.
        * 
        * @see mutation_method
        */
        using MutationCallable = std::function<void(const GA&, Candidate<T>&)>;

        /**
        * The general callable type that can be used as a repair function in the %GA.
        * The function takes a candidate solution, and performs some operation on it.
        * 
        * @see repair_function
        */
        using RepairCallable = std::function<Chromosome<T>(const GA&, const Chromosome<T>&)>;


        /**
        * Create a genetic algorithm using the default genetic operators.
        * 
        * The algorithm used will be deduced from the number of objectives of the fitness function
        * (single- or multi-objective), along with the mutation probability used, which will be
        * deduced from the chromosome length.
        *
        * @param population_size The number of candidates in the population. Must be at least 1.
        */
        GA(Positive<size_t> population_size = DEFAULT_POPSIZE);

        /**
        * Create a genetic algorithm using the default genetic operators.
        * The mutation probability used will be deduced from the chromosome length.
        *
        * @param population_size The number of candidates in the population. Must be at least 1.
        * @param algorithm The algorithm to use. Can't be a nullptr.
        */
        GA(Positive<size_t> population_size, std::unique_ptr<algorithm::Algorithm> algorithm);

        /**
        * Create a genetic algorithm using the specified operators.
        * The algorithm used will be deduced from the number of objectives of the fitness function
        * (single- or multi-objective).
        *
        * @param population_size The number of candidates in the population. Must be at least 1.
        * @param crossover The crossover operator to use. Can't be a nullptr.
        * @param mutation The mutation operator to use. Can't be a nullptr.
        * @param stop_condition The early-stop condition to use. Can't be a nullptr.
        */
        GA(Positive<size_t> population_size,
           std::unique_ptr<crossover::Crossover<T>> crossover,
           std::unique_ptr<mutation::Mutation<T>> mutation,
           std::unique_ptr<stopping::StopCondition> stop_condition = std::make_unique<stopping::NoEarlyStop>());

        /**
        * Create a genetic algorithm using the specified algorithm and operators.
        *
        * @param population_size The number of candidates in the population. Must be at least 1.
        * @param algorithm The algorithm to use. Can't be a nullptr.
        * @param crossover The crossover operator to use. Can't be a nullptr.
        * @param mutation The mutation operator to use. Can't be a nullptr.
        * @param stop_condition The early-stop condition to use. Can't be a nullptr.
        */
        GA(Positive<size_t> population_size,
           std::unique_ptr<algorithm::Algorithm> algorithm,
           std::unique_ptr<crossover::Crossover<T>> crossover,
           std::unique_ptr<mutation::Mutation<T>> mutation,
           std::unique_ptr<stopping::StopCondition> stop_condition = std::make_unique<stopping::NoEarlyStop>());

        /**
        * Create a genetic algorithm using the default genetic operators.
        * The mutation probability will be deduced from the chromosome length.
        *
        * @param population_size The number of candidates in the population. Must be at least 1.
        * @param algorithm The algorithm to use.
        */
        template<typename AlgorithmType>
        requires std::derived_from<AlgorithmType, algorithm::Algorithm> && std::is_final_v<AlgorithmType>
        GA(Positive<size_t> population_size, AlgorithmType algorithm);

        /**
        * Create a genetic algorithm using the specified operators.
        * The algorithm used will be deduced from the number of objectives of the fitness function
        * (single- or multi-objective).
        *
        * @param population_size The number of candidates in the population. Must be at least 1.
        * @param crossover The crossover operator to use.
        * @param mutation The mutation operator to use.
        * @param stop_condition The early-stop condition to use.
        */
        template<typename CrossoverType, typename MutationType, typename StoppingType = stopping::NoEarlyStop>
        requires std::derived_from<CrossoverType, crossover::Crossover<T>> && std::is_final_v<CrossoverType> &&
                 std::derived_from<MutationType, mutation::Mutation<T>> && std::is_final_v<MutationType> &&
                 std::derived_from<StoppingType, stopping::StopCondition> && std::is_final_v<StoppingType>
        GA(Positive<size_t> population_size, CrossoverType crossover, MutationType mutation, StoppingType stop_condition = {});

        /**
        * Create a genetic algorithm using the specified algorithm and operators.
        *
        * @param population_size The number of candidates in the population. Must be at least 1.
        * @param algorithm The algorithm to use.
        * @param crossover The crossover operator to use.
        * @param mutation The mutation operator to use.
        * @param stop_condition The early-stop condition to use.
        */
        template<typename AlgorithmType, typename CrossoverType, typename MutationType, typename StoppingType = stopping::NoEarlyStop>
        requires std::derived_from<AlgorithmType, algorithm::Algorithm> && std::is_final_v<AlgorithmType> &&
                 std::derived_from<CrossoverType, crossover::Crossover<T>> && std::is_final_v<CrossoverType> &&
                 std::derived_from<MutationType, mutation::Mutation<T>> && std::is_final_v<MutationType> &&
                 std::derived_from<StoppingType, stopping::StopCondition> && std::is_final_v<StoppingType>
        GA(Positive<size_t> population_size, AlgorithmType algorithm, CrossoverType crossover, MutationType mutation, StoppingType stop_condition = {});


        /** @returns The chromosome length of the candidates of the population. */
        [[nodiscard]]
        size_t chrom_len() const noexcept final;

        /** @returns True if variable chromosome lengths are allowed and used. */
        [[nodiscard]]
        bool variable_chrom_len() const noexcept final;

        /** @returns True if a dynamic fitness function is used. */
        [[nodiscard]]
        bool dynamic_fitness() const noexcept final;


        /** @returns The lower and upper bounds for each of the chromosomes' genes (the ranges are inclusive). @see Bounds */
        [[nodiscard]]
        const BoundsVector<T>& gene_bounds() const noexcept requires (is_bounded<T>);

        /**
        * Set the crossover method the %GA will use.
        * The crossover method should be thread-safe if parallel execution is enabled (true by default).
        * 
        * @note As the crossover probability is associated with a crossover method, setting a new
        *   crossover method will override any previous crossover probability set using crossover_rate().
        * 
        * @param f The crossover method to use.
        */
        template<typename F>
        requires std::derived_from<F, crossover::Crossover<T>> && std::is_final_v<F>
        void crossover_method(F f);

        /**
        * Set the crossover method the %GA will use.
        * The crossover method should be thread-safe if parallel execution is enabled (true by default).
        *
        * @note As the crossover probability is associated with a crossover method, setting a new
        *   crossover method will override any previous crossover probability set using crossover_rate().
        *
        * @param f The crossover method to use. Can't be a nullptr.
        */
        void crossover_method(std::unique_ptr<crossover::Crossover<T>> f);

        /**
        * Set the crossover method the %GA will use.
        * The crossover method should be thread-safe if parallel execution is enabled (true by default).
        * 
        * @note As the crossover probability is associated with a crossover method, setting a new
        *   crossover method will override any previous crossover probability set using crossover_rate().
        *   The crossover probability that will be used for simple callables will be the default
        *   crossover rate by default, but it can be changed using crossover_rate().
        * 
        * @see CrossoverCallable
        *
        * @param f The crossover method to use.
        */
        void crossover_method(CrossoverCallable f);

        /** @returns The crossover operator used by the %GA. */
        [[nodiscard]]
        const crossover::Crossover<T>& crossover_method() const& noexcept;

        /**
        * Set the crossover rate of the crossover operator used.
        * 
        * @note This crossover rate will be set for the current crossover method, and it
        *   will be changed if a new crossover method is set.
        * 
        * @see crossover_method
        *
        * @param pc The crossover probability. Must be in the closed interval [0.0, 1.0].
        */
        void crossover_rate(Probability pc) noexcept final;

        /** @returns The crossover rate set for the crossover operator. */
        [[nodiscard]]
        Probability crossover_rate() const noexcept final;


        /**
        * Set the mutation method the %GA will use.
        * The mutation method should be thread-safe if parallel execution is enabled (true by default).
        * 
        * @note As the mutation probability is associated with a mutation method, setting a new
        *   mutation method will override any previous mutation rate set using mutation_rate().
        *
        * @param f The mutation method to use.
        */
        template<typename F>
        requires std::derived_from<F, mutation::Mutation<T>> && std::is_final_v<F>
        void mutation_method(F f);

        /**
        * Set the mutation method the %GA will use.
        * The mutation method should be thread-safe if parallel execution is enabled (true by default).
        *
        * @note As the mutation probability is associated with a mutation method, setting a new
        *   mutation method will override any previous mutation rate set using mutation_rate().
        *
        * @param f The mutation method to use. Can't be a nullptr.
        */
        void mutation_method(std::unique_ptr<mutation::Mutation<T>> f);

        /**
        * Set the mutation method the %GA will use.
        * The mutation method should be thread-safe if parallel execution is enabled (true by default).
        *
        * @note As the mutation probability is associated with a mutation method, setting a new
        *   mutation method will override any previous mutation rate set using mutation_rate().
        *   The mutation probability that will be used for simple callables will be the default
        *   (deduced) mutation rate by default, but it can be changed using mutation_rate().
        * 
        * @see MutationCallable
        *
        * @param f The mutation method to use.
        */
        void mutation_method(MutationCallable f);

        /** @returns The mutation operator used by the %GA. */
        [[nodiscard]]
        const mutation::Mutation<T>& mutation_method() const& noexcept;

        /**
        * Set the mutation rate of the mutation operator.
        * 
        * @note This mutation rate will be set for the current mutation method, and it
        *   will be changed if a new mutation method is set.
        *
        * @see mutation_method
        *
        * @param pm The mutation probability. Must be in the closed interval [0.0, 1.0].
        */
        void mutation_rate(Probability pm) noexcept final;

        /** @returns The mutation probability set for the mutation operator. */
        [[nodiscard]]
        Probability mutation_rate() const noexcept final;

        /**
        * Set a repair function for the %GA.
        * This function will be invoked for each chromosome in every generation,
        * after the chromosomes were mutated.
        * 
        * Using a repair function is optional, and no repair function is used
        * in the %GA by default. A nullptr can be used as the argument to this function
        * to disable the use of a repair function.
        * 
        * The repair function should be thread-safe if parallel execution is enabled (true by default).
        * 
        * @param f The repair function to use. Allowed to be a nullptr.
        */
        void repair_function(RepairCallable f);


        /**
        * @returns The pareto-optimal solutions found by the %GA.
        *   These are the optimal solutions of the last generation's population if
        *   keep_all_optimal_solutions() is not set, otherwise it contains every optimal solution
        *   found during the run (the solution set is updated in every generation in this case).
        */
        [[nodiscard]]
        const Candidates<T>& solutions() const& noexcept { return solutions_; }

        /**
        * @returns The current population of the algorithm. This is the entire population,
        *   which may contain non optimal solutions, unlike the population returned by solutions().
        */
        [[nodiscard]]
        const Population<T>& population() const& noexcept { return population_; }


        /********************* NON-BOUNDED ALGORITHM SOLVE METHODS *********************/

        /**
        * Find the maximum of a fitness function using the genetic algorithm.
        * 
        * Specifying an initial population is optional and it may have any number of solutions,
        * but only the first population_size() elements will be used if it's greater than the
        * population size set for the %GA. If it has less candidates than the population size,
        * the rest of the candidates will be generated using the generateCandidate() method.
        *
        * @param fitness_function The fitness function to find the maximum of. Can't be a nullptr.
        * @param generations The maximum number of generations the run can last for.
        * @param initial_population The starting population to use as the first generation of the %GA.
        * @returns The pareto-optimal solutions found (this is not the final population).
        */
        Candidates<T> solve(std::unique_ptr<FitnessFunctionBase<T>> fitness_function, size_t generations, Population<T> initial_population = {}) requires (!is_bounded<T>);

        /**
        * Find the maximum of a fitness function using the genetic algorithm.
        * 
        * Specifying an initial population is optional and it may have any number of solutions,
        * but only the first population_size() elements will be used if it's greater than the
        * population size set for the %GA. If it has less candidates than the population size,
        * the rest of the candidates will be generated using the generateCandidate() method.
        *
        * @param fitness_function The fitness function to find the maximum of. Can't be a nullptr.
        * @param initial_population The starting population to use as the first generation of the %GA.
        * @returns The pareto-optimal solutions found (this is not the final population).
        */
        Candidates<T> solve(std::unique_ptr<FitnessFunctionBase<T>> fitness_function, Population<T> initial_population = {}) requires (!is_bounded<T>);

        /**
        * Find the maximum of a fitness function using the genetic algorithm.
        * 
        * Specifying an initial population is optional and it may have any number of solutions,
        * but only the first population_size() elements will be used if it's greater than the
        * population size set for the %GA. If it has less candidates than the population size,
        * the rest of the candidates will be generated using the generateCandidate() method.
        *
        * @param fitness_function The fitness function to find the maximum of.
        * @param initial_population The starting population to use as the first generation of the %GA.
        * @returns The pareto-optimal solutions found (this is not the final population).
        */
        template<typename F>
        requires (!is_bounded<T> && std::derived_from<F, FitnessFunctionBase<T>> && std::is_final_v<F>)
        Candidates<T> solve(F fitness_function, Population<T> initial_population = {});

        /**
        * Find the maximum of a fitness function using the genetic algorithm.
        * 
        * Specifying an initial population is optional and it may have any number of solutions,
        * but only the first population_size() elements will be used if it's greater than the
        * population size set for the %GA. If it has less candidates than the population size,
        * the rest of the candidates will be generated using the generateCandidate() method.
        *
        * @param fitness_function The fitness function to find the maximum of.
        * @param generations The maximum number of generations the run can last for.
        * @param initial_population The starting population to use as the first generation of the %GA.
        * @returns The pareto-optimal solutions found (this is not the final population).
        */
        template<typename F>
        requires (!is_bounded<T> && std::derived_from<F, FitnessFunctionBase<T>> && std::is_final_v<F>)
        Candidates<T> solve(F fitness_function, size_t generations, Population<T> initial_population = {});


        /********************* BOUNDED ALGORITHM SOLVE METHODS *********************/

        /**
        * Find the maximum of a fitness function using the genetic algorithm.
        * 
        * Specifying an initial population is optional and it may have any number of solutions,
        * but only the first population_size() elements will be used if it's greater than the
        * population size set for the %GA. If it has less candidates than the population size,
        * the rest of the candidates will be generated using the generateCandidate() method.
        *
        * @param fitness_function The fitness function to find the maximum of. Can't be a nullptr.
        * @param bounds The lower and upper bounds of each of the chromosomes' genes.
        * @param generations The maximum number of generations the run can last for.
        * @param initial_population The starting population to use as the first generation of the %GA.
        * @returns The pareto-optimal solutions found (this is not the final population).
        */
        Candidates<T> solve(std::unique_ptr<FitnessFunctionBase<T>> fitness_function, BoundsVector<T> bounds, size_t generations, Population<T> initial_population = {}) requires (is_bounded<T>);

        /**
        * Find the maximum of a fitness function using the genetic algorithm.
        *
        * Specifying an initial population is optional and it may have any number of solutions,
        * but only the first population_size() elements will be used if it's greater than the
        * population size set for the %GA. If it has less candidates than the population size,
        * the rest of the candidates will be generated using the generateCandidate() method.
        *
        * @param fitness_function The fitness function to find the maximum of. Can't be a nullptr.
        * @param bounds The lower and upper bounds of a gene. The same bounds will be used for every gene of the chromosomes.
        * @param generations The maximum number of generations the run can last for.
        * @param initial_population The starting population to use as the first generation of the %GA.
        * @returns The pareto-optimal solutions found (this is not the final population).
        */
        Candidates<T> solve(std::unique_ptr<FitnessFunctionBase<T>> fitness_function, Bounds<T> bounds, size_t generations, Population<T> initial_population = {}) requires (is_bounded<T>);

        /**
        * Find the maximum of a fitness function using the genetic algorithm.
        * 
        * Specifying an initial population is optional and it may have any number of solutions,
        * but only the first population_size() elements will be used if it's greater than the
        * population size set for the %GA. If it has less candidates than the population size,
        * the rest of the candidates will be generated using the generateCandidate() method.
        *
        * @param fitness_function The fitness function to find the maximum of. Can't be a nullptr.
        * @param bounds The lower and upper bounds of each of the chromosomes' genes.
        * @param initial_population The starting population to use as the first generation of the %GA.
        * @returns The pareto-optimal solutions found (this is not the final population).
        */
        Candidates<T> solve(std::unique_ptr<FitnessFunctionBase<T>> fitness_function, BoundsVector<T> bounds, Population<T> initial_population = {}) requires (is_bounded<T>);

        /**
        * Find the maximum of a fitness function using the genetic algorithm.
        * 
        * Specifying an initial population is optional and it may have any number of solutions,
        * but only the first population_size() elements will be used if it's greater than the
        * population size set for the %GA. If it has less candidates than the population size,
        * the rest of the candidates will be generated using the generateCandidate() method.
        *
        * @param fitness_function The fitness function to find the maximum of. Can't be a nullptr.
        * @param bounds The lower and upper bounds of a gene. The same bounds will be used for every gene of the chromosomes.
        * @param initial_population The starting population to use as the first generation of the %GA.
        * @returns The pareto-optimal solutions found (this is not the final population).
        */
        Candidates<T> solve(std::unique_ptr<FitnessFunctionBase<T>> fitness_function, Bounds<T> bounds, Population<T> initial_population = {}) requires (is_bounded<T>);

        /**
        * Find the maximum of a fitness function using the genetic algorithm.
        * 
        * Specifying an initial population is optional and it may have any number of solutions,
        * but only the first population_size() elements will be used if it's greater than the
        * population size set for the %GA. If it has less candidates than the population size,
        * the rest of the candidates will be generated using the generateCandidate() method.
        *
        * @param fitness_function The fitness function to find the maximum of.
        * @param bounds The lower and upper bounds of each of the chromosomes' genes.
        * @param initial_population The starting population to use as the first generation of the %GA.
        * @returns The pareto-optimal solutions found (this is not the final population).
        */
        template<typename F>
        requires (is_bounded<T> && std::derived_from<F, FitnessFunctionBase<T>> && std::is_final_v<F>)
        Candidates<T> solve(F fitness_function, BoundsVector<T> bounds, Population<T> initial_population = {});

        /**
        * Find the maximum of a fitness function using the genetic algorithm.
        * 
        * Specifying an initial population is optional and it may have any number of solutions,
        * but only the first population_size() elements will be used if it's greater than the
        * population size set for the %GA. If it has less candidates than the population size,
        * the rest of the candidates will be generated using the generateCandidate() method.
        *
        * @param fitness_function The fitness function to find the maximum of.
        * @param bounds The lower and upper bounds of a gene. The same bounds will be used for every gene of the chromosomes.
        * @param initial_population The starting population to use as the first generation of the %GA.
        * @returns The pareto-optimal solutions found (this is not the final population).
        */
        template<typename F>
        requires (is_bounded<T> && std::derived_from<F, FitnessFunctionBase<T>> && std::is_final_v<F>)
        Candidates<T> solve(F fitness_function, Bounds<T> bounds, Population<T> initial_population = {});

        /**
        * Find the maximum of a fitness function using the genetic algorithm.
        * 
        * Specifying an initial population is optional and it may have any number of solutions,
        * but only the first population_size() elements will be used if it's greater than the
        * population size set for the %GA. If it has less candidates than the population size,
        * the rest of the candidates will be generated using the generateCandidate() method.
        *
        * @param fitness_function The fitness function to find the maximum of.
        * @param bounds The lower and upper bounds of each of the chromosomes' genes.
        * @param generations The maximum number of generations the run can last for.
        * @param initial_population The starting population to use as the first generation of the %GA.
        * @returns The pareto-optimal solutions found (this is not the final population).
        */
        template<typename F>
        requires (is_bounded<T> && std::derived_from<F, FitnessFunctionBase<T>> && std::is_final_v<F>)
        Candidates<T> solve(F fitness_function, BoundsVector<T> bounds, size_t generations, Population<T> initial_population = {});

        /**
        * Find the maximum of a fitness function using the genetic algorithm.
        * 
        * Specifying an initial population is optional and it may have any number of solutions,
        * but only the first population_size() elements will be used if it's greater than the
        * population size set for the %GA. If it has less candidates than the population size,
        * the rest of the candidates will be generated using the generateCandidate() method.
        *
        * @param fitness_function The fitness function to find the maximum of.
        * @param bounds The lower and upper bounds of a gene. The same bounds will be used for every gene of the chromosomes.
        * @param generations The maximum number of generations the run can last for.
        * @param initial_population The starting population to use as the first generation of the %GA.
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

        /**
        * Initialize the derived genetic algorithm. This method will be called exactly once
        * at the start of each run.
        * 
        * The implementation of this method should only contain the initialization logic
        * specific to the derived class, as the initialization of the %GA class is separate.
        * 
        * Implementing this in the derived classes is optional,
        * the default implementation does nothing.
        */
        virtual void initialize() {};

        /**
        * Generate a candidate solution. This method will be used to generate the initial
        * population of the %GA if the initial population is not fully specified in solve().
        * 
        * Has to be implemented in the derived classes, and every generated candidate must be
        * valid (eg. the chromosome sizes must be correct, its genes must be within the specified
        * bounds etc.).
        */
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