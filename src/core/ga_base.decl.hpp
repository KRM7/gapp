/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CORE_GA_BASE_DECL_HPP
#define GA_CORE_GA_BASE_DECL_HPP

#include "ga_info.hpp"
#include "../population/candidate.hpp"
#include "../mutation/mutation_base.fwd.hpp"
#include "../crossover/crossover_base.fwd.hpp"
#include <algorithm>
#include <vector>
#include <utility>
#include <functional>
#include <type_traits>
#include <atomic>
#include <cstddef>
#include <memory>
#include <concepts>

namespace genetic_algorithm
{
    /**
    * Base GA class. \n
    * Contains everything of the genetic algorithm, except for the
    * generateCandidate() function, which should be implemented in the derived classes.
    *
    * @tparam Gene The type of the genes in the candidates' chromosomes.
    * @tparam Derived The type of the derived class.
    */
    template<Gene T, typename Derived>
    class GA : public GaInfo
    {
    public:

        using GeneType = T;                             /**< The gene type used in the chromosomes. */
        using Candidate = Candidate<GeneType>;          /**< The type used for the candidates in the algorithm. */
        using Chromosome = std::vector<GeneType>;       /**< The type of the chromosomes of the Candidates, representing a solution. */
        using CandidatePair = CandidatePair<GeneType>;  /**< A pair of Candidates. */
        using Candidates = std::vector<Candidate>;      /**< A vector of Candidates. */
        using Population = std::vector<Candidate>;      /**< The type of the population of the algorithm. */

        /**
        * The type of the fitness function. \n
        * Takes a vector of genes (chromosome) and returns the fitness vector of the chromosome. \n
        * The returned fitness vector should only contain finite values, and always be the same length
        * during a run. \n
        * The fitness function is allowed to return different fitness vectors for the
        * same chromosome if dynamic_fitness is set to true.
        */
        using FitnessFunction = std::function<std::vector<double>(const Chromosome&)>;

        /**
        * The type of the crossover function. \n
        * Takes two candidate solutions (parents), and returns two candidates generated from these
        * parent solutions (children). \n
        * The returned children's chromosomes don't have to be the same length as the parents, but
        * all of the other operators used should be able to handle different chromosome lengths in this case
        * (fitness function, mutation, repair).
        */
        using CrossoverFunction = std::function<CandidatePair(const GaInfo&, const Candidate&, const Candidate&)>;

        /**
        * The type of the mutation function. \n
        * Takes a candidate solution, and mutates this solution. \n
        * The function is allowed to change the candidate's chromosome's length, but
        * all of the other operators used should be able to handle different chromosome lengths in this
        * case (fitness function, crossover, repair).
        */
        using MutationFunction = std::function<void(const GaInfo&, Candidate&)>;

        /**
        * The type of the repair function. \n
        * Takes a candidate solution, and performs some operation on it. \n
        * This function is allowed the change the candidate's chromosome's length, but
        * all of the other operators used should be able to handle different chromosome lengths in this
        * case (fitness function, crossover, mutation).
        */
        using RepairFunction = std::function<Chromosome(const Chromosome&)>;

        /** The type of the callback functions used in the algorithm. */
        using CallbackFunction = std::function<void(const GA<T, Derived>&)>;

        /**
        * Create a genetic algorithm.
        *
        * @param chrom_len The length of the chromosomes (number of genes).
        * @param fitness_function The function to find the maximum of. @see FitnessFunction
        */
        GA(size_t chrom_len, FitnessFunction fitness_function);

        /**
        * Create a genetic algorithm.
        *
        * @param population_size The number of candidates in the population.
        * @param chrom_len The length of the chromosomes (number of genes).
        * @param fitness_function The function to find the maximum of. @see FitnessFunction
        */
        GA(size_t population_size, size_t chrom_len, FitnessFunction fitness_function);

        /**
        * Set the fitness function used by the algorithm. \n
        * The fitness function should return a vector with a size equal to the number of objectives. \n
        * The fitness function should be thread-safe if parallel execution is enabled (it is enabled by default).
        * @see FitnessFunction
        *
        * @param f The function the algorithm should find the maximum of.
        */
        void fitness_function(FitnessFunction f);

        /**
        * Set the crossover method the algorithm will use. \n
        * The crossover function should be thread-safe if parallel execution is enabled (it is enabled by default).
        * @see CrossoverFunction
        * @see Crossover
        * 
        * @param f The crossover method used by the algorithm.
        */
        template<typename F>
        requires crossover::CrossoverType<F, T> && std::is_final_v<F>
        void crossover_method(F&& f);

        /**
        * Set the crossover method the algorithm will use. \n
        * The crossover function should be thread-safe if parallel execution is enabled (it is enabled by default).
        * @see CrossoverFunction
        * @see Crossover
        *
        * @param f The crossover method used by the algorithm.
        */
        template<crossover::CrossoverType<T> F>
        void crossover_method(std::unique_ptr<F>&& f);

        /**
        * Set the crossover method the algorithm will use. \n
        * The crossover function should be thread-safe if parallel execution is enabled (it is enabled by default).
        * @see CrossoverFunction
        * @see Crossover
        *
        * @param f The crossover method used by the algorithm.
        */
        void crossover_method(CrossoverFunction f);

        /** @returns The crossover operator used by the algorithm, cast to type @p F. */
        template<crossover::CrossoverType<T> F = crossover::Crossover<GeneType>>
        F& crossover_method() &;

        /**
        * Set the crossover rate of the crossover operator used by the algorithm.
        *
        * @param pc The crossover probability. Must be in the closed interval [0.0, 1.0].
        */
        void crossover_rate(double pc) final;

        /** @returns The crossover rate set for the crossover operator. */
        [[nodiscard]]
        double crossover_rate() const noexcept final;

        /**
        * Set the mutation method the algorithm will use. \n
        * The mutation function should be thread-safe if parallel execution is enabled (it is enabled by default).
        * @see MutationFunction
        * @see Mutation
        *
        * @param f The crossover function that will be used by the algorithm.
        */
        template<typename F>
        requires mutation::MutationType<F, T>&& std::is_final_v<F>
        void mutation_method(F&& f);

        /**
        * Set the mutation method the algorithm will use. \n
        * The mutation function should be thread-safe if parallel execution is enabled (it is enabled by default).
        * @see MutationFunction
        * @see Mutation
        *
        * @param f The crossover function that will be used by the algorithm.
        */
        template<mutation::MutationType<T> F>
        void mutation_method(std::unique_ptr<F>&& f);

        /**
        * Set the mutation method the algorithm will use. \n
        * The mutation function should be thread-safe if parallel execution is enabled (it is enabled by default).
        * @see MutationFunction
        * @see Mutation
        *
        * @param f The crossover function that will be used by the algorithm.
        */
        void mutation_method(MutationFunction f);

        /** @returns The mutation operator used by the algorithm, cast to type @p F. */
        template<mutation::MutationType<T> F = mutation::Mutation<GeneType>>
        F& mutation_method() &;

        /**
        * Set the mutation rate of the mutation operator.
        *
        * @param pm The mutation probability. Must be in the closed interval [0.0, 1.0].
        */
        void mutation_rate(double pm) final;

        /** @returns The mutation rate of the mutation operator. */
        [[nodiscard]]
        double mutation_rate() const noexcept final;

        /**
        * Set a repair function for the algorithm. \n
        * This function will be applied to each chromosome after the mutations have been performed. \n
        * This can be used to perform a local search for example. \n
        * 
        * No repair function is used in the algorithm by default, set to nullptr if no repair function should be used.
        * 
        * The repair function should be thread-safe if parallel execution is enabled (it is enabled by default).
        * 
        * @param f The repair function the algorithm will use.
        */
        void repair_function(const RepairFunction& f);

        /** 
        * @returns The pareto optimal solutions found by the algorithm.
        * These are just the optimal solutions of the last generation's population if
        * keep_all_optimal_sols is not set, otherwise it contains every optimal solution
        * found during the run, updated in every generation.
        */
        [[nodiscard]]
        const Candidates& solutions() const noexcept { return solutions_; }

        /**
        * @returns The current population of the algorithm.
        * This is not the same as the solutions (this will always the same size as the population_size set,
        * and may also contain non pareto-optimal solutions).
        */
        [[nodiscard]]
        const Population& population() const noexcept { return population_; }

        /**
        * Run the genetic algorithm for a set number of generations. \n
        * The algorithm will always stop when reaching this number of generations.
        *
        * @param num_generations The maximum number of generations.
        * @returns The pareto-optimal solutions found. @see solutions
        */
        Candidates run(size_t num_generations = 1000);

        /**
        * Continue running the genetic algorithm for the set number of generations. \n
        * The algorithm can only continue if run has already been called at least once,
        * and the selection_method hasn't been changed since then.
        * (Otherwise this function will behave the same as calling run().)
        *
        * @param num_generations The number of generations to run the algorithm for.
        * @returns The optimal solutions.
        */
        Candidates continueFor(size_t num_generations);

        /** This function will be called exactly once at the end of each generation. */
        CallbackFunction endOfGenerationCallback = nullptr;

        /** This function will be called exactly once at the end of a run. */
        CallbackFunction endOfRunCallback = nullptr;

    protected:

        void setDefaultAlgorithm();

    private:

        Population population_;
        Candidates solutions_;

        FitnessFunction fitness_function_;
        std::unique_ptr<crossover::Crossover<GeneType>> crossover_;
        std::unique_ptr<mutation::Mutation<GeneType>> mutation_;
        RepairFunction repair_ = nullptr;

        void initializeAlgorithm();
        [[nodiscard]] size_t findNumObjectives(const FitnessFunction& f) const;
        [[nodiscard]] Candidate generateCandidate() const;
        [[nodiscard]] Population generatePopulation(size_t pop_size) const;
        void prepareSelections() const;
        [[nodiscard]] const Candidate& select() const;
        [[nodiscard]] CandidatePair crossover(const Candidate& parent1, const Candidate& parent2) const;
        void mutate(Candidate& sol) const;
        void repair(Candidate& sol) const;
        void updatePopulation(Population& pop, Population&& children);
        [[nodiscard]] bool stopCondition() const;

        void evaluateSolution(Candidate& sol);
        [[nodiscard]] FitnessMatrix evaluatePopulation(Population& pop);
        void updateOptimalSolutions(Candidates& optimal_sols, const Population& pop) const;

        void advance();

        Derived& derived() noexcept;
        const Derived& derived() const noexcept;

        bool fitnessMatrixIsValid() const noexcept;

        /* Make the protected fields of GaInfo private. */
        using GaInfo::fitness_matrix_;
        using GaInfo::algorithm_;
        using GaInfo::stop_condition_;
        using GaInfo::num_fitness_evals_;
        using GaInfo::generation_cntr_;
        using GaInfo::num_objectives_;
        using GaInfo::chrom_len_;
        using GaInfo::population_size_;
        using GaInfo::max_gen_;
        using GaInfo::can_continue_;
    };

    /** Genetic algorithm types. */
    template<typename T>
    concept GeneticAlgorithmType = detail::DerivedFromSpecializationOf<T, GA>;

} // namespace genetic_algorithm

#endif // !GA_CORE_GA_BASE_DECL_HPP