/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CORE_GA_INFO_HPP
#define GA_CORE_GA_INFO_HPP

#include "../population/population.hpp"
#include "../utility/bounded_value.hpp"
#include "../utility/utility.hpp"
#include <functional>
#include <type_traits>
#include <concepts>
#include <memory>
#include <cstddef>

namespace genetic_algorithm::algorithm
{
    class Algorithm;

} // namespace genetic_algorithm::algorithm

namespace genetic_algorithm::stopping
{
    class StopCondition;

} // namespace genetic_algorithm::stopping


namespace genetic_algorithm
{
    /**
    * Base class that all GAs are derived from. \n
    * Contains all of the general information about a GA that does not depend on the encoding (gene) type. \n
    * Move-only.
    */
    class GaInfo
    {
    public:
        /**
        * Type returned by fitness_matrix(), used to represent the fitness matrix of the population. \n
        * Each element of the matrix is the fitness vector of the corresponding solution of the population. \n
        * Eg. fmat[0] is the fitness vector of the first member of the population.
        */
        using FitnessMatrix = detail::FitnessMatrix;

        /**
        * The type of the candidate solutions' fitness vectors. \n
        * Contains fitness values along each objective axis.
        */
        using FitnessVector = detail::FitnessVector;

        /**
        * The general callable type that can be used as a stop condition in the algorithm. \n
        * It should return true if the algorithm should stop.
        */
        using StopConditionCallable = std::function<bool(const GaInfo&)>;

        /** The type of the callback functions that can be used in the algorithm. */
        using GaCallback = std::function<void(const GaInfo&)>;

        /** 
        * Create a genetic algorithm.
        *
        * @param population_size The size of the population. Must be at least 1.
        * @param algorithm The algorithm to use. Should be consistent with the objective function (single- or multi-objective).
        * @param stop_condition The stop condition to use for the algorithm.
        */
        GaInfo(Positive<size_t> population_size, std::unique_ptr<algorithm::Algorithm> algorithm, std::unique_ptr<stopping::StopCondition> stop_condition) noexcept;

        /**
        * Set the number of candidates used in the population.
        *
        * @param size The number of candidates in a population. Must be at least 1.
        */
        void population_size(Positive<size_t> size) noexcept { population_size_ = size; }

        /** @returns The number of candidates in the population. */
        [[nodiscard]]
        size_t population_size() const noexcept { return population_size_; }

        /**
        * Set the maximum number of generations a run can last for. \n
        * The algorthm will always stop when reaching this generation, even if another stop condition was set.
        * 
        * @param max_gen The maximum number of generations. Must be at least 1.
        */
        void max_gen(Positive<size_t> max_gen) noexcept { max_gen_ = max_gen; }

        /** @returns The maximum number of generations set for the algorithm. */
        [[nodiscard]]
        size_t max_gen() const noexcept { return max_gen_; }

        /** @returns The chromosome length used for the candidates of the population. */
        [[nodiscard]]
        virtual size_t chrom_len() const noexcept = 0;

        /** @returns True if variable chromosome lengths are allowed and used. */
        [[nodiscard]]
        virtual bool variable_chrom_len() const noexcept = 0;

        /** @returns The number of objectives of the fitness function. */
        [[nodiscard]]
        virtual size_t num_objectives() const noexcept = 0;

        /** @returns True if a dynamic fitness function is used. */
        [[nodiscard]]
        virtual bool dynamic_fitness() const noexcept = 0;

        /**
        * @returns The fitness matrix of the population.
        * Each element of the matrix is the fitness vector of the corresponding solution in the population. \n
        * Eg. fmat[0] is the fitness vector of the first member of the population.
        */
        [[nodiscard]]
        const FitnessMatrix& fitness_matrix() const noexcept { return fitness_matrix_; }

        /** @returns The number of fitness evaluations performed during the run. Updated after every objective function evaluation. */
        [[nodiscard]]
        size_t num_fitness_evals() const noexcept;

        /** @returns The current generation's number. */
        [[nodiscard]]
        size_t generation_cntr() const noexcept { return generation_cntr_; }
        

        /**
        * Set the crossover rate of the crossover operator used by the algorithm.
        *
        * @param pc The crossover probability. Must be in the closed interval [0.0, 1.0].
        */
        virtual void crossover_rate(Probability pc) noexcept = 0;

        /** @returns The crossover rate set for the crossover operator. */
        [[nodiscard]]
        virtual Probability crossover_rate() const noexcept = 0;

        /**
        * Set the mutation rate of the mutation operator used by the algorithm.
        *
        * @param pm The mutation probability. Must be in the closed interval [0.0, 1.0].
        */
        virtual void mutation_rate(Probability pm) noexcept = 0;

        /** @returns The mutation rate set for the mutation operator. */
        [[nodiscard]]
        virtual Probability mutation_rate() const noexcept = 0;


        /**
        * Set the algorithm used by the GA. \n
        * The algorithm type should be consistent with the size of the fitness vectors returned by
        * the fitness function (single- or multi-objective).
        * 
        * @see Algorithm
        *
        * @param f The algorithm used by the GA.
        */
        template<typename F>
        requires std::derived_from<F, algorithm::Algorithm> && std::is_final_v<F>
        void algorithm(F f);

        /**
        * Set the algorithm used by the GA. \n
        * The algorithm type should be consistent with the size of the fitness vectors returned by
        * the fitness function (single- or multi-objective).
        * 
        * @see Algorithm
        *
        * @param f The algorithm used by the GA. Can't be a nullptr.
        */
        void algorithm(std::unique_ptr<algorithm::Algorithm> f);

        /** @returns The algorithm used by the GA. */
        [[nodiscard]]
        const algorithm::Algorithm& algorithm() const& noexcept { GA_ASSERT(algorithm_); return *algorithm_; }


        /**
        * Set an early-stop condition for the genetic algorithm. \n
        * The algorithm will always stop when reaching the maximum generations set,
        * regardless of the stop condition set here.
        * 
        * @see StopCondition
        *
        * @param f The StopCondition the algorithm should use.
        */
        template<typename F>
        requires std::derived_from<F, stopping::StopCondition> && std::is_final_v<F>
        void stop_condition(F f);

        /**
        * Set an early-stop condition for the genetic algorithm. \n
        * The algorithm will always stop when reaching the maximum generations set, regardless of the stop
        * condition set here.
        * 
        * @see StopCondition
        *
        * @param f The StopCondition the algorithm should use. Can't be a nullptr.
        */
        void stop_condition(std::unique_ptr<stopping::StopCondition> f);

        /**
        * Set an early-stop condition for the genetic algorithm. \n
        * The algorithm will always stop when reaching the maximum generations set,
        * regardless of the stop condition set here.
        * 
        * @see StopCondition
        * @see StopConditionCallable
        *
        * @param f The function used to check for the early-stop condition.
        */
        void stop_condition(StopConditionCallable f);

        /** @returns The stop condition used by the algorithm. */
        [[nodiscard]]
        const stopping::StopCondition& stop_condition() const& noexcept { GA_ASSERT(stop_condition_); return *stop_condition_; }


        /**
        * When set to true, all pareto optimal Candidates found by the algorithm during a run
        * will be kept and stored in the solutions set, regardless of which generation they were
        * found in. \n
        * When set to false, the optimal solutions returned by the algorithm at the end of a run
        * are the optimal solutions of the last generation. \n
        *
        * Disabled by default. \n
        *
        * @param enable Whether all pareto optimal solutions should be kept.
        */
        void keep_all_optimal_solutions(bool enable) noexcept { keep_all_optimal_sols_ = enable; }

        /** @returns true if all pareto optimal solutions are kept during a run. */
        [[nodiscard]]
        bool keep_all_optimal_solutions() const noexcept { return keep_all_optimal_sols_; }

        /** This function will be called once at the end of each generation. */
        GaCallback endOfGenerationCallback = nullptr;


        GaInfo(const GaInfo&)            = delete;
        GaInfo& operator=(const GaInfo&) = delete;

        /** Destructor. */
        virtual ~GaInfo();

    protected:

        GaInfo(GaInfo&&) noexcept;
        GaInfo& operator=(GaInfo&&) noexcept;

        FitnessMatrix fitness_matrix_;

        std::unique_ptr<algorithm::Algorithm> algorithm_;
        std::unique_ptr<stopping::StopCondition> stop_condition_;

        Positive<size_t> population_size_ = DEFAULT_POPSIZE;
        Positive<size_t> max_gen_ = 500;
        size_t generation_cntr_ = 0;
        size_t num_fitness_evals_ = 0;

        bool keep_all_optimal_sols_ = false;
        bool use_default_algorithm_ = false;

        static constexpr size_t DEFAULT_POPSIZE = 100;
    };

} // namespace genetic_algorithm


/* IMPLEMENTATION */

#include <utility>

namespace genetic_algorithm
{
    template<typename F>
    requires std::derived_from<F, algorithm::Algorithm> && std::is_final_v<F>
    inline void GaInfo::algorithm(F f)
    {
        algorithm_ = std::make_unique<F>(std::move(f));
        use_default_algorithm_ = false;
    }

    template<typename F>
    requires std::derived_from<F, stopping::StopCondition> && std::is_final_v<F>
    inline void GaInfo::stop_condition(F f)
    {
        stop_condition_ = std::make_unique<F>(std::move(f));
    }

} // namespace genetic_algorithm

#endif // !GA_CORE_GA_INFO_HPP