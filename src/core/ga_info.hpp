/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_CORE_GA_INFO_HPP
#define GAPP_CORE_GA_INFO_HPP

#include "population.hpp"
#include "fitness_function.hpp"
#include "../metrics/metric_set.hpp"
#include "../utility/small_vector.hpp"
#include "../utility/bounded_value.hpp"
#include "../utility/type_id.hpp"
#include "../utility/utility.hpp"
#include <functional>
#include <type_traits>
#include <concepts>
#include <memory>
#include <utility>
#include <cstddef>

namespace gapp::algorithm
{
    class Algorithm;

} // namespace gapp::algorithm

namespace gapp::stopping
{
    class StopCondition;

} // namespace gapp::stopping

namespace gapp::crossover
{
    template<typename T>
    class Crossover;

} // namespace gapp::crossover

namespace gapp::mutation
{
    template<typename T>
    class Mutation;

} // namespace gapp::mutation

namespace gapp
{
    /**
    * The base class that all GAs are derived from.
    * Contains all of the general properties of a genetic algorithm that 
    * do not depend on the encoding type.
    * 
    * %GA implementations should use the GA class as their base class instead
    * of inheriting from this class directly.
    * 
    * Move-only.
    */
    class GaInfo
    {
    public:
        /**
        * The general callable type that can be used as a stop condition in the algorithm (when
        * not using a stop condition that is derived from StopCondition).
        * The function should return true when the algorithm should be stopped.
        * @see stop_condition
        */
        using StopConditionCallable = std::function<bool(const GaInfo&)>;

        /**
        * The type of a generic callback function that can be provided for the algorithm.
        * @see on_generation_end
        */
        using GaInfoCallback = std::function<void(const GaInfo&)>;

        /** 
        * Create a genetic algorithm.
        *
        * @param population_size The size of the population. Must be at least 1.
        * @param algorithm The algorithm to use. Should be consistent with the objective function (ie. a single-objective
        *       algorithm should be used for single-objective fitness functions, and a multi-objective algorithm should be
        *       used for multi-objective problems).
        * @param stop_condition The early-stop condition to use for the algorithm. The algorithm will stop when
        *       reaching the maximum number of generations set regardless of the stop condition specified here.
        */
        GaInfo(Positive<size_t> population_size, std::unique_ptr<algorithm::Algorithm> algorithm, std::unique_ptr<stopping::StopCondition> stop_condition);


        /**
        * Set the number of candidate solutions used in the population.
        *
        * @param size The number of candidates in the population. Must be at least 1.
        */
        void population_size(Positive<size_t> size) noexcept { population_size_ = size; }

        /** @returns The number of candidate solutions in the population. */
        [[nodiscard]]
        size_t population_size() const noexcept { return population_size_; }

        /**
        * Set the maximum number of generations a run can last for.
        * The algorthm will always stop when reaching this generation, even if there is
        * another stop-condition being used.
        * 
        * @see stop_condition
        * 
        * @param max_gen The maximum number of generations. Must be at least 1.
        */
        void max_gen(Positive<size_t> max_gen) noexcept { max_gen_ = max_gen; }

        /** @returns The maximum number of generations set for the algorithm. */
        [[nodiscard]]
        size_t max_gen() const noexcept { return max_gen_; }


        /** @returns The fitness function used. A nullptr is returned if no fitness function is set. */
        [[nodiscard]]
        virtual const FitnessFunctionInfo* fitness_function() const& noexcept;

        /**
        * @returns The chromosome length used for each chromosome of the encoding.
        * The size of the returned vector will be equal to the number of gene types in the encoding.
        * This means a size of 1 for simple encodings, and a size equal to the number of component
        * genes for mixed encodings.
        * The order of the values in the vector will match the order of the genes as specified in
        * the encoding.
        * An empty vector is returned if no fitness function is set.
        */
        [[nodiscard]]
        const small_vector<size_t>& chrom_lens() const noexcept;

        /**
        * @returns The chromosome length used for the specified gene type of the encoding.
        * Returns 0 if the candidates don't have a chromosome for the specified gene type or
        * if no fitness function is set.
        */
        template<typename GeneType>
        [[nodiscard]] size_t chrom_len() const noexcept;

        /** @returns The number of objectives of the fitness function, or 0 if no fitness function is set. */
        [[nodiscard]]
        size_t num_objectives() const noexcept { return num_objectives_; }

        /** @returns The number of constraints associated with the fitness function, or 0 if no fitness function is set. */
        [[nodiscard]]
        size_t num_constraints() const noexcept { return num_constraints_; }


        /**
        * @returns A view of the population, without the encoding dependent parts of the candidates
        * (i.e. without the chromosomes). Each element of this view is a reference to an actual
        * candidate of the population, not a separate data structure. They may be cast to the concrete
        * candidate type if the encoding type is known.
        */
        [[nodiscard]]
        virtual PopulationView population_view() const& noexcept = 0;

        /**
        * @returns The fitness matrix of the population.
        * Each row of the matrix is the fitness vector of the corresponding solution in the population. \n
        * For example, fmat[0] is the fitness vector of the first member of the population.
        */
        [[nodiscard]]
        const FitnessMatrix& fitness_matrix() const& noexcept { return fitness_matrix_; }

        /**
        * @returns The number of fitness evaluations performed during the run so far.
        * This value is updated after every objective function evaluation.
        */
        [[nodiscard]]
        size_t num_fitness_evals() const noexcept;

        /** 
        * @returns The current generation's number. This value will be in the range [0, max_gen),
        * where 0 corresponds to the initial/first generation.
        */
        [[nodiscard]]
        size_t generation_cntr() const noexcept { return generation_cntr_; }


        /**
        * @returns The crossover method associated with the specified gene type.
        * The gene type may either be the mixed gene or one of the component genes for mixed
        * encodings. A nullptr is returned if the specified gene type is not part of the
        * encoding.
        */
        template<typename GeneType>
        [[nodiscard]] crossover::Crossover<GeneType>* crossover_method() const& noexcept;

        /**
        * @returns The mutation method associated with the specified gene type.
        * The gene type may either be the mixed gene or one of the component genes for mixed
        * encodings. A nullptr is returned if the specified gene type is not part of the
        * encoding.
        */
        template<typename GeneType>
        [[nodiscard]] mutation::Mutation<GeneType>* mutation_method() const& noexcept;


        /**
        * Set the algorithm used by the %GA. \n
        * The algorithm type should be consistent with the fitness function (ie. it should
        * be a single-objective algorithm for single-objective problems, and a multi-objective
        * algorithm for multi-objective problems).
        *
        * @param f The algorithm used by the %GA.
        */
        template<typename F>
        requires std::derived_from<F, algorithm::Algorithm>
        void algorithm(F f);

        /**
        * Set the algorithm used by the %GA. \n
        * The algorithm type should be consistent with the fitness function (ie. it should
        * be a single-objective algorithm for single-objective problems, and a multi-objective
        * algorithm for multi-objective problems).
        *
        * @param f The algorithm used by the %GA. The default algorithm will be used if it's a
        *   nullptr.
        */
        void algorithm(std::unique_ptr<algorithm::Algorithm> f);

        /**
        * Clear the algorithm currently set for the GA. \n
        * The GA will use the default algorithm that is selected based on the number of
        * objectives of the fitness functions.
        */
        void algorithm(std::nullptr_t);

        /** @returns The algorithm used by the %GA. */
        [[nodiscard]]
        const algorithm::Algorithm& algorithm() const& noexcept { GAPP_ASSERT(algorithm_); return *algorithm_; }


        /**
        * Set an early-stop condition for the algorithm. \n
        * This is an optional early-stop condition that can be used to stop the run
        * before the maximum number of generations is reached, but the algorithm will
        * always stop when reaching the maximum generations regardless of the stop condition set here.
        *
        * @param f The early stop condition the %GA should use.
        */
        template<typename F>
        requires std::derived_from<F, stopping::StopCondition>
        void stop_condition(F f);

        /**
        * Set an early-stop condition for the algorithm. \n
        * This is an optional early-stop condition that can be used to stop the run
        * before the maximum number of generations is reached, but the algorithm will
        * always stop when reaching the maximum generations regardless of the stop condition set here.
        *
        * @param f The stop condition the %GA should use. No early-stopping will be used if it's a nullptr.
        */
        void stop_condition(std::unique_ptr<stopping::StopCondition> f);

        /**
        * Clear the early-stop condition currently set for the GA. \n
        * The GA will run for the maximum number of generations set without
        * the possibility of stopping earlier.
        */
        void stop_condition(std::nullptr_t);

        /**
        * Set an early-stop condition for the algorithm. \n
        * This is an optional early-stop condition that can be used to stop the run
        * before the maximum number of generations is reached, but the algorithm will
        * always stop when reaching the maximum generations regardless of the stop condition set here.
        * 
        * @see StopConditionCallable
        *
        * @param f The function used to check for the early-stop condition.
        */
        void stop_condition(StopConditionCallable f);

        /** @returns The stop condition used by the %GA. */
        [[nodiscard]]
        const stopping::StopCondition& stop_condition() const& noexcept { GAPP_ASSERT(stop_condition_); return *stop_condition_; }


        /**
        * Set a number of metrics to track throughout the run of the algorithm.
        * The values of these metrics will be computed and saved for each generation of the run.
        * 
        * Example:
        * ```
        *   GA.track(metrics::FitnessMin{}, metrics::FitnessMax{});
        *   GA.solve(...);
        *   const auto& fmin = GA.get_metric<metrics::FitnessMin>();
        * ```
        * 
        * The same metric type may be present multiple times in the parameters, but it is
        * unspecified which one get_metric() will return in this case.
        *
        * @warning
        * Setting a new set of metrics to be tracked will invalidate all references returned by get_metric().
        * 
        * @param metrics The metrics to track during the runs of the algorithm.
        */
        template<typename... Metrics>
        requires (std::derived_from<Metrics, metrics::MonitorBase> && ...)
        void track(Metrics... metrics);

        /**
        * Get one of the metrics tracked by the algorithm.
        * 
        * Example:
        * ```
        *   GA.track(metrics::FitnessMin{}, metrics::FitnessMax{});
        *   GA.solve(...);
        *   const auto& fmin = GA.get_metric<metrics::FitnessMin>();
        * ```
        * 
        * If the set of tracked metrics contains the same metric multiple times,
        * one of them will be returned, but it is unspecified which one.
        * 
        * @warning
        * The type parameter Metric must be the type of one of the metrics being tracked, otherwise
        * the behaviour of this method is undefined.
        * 
        * @warning
        * Setting a new set of metrics to be tracked via track() will invalidate all references
        * returned by this function.
        * 
        * @tparam Metric The tracked metric to get.
        * @returns The metric of the specified type if it is tracked, otherwise undefined.
        */
        template<typename Metric>
        requires std::derived_from<Metric, metrics::MonitorBase>
        const Metric& get_metric() const noexcept;

        /**
        * Get one of the metrics tracked by the algorithm.
        *
        * This function is similar to the get_metric() method, but it returns a pointer
        * to the metric, which will be a nullptr if the given metric is not being tracked
        * by the %GA.
        *
        * @tparam Metric The tracked metric to get.
        * @returns The metric of the specified type if it is tracked, or nullptr otherwise.
        */
        template<typename Metric>
        requires std::derived_from<Metric, metrics::MonitorBase>
        const Metric* get_metric_if() const noexcept;


        /**
        * When set to true, all pareto-optimal candidates found by the %GA during a run
        * will be kept and stored in the solutions set, regardless of which generation they were
        * found in. \n
        * When set to false, the optimal solutions returned by the %GA at the end of a run
        * are the optimal solutions of the last generation only, optimal solutions found
        * in earlier generations are discarded.
        *
        * Disabled by default. \n
        * 
        * @note This will most likely have no effect for single-objective problems, as they likely
        *   only have a single optimal solution.
        * 
        * @warning Setting this to True can significantly increase the run-time of the GA.
        *
        * @param enable Whether all pareto optimal solutions should be kept.
        */
        void keep_all_optimal_solutions(bool enable) noexcept { keep_all_optimal_sols_ = enable; }

        /** @returns True if all pareto-optimal solutions are kept during a run. */
        [[nodiscard]]
        bool keep_all_optimal_solutions() const noexcept { return keep_all_optimal_sols_; }

        /**
        * Set a generic callback function that will be called exactly once at the end
        * of each generation of a run.
        * 
        * @param f The function that will be invoked at the end of each generation.
        */
        void on_generation_end(GaInfoCallback f) noexcept { on_generation_end_ = std::move(f); }


        GaInfo(const GaInfo&)            = delete;
        GaInfo& operator=(const GaInfo&) = delete;

        /** Destructor. */
        virtual ~GaInfo() noexcept;

    private:

        GaInfo(GaInfo&&) noexcept;
        GaInfo& operator=(GaInfo&&) noexcept;

        template<typename T>
        friend class GA;

        FitnessMatrix fitness_matrix_;

        std::unique_ptr<FitnessFunctionInfo> fitness_function_;
        std::unique_ptr<algorithm::Algorithm> algorithm_;
        std::unique_ptr<stopping::StopCondition> stop_condition_;
        detail::MetricSet metrics_;
        GaInfoCallback on_generation_end_ = nullptr;

        Positive<size_t> population_size_ = DEFAULT_POPSIZE;
        Positive<size_t> max_gen_ = 500;
        size_t num_objectives_ = 0;
        size_t num_constraints_ = 0;
        size_t generation_cntr_ = 0;
        size_t num_fitness_evals_ = 0;

        bool keep_all_optimal_sols_ = false;

        /** The default population size used in the %GA if none is specified. */
        static constexpr size_t DEFAULT_POPSIZE = 100;

        virtual size_t index_of_gene(size_t type_id) const noexcept = 0;

        virtual void* crossover_method_impl(size_t type_id) const noexcept = 0;
        virtual void* mutation_method_impl(size_t type_id) const noexcept = 0;
    };

} // namespace gapp

namespace gapp
{
    template<typename GeneType>
    size_t GaInfo::chrom_len() const noexcept
    {
        if (!fitness_function_) return 0;

        return chrom_lens()[index_of_gene(detail::type_id<GeneType>())];
    }

    template<typename GeneType>
    crossover::Crossover<GeneType>* GaInfo::crossover_method() const& noexcept
    {
        return static_cast<crossover::Crossover<GeneType>*>(crossover_method_impl(detail::type_id<GeneType>()));
    }

    template<typename GeneType>
    mutation::Mutation<GeneType>* GaInfo::mutation_method() const& noexcept
    {
        return static_cast<mutation::Mutation<GeneType>*>(mutation_method_impl(detail::type_id<GeneType>()));
    }

    template<typename F>
    requires std::derived_from<F, algorithm::Algorithm>
    inline void GaInfo::algorithm(F f)
    {
        algorithm_ = std::make_unique<F>(std::move(f));
    }

    template<typename F>
    requires std::derived_from<F, stopping::StopCondition>
    inline void GaInfo::stop_condition(F f)
    {
        stop_condition_ = std::make_unique<F>(std::move(f));
    }

    template<typename... Metrics>
    requires (std::derived_from<Metrics, metrics::MonitorBase> && ...)
    inline void GaInfo::track(Metrics... metrics)
    {
        metrics_ = detail::MetricSet{ std::move(metrics)... };
    }

    template<typename Metric>
    requires std::derived_from<Metric, metrics::MonitorBase>
    inline const Metric& GaInfo::get_metric() const noexcept
    {
        GAPP_ASSERT(metrics_.get<Metric>(), "Attempting to get an untracked metric type is invalid.");

        return *metrics_.get<Metric>();
    }

    template<typename Metric>
    requires std::derived_from<Metric, metrics::MonitorBase>
    inline const Metric* GaInfo::get_metric_if() const noexcept
    {
        return metrics_.get<Metric>();
    }

} // namespace gapp

#endif // !GAPP_CORE_GA_INFO_HPP
