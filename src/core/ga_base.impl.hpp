/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_CORE_GA_BASE_IMPL_HPP
#define GAPP_CORE_GA_BASE_IMPL_HPP

#include "ga_base.decl.hpp"
#include "ga_traits.hpp"
#include "candidate.hpp"
#include "population.hpp"
#include "fitness_function.hpp"
#include "../encoding/gene_types.hpp"
#include "../algorithm/algorithm_base.hpp"
#include "../algorithm/any_objective.hpp"
#include "../crossover/crossover_base.hpp"
#include "../crossover/lambda.hpp"
#include "../mutation/mutation_base.hpp"
#include "../mutation/lambda.hpp"
#include "../stop_condition/stop_condition_base.hpp"
#include "../utility/algorithm.hpp"
#include "../utility/functional.hpp"
#include "../utility/thread_pool.hpp"
#include "../utility/scope_exit.hpp"
#include "../utility/type_id.hpp"
#include "../utility/type_list.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <functional>
#include <numeric>
#include <type_traits>
#include <tuple>
#include <memory>
#include <atomic>
#include <utility>
#include <cstddef>

namespace gapp
{
    template<typename T>
    GA<T>::GA(Positive<size_t> population_size,
              std::unique_ptr<typename algorithm::Algorithm> algorithm,
              std::unique_ptr<typename crossover::Crossover<T>> crossover,
              std::unique_ptr<typename mutation::Mutation<T>> mutation,
              std::unique_ptr<typename stopping::StopCondition> stop_condition) :
        GaInfo(population_size, std::move(algorithm), std::move(stop_condition)), crossover_(std::move(crossover)), mutation_(std::move(mutation))
    {
        GAPP_ASSERT(crossover_, "The crossover method can't be a nullptr.");
        GAPP_ASSERT(mutation_, "The mutation method can't be a nullptr.");
    }

    template<typename T>
    GA<T>::GA(Positive<size_t> population_size) :
        GA(population_size, std::make_unique<algorithm::AnyObjective<>>(), std::make_unique<typename GaTraits<T>::DefaultCrossover>(), std::make_unique<typename GaTraits<T>::DefaultMutation>())
    {}

    template<typename T>
    GA<T>::GA(Positive<size_t> population_size, std::unique_ptr<typename algorithm::Algorithm> algorithm) :
        GA(population_size, std::move(algorithm), std::make_unique<typename GaTraits<T>::DefaultCrossover>(), std::make_unique<typename GaTraits<T>::DefaultMutation>())
    {}

    template<typename T>
    GA<T>::GA(Positive<size_t> population_size,
              std::unique_ptr<typename crossover::Crossover<T>> crossover,
              std::unique_ptr<typename mutation::Mutation<T>> mutation,
              std::unique_ptr<typename stopping::StopCondition> stop_condition) :
        GA(population_size, std::make_unique<algorithm::AnyObjective<>>(), std::move(crossover), std::move(mutation), std::move(stop_condition))
    {}

    template<typename T>
    template<typename AlgorithmType>
    requires std::derived_from<AlgorithmType, algorithm::Algorithm>
    GA<T>::GA(Positive<size_t> population_size, AlgorithmType algorithm) :
        GA(population_size, std::make_unique<AlgorithmType>(std::move(algorithm)))
    {}

    template<typename T>
    template<typename CrossoverType, typename MutationType, typename StoppingType>
    requires std::derived_from<CrossoverType, crossover::Crossover<T>> &&
             std::derived_from<MutationType, mutation::Mutation<T>> &&
             std::derived_from<StoppingType, stopping::StopCondition>
    GA<T>::GA(Positive<size_t> population_size, CrossoverType crossover, MutationType mutation, StoppingType stop_condition) :
        GA(population_size,
           std::make_unique<algorithm::AnyObjective<>>(),
           std::make_unique<CrossoverType>(std::move(crossover)),
           std::make_unique<MutationType>(std::move(mutation)),
           std::make_unique<StoppingType>(std::move(stop_condition)))
    {}

    template<typename T>
    template<typename AlgorithmType, typename CrossoverType, typename MutationType, typename StoppingType>
    requires std::derived_from<AlgorithmType, algorithm::Algorithm> &&
             std::derived_from<CrossoverType, crossover::Crossover<T>> &&
             std::derived_from<MutationType, mutation::Mutation<T>> &&
             std::derived_from<StoppingType, stopping::StopCondition>
    GA<T>::GA(Positive<size_t> population_size, AlgorithmType algorithm, CrossoverType crossover, MutationType mutation, StoppingType stop_condition) :
        GA(population_size,
           std::make_unique<AlgorithmType>(std::move(algorithm)),
           std::make_unique<CrossoverType>(std::move(crossover)),
           std::make_unique<MutationType>(std::move(mutation)),
           std::make_unique<StoppingType>(std::move(stop_condition)))
    {}

    
    template<typename T>
    size_t GA<T>::index_of_gene(size_t type_id) const noexcept
    {
        return component_genes_t<T>::index_of_id(type_id).value_or(0);
    }

    template<typename T>
    template<typename U>
    BoundsVector<U> GA<T>::boundsToUniformBoundsVector(const FitnessFunctionBase<T>& fitness_function, Bounds<U> bounds) const
    {
        return BoundsVector<U>(fitness_function.template chrom_len<U>(), bounds);
    }

    template<typename T>
    template<typename... Us>
    auto GA<T>::boundsToUniformBoundsVector(const FitnessFunctionBase<T>& fitness_function, std::tuple<Bounds<Us>...> bounds) const
    {
        return std::tuple{ boundsToUniformBoundsVector(fitness_function, std::get<Bounds<Us>>(bounds))... };
    }
    
    template<typename T>
    inline const FitnessFunctionBase<T>* GA<T>::fitness_function() const& noexcept
    {
        return static_cast<const FitnessFunctionBase<T>*>(fitness_function_.get());
    }

    template<typename T>
    template<typename GeneType>
    inline size_t GA<T>::chrom_len() const noexcept
    {
        if (!fitness_function_) return 0;

        return fitness_function()->template chrom_len<GeneType>();
    }

    template<typename T>
    template<typename GeneType>
    inline const BoundsVector<GeneType>& GA<T>::gene_bounds() const noexcept
    {
        static_assert(is_bounded_gene_v<GeneType>, "The gene type must be a bounded gene type.");

        return std::get<BoundsVector<GeneType>>(bounds_);
    }


    template<typename T>
    template<typename F>
    requires std::derived_from<F, crossover::Crossover<T>>
    inline void GA<T>::crossover_method(F f)
    {
        crossover_ = std::make_unique<F>(std::move(f));
    }

    template<typename T>
    inline void GA<T>::crossover_method(std::unique_ptr<typename crossover::Crossover<T>> f)
    {
        GAPP_ASSERT(f, "The crossover method can't be a nullptr.");

        crossover_ = std::move(f);
    }

    template<typename T>
    inline void GA<T>::crossover_method(CrossoverCallable f)
    {
        crossover_ = std::make_unique<crossover::Lambda<T>>(std::move(f));
    }

    template<typename T>
    template<typename GeneType>
    inline crossover::Crossover<GeneType>& GA<T>::crossover_method() const& noexcept
    {
        GAPP_ASSERT(crossover_);

        if constexpr (std::is_same_v<GeneType, T>)
        {
            return *crossover_;
        }
        else
        {
            return crossover_->template component<GeneType>();
        }
    }

    template<typename T>
    void* GA<T>::crossover_method_impl(size_t type_id) const noexcept
    {
        GAPP_ASSERT(crossover_);

        if (detail::type_id<T>() == type_id)
        {
            return static_cast<void*>(crossover_.get());
        }

        if constexpr (is_mixed_gene_v<T>)
        {
            return component_genes_t<T>::apply([&]<typename... Ts>()
            {
                void* crossover_method = nullptr;
                ((detail::type_id<Ts>() == type_id ? crossover_method = (void*)std::addressof(crossover_->template component<Ts>())
                                                   : crossover_method), ...);
                return crossover_method;
            });
        }
        else return nullptr;
    }

    template<typename T>
    template<typename F>
    requires std::derived_from<F, mutation::Mutation<T>>
    inline void GA<T>::mutation_method(F f)
    {
        mutation_ = std::make_unique<F>(std::move(f));
    }

    template<typename T>
    inline void GA<T>::mutation_method(std::unique_ptr<typename mutation::Mutation<T>> f)
    {
        GAPP_ASSERT(f, "The mutation method can't be a nullptr.");

        mutation_ = std::move(f);
    }

    template<typename T>
    inline void GA<T>::mutation_method(MutationCallable f)
    {
        mutation_ = std::make_unique<mutation::Lambda<T>>(std::move(f));
    }

    template<typename T>
    template<typename GeneType>
    inline mutation::Mutation<GeneType>& GA<T>::mutation_method() const& noexcept
    {
        GAPP_ASSERT(mutation_);

        if constexpr (std::is_same_v<GeneType, T>)
        {
            return *mutation_;
        }
        else
        {
            return mutation_->template component<GeneType>();
        }
    }

    template<typename T>
    void* GA<T>::mutation_method_impl(size_t type_id) const noexcept
    {
        GAPP_ASSERT(mutation_);

        if (detail::type_id<T>() == type_id)
        {
            return static_cast<void*>(mutation_.get());
        }

        if constexpr (is_mixed_gene_v<T>)
        {
            return component_genes_t<T>::apply([&]<typename... Ts>()
            {
                void* mutation_method = nullptr;
                ((detail::type_id<Ts>() == type_id ? mutation_method = (void*)std::addressof(mutation_->template component<Ts>())
                                                   : mutation_method), ...);
                return mutation_method;
            });
        }
        else return nullptr;
    }

    template<typename T>
    inline void GA<T>::constraints_function(ConstraintsFunction f) noexcept
    {
        /* Nullptr is fine here, it just won't be called */
        constraints_function_ = std::move(f);
    }

    template<typename T>
    inline void GA<T>::repair_function(RepairCallable f) noexcept
    {
        /* Nullptr is fine here, it just won't be called */
        repair_ = std::move(f);
    }

    template<typename T>
    inline void GA<T>::cache_size(size_t generations) noexcept
    {
        cached_generations_ = generations;
    }

    template<typename T>
    inline std::pair<size_t, size_t> GA<T>::findObjectiveProperties() const
    {
        GAPP_ASSERT(fitness_function_);

        Candidate<T> candidate = generateCandidate();

        if (constraints_function_)
        {
            candidate.constraint_violation = constraints_function_(*this, candidate);
        }

        const FitnessVector fitness = std::invoke(*fitness_function(), candidate);
        GAPP_ASSERT(!fitness.empty(), "The number of objectives must be greater than 0.");

        return { fitness.size(), candidate.num_constraints() };
    }

    template<typename T>
    template<typename GeneType>
    inline void GA<T>::setDefaultMutationRate() const
    {
        GAPP_ASSERT(mutation_);
        GAPP_ASSERT(fitness_function_);

        if (!mutation_method<GeneType>().use_default_mutation_rate()) return;
        mutation_method<GeneType>().mutation_rate(GaTraits<GeneType>::defaultMutationRate(chrom_len<GeneType>()));
    }

    template<typename T>
    inline void GA<T>::setDefaultMutationRates() const
    {
        GAPP_ASSERT(mutation_);
        GAPP_ASSERT(fitness_function_);

        component_genes_t<T>::apply([&]<typename... Ts>()
        {
            (setDefaultMutationRate<Ts>(), ...);
        });
    }

    template<typename T>
    inline bool GA<T>::hasValidFitness(const Candidate<T>& sol) const noexcept
    {
        GAPP_ASSERT(fitness_function_);

        return sol.is_evaluated() && (sol.fitness.size() == num_objectives());
    }

    template<typename T>
    inline bool GA<T>::hasValidConstraints(const Candidate<T>& sol) const noexcept
    {
        GAPP_ASSERT(fitness_function_);

        return sol.num_constraints() == num_constraints();
    }

    template<typename T>
    inline bool GA<T>::hasValidChromosome(const Candidate<T>& sol) const noexcept
    {
        return component_genes_t<T>::apply([&]<typename... Ts>()
        {
            return ((sol.template chrom_len<Ts>() == chrom_len<Ts>() ||
                    (crossover_method<Ts>().allow_variable_chrom_length() && mutation_method<Ts>().allow_variable_chrom_length())) && ...);
        });
    }

    template<typename T>
    inline bool GA<T>::isValidEvaluatedPopulation(const Population<T>& pop) const
    {
        return std::all_of(pop.begin(), pop.end(), [this](const Candidate<T>& sol)
        {
            return hasValidChromosome(sol) && hasValidFitness(sol) && hasValidConstraints(sol);
        });
    }

    template<typename T>
    inline bool GA<T>::isValidUnevaluatedPopulation(const Population<T>& pop) const
    {
        return std::all_of(pop.begin(), pop.end(), [this](const Candidate<T>& sol)
        {
            return hasValidChromosome(sol) && (!sol.is_evaluated() || hasValidFitness(sol));
        });
    }

    template<typename T>
    inline bool GA<T>::isValidBoundsVectors(const BoundsVectors& bounds) const
    {
        return bounded_component_genes_t<T>::apply([&]<typename... Us>()
        {
            return ((std::get<BoundsVector<Us>>(bounds).size() == chrom_len<Us>()) && ...);
        });
    }

    template<typename T>
    inline bool GA<T>::fitnessMatrixIsSynced() const
    {
        return std::equal(fitness_matrix_.begin(), fitness_matrix_.end(), population_.begin(), population_.end(),
        [](const auto& fvec, const auto& sol)
        {
            return fvec == sol.fitness;
        });
    }


    template<typename T>
    void GA<T>::initializeAlgorithm(BoundsVectors bounds, Population<T> initial_population)
    {
        GAPP_ASSERT(fitness_function_);
        GAPP_ASSERT(algorithm_ && stop_condition_);
        GAPP_ASSERT(crossover_ && mutation_);
        GAPP_ASSERT(isValidBoundsVectors(bounds));

        detail::execution_context::global_thread_pool.reset_scheduler();

        /* Reset state in case solve() has already been called before. */
        generation_cntr_ = 0;
        num_fitness_evals_ = 0;
        solutions_.clear();
        population_.clear();

        fitness_cache_.reset(size_t(!fitness_function_->is_dynamic()) * cached_generations_ * population_size_);

        bounds_ = std::move(bounds);

        /* Derived GA. */
        initialize();

        /* Create and evaluate the initial population of the algorithm. */
        std::tie(num_objectives_, num_constraints_) = findObjectiveProperties();
        population_ = generatePopulation(population_size_, std::move(initial_population));
        detail::parallel_for(population_.begin(), population_.end(), [this](Candidate<T>& sol) { validate(sol); repair(sol); evaluate(sol); });
        fitness_matrix_ = detail::toFitnessMatrix(population_);
        if (keep_all_optimal_sols_) solutions_ = detail::findParetoFront(population_);

        GAPP_ASSERT(isValidEvaluatedPopulation(population_));
        GAPP_ASSERT(fitnessMatrixIsSynced());

        /* Initialize the algorithm used.
         * This must be done after the initial population has been created and evaluted,
         * as it might want to use the population's fitness values (fitness_matrix_). */
        algorithm_->initialize(*this);

        setDefaultMutationRates();

        stop_condition_->initialize(*this);

        metrics_.initialize(*this);
        metrics_.update(*this);

        if (on_generation_end_) on_generation_end_(*this);
    }

    template<typename T>
    template<typename GeneType>
    Candidate<GeneType> GA<T>::generateCandidate() const
    {
        static_assert(!is_mixed_gene_v<GeneType>);

        const size_t chrom_size = chrom_len<GeneType>();

        if constexpr (is_bounded_gene_v<GeneType>)
        {
            const BoundsVector<GeneType>& bounds = gene_bounds<GeneType>();
            return Candidate<GeneType>(GaTraits<GeneType>::randomChromosome(chrom_size, bounds), bounds);
        }
        else
        {
            return Candidate<GeneType>(GaTraits<GeneType>::randomChromosome(chrom_size));
        }
    }

    template<typename T>
    Candidate<T> GA<T>::generateCandidate() const
    {
        return component_genes_t<T>::apply([&]<typename... Ts>()
        {
            return Candidate<T>{{ generateCandidate<Ts>()... }};
        });
    }

    template<typename T>
    Population<T> GA<T>::generatePopulation(Positive<size_t> pop_size, Population<T> initial_population) const
    {
        GAPP_ASSERT(isValidUnevaluatedPopulation(initial_population), "An invalid initial population was specified for the GA.");

        Population<T> population;
        population.reserve(pop_size);

        const size_t npreset = std::min(size_t(pop_size), initial_population.size());
        std::move(initial_population.begin(), initial_population.begin() + npreset, std::back_inserter(population));

        while (population.size() < pop_size)
        {
            population.push_back(generateCandidate());
            GAPP_ASSERT(hasValidChromosome(population.back()), "An invalid solution was returned by generateCandidate().");
        }

        return population;
    }

    template<typename T>
    inline void GA<T>::prepareSelections() const
    {
        GAPP_ASSERT(algorithm_);
        GAPP_ASSERT(isValidEvaluatedPopulation(population_));
        GAPP_ASSERT(fitnessMatrixIsSynced());

        algorithm_->prepareSelections(*this, population_);
    }

    template<typename T>
    inline const Candidate<T>& GA<T>::select() const
    {
        GAPP_ASSERT(algorithm_);

        return algorithm_->select(*this, population_);
    }

    template<typename T>
    inline CandidatePair<T> GA<T>::crossover(const Candidate<T>& parent1, const Candidate<T>& parent2) const
    {
        GAPP_ASSERT(crossover_);

        auto [child1, child2] = (*crossover_)(*this, parent1, parent2);

        child1.fitness.clear();
        child2.fitness.clear();

        /*
        * Check if either of the children is the same as one of the parents.
        * This can happen in edge cases even if the parents are different.
        * If one of the children is the same as one of the parents, then the fitness function
        * evaluation for that child can be skipped by assigning it the same fitness as the parent.
        */
        if (child1 == parent1)
        {
            child1.fitness = parent1.fitness;
        }
        else if (child1 == parent2)
        {
            child1.fitness = parent2.fitness;
        }
        if (child2 == parent1)
        {
            child2.fitness = parent1.fitness;
        }
        else if (child2 == parent2)
        {
            child2.fitness = parent2.fitness;
        }

        return { std::move(child1), std::move(child2) };
    }

    template<typename T>
    inline void GA<T>::mutate(Candidate<T>& candidate) const
    {
        GAPP_ASSERT(mutation_);

        if (!candidate.is_evaluated())
        {
            (*mutation_)(*this, candidate);
            return;
        }

        /* 
        * If the candidate has a valid fitness vector, which can happen when the crossover
        * didn't change the candidate, we can save a fitness function call in the cases 
        * the mutation also doesn't change the candidate. 
        */
        thread_local Candidate<T> old_candidate;
        old_candidate = candidate;

        (*mutation_)(*this, candidate);

        if (candidate != old_candidate)
        {
            candidate.fitness.clear();
        }
    }

    template<typename T>
    void GA<T>::validate(Candidate<T>& candidate) const
    {
        GAPP_ASSERT(hasValidChromosome(candidate));

        /* Don't do anything for unconstrained optimization problems. */
        if (!constraints_function_) return;
        candidate.constraint_violation = constraints_function_(*this, candidate);

        GAPP_ASSERT(hasValidConstraints(candidate), "An invalid constraints vector was returned by the constraints function.");
    }

    template<typename T>
    void GA<T>::repair(Candidate<T>& candidate) const
    {
        GAPP_ASSERT(hasValidChromosome(candidate));

        /* Don't try to do anything unless a repair function is set. */
        if (!repair_) return;

        if (repair_(*this, candidate))
        {
            candidate.fitness.clear();
            validate(candidate);
        }

        GAPP_ASSERT(hasValidChromosome(candidate), "Invalid chromosome returned by the repair function.");
    }

    template<typename T>
    void GA<T>::updatePopulation(Population<T>&& children)
    {
        GAPP_ASSERT(algorithm_);
        GAPP_ASSERT(isValidEvaluatedPopulation(population_));
        GAPP_ASSERT(fitnessMatrixIsSynced());

        fitness_cache_.insert(population_.begin(), population_.end(), &Candidate<T>::fitness);
        population_ = algorithm_->nextPopulation(*this, std::move(population_), std::move(children));
        fitness_matrix_ = detail::toFitnessMatrix(population_);
    }

    template<typename T>
    inline bool GA<T>::stopCondition() const
    {
        GAPP_ASSERT(stop_condition_);

        return (*stop_condition_)(*this);
    }

    template<typename T>
    inline void GA<T>::evaluate(Candidate<T>& candidate)
    {
        GAPP_ASSERT(fitness_function_);
        GAPP_ASSERT(hasValidChromosome(candidate));

        /* If the fitness function is static, and the solution has already
         * been evaluted sometime earlier (in an earlier generation), there
         * is no point doing it again. */
        if (!fitness_function_->is_dynamic() && candidate.is_evaluated()) return;
        
        if (cached_generations_)
        {
            GAPP_ASSERT(!fitness_function_->is_dynamic());

            if (const FitnessVector* fitness = fitness_cache_.get(candidate))
            {
                candidate.fitness = *fitness;
                return;
            }
        }

        std::atomic_ref{ num_fitness_evals_ }.fetch_add(1, std::memory_order_release);
        candidate.fitness = (*fitness_function())(candidate);

        GAPP_ASSERT(hasValidFitness(candidate), "Invalid fitness vector returned by the fitness function.");
    }

    template<typename T>
    void GA<T>::updateOptimalSolutions(Candidates<T>& optimal_sols, const Population<T>& pop) const
    {
        GAPP_ASSERT(algorithm_);

        auto optimal_pop = algorithm_->optimalSolutions(*this, pop);

        optimal_sols = detail::mergeParetoSets(std::move(optimal_sols), std::move(optimal_pop));

        /* Duplicate elements are removed from optimal_sols using exact comparison
        *  of the chromosomes in order to avoid issues with using a non-transitive
        *  comparison function for std::sort and std::unique. */
        std::sort(optimal_sols.begin(), optimal_sols.end(), detail::CandidateLessThanExact{});
        auto last = std::unique(optimal_sols.begin(), optimal_sols.end(), detail::CandidateEqualExact{});
        optimal_sols.erase(last, optimal_sols.end());
    }

    template<typename T>
    void GA<T>::advance()
    {
        GAPP_ASSERT(population_.size() == population_size_);

        Population<T> children(population_size_);

        prepareSelections();

        detail::parallel_for(detail::iota_iterator(0_sz), detail::iota_iterator(population_size_ / 2), [&](size_t i)
        {
            CandidatePair child_pair = crossover(select(), select());
            children[2 * i]     = std::move(child_pair.first);
            children[2 * i + 1] = std::move(child_pair.second);
        });

        if (population_size_ % 2) children.back() = crossover(select(), select()).first;

        detail::parallel_for(children.begin(), children.end(), [this](Candidate<T>& child)
        {
            mutate(child);
            validate(child);
            repair(child);
            evaluate(child);
        });

        updatePopulation(std::move(children));

        if (keep_all_optimal_sols_) updateOptimalSolutions(solutions_, population_);
        metrics_.update(*this);

        if (on_generation_end_) on_generation_end_(*this);
        generation_cntr_++;
    }

    template<typename T>
    Candidates<T> GA<T>::solve(std::unique_ptr<FitnessFunctionBase<T>> fitness_function, size_t generations, Population<T> initial_population) requires(!is_partially_bounded_gene_v<T>)
    {
        GAPP_ASSERT(fitness_function, "The fitness function can't be a nullptr.");

        detail::restore_on_exit _{ max_gen_ };

        fitness_function_ = std::move(fitness_function);
        max_gen(generations);

        initializeAlgorithm({ /* no bounds */ }, std::move(initial_population));
        while (!stopCondition())
        {
            advance();
        }
        if (!keep_all_optimal_sols_) updateOptimalSolutions(solutions_, population_);

        return solutions_;
    }

    template<typename T>
    Candidates<T> GA<T>::solve(std::unique_ptr<FitnessFunctionBase<T>> fitness_function, BoundsVectors bounds, size_t generations, Population<T> initial_population) requires(is_partially_bounded_gene_v<T>)
    {
        GAPP_ASSERT(fitness_function, "The fitness function can't be a nullptr.");

        detail::restore_on_exit _{ max_gen_ };

        fitness_function_ = std::move(fitness_function);
        max_gen(generations);

        initializeAlgorithm(std::move(bounds), std::move(initial_population));
        while (!stopCondition())
        {
            advance();
        }
        if (!keep_all_optimal_sols_) updateOptimalSolutions(solutions_, population_);

        return solutions_;
    }

    template<typename T>
    template<typename F>
    requires (!is_partially_bounded_gene_v<T> && std::derived_from<F, FitnessFunctionBase<T>>)
    Candidates<T> GA<T>::solve(F fitness_function, Population<T> initial_population)
    {
        return solve(std::make_unique<F>(std::move(fitness_function)), max_gen(), std::move(initial_population));
    }

    template<typename T>
    template<typename F>
    requires (!is_partially_bounded_gene_v<T> && std::derived_from<F, FitnessFunctionBase<T>>)
    Candidates<T> GA<T>::solve(F fitness_function, size_t generations, Population<T> initial_population)
    {
        return solve(std::make_unique<F>(std::move(fitness_function)), generations, std::move(initial_population));
    }

    template<typename T>
    template<typename F>
    requires (is_partially_bounded_gene_v<T> && std::derived_from<F, FitnessFunctionBase<T>>)
    Candidates<T> GA<T>::solve(F fitness_function, BoundsVectors bounds, Population<T> initial_population)
    {
        return solve(std::make_unique<F>(std::move(fitness_function)), std::move(bounds), max_gen(), std::move(initial_population));
    }

    template<typename T>
    template<typename F>
    requires (is_partially_bounded_gene_v<T> && std::derived_from<F, FitnessFunctionBase<T>>)
    Candidates<T> GA<T>::solve(F fitness_function, BoundsList bounds, Population<T> initial_population)
    {
        BoundsVectors bounds_vectors = boundsToUniformBoundsVector(fitness_function, bounds);
        return solve(std::make_unique<F>(std::move(fitness_function)), std::move(bounds_vectors), max_gen(), std::move(initial_population));
    }

    template<typename T>
    template<typename F>
    requires (is_partially_bounded_gene_v<T> && std::derived_from<F, FitnessFunctionBase<T>>)
    Candidates<T> GA<T>::solve(F fitness_function, BoundsVectors bounds, size_t generations, Population<T> initial_population)
    {
        return solve(std::make_unique<F>(std::move(fitness_function)), std::move(bounds), generations, std::move(initial_population));
    }

    template<typename T>
    template<typename F>
    requires (is_partially_bounded_gene_v<T> && std::derived_from<F, FitnessFunctionBase<T>>)
    Candidates<T> GA<T>::solve(F fitness_function, BoundsList bounds, size_t generations, Population<T> initial_population)
    {
        BoundsVectors bounds_vectors = boundsToUniformBoundsVector(fitness_function, bounds);
        return solve(std::make_unique<F>(std::move(fitness_function)), std::move(bounds_vectors), generations, std::move(initial_population));
    }

} // namespace gapp

#endif // !GAPP_CORE_GA_BASE_IMPL_HPP
