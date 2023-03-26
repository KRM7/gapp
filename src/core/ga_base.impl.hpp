/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CORE_GA_BASE_IMPL_HPP
#define GA_CORE_GA_BASE_IMPL_HPP

#include "ga_base.decl.hpp"
#include "fitness_function.hpp"
#include "../population/population.hpp"
#include "../algorithm/algorithm_base.hpp"
#include "../algorithm/single_objective.hpp"
#include "../crossover/crossover_base.hpp"
#include "../crossover/lambda.hpp"
#include "../mutation/mutation_base.hpp"
#include "../mutation/lambda.hpp"
#include "../stop_condition/stop_condition_base.hpp"
#include "../utility/algorithm.hpp"
#include "../utility/functional.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <functional>
#include <numeric>
#include <execution>
#include <type_traits>
#include <atomic>
#include <utility>

namespace genetic_algorithm
{
    template<Gene T>
    GA<T>::GA(std::unique_ptr<FitnessFunction<T>> f, Positive<size_t> population_size) :
        GaInfo(population_size, f->num_objectives())
    {
        GA_ASSERT(f, "The fitness function can't be a nullptr.");

        fitness_function_ = std::move(f);
    }

    template<Gene T>
    template<typename F>
    requires std::derived_from<F, FitnessFunction<T>> && std::is_final_v<F>
    void GA<T>::fitness_function(F f)
    {
        fitness_function_ = std::make_unique<F>(std::move(f));
        is_initialized_ = false;
    }

    template<Gene T>
    void GA<T>::fitness_function(std::unique_ptr<FitnessFunction<T>> f)
    {
        GA_ASSERT(f, "The fitness function can't be a nullptr.");

        fitness_function_ = std::move(f);
        is_initialized_ = false;
    }

    template<Gene T>
    inline FitnessFunction<T>& GA<T>::fitness_function() & noexcept
    {
        GA_ASSERT(fitness_function_);

        return *fitness_function_;
    }

    template<Gene T>
    inline const FitnessFunction<T>& GA<T>::fitness_function() const& noexcept
    {
        GA_ASSERT(fitness_function_);

        return *fitness_function_;
    }

    template<Gene T>
    inline size_t GA<T>::chrom_len() const noexcept
    {
        GA_ASSERT(fitness_function_);

        return fitness_function_->chrom_len();
    }

    template<Gene T>
    inline bool GA<T>::variable_chrom_len() const noexcept
    {
        GA_ASSERT(fitness_function_);

        return fitness_function_->variable_chrom_len();
    }

    template<Gene T>
    inline size_t GA<T>::num_objectives() const noexcept
    {
        GA_ASSERT(fitness_function_);

        return fitness_function_->num_objectives();
    }

    template<Gene T>
    inline bool GA<T>::dynamic_fitness() const noexcept
    {
        GA_ASSERT(fitness_function_);

        return fitness_function_->dynamic();
    }


    template<Gene T>
    inline auto GA<T>::gene_bounds() const noexcept -> const BoundsVector<T>&
    {
        return bounds_;
    }


    template<Gene T>
    template<typename F>
    requires std::derived_from<F, crossover::Crossover<T>> && std::is_final_v<F>
    inline void GA<T>::crossover_method(F f)
    {
        crossover_ = std::make_unique<F>(std::move(f));
    }

    template<Gene T>
    inline void GA<T>::crossover_method(std::unique_ptr<typename crossover::Crossover<T>> f)
    {
        GA_ASSERT(f, "The crossover method can't be a nullptr.");

        crossover_ = std::move(f);
    }

    template<Gene T>
    inline void GA<T>::crossover_method(CrossoverCallable f)
    {
        crossover_ = std::make_unique<crossover::Lambda<T>>(std::move(f));
    }

    template<Gene T>
    inline crossover::Crossover<T>& GA<T>::crossover_method() & noexcept
    {
        GA_ASSERT(crossover_);

        return *crossover_;
    }

    template<Gene T>
    inline const crossover::Crossover<T>& GA<T>::crossover_method() const& noexcept
    {
        GA_ASSERT(crossover_);

        return *crossover_;
    }

    template<Gene T>
    inline void GA<T>::crossover_rate(Probability pc) noexcept
    {
        GA_ASSERT(crossover_);

        crossover_->crossover_rate(pc);
    }

    template<Gene T>
    inline Probability GA<T>::crossover_rate() const noexcept
    {
        GA_ASSERT(crossover_);

        return crossover_->crossover_rate();
    }

    template<Gene T>
    template<typename F>
    requires std::derived_from<F, mutation::Mutation<T>> && std::is_final_v<F>
    inline void GA<T>::mutation_method(F f)
    {
        mutation_ = std::make_unique<F>(std::move(f));
    }

    template<Gene T>
    inline void GA<T>::mutation_method(std::unique_ptr<typename mutation::Mutation<T>> f)
    {
        GA_ASSERT(f, "The mutation method can't be a nullptr.");

        mutation_ = std::move(f);
    }

    template<Gene T>
    inline void GA<T>::mutation_method(MutationCallable f)
    {
        mutation_ = std::make_unique<mutation::Lambda<T>>(std::move(f));
    }

    template<Gene T>
    inline mutation::Mutation<T>& GA<T>::mutation_method() & noexcept
    {
        GA_ASSERT(mutation_);

        return *mutation_;
    }

    template<Gene T>
    inline const mutation::Mutation<T>& GA<T>::mutation_method() const& noexcept
    {
        GA_ASSERT(mutation_);

        return *mutation_;
    }

    template<Gene T>
    inline void GA<T>::mutation_rate(Probability pm) noexcept
    {
        GA_ASSERT(mutation_);

        mutation_->mutation_rate(pm);
    }

    template<Gene T>
    inline Probability GA<T>::mutation_rate() const noexcept
    {
        GA_ASSERT(mutation_);

        return mutation_->mutation_rate();
    }

    template<Gene T>
    inline void GA<T>::repair_function(RepairCallable f)
    {
        /* Nullptr is fine here, it just won't be called */
        repair_ = std::move(f);
    }


    template<Gene T>
    inline bool GA<T>::hasValidFitness(const Candidate<T>& sol) const noexcept
    {
        GA_ASSERT(fitness_function_);

        return sol.is_evaluated && (sol.fitness.size() == fitness_function_->num_objectives());
    }

    template<Gene T>
    inline bool GA<T>::hasValidChromosome(const Candidate<T>& sol) const noexcept
    {
        return variable_chrom_len() || (sol.chromosome.size() == chrom_len());
    }

    template<Gene T>
    inline bool GA<T>::isValidEvaluatedPopulation(const Population<T>& pop) const
    {
        return std::all_of(pop.begin(), pop.end(), [this](const Candidate<T>& sol)
        {
            return hasValidChromosome(sol) && hasValidFitness(sol);
        });
    }

    template<Gene T>
    inline bool GA<T>::isValidUnevaluatedPopulation(const Population<T>& pop) const
    {
        return std::all_of(pop.begin(), pop.end(), [this](const Candidate<T>& sol)
        {
            return hasValidChromosome(sol) && (!sol.is_evaluated || hasValidFitness(sol));
        });
    }

    template<Gene T>
    inline bool GA<T>::fitnessMatrixIsSynced() const
    {
        return std::equal(fitness_matrix_.begin(), fitness_matrix_.end(), population_.begin(), population_.end(),
        [this](const auto& fvec, const auto& sol)
        {
            return fvec == sol.fitness;
        });
    }


    template<Gene T>
    void GA<T>::initializeAlgorithm(Population<T> initial_population)
    {
        GA_ASSERT(fitness_function_);
        GA_ASSERT(algorithm_ && stop_condition_);
        GA_ASSERT(crossover_ && mutation_);

        /* Derived GA. */
        initialize();

        /* Reset state variables in case run() has already been called before. */
        generation_cntr_ = 0;
        num_fitness_evals_ = 0;
        solutions_.clear();
        population_.clear();

        /* Create and evaluate the initial population of the algorithm. */
        population_ = generatePopulation(population_size_, std::move(initial_population));
        std::for_each(GA_EXECUTION_UNSEQ, population_.begin(), population_.end(), [this](Candidate<T>& sol) { evaluate(sol); });
        fitness_matrix_ = detail::toFitnessMatrix(population_);

        GA_ASSERT(isValidEvaluatedPopulation(population_));
        GA_ASSERT(fitnessMatrixIsSynced());

        /* Initialize the algorithm used.
         * This must be done after the initial population has been created and evaluted,
         * as it might want to use the population's fitness values (fitness_matrix_). */
        algorithm_->initialize(*this);

        is_initialized_ = true;
    }

    template<Gene T>
    Population<T> GA<T>::generatePopulation(Positive<size_t> pop_size, Population<T> initial_population) const
    {
        GA_ASSERT(isValidUnevaluatedPopulation(initial_population), "An invalid initial population was specified for the GA.");

        Population<T> population;
        population.reserve(pop_size);

        const size_t npreset = std::min(size_t(pop_size), initial_population.size());
        std::move(initial_population.begin(), initial_population.begin() + npreset, std::back_inserter(population));

        while (population.size() < pop_size)
        {
            population.push_back(generateCandidate());
            GA_ASSERT(hasValidChromosome(population.back()), "An invalid solution was returned by generateCandidate().");
        }

        return population;
    }

    template<Gene T>
    inline void GA<T>::prepareSelections() const
    {
        GA_ASSERT(algorithm_);
        GA_ASSERT(isValidEvaluatedPopulation(population_));
        GA_ASSERT(fitnessMatrixIsSynced());

        algorithm_->prepareSelections(*this, fitness_matrix());
    }

    template<Gene T>
    inline const Candidate<T>& GA<T>::select() const
    {
        GA_ASSERT(algorithm_);

        return algorithm_->select(*this, population(), fitness_matrix());
    }

    template<Gene T>
    inline CandidatePair<T> GA<T>::crossover(const Candidate<T>& parent1, const Candidate<T>& parent2) const
    {
        GA_ASSERT(crossover_);

        return (*crossover_)(*this, parent1, parent2);
    }

    template<Gene T>
    inline void GA<T>::mutate(Candidate<T>& sol) const
    {
        GA_ASSERT(mutation_);

        (*mutation_)(*this, sol);
    }

    template<Gene T>
    void GA<T>::repair(Candidate<T>& sol) const
    {
        GA_ASSERT(hasValidChromosome(sol));

        /* Don't try to do anything unless a repair function is set. */
        if (!repair_) return;

        auto improved_chrom = repair_(*this, sol.chromosome);

        if (improved_chrom != sol.chromosome)
        {
            sol.is_evaluated = false;
            sol.chromosome = std::move(improved_chrom);
        }

        GA_ASSERT(hasValidChromosome(sol), "Invalid chromosome returned by the repair function.");
    }

    template<Gene T>
    void GA<T>::updatePopulation(Population<T>&& children)
    {
        GA_ASSERT(algorithm_);
        GA_ASSERT(isValidEvaluatedPopulation(population_));
        GA_ASSERT(fitnessMatrixIsSynced());

        population_ = algorithm_->nextPopulation(*this, std::move(population_), std::move(children));
        fitness_matrix_ = detail::toFitnessMatrix(population_);
    }

    template<Gene T>
    inline bool GA<T>::stopCondition() const
    {
        GA_ASSERT(stop_condition_);

        return (*stop_condition_)(*this);
    }

    template<Gene T>
    inline void GA<T>::evaluate(Candidate<T>& sol)
    {
        GA_ASSERT(fitness_function_);
        GA_ASSERT(hasValidChromosome(sol));

        /* If the fitness function is static, and the solution has already
         * been evaluted sometime earlier (in an earlier generation), there
         * is no point doing it again. */
        if (!sol.is_evaluated || fitness_function_->dynamic())
        {
            sol.fitness = (*fitness_function_)(sol.chromosome);
            sol.is_evaluated = true;

            std::atomic_ref num_evals{ num_fitness_evals_ };
            num_evals.fetch_add(1_sz, std::memory_order_acq_rel);
        }

        GA_ASSERT(hasValidFitness(sol));
    }

    template<Gene T>
    void GA<T>::updateOptimalSolutions(Candidates<T>& optimal_sols, const Population<T>& pop) const
    {
        GA_ASSERT(algorithm_);

        auto optimal_pop = algorithm_->optimalSolutions(*this, pop);

        optimal_sols = detail::mergeParetoSets(std::move(optimal_sols), std::move(optimal_pop));
        detail::erase_duplicates(optimal_sols);
    }

    template<Gene T>
    void GA<T>::advance()
    {
        GA_ASSERT(population_.size() == population_size_);

        if (keep_all_optimal_sols_) updateOptimalSolutions(solutions_, population_);

        const size_t num_children = population_size_ + population_size_ % 2;
        std::vector<CandidatePair<T>> child_pairs(num_children / 2);

        prepareSelections();
        std::generate(GA_EXECUTION_UNSEQ, child_pairs.begin(), child_pairs.end(),
        [this]
        {
            const auto& parent1 = select();
            const auto& parent2 = select();

            return crossover(parent1, parent2);
        });

        auto children = detail::flatten(std::move(child_pairs));

        /* If the population size is odd, one too many child candidates were generated by the crossovers. */
        if (children.size() > population_size_) children.pop_back();

        std::for_each(GA_EXECUTION_UNSEQ, children.begin(), children.end(),
        [this](Candidate<T>& child)
        {
            mutate(child);
            repair(child);
            evaluate(child);
        });

        updatePopulation(std::move(children));

        if (endOfGenerationCallback) endOfGenerationCallback(*this);
        generation_cntr_++;
    }

    template<Gene T>
    inline Candidates<T> GA<T>::run(Population<T> initial_population)
    {
        return run(std::move(initial_population), max_gen());
    }

    template<Gene T>
    inline Candidates<T> GA<T>::run(size_t num_generations)
    {
        return run({}, num_generations);
    }

    template<Gene T>
    Candidates<T> GA<T>::run(Population<T> initial_population, size_t num_generations)
    {
        max_gen(num_generations);

        initializeAlgorithm(std::move(initial_population));
        while (!stopCondition())
        {
            advance();
        }
        updateOptimalSolutions(solutions_, population_);

        return solutions_;
    }

    template<Gene T>
    Candidates<T> GA<T>::continueFor(size_t num_generations)
    {
        max_gen(max_gen_ + num_generations);

        if (!is_initialized_) initializeAlgorithm(); // fallback to run() if called on uninitialized GA
        while (!stopCondition())
        {
            advance();
        }
        updateOptimalSolutions(solutions_, population_);

        return solutions_;
    }

} // namespace genetic_algorithm

#endif // !GA_CORE_GA_BASE_IMPL_HPP