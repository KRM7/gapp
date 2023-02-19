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
#include <stdexcept>
#include <cassert>

namespace genetic_algorithm
{
    template<Gene T>
    GA<T>::GA(std::unique_ptr<FitnessFunction<T>> f, size_t population_size) :
        GaInfo(population_size, f->num_objectives())
    {
        if (!f) GA_THROW(std::invalid_argument, "The fitness function can't be a nullptr.");

        fitness_function_ = std::move(f);
    }

    template<Gene T>
    template<typename F>
    requires FitnessFunctionType<F, T> && std::is_final_v<F>
    void GA<T>::fitness_function(F f)
    {
        fitness_function_ = std::make_unique<F>(std::move(f));
        is_initialized_ = false;
    }

    template<Gene T>
    void GA<T>::fitness_function(std::unique_ptr<FitnessFunction<T>> f)
    {
        if (!f) GA_THROW(std::invalid_argument, "The fitness function can't be a nullptr.");

        fitness_function_ = std::move(f);
        is_initialized_ = false;
    }

    template<Gene T>
    inline FitnessFunction<T>& GA<T>::fitness_function() &
    {
        assert(fitness_function_ != nullptr);

        return *fitness_function_;
    }

    template<Gene T>
    inline const FitnessFunction<T>& GA<T>::fitness_function() const&
    {
        assert(fitness_function_ != nullptr);

        return *fitness_function_;
    }

    template<Gene T>
    inline size_t GA<T>::chrom_len() const noexcept
    {
        assert(fitness_function_ != nullptr);

        return fitness_function_->chrom_len();
    }

    template<Gene T>
    inline bool GA<T>::variable_chrom_len() const noexcept
    {
        assert(fitness_function_ != nullptr);

        return fitness_function_->variable_chrom_len();
    }

    template<Gene T>
    inline size_t GA<T>::num_objectives() const noexcept
    {
        assert(fitness_function_ != nullptr);

        return fitness_function_->num_objectives();
    }

    template<Gene T>
    inline bool GA<T>::dynamic_fitness() const noexcept
    {
        assert(fitness_function_ != nullptr);

        return fitness_function_->dynamic_fitness();
    }


    template<Gene T>
    inline auto GA<T>::gene_bounds() const noexcept -> const BoundsVector&
    {
        return bounds_;
    }


    template<Gene T>
    template<typename F>
    requires crossover::CrossoverType<F, T> && std::is_final_v<F>
    inline void GA<T>::crossover_method(F f)
    {
        crossover_ = std::make_unique<F>(std::move(f));
    }

    template<Gene T>
    template<crossover::CrossoverType<T> F>
    inline void GA<T>::crossover_method(std::unique_ptr<F> f)
    {
        if (!f) GA_THROW(std::invalid_argument, "The crossover method can't be a nullptr.");

        crossover_ = std::move(f);
    }

    template<Gene T>
    inline void GA<T>::crossover_method(CrossoverCallable f)
    {
        crossover_ = std::make_unique<crossover::Lambda<T>>(std::move(f));
    }

    template<Gene T>
    inline crossover::Crossover<T>& GA<T>::crossover_method() &
    {
        assert(crossover_ != nullptr);

        return *crossover_;
    }

    template<Gene T>
    inline const crossover::Crossover<T>& GA<T>::crossover_method() const&
    {
        assert(crossover_ != nullptr);

        return *crossover_;
    }

    template<Gene T>
    inline void GA<T>::crossover_rate(Probability pc)
    {
        assert(crossover_ != nullptr);

        crossover_->crossover_rate(pc);
    }

    template<Gene T>
    inline Probability GA<T>::crossover_rate() const noexcept
    {
        assert(crossover_ != nullptr);

        return crossover_->crossover_rate();
    }

    template<Gene T>
    template<typename F>
    requires mutation::MutationType<F, T> && std::is_final_v<F>
    inline void GA<T>::mutation_method(F f)
    {
        mutation_ = std::make_unique<F>(std::move(f));
    }

    template<Gene T>
    template<mutation::MutationType<T> F>
    inline void GA<T>::mutation_method(std::unique_ptr<F> f)
    {
        if (!f) GA_THROW(std::invalid_argument, "The mutation method can't be a nullptr.");

        mutation_ = std::move(f);
    }

    template<Gene T>
    inline void GA<T>::mutation_method(MutationCallable f)
    {
        mutation_ = std::make_unique<mutation::Lambda<T>>(std::move(f));
    }

    template<Gene T>
    inline mutation::Mutation<T>& GA<T>::mutation_method() &
    {
        assert(mutation_ != nullptr);

        return *mutation_;
    }

    template<Gene T>
    inline const mutation::Mutation<T>& GA<T>::mutation_method() const&
    {
        assert(mutation_ != nullptr);

        return *mutation_;
    }

    template<Gene T>
    inline void GA<T>::mutation_rate(Probability pm)
    {
        assert(mutation_ != nullptr);

        mutation_->mutation_rate(pm);
    }

    template<Gene T>
    inline Probability GA<T>::mutation_rate() const noexcept
    {
        assert(mutation_ != nullptr);

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
        return sol.fitness.size() == fitness_function_->num_objectives();
    }

    template<Gene T>
    inline bool GA<T>::hasValidChromosome(const Candidate<T>& sol) const noexcept
    {
        return variable_chrom_len() || (sol.chromosome.size() == chrom_len());
    }

    template<Gene T>
    inline bool GA<T>::fitnessMatrixIsSynced() const
    {
        return std::equal(fitness_matrix_.begin(), fitness_matrix_.end(), population_.begin(), population_.end(),
        [this](const FitnessVector& fvec, const Candidate<T>& sol)
        {
            return fvec == sol.fitness;
        });
    }

    template<Gene T>
    inline bool GA<T>::populationIsValid(const Population<T>& pop) const
    {
        return std::all_of(pop.begin(), pop.end(), [this](const Candidate<T>& sol)
        {
            return hasValidChromosome(sol) && hasValidFitness(sol);
        });
    }

    template<Gene T>
    void GA<T>::initializeAlgorithm()
    {
        assert(fitness_function_);
        assert(algorithm_ && crossover_ && mutation_ && stop_condition_);

        /* Derived GA. */
        initialize();

        /* Reset state variables in case run() has already been called before. */
        generation_cntr_ = 0;
        num_fitness_evals_ = 0;
        solutions_.clear();
        population_.clear();

        /* Create and evaluate the initial population of the algorithm. */
        population_ = generatePopulation(population_size_);
        std::for_each(GA_EXECUTION_UNSEQ, population_.begin(), population_.end(), [this](Candidate<T>& sol) { evaluate(sol); });
        fitness_matrix_ = detail::toFitnessMatrix(population_);

        assert(populationIsValid(population_));
        assert(fitnessMatrixIsSynced());

        /* Initialize the algorithm used.
         * This must be done after the initial population has been created and evaluted,
         * as it might want to use the population's fitness values (fitness_matrix_). */
        algorithm_->initialize(*this);

        is_initialized_ = true;
    }

    template<Gene T>
    Population<T> GA<T>::generatePopulation(size_t pop_size) const
    {
        assert(chrom_len() > 0);
        assert(pop_size > 0);

        Population<T> population;
        population.reserve(pop_size);

        while (pop_size--)
        {
            auto sol = generateCandidate();

            if (!hasValidChromosome(sol))
            {
                GA_THROW(std::logic_error, "A candidate with incorrect chromosome length was generated by generateCandidate().");
            }
            population.push_back(std::move(sol));
        }

        return population;
    }

    template<Gene T>
    inline void GA<T>::prepareSelections() const
    {
        assert(std::all_of(population_.begin(), population_.end(), [this](const Candidate<T>& sol) { return sol.is_evaluated && hasValidFitness(sol); }));
        assert(fitnessMatrixIsSynced());

        algorithm_->prepareSelections(*this, fitness_matrix());
    }

    template<Gene T>
    inline const Candidate<T>& GA<T>::select() const
    {
        return algorithm_->select(*this, population(), fitness_matrix());
    }

    template<Gene T>
    inline CandidatePair<T> GA<T>::crossover(const Candidate<T>& parent1, const Candidate<T>& parent2) const
    {
        return (*crossover_)(*this, parent1, parent2);
    }

    template<Gene T>
    inline void GA<T>::mutate(Candidate<T>& sol) const
    {
        (*mutation_)(*this, sol);
    }

    template<Gene T>
    void GA<T>::repair(Candidate<T>& sol) const
    {
        assert(hasValidChromosome(sol));

        /* Don't try to do anything unless a repair function is set. */
        if (!repair_) return;

        auto improved_chrom = repair_(*this, sol.chromosome);

        if (!variable_chrom_len() && improved_chrom.size() != chrom_len())
        {
            GA_THROW(std::logic_error, "The repair function returned a chromosome with incorrect length.");
        }

        /* If the repair function did something. */
        if (improved_chrom != sol.chromosome)
        {
            sol.is_evaluated = false;
            sol.chromosome = std::move(improved_chrom);
        }
    }

    template<Gene T>
    void GA<T>::updatePopulation(Population<T>&& children)
    {
        assert(fitnessMatrixIsSynced());
        assert(populationIsValid(population_));
        assert(std::all_of(population_.begin(), population_.end(), std::mem_fn(&Candidate<T>::is_evaluated)));

        population_ = algorithm_->nextPopulation(*this, std::move(population_), std::move(children));
        fitness_matrix_ = detail::toFitnessMatrix(population_);
    }

    template<Gene T>
    inline bool GA<T>::stopCondition() const
    {
        return (*stop_condition_)(*this);
    }

    template<Gene T>
    inline void GA<T>::evaluate(Candidate<T>& sol)
    {
        assert(fitness_function_);
        assert(hasValidChromosome(sol));

        /* If the fitness function is static, and the solution has already
         * been evaluted sometime earlier (in an earlier generation), there
         * is no point doing it again. */
        if (!sol.is_evaluated || fitness_function_->dynamic_fitness())
        {
            sol.fitness = (*fitness_function_)(sol.chromosome);
            sol.is_evaluated = true;

            std::atomic_ref num_evals{ num_fitness_evals_ };
            num_evals.fetch_add(1, std::memory_order_acq_rel);
        }

        assert(hasValidFitness(sol));
    }

    template<Gene T>
    void GA<T>::updateOptimalSolutions(Candidates<T>& optimal_sols, const Population<T>& pop) const
    {
        auto optimal_pop = algorithm_->optimalSolutions(*this, pop);

        optimal_sols = detail::mergeParetoSets(std::move(optimal_sols), std::move(optimal_pop));
        detail::erase_duplicates(optimal_sols);
    }

    template<Gene T>
    void GA<T>::advance()
    {
        assert(population_.size() == population_size_);

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
    inline Candidates<T> GA<T>::run()
    {
        return run(max_gen());
    }

    template<Gene T>
    Candidates<T> GA<T>::run(size_t num_generations)
    {
        max_gen(num_generations);

        initializeAlgorithm();
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