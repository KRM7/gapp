/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CORE_GA_BASE_IMPL_HPP
#define GA_CORE_GA_BASE_IMPL_HPP

#include "ga_base.decl.hpp"
#include "../population/population.hpp"
#include "../selection/selection.hpp"
#include "../crossover/crossover_base.hpp"
#include "../crossover/lambda.hpp"
#include "../mutation/mutation_base.hpp"
#include "../mutation/lambda.hpp"
#include "../utility/rng.hpp"
#include "../utility/math.hpp"
#include "../utility/utility.hpp"
#include "../utility/algorithm.hpp"
#include <execution>
#include <numeric>
#include <limits>
#include <stdexcept>
#include <cstdlib>
#include <cassert>
#include <cmath>

namespace genetic_algorithm
{
    template<Gene T, typename D>
    GA<T, D>::GA(size_t chrom_len, FitnessFunction fitness_function)
        : GaInfo(chrom_len)
    {
        if (!fitness_function)
        {
            throw std::invalid_argument("The fitness function is requires for the GA.");
        }

        fitness_function_ = std::move(fitness_function);
    }

    template<Gene T, typename D>
    GA<T, D>::GA(size_t population_size, size_t chrom_len, FitnessFunction fitness_function)
        : GaInfo(population_size, chrom_len)
    {
        if (!fitness_function)
        {
            throw std::invalid_argument("The fitness function is requires for the GA.");
        }

        fitness_function_ = std::move(fitness_function);
    }

    template<Gene T, typename D>
    D& GA<T, D>::derived() noexcept
    {
        return static_cast<D&>(*this);
    }

    template<Gene T, typename D>
    const D& GA<T, D>::derived() const noexcept
    {
        return static_cast<const D&>(*this);
    }

    template<Gene T, typename D>
    auto GA<T, D>::solutions() const noexcept -> Population
    {
        return solutions_;
    }

    template<Gene T, typename D>
    auto GA<T, D>::population() const noexcept -> Population
    {
        return population_;
    }

    template<Gene T, typename D>
    void GA<T, D>::fitness_function(FitnessFunction f)
    {
        if (!f)
        {
            throw std::invalid_argument("The fitness function is requires for the GA.");
        }

        fitness_function_ = std::move(f);

        num_objectives(findNumObjectives(fitness_function_));
    }

    template<Gene T, typename D>
    template<crossover::CrossoverMethod<T> F>
    void GA<T, D>::crossover_method(const F& f)
    {
        crossover_ = std::make_unique<F>(f);
    }

    template<Gene T, typename D>
    template<crossover::CrossoverMethod<T> F>
    void GA<T, D>::crossover_method(std::unique_ptr<F>&& f)
    {
        if (!f)
        {
            throw std::invalid_argument("The crossover method can't be a nullptr.");
        }

        crossover_ = std::move(f);
    }

    template<Gene T, typename D>
    void GA<T, D>::crossover_method(CrossoverFunction f)
    {
        if (!f)
        {
            throw std::invalid_argument("The crossover function can't be empty.");
        }

        crossover_ = std::make_unique<crossover::dtl::Lambda<T>>(std::move(f));
    }

    template<Gene T, typename D>
    template<crossover::CrossoverMethod<T> F>
    F& GA<T, D>::crossover_method() &
    {
        return dynamic_cast<F&>(*crossover_);
    }

    template<Gene T, typename D>
    void GA<T, D>::crossover_rate(double pc)
    {
        crossover_->crossover_rate(pc);
    }

    template<Gene T, typename D>
    double GA<T, D>::crossover_rate() const noexcept
    {
        return crossover_->crossover_rate();
    }

    template<Gene T, typename D>
    template<mutation::MutationMethod<T> F>
    void GA<T, D>::mutation_method(const F& f)
    {
        mutation_ = std::make_unique<F>(f);
    }

    template<Gene T, typename D>
    template<mutation::MutationMethod<T> F>
    void GA<T, D>::mutation_method(std::unique_ptr<F>&& f)
    {
        if (!f)
        {
            throw std::invalid_argument("The mutation method can't be a nullptr.");
        }

        mutation_ = std::move(f);
    }

    template<Gene T, typename D>
    void GA<T, D>::mutation_method(MutationFunction f)
    {
        if (!f)
        {
            throw std::invalid_argument("Thhe mutation method can't be empty.");
        }

        mutation_ = std::make_unique<mutation::dtl::Lambda<T>>(std::move(f));
    }

    template<Gene T, typename D>
    template<mutation::MutationMethod<T> F>
    F& GA<T, D>::mutation_method() &
    {
        return dynamic_cast<F&>(*mutation_);
    }

    template<Gene T, typename D>
    void GA<T, D>::mutation_rate(double pm)
    {
        mutation_->mutation_rate(pm);
    }

    template<Gene T, typename D>
    double GA<T, D>::mutation_rate() const noexcept
    {
        return mutation_->mutation_rate();
    }

    template<Gene T, typename D>
    void GA<T, D>::repair_function(const RepairFunction& f)
    {
        repair_ = f;
    }

    template<Gene T, typename D>
    void GA<T, D>::advance()
    {
        if (keep_all_optimal_solutions)
        {
            updateOptimalSolutions(solutions_, population_);
        }

        size_t num_children = population_size_ + population_size_ % 2;
        std::vector<CandidatePair> child_pairs(num_children / 2);

        (*selection_).prepare(*this, fitness_matrix());
        std::generate(GA_EXECUTION_UNSEQ, child_pairs.begin(), child_pairs.end(),
        [this]
        {
            const auto& parent1 = selectCandidate();
            const auto& parent2 = selectCandidate();

            return crossover(parent1, parent2);
        });

        auto children = detail::flatten(std::move(child_pairs));

        std::for_each(GA_EXECUTION_UNSEQ, children.begin(), children.end(),
        [this](Candidate& child)
        {
            mutate(child);
            repair(child);
        });

        updatePopulation(population_, std::move(children));

        if (endOfGenerationCallback != nullptr) endOfGenerationCallback(*this);
        generation_cntr_++;
    }

    template<Gene T, typename D>
    auto GA<T, D>::run(size_t num_generations) -> Candidates
    {
        max_gen(num_generations);

        initializeAlgorithm();
        while (!stopCondition())
        {
            advance();
        }
        updateOptimalSolutions(solutions_, population_);

        can_continue_ = true;

        return solutions_;
    }

    template<Gene T, typename D>
    auto GA<T, D>::continueFor(size_t num_generations) -> Candidates
    {
        if (!can_continue_) { return run(num_generations); }

        max_gen(max_gen_ + num_generations);

        while (!stopCondition())
        {
            advance();
        }
        updateOptimalSolutions(solutions_, population_);

        return solutions_;
    }

    template<Gene T, typename D>
    void GA<T, D>::initializeAlgorithm()
    {
        /* The number of objectives is determined from the return value of the fitness func,
        *  this assumes that the returned vector will always be the same size during a run. 
        *
        *  This also needs generateCandidate() to create a dummy solution to pass to the fitness
        *  function, which might only return valid values after the derived class contructor set
        *  up some stuff for the candidate generation (eg. bounds), so this can't be called earlier,
        *  eg. in the base ctor. */
        num_objectives(findNumObjectives(fitness_function_));

        /* Reset state variables just in case GA.run() has already been called before. */
        can_continue_ = false;
        generation_cntr_ = 0;
        num_fitness_evals_ = 0;
        solutions_.clear();
        population_.clear();

        /* Create and evaluate the initial population of the algorithm. */
        population_ = generatePopulation(population_size_);
        fitness_matrix_ = evaluatePopulation(population_);

        /* Initialize the selection method. 
         * This must be done after the initial population has been created and evaluted,
         * as it might rely on the population's fitness values (fitness_matrix). */
        (*selection_).init(*this);
    }
    
    template<Gene T, typename D>
    size_t GA<T, D>::findNumObjectives(const FitnessFunction& f) const
    {
        Candidate c = generateCandidate();
        return f(c.chromosome).size();
    }

    template<Gene T, typename D>
    void GA<T, D>::setDefaultAlgorithm()
    {
        num_objectives(findNumObjectives(fitness_function_));

        if (num_objectives_ == 1)
        {
            selection_method(std::make_unique<selection::single_objective::Tournament>());
        }
        else
        {
            selection_method(std::make_unique<selection::multi_objective::NSGA3>());
        }
    }

    template<Gene T, typename D>
    auto GA<T, D>::generateCandidate() const -> Candidate
    {
        return derived().generateCandidate();
    }

    template<Gene T, typename D>
    auto GA<T, D>::generatePopulation(size_t pop_size) const -> Population
    {
        assert(pop_size > 0);

        Population pop;
        pop.reserve(pop_size);

        while (pop_size--)
        {
            pop.push_back(generateCandidate());
        }

        return pop;
    }

    template<Gene T, typename D>
    auto GA<T, D>::selectCandidate() const -> const Candidate&
    {
        size_t idx = (*selection_).select(*this, fitness_matrix());

        return population_[idx];
    }

    template<Gene T, typename D>
    auto GA<T, D>::crossover(const Candidate& parent1, const Candidate& parent2) const -> CandidatePair
    {
        return (*crossover_)(*this, parent1, parent2);
    }

    template<Gene T, typename D>
    void GA<T, D>::mutate(Candidate& sol) const
    {
        return (*mutation_)(*this, sol);
    }

    template<Gene T, typename D>
    auto GA<T, D>::evaluatePopulation(Population& pop) -> FitnessMatrix
    {
        assert(fitness_function_ != nullptr);

        std::for_each(GA_EXECUTION_UNSEQ, pop.begin(), pop.end(),
        [this](Candidate& sol)
        {
            if (dynamic_fitness || !sol.is_evaluated)
            {
                sol.fitness = fitness_function_(sol.chromosome);
                sol.is_evaluated = true;

                num_fitness_evals_++;
            }
        });

        for (const auto& sol : pop)
        {
            if (sol.fitness.size() != num_objectives_)
            {
                throw std::domain_error("A fitness vector returned by the fitness function has incorrect size.");
            }
            if (!std::all_of(sol.fitness.begin(), sol.fitness.end(), std::isfinite<double>))
            {
                throw std::domain_error("A non-finite fitness value was returned by the fitness function.");
            }
        }

        return detail::toFitnessMatrix(pop);
    }

    template<Gene T, typename D>
    void GA<T, D>::updateOptimalSolutions(Candidates& optimal_sols, const Population& pop) const
    {
        assert(std::all_of(pop.begin(), pop.end(), [](const Candidate& sol) { return sol.is_evaluated; }));

        optimal_sols.insert(optimal_sols.end(), pop.begin(), pop.end());

        optimal_sols = detail::findParetoFront(optimal_sols);

        detail::erase_duplicates(optimal_sols);
    }

    template<Gene T, typename D>
    void GA<T, D>::repair(Candidate& sol) const
    {
        /* Don't do anything unless a repair function is specified. */
        if (!repair_) return;

        auto improved_chrom = repair_(sol.chromosome); // assert=chrom_len
        if (improved_chrom != sol.chromosome)
        {
            sol.is_evaluated = false;
            sol.chromosome = std::move(improved_chrom);
        }
    }

    template<Gene T, typename D>
    void GA<T, D>::updatePopulation(Population& current_pop, Population&& children)
    {
        auto child_fmat = evaluatePopulation(children);

        current_pop.insert(current_pop.end(),
                           std::make_move_iterator(children.begin()),
                           std::make_move_iterator(children.end()));
        fitness_matrix_.insert(fitness_matrix_.end(),
                               std::make_move_iterator(child_fmat.begin()),
                               std::make_move_iterator(child_fmat.end()));

        auto next_indices = (*selection_).nextPopulation(*this, fitness_matrix_);

        current_pop = detail::map(next_indices, [&current_pop](size_t idx) { return current_pop[idx]; });
        fitness_matrix_ = detail::map(next_indices, [this](size_t idx) { return fitness_matrix_[idx]; });
    }

    template<Gene T, typename D>
    bool GA<T, D>::stopCondition() const
    {
        return (*stop_condition_)(*this);
    }

} // namespace genetic_algorithm

#endif // !GA_CORE_GA_BASE_IMPL_HPP