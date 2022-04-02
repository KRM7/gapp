/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_GA_BASE_IMPL_HPP
#define GA_GA_BASE_IMPL_HPP

#include "ga_base.decl.hpp"
#include "../selection/selection_base.hpp"
#include "../crossover/crossover_base.hpp"
#include "../mutation/mutation_base.hpp"
#include "../stop_condition/stop_condition_base.hpp"
#include "../utility/rng.hpp"
#include "../utility/math.hpp"
#include "../population/population.hpp"
#include "../utility/utils.hpp"

#include <execution>
#include <numeric>
#include <limits>
#include <stdexcept>
#include <cstdlib>
#include <cassert>
#include <cmath>

namespace genetic_algorithm
{
    template<typename T>
    GA<T>::GA(size_t chrom_len, FitnessFunction fitness_function)
        : GaBase(chrom_len)
    {
        if (!fitness_function)
        {
            throw std::invalid_argument("The fitness function is required for the GA.");
        }
        fitness_function_ = fitness_function;
    }

    template<typename T>
    typename const GA<T>::CandidateVec& GA<T>::solutions() const
    {
        return solutions_;
    }

    template<typename T>
    typename const GA<T>::Population& GA<T>::population() const
    {
        return population_;
    }

    template<typename T>
    std::vector<double> GA<T>::fitness_vec() const
    {
        return detail::fitnessVector(population_);
    }

    template<typename T>
    std::vector<std::vector<double>> GA<T>::fitness_matrix() const
    {
        return detail::fitnessMatrix(population_);
    }

    template<typename geneType>
    void GA<geneType>::initial_population(const Population& pop)
    {
        if (!std::all_of(pop.begin(), pop.end(), [this](const Candidate& c) { return c.chromosome.size() == chrom_len_; }))
        {
            throw std::invalid_argument("The length of each chromosome in the preset pop must be equal to chrom_len.");
        }

        initial_population_preset_ = pop;
    }

    template<typename geneType>
    void GA<geneType>::fitness_function(FitnessFunction f)
    {
        if (!f)
        {
            throw std::invalid_argument("The fitness function is requires for the GA.");
        }

        fitness_function_ = f;
    }

    /* SELECTION METHOD */
    template<typename GeneType>
    template<typename SelectionType>
    //requires std::derived_from<SelectionType, selection::Selection<GeneType>> && std::copy_constructible<SelectionType>
    void GA<GeneType>::selection_method(const SelectionType& f)
    {
        selection_ = std::make_unique<SelectionType>(f);
    }

    template<typename GeneType>
    void GA<GeneType>::selection_method(std::unique_ptr<selection::Selection<GeneType>>&& f)
    {
        if (f == nullptr) throw std::invalid_argument("The selection method can't be a nullptr.");

        selection_ = std::move(f);
    }

    template<typename GeneType>
    template<typename SelectionType>
    SelectionType& GA<GeneType>::selection_method() const
    {
        return dynamic_cast<SelectionType&>(*selection_);
    }

    template<typename GeneType>
    template<typename CrossoverType>
    //requires std::derived_from<CrossoverType, crossover::Crossover<GeneType>> && std::copy_constructible<CrossoverType>
    void GA<GeneType>::crossover_method(const CrossoverType& f)
    {
        crossover_ = std::make_unique<CrossoverType>(f);
    }

    template<typename GeneType>
    void GA<GeneType>::crossover_method(std::unique_ptr<crossover::Crossover<GeneType>>&& f)
    {
        if (f == nullptr) throw std::invalid_argument("The crossover method can't be a nullptr.");

        crossover_ = std::move(f);
    }

    template<typename GeneType>
    template<typename CrossoverType>
    //requires std::derived_from<CrossoverType, crossover::Crossover<GeneType>>
    CrossoverType& GA<GeneType>::crossover_method() const
    {
        return dynamic_cast<CrossoverType&>(*crossover_);
    }

    template<typename GeneType>
    template<typename MutationType>
    //requires std::derived_from<MutationType, mutation::Mutation<GeneType>> && std::copy_constructible<MutationType>
    void GA<GeneType>::mutation_method(const MutationType& f)
    {
        mutation_ = std::make_unique<MutationType>(f);
    }

    template<typename GeneType>
    void GA<GeneType>::mutation_method(std::unique_ptr<mutation::Mutation<GeneType>>&& f)
    {
        if (f == nullptr) throw std::invalid_argument("The mutation method can't be a nullptr.");

        mutation_ = std::move(f);
    }

    template<typename GeneType>
    template<typename MutationType>
    //requires std::derived_from<MutationType, mutation::Crossover<GeneType>>
    MutationType& GA<GeneType>::mutation_method() const
    {
        return dynamic_cast<MutationType&>(*mutation_);
    }

    template<typename GeneType>
    template<typename StopType>
    //requires std::derived_from<StopType, stopping::StopCondition<GA<GeneType>>> && std::copy_constructible<StopType>
    void GA<GeneType>::stop_condition(const StopType& f)
    {
        stop_condition_ = std::make_unique<StopType>(f);
    }

    template<typename GeneType>
    void GA<GeneType>::stop_condition(std::unique_ptr<stopping::StopCondition<GeneType>>&& f)
    {
        stop_condition_ = std::move(f);
    }

    template<typename GeneType>
    template<typename StopType>
    //requires std::derived_from<StopType, stopping::StopCondition<GeneType>>
    StopType& GA<GeneType>::stop_condition() const
    {
        return dynamic_cast<StopType&>(*stop_condition_);
    }

    template<typename GeneType>
    bool GA<GeneType>::stopCondition()
    {
        if (stop_condition_)
        {
            return (generation_cntr_ >= (max_gen_ - 1)) || std::invoke(*stop_condition_, *this);;
        }
        else
        {
            return (generation_cntr_ >= (max_gen_ - 1));
        }
    }

    template<typename geneType>
    typename const GA<geneType>::CandidateVec& GA<geneType>::run()
    {
        using namespace std;

        init();

        /* Create and evaluate the initial population. */
        population_ = generateInitialPopulation();
        evaluate(population_);

        (*selection_).init(*this);

        /* Other generations. */
        size_t num_children = population_size_ + population_size_ % 2;
        while (!stopCondition())
        {
            vector<CandidatePair> parent_pairs(num_children / 2);

            (*selection_).prepare(*this, population_);
            if (archive_optimal_solutions) updateOptimalSolutions(solutions_, population_);

            /* Selections. */
            generate(GA_EXECUTION_UNSEQ, parent_pairs.begin(), parent_pairs.end(),
            [this]() -> CandidatePair
            {
                return make_pair((*selection_).select(*this, population_),
                                 (*selection_).select(*this, population_));
            });

            /* Crossovers. */
            for_each(GA_EXECUTION_UNSEQ, parent_pairs.begin(), parent_pairs.end(),
            [this](CandidatePair& p) -> void
            {
                p = (*crossover_)(*this, p.first, p.second);
            });

            vector<Candidate> children;
            children.reserve(num_children);
            for (size_t i = 0; i < parent_pairs.size(); i++)
            {
                children.push_back(move(parent_pairs[i].first));
                children.push_back(move(parent_pairs[i].second));
            }

            /* Mutations. */
            for_each(GA_EXECUTION_UNSEQ, children.begin(), children.end(),
            [this](Candidate& c) -> void
            {
                (*mutation_)(*this, c);
            });

            /* Apply repair function to the children if set. */
            repair(children);

            /* Overwrite the current population with the children. */
            evaluate(children);
            population_ = (*selection_).nextPopulation(*this, population_, children);

            if (endOfGenerationCallback != nullptr) endOfGenerationCallback(*this);
            generation_cntr_++;
        }
        updateOptimalSolutions(solutions_, population_);

        return solutions_;
    }

    template<typename geneType>
    void GA<geneType>::init()
    {
        Candidate temp = generateCandidate();
        temp.fitness = fitness_function_(temp.chromosome);
        num_objectives_ = temp.fitness.size();

        generation_cntr_ = 0;
        num_fitness_evals_ = 0;
        solutions_.clear();
        population_.clear();
    }

    template<typename geneType>
    typename GA<geneType>::Population GA<geneType>::generateInitialPopulation() const
    {
        assert(population_size_ > 0);

        if (!std::all_of(initial_population_preset_.begin(), initial_population_preset_.end(),
            [this](const Candidate& sol)
            {
                return sol.chromosome.size() == chrom_len_;
            }))
        {
            throw std::domain_error("The chromosome lengths in the preset initial population must be equal to the chrom_len set.");
        }

        Population pop;
        pop.reserve(population_size_);

        for (size_t i = 0; i < std::min(population_size_, initial_population_preset_.size()); i++)
        {
            pop.push_back(initial_population_preset_[i]);
        }
        while (pop.size() < population_size_)
        {
            pop.push_back(generateCandidate());
        }

        return pop;
    }

    template<typename geneType>
    void GA<geneType>::evaluate(Population& pop)
    {
        assert(fitness_function_ != nullptr);

        std::for_each(GA_EXECUTION_UNSEQ, pop.begin(), pop.end(),
        [this](Candidate& sol)
        {
            if (changing_fitness_func || !sol.is_evaluated)
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
            if (!std::all_of(sol.fitness.begin(), sol.fitness.end(), [](double val) { return std::isfinite(val); }))
            {
                throw std::domain_error("A non-finite fitness value was returned by the fitness function.");
            }
        }
    }

    template<typename geneType>
    void GA<geneType>::updateOptimalSolutions(CandidateVec& optimal_sols, const Population& pop) const
    {
        assert(std::all_of(pop.begin(), pop.end(), [](const Candidate& sol) { return sol.is_evaluated; }));

        optimal_sols.insert(optimal_sols.end(), pop.begin(), pop.end());
        if (num_objectives_ == 1)
        {
            optimal_sols = detail::findParetoFront1D(optimal_sols);
        }
        else
        {
            optimal_sols = detail::findParetoFrontKung(optimal_sols);
        }

        /* Remove duplicate solutions. */
        CandidateSet unique_sols;
        for (const auto& sol : optimal_sols) unique_sols.insert(std::move(sol));    /* Insertion is faster than ctor. */
        optimal_sols.assign(unique_sols.begin(), unique_sols.end());
    }

    template<typename geneType>
    void GA<geneType>::repair(Population& pop) const
    {
        /* This function doesn't do anything unless a repair function is specified. */
        if (repairFunction == nullptr) return;

        std::for_each(GA_EXECUTION_UNSEQ, pop.begin(), pop.end(),
        [this](Candidate& sol)
        {
            Chromosome improved_chrom = repairFunction(sol.chromosome);
            if (improved_chrom != sol.chromosome)
            {
                sol.is_evaluated = false;
                sol.chromosome = std::move(improved_chrom);
            }
        });

        for (const auto& sol : pop)
        {
            if (sol.chromosome.size() != chrom_len_)
            {
                throw std::domain_error("The repair function must return chromosomes of chrom_len length.");
            }
        }
    }

} // namespace genetic_algorithm

#endif // !GA_GA_BASE_IMPL_HPP