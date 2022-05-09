/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_POPULATION_HPP
#define GA_POPULATION_HPP

#include "candidate.hpp"
#include <vector>

namespace genetic_algorithm
{
    /** The Population type used in all of the genetic algorithms. */
    template<Gene T>
    using Population = std::vector<Candidate<T>>;

    /** A vector of Candidates. */
    template<Gene T>
    using Candidates = std::vector<Candidate<T>>;

} // namespace genetic_algorithm

namespace genetic_algorithm::detail
{
    using FitnessVector = std::vector<double>;
    using FitnessMatrix = std::vector<std::vector<double>>;

    /* Get the fitness vector of the population (single-objective). */
    template<Gene T>
    FitnessVector toFitnessVector(const Population<T>& pop);

    /* Get the fitness matrix of the population (multi-objective). */
    template<Gene T>
    FitnessMatrix toFitnessMatrix(const Population<T>& pop);

    /* Get the fitness vector of the population along the first objective axis from the fitness matrix. */
    FitnessVector toFitnessVector(const FitnessMatrix& pop);

    /* Return the minimum fitness value of the population along each objective axis. */
    FitnessVector populationFitnessMin(const FitnessMatrix& pop);

    /* Return the minimum fitness value of the population along each objective axis. */
    template<Gene T>
    FitnessVector populationFitnessMin(const Population<T>& pop);

    /* Return the maximum fitness value of the population along each objective axis. */
    FitnessVector populationFitnessMax(const FitnessMatrix& pop);

    /* Return the maximum fitness value of the population along each objective axis. */
    template<Gene T>
    FitnessVector populationFitnessMax(const Population<T>& pop);

    /* Return the mean fitness value of the population along each objective axis. */
    FitnessVector populationFitnessMean(const FitnessMatrix& pop);

    /* Return the mean fitness value of the population along each objective axis. */
    template<Gene T>
    FitnessVector populationFitnessMean(const Population<T>& pop);

    /* Return the standard deviation of the fitness values of the population along each objective axis. */
    FitnessVector populationFitnessSD(const FitnessMatrix& pop);

    /* Return the standard deviation of the fitness values of the population along each objective axis. */
    FitnessVector populationFitnessSD(const FitnessMatrix& pop, const FitnessVector& mean);

    /* Return the standard deviation of the fitness values of the population along each objective axis. */
    template<Gene T>
    FitnessVector populationFitnessSD(const Population<T>& pop);

    /* Return the standard deviation of the fitness values of the population along each objective axis. */
    template<Gene T>
    FitnessVector populationFitnessSD(const Population<T>& pop, const FitnessVector& mean);

    /* Finds the pareto-optimal solutions in a population. */
    template<Gene T>
    Candidates<T> findParetoFront(const Population<T>& pop);

} // namespace genetic_algorithm::detail

namespace genetic_algorithm::detail::_
{
    template<Gene T>
    Candidates<T> findParetoFront1D(const Population<T>& pop);

    template<Gene T>
    Candidates<T> findParetoFrontKung(const Population<T>& pop);

    template<Gene T>
    std::vector<size_t> findParetoFrontKungImpl(const Population<T>& pop, std::vector<size_t>::iterator first, std::vector<size_t>::iterator last);

    template<Gene T>
    Candidates<T> findParetoFrontSimple(const Population<T>& pop);

    inline bool kungCompareLess(const std::vector<double>& lhs, const std::vector<double>& rhs);

} // namespace genetic_algorithm::detail::_


/* IMPLEMENTATION */

#include "../utility/math.hpp"
#include "../utility/algorithm.hpp"
#include "../utility/functional.hpp"
#include <algorithm>
#include <iterator>
#include <cstddef>
#include <cassert>

namespace genetic_algorithm::detail
{
    template<Gene T>
    FitnessVector toFitnessVector(const Population<T>& pop)
    {
        assert(std::all_of(pop.begin(), pop.end(), [](const Candidate<T>& sol) { return sol.fitness.size() >= 1; }));

        return detail::map(pop, [](const Candidate<T>& sol) { return sol.fitness[0]; });
    }

    template<Gene T>
    FitnessMatrix toFitnessMatrix(const Population<T>& pop)
    {
        return detail::map(pop, [](const Candidate<T>& sol) { return sol.fitness; });
    }

    template<Gene T>
    FitnessVector populationFitnessMin(const Population<T>& pop)
    {
        return populationFitnessMin(toFitnessMatrix(pop));
    }

    template<Gene T>
    FitnessVector populationFitnessMax(const Population<T>& pop)
    {
        return populationFitnessMax(toFitnessMatrix(pop));
    }

    template<Gene T>
    FitnessVector populationFitnessMean(const Population<T>& pop)
    {
        return populationFitnessMean(toFitnessMatrix(pop));
    }

    template<Gene T>
    FitnessVector populationFitnessSD(const Population<T>& pop)
    {
        return populationFitnessSD(toFitnessMatrix(pop));
    }

    template<Gene T>
    FitnessVector populationFitnessSD(const Population<T>& pop, const FitnessVector& mean)
    {
        return populationFitnessSD(toFitnessMatrix(pop), mean);
    }

    template<Gene T>
    Candidates<T> findParetoFront(const Population<T>& pop)
    {
        assert(!pop.empty());
        assert(pop[0].fitness.size() > 0);
        assert(std::all_of(pop.begin(), pop.end(), [](const auto& sol) { sol.fitness.size() == pop[0].fitness.size(); }));

        size_t dim = pop[0].fitness.size();

        if (dim == 1)
        {
            return _::findParetoFront1D(pop);
        }
        else
        {
            return _::findParetoFrontSimple(pop);
        }
    }

} // namespace genetic_algorithm::detail

namespace genetic_algorithm::detail::_
{
    template<Gene T>
    Candidates<T> findParetoFront1D(const Population<T>& pop)
    {
        auto best = std::max_element(pop.begin(), pop.end(),
        [](const Candidate<T>& lhs, const Candidate<T>& rhs)
        {
            return lhs.fitness[0] < rhs.fitness[0];
        });

        return detail::find_all_v(pop.begin(), pop.end(),
        [&best](const Candidate<T>& sol)
        {
            return detail::floatIsEqual(best->fitness[0], sol.fitness[0]);
        });
    }

    template<Gene T>
    Candidates<T> findParetoFrontKung(const Population<T>& pop)
    {
        /* See: Kung et al. "On finding the maxima of a set of vectors." Journal of the ACM (JACM) 22.4 (1975): 469-476. */
        /* Doesn't work for d = 1 (single-objective optimization). */

        auto indices = detail::argsort(pop.begin(), pop.end(),
        [](const Candidate<T>& lhs, const Candidate<T>& rhs)
        {
            for (size_t i = 0; i < lhs.fitness.size(); i++)
            {
                if (lhs.fitness[i] != rhs.fitness[i])
                {
                    return lhs.fitness[i] > rhs.fitness[i];
                }
            }
            return false;
        });

        indices = _::findParetoFrontKungImpl(pop, indices.begin(), indices.end());

        return detail::map(indices, [&pop](size_t idx) { return pop[idx]; });
    }

    template<Gene T>
    std::vector<size_t> findParetoFrontKungImpl(const Population<T>& pop, std::vector<size_t>::iterator first, std::vector<size_t>::iterator last)
    {
        if (std::distance(first, last) == 1) return { *first };

        auto middle = first + std::distance(first, last) / 2;

        auto top_half    = _::findParetoFrontKungImpl(pop, first, middle);
        auto bottom_half = _::findParetoFrontKungImpl(pop, middle, last);

        for (const auto& bad : bottom_half)
        {
            bool is_dominated = false;
            for (const auto& good : top_half)
            {
                if (_::kungCompareLess(pop[bad].fitness, pop[good].fitness))
                {
                    is_dominated = true;
                    break;
                }
            }
            if (!is_dominated) top_half.push_back(bad);
        }

        return top_half;
    }

    template<Gene T>
    Candidates<T> findParetoFrontSimple(const Population<T>& pop)
    {
        auto indices = detail::argsort(pop.begin(), pop.end(),
        [](const Candidate<T>& lhs, const Candidate<T>& rhs)
        {
            for (size_t i = 0; i < lhs.fitness.size(); i++)
            {
                if (lhs.fitness[i] != rhs.fitness[i])
                {
                    return lhs.fitness[i] > rhs.fitness[i];
                }
            }
            return false;
        });

        std::vector<size_t> optimal_indices;

        for (const auto& idx : indices)
        {
            bool dominated = std::any_of(optimal_indices.begin(), optimal_indices.end(),
            [&pop, idx](size_t optimal_idx)
            {
                return detail::paretoCompareLess(pop[idx].fitness, pop[optimal_idx].fitness);
            });
            if (!dominated) optimal_indices.push_back(idx);
        }

        return detail::map(optimal_indices, [&pop](size_t idx) { return pop[idx]; });
    }

    inline bool kungCompareLess(const std::vector<double>& lhs, const std::vector<double>& rhs)
    {
        bool is_dominated = detail::paretoCompareLess(lhs, rhs, 1);
        bool is_equal = !detail::floatIsEqual(lhs[0], rhs[0]) &&
                        std::equal(lhs.begin() + 1, lhs.end(),
                                   rhs.begin() + 1,
                                   detail::floatIsEqual<double>);

        return is_dominated || is_equal;
    }

} // namespace genetic_algorithm::detail::_

#endif // !GA_POPULATION_HPP