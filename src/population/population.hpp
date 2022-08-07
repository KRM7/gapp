/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_POPULATION_HPP
#define GA_POPULATION_HPP

#include "candidate.hpp"
#include "../utility/utility.hpp"
#include <vector>
#include <atomic>

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

    /* Return the fitness vector of the population (single-objective). */
    template<Gene T>
    FitnessVector toFitnessVector(const Population<T>& pop);

    /* Return the fitness matrix of the population (multi-objective). */
    template<Gene T>
    FitnessMatrix toFitnessMatrix(const Population<T>& pop);

    /* Return the fitness vector of a fitness matrix along the first objective axis. */
    FitnessVector toFitnessVector(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);

    /* Return the minimum fitness values of a fitness matrix along each objective axis. */
    FitnessVector minFitness(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);

    /* Return the maximum fitness values of a fitness matrix along each objective axis. */
    FitnessVector maxFitness(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);

    /* Return the mean fitness values of a fitness matrix along each objective axis. */
    FitnessVector fitnessMean(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);

    /* Return the standard deviation of the fitness values of a fitness matrix along each objective axis. */
    FitnessVector fitnessStdDev(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);

    /* Return the standard deviation of the fitness values of a fitness matrix along each objective axis. */
    FitnessVector fitnessStdDev(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last, const FitnessVector& mean);

    /* Find the pareto-optimal solutions in a population. */
    template<Gene T>
    Candidates<T> findParetoFront(const Population<T>& pop);

    /* Find the pareto-optimal solutions in the set (lhs U rhs), assuming both lhs and rhs are pareto sets. */
    template<Gene T>
    Candidates<T> mergeParetoSets(Candidates<T>&& lhs, Candidates<T>&& rhs);

} // namespace genetic_algorithm::detail

namespace genetic_algorithm::detail::_
{
    std::vector<size_t> findParetoFront1D(const FitnessMatrix& fmat);

    std::vector<size_t> findParetoFrontSort(const FitnessMatrix& fmat);

    std::vector<size_t> findParetoFrontBest(const FitnessMatrix& fmat);

    std::vector<size_t> findParetoFrontKung(const FitnessMatrix& fmat);

} // namespace genetic_algorithm::detail::_


/* IMPLEMENTATION */

#include "../utility/functional.hpp"
#include <algorithm>
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
        return detail::map(pop, &Candidate<T>::fitness);
    }

    template<Gene T>
    Candidates<T> findParetoFront(const Population<T>& pop)
    {
        assert(!pop.empty());
        assert(pop[0].fitness.size() > 0);
        assert(std::all_of(pop.begin(), pop.end(), [&pop](const auto& sol) { return sol.fitness.size() == pop[0].fitness.size(); }));

        auto fitness_matrix = detail::toFitnessMatrix(pop);

        auto optimal_indices = fitness_matrix[0].size() == 1 ?
            _::findParetoFront1D(fitness_matrix) :
            _::findParetoFrontSort(fitness_matrix);

        return detail::map(optimal_indices, [&pop](size_t idx) { return pop[idx]; });
    }
    
    template<Gene T>
    Candidates<T> mergeParetoSets(Candidates<T>&& A, Candidates<T>&& B)
    {
        if (A.empty()) return B;
        if (B.empty()) return A;

        auto [lhs, rhs] = (A.size() > B.size()) ? std::tie(A, B) : std::tie(B, A);

        enum Dominance : unsigned char { UNKNOWN = 0, OPTIMAL = 1, DOMINATED = 2 };

        std::vector<Dominance> lhs_state(lhs.size());
        std::vector<std::atomic<Dominance>> rhs_state(rhs.size());

        std::vector<size_t> lhs_indices(lhs.size());
        std::iota(lhs_indices.begin(), lhs_indices.end(), 0_sz);

        std::for_each(GA_EXECUTION_UNSEQ, lhs_indices.begin(), lhs_indices.end(), [&](size_t i) noexcept
        {
            for (size_t j = 0; j < rhs.size(); j++)
            {
                if (rhs_state[j] == DOMINATED) continue;

                if (rhs_state[j] == OPTIMAL && detail::paretoCompareLess(lhs[i].fitness, rhs[j].fitness))
                {
                    lhs_state[i] = DOMINATED;
                    continue;
                }

                if (lhs_state[i] == DOMINATED) continue;

                if (lhs_state[i] == OPTIMAL && detail::paretoCompareLess(rhs[j].fitness, lhs[i].fitness))
                {
                    rhs_state[j] = DOMINATED;
                    continue;
                }

                auto comp = detail::paretoCompare(lhs[i].fitness, rhs[j].fitness);
                if (comp < 0)
                {
                    lhs_state[i] = DOMINATED;
                    rhs_state[j] = OPTIMAL;
                }
                else if (comp > 0)
                {
                    lhs_state[i] = OPTIMAL;
                    rhs_state[j] = DOMINATED;
                }
            }
        });

        Candidates<T> optimal_solutions;
        optimal_solutions.reserve(lhs.size() + rhs.size());

        for (size_t i = 0; i < lhs.size(); i++)
        {
            if (lhs_state[i] != DOMINATED) optimal_solutions.emplace_back(std::move(lhs[i]));
        }
        for (size_t i = 0; i < rhs.size(); i++)
        {
            if (rhs_state[i] != DOMINATED) optimal_solutions.emplace_back(std::move(rhs[i]));
        }

        return optimal_solutions;
    }

} // namespace genetic_algorithm::detail

#endif // !GA_POPULATION_HPP