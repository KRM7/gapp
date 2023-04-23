/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_POPULATION_HPP
#define GA_POPULATION_HPP

#include "candidate.hpp"
#include "../utility/matrix.hpp"
#include <vector>
#include <cstddef>

namespace genetic_algorithm
{
    /** The Population type used in all of the genetic algorithms. */
    template<typename T>
    using Population = std::vector<Candidate<T>>;

    /** A vector of Candidates. */
    template<typename T>
    using Candidates = std::vector<Candidate<T>>;

} // namespace genetic_algorithm

namespace genetic_algorithm::detail
{
    /* Return the fitness matrix of the population (multi-objective). */
    template<typename T>
    FitnessMatrix toFitnessMatrix(const Population<T>& pop);

    /* Return the fitness vector of a fitness matrix along the first objective axis. */
    FitnessVector toFitnessVector(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);

    /* Return the minimum fitness values of a fitness matrix along each objective axis. */
    FitnessVector minFitness(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);

    /* Return the maximum fitness values of a fitness matrix along each objective axis. */
    FitnessVector maxFitness(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);

    /* Return the mean fitness values of a fitness matrix along each objective axis. */
    FitnessVector fitnessMean(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);

    /* Return the variance of the fitness values of a fitness matrix along each objective axis. */
    FitnessVector fitnessVariance(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);

    /* Return the variance of the fitness values of a fitness matrix along each objective axis. */
    FitnessVector fitnessVariance(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last, const FitnessVector& fitness_mean);

    /* Return the standard deviation of the fitness values of a fitness matrix along each objective axis. */
    FitnessVector fitnessStdDev(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);

    /* Return the standard deviation of the fitness values of a fitness matrix along each objective axis. */
    FitnessVector fitnessStdDev(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last, const FitnessVector& fitness_mean);


    /* Find the pareto-optimal solutions in a population. */
    template<typename T>
    Candidates<T> findParetoFront(const Population<T>& pop);

    std::vector<size_t> findParetoFront1D(const FitnessMatrix& fmat);

    std::vector<size_t> findParetoFrontSort(const FitnessMatrix& fmat);

    std::vector<size_t> findParetoFrontBest(const FitnessMatrix& fmat);

    std::vector<size_t> findParetoFrontKung(const FitnessMatrix& fmat);


    /* Find the pareto-optimal solutions in the set (lhs U rhs), assuming both lhs and rhs are pareto sets. */
    template<typename T>
    Candidates<T> mergeParetoSets(Candidates<T> lhs, Candidates<T> rhs);

} // namespace genetic_algorithm::detail


/* IMPLEMENTATION */

#include "../utility/algorithm.hpp"
#include "../utility/functional.hpp"
#include "../utility/iterators.hpp"
#include "../utility/utility.hpp"
#include "../utility/math.hpp"
#include <algorithm>
#include <functional>
#include <atomic>

namespace genetic_algorithm::detail
{
    template<typename T>
    FitnessMatrix toFitnessMatrix(const Population<T>& pop)
    {
        if (pop.empty()) return {};

        FitnessMatrix fitness_matrix;
        fitness_matrix.reserve(pop.size(), pop[0].fitness.size());

        for (const Candidate<T>& sol : pop)
        {
            fitness_matrix.append_row(sol.fitness);
        }

        return fitness_matrix;
    }

    template<typename T>
    Candidates<T> findParetoFront(const Population<T>& pop)
    {
        if (pop.empty()) return {};

        GA_ASSERT(!pop[0].fitness.empty());
        GA_ASSERT(std::all_of(pop.begin(), pop.end(), [&](const Candidate<T>& sol) { return sol.fitness.size() == pop[0].fitness.size(); }));

        auto fitness_matrix = detail::toFitnessMatrix(pop);

        auto optimal_indices = fitness_matrix.ncols() == 1 ?
            findParetoFront1D(fitness_matrix) :
            findParetoFrontSort(fitness_matrix);

        return detail::select(pop, optimal_indices);
    }
    
    template<typename T>
    Candidates<T> mergeParetoSets(Candidates<T> lhs, Candidates<T> rhs)
    {
        if (lhs.empty()) return rhs;
        if (rhs.empty()) return lhs;

        if (rhs.size() > lhs.size()) std::swap(lhs, rhs);

        enum Dominance : unsigned char { UNKNOWN = 0, OPTIMAL = 1, DOMINATED = 2 };

        std::vector<Dominance> lhs_state(lhs.size());
        std::vector<std::atomic<Dominance>> rhs_state(rhs.size());

        std::for_each(GA_EXECUTION_UNSEQ, iota_iterator(0_sz), iota_iterator(lhs.size()), [&](size_t i) noexcept
        {
            for (size_t j = 0; j < rhs.size(); j++)
            {
                if (lhs_state[i] == DOMINATED) continue;
                if (rhs_state[j] == DOMINATED) continue;

                if (lhs_state[i] == OPTIMAL)
                {
                    if (math::paretoCompareLess(rhs[j].fitness, lhs[i].fitness))
                    {
                        rhs_state[j] = DOMINATED;
                    }
                    continue;
                }
                if (rhs_state[j] == OPTIMAL)
                {
                    if (math::paretoCompareLess(lhs[i].fitness, rhs[j].fitness))
                    {
                        lhs_state[i] = DOMINATED;
                    }
                    continue;
                }

                const auto comp = math::paretoCompare(lhs[i].fitness, rhs[j].fitness);
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
                // comp == 0 --> both are OPTIMAL or DOMINATED, can't know
            }
        });

        Candidates<T> optimal_solutions;
        optimal_solutions.reserve(lhs.size() + rhs.size());

        for (size_t i = 0; i < rhs.size(); i++)
        {
            if (rhs_state[i] != DOMINATED) optimal_solutions.push_back(std::move(rhs[i]));
        }
        for (size_t i = 0; i < lhs.size(); i++)
        {
            if (lhs_state[i] != DOMINATED) optimal_solutions.push_back(std::move(lhs[i]));
        }

        return optimal_solutions;
    }

} // namespace genetic_algorithm::detail

#endif // !GA_POPULATION_HPP