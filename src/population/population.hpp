/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_POPULATION_HPP
#define GA_POPULATION_HPP

#include "candidate.hpp"
#include <vector>
#include <cstddef>

namespace gapp
{
    /** The population type used in all of the algorithms. */
    template<typename Gene>
    using Population = std::vector<Candidate<Gene>>;

    /** A vector of candidates, same as the population type. */
    template<typename Gene>
    using Candidates = std::vector<Candidate<Gene>>;

} // namespace gapp

namespace gapp::detail
{
    /* Return the fitness matrix of the population (multi-objective). */
    template<typename T>
    FitnessMatrix toFitnessMatrix(const Population<T>& pop);

    /* Return the fitness vector of a fitness matrix along the first objective axis. */
    FitnessVector toFitnessVector(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);

    /* Find the pareto-optimal solutions in a population. Assumes fitness maximization, duplicates are not eliminated. */
    template<typename T>
    Candidates<T> findParetoFront(const Population<T>& pop);

    std::vector<size_t> findParetoFront(const FitnessMatrix& fmat);

    std::vector<size_t> findParetoFront1D(const FitnessMatrix& fmat);

    std::vector<size_t> findParetoFrontSort(const FitnessMatrix& fmat);

    std::vector<size_t> findParetoFrontBest(const FitnessMatrix& fmat);

    std::vector<size_t> findParetoFrontKung(const FitnessMatrix& fmat);

    /* Find the pareto-optimal solutions in the set (lhs U rhs), assuming both lhs and rhs are pareto sets. */
    template<typename T>
    Candidates<T> mergeParetoSets(Candidates<T> lhs, Candidates<T> rhs);

    /* Find the nadir point of a fitness matrix assuming fitness maximization. */
    FitnessVector findNadirPoint(const FitnessMatrix& fitness_matrix);

} // namespace gapp::detail


/* IMPLEMENTATION */

#include "../utility/algorithm.hpp"
#include "../utility/functional.hpp"
#include "../utility/iterators.hpp"
#include "../utility/utility.hpp"
#include "../utility/math.hpp"
#include <algorithm>
#include <functional>
#include <atomic>

namespace gapp::detail
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

        auto optimal_indices = detail::findParetoFront(fitness_matrix);

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

} // namespace gapp::detail

#endif // !GA_POPULATION_HPP