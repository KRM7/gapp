/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_CORE_POPULATION_HPP
#define GAPP_CORE_POPULATION_HPP

#include "candidate.hpp"
#include "../utility/small_vector.hpp"
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

    small_vector<size_t> findParetoFront(const FitnessMatrix& fmat);

    small_vector<size_t> findParetoFront1D(const FitnessMatrix& fmat);
    small_vector<size_t> findParetoFrontSort(const FitnessMatrix& fmat);
    small_vector<size_t> findParetoFrontBest(const FitnessMatrix& fmat);
    small_vector<size_t> findParetoFrontKung(const FitnessMatrix& fmat);

    /* Find the pareto-optimal solutions in the set (lhs U rhs), assuming both lhs and rhs are pareto sets. */
    template<typename T>
    Candidates<T> mergeParetoSets(Candidates<T> lhs, Candidates<T> rhs);

    /* Find the nadir point of a fitness matrix assuming fitness maximization. */
    FitnessVector findNadirPoint(const FitnessMatrix& fitness_matrix);

    /* Find the nadir point of a pareto front assuming fitness maximization. */
    FitnessVector findFrontNadirPoint(const FitnessMatrix& optimal_points);

} // namespace gapp::detail


/* IMPLEMENTATION */

#include "../utility/algorithm.hpp"
#include "../utility/functional.hpp"
#include "../utility/iterators.hpp"
#include "../utility/thread_pool.hpp"
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

        GAPP_ASSERT(!pop[0].fitness.empty());
        GAPP_ASSERT(std::all_of(pop.begin(), pop.end(), [&](const Candidate<T>& sol) { return sol.fitness.size() == pop[0].fitness.size(); }));

        auto fitness_matrix = detail::toFitnessMatrix(pop);
        auto optimal_indices = detail::findParetoFront(fitness_matrix);

        return detail::select(pop, optimal_indices);
    }

    enum class ParetoDominance : char { UNKNOWN = 0, OPTIMAL = 1, DOMINATED = 2 };

    template<typename T>
    Candidates<T> mergeParetoSets(Candidates<T> lhs, Candidates<T> rhs)
    {
        using enum ParetoDominance;

        if (lhs.empty()) return rhs;
        if (rhs.empty()) return lhs;

        if (rhs.size() > lhs.size()) std::swap(lhs, rhs);

        std::vector<ParetoDominance> lhs_state(lhs.size());
        std::vector<ParetoDominance> rhs_state(rhs.size());

        detail::parallel_for(iota_iterator(0_sz), iota_iterator(lhs.size()), [&](size_t i) noexcept
        {
            for (size_t j = 0; j < rhs.size(); j++)
            {
                const ParetoDominance rhs_state_j = std::atomic_ref{ rhs_state[j] }.load(std::memory_order_relaxed);

                if (lhs_state[i] == DOMINATED) continue;
                if (rhs_state_j  == DOMINATED) continue;

                if (lhs_state[i] == OPTIMAL)
                {
                    if (rhs_state_j == UNKNOWN && math::paretoCompareLess(rhs[j].fitness, lhs[i].fitness))
                    {
                        std::atomic_ref{ rhs_state[j] }.store(DOMINATED, std::memory_order_relaxed);
                    }
                    continue;
                }
                if (rhs_state_j == OPTIMAL)
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
                    std::atomic_ref{ rhs_state[j] }.store(OPTIMAL, std::memory_order_relaxed);
                }
                else if (comp > 0)
                {
                    lhs_state[i] = OPTIMAL;
                    std::atomic_ref{ rhs_state[j] }.store(DOMINATED, std::memory_order_relaxed);
                }
                // comp == 0 --> both are OPTIMAL or DOMINATED, can't know
            }
        });

        Candidates<T> optimal_solutions;
        optimal_solutions.reserve(lhs.size() + rhs.size());

        for (size_t i = 0; i < lhs.size(); i++)
        {
            if (lhs_state[i] != DOMINATED) optimal_solutions.push_back(std::move(lhs[i]));
        }
        for (size_t i = 0; i < rhs.size(); i++)
        {
            if (rhs_state[i] != DOMINATED) optimal_solutions.push_back(std::move(rhs[i]));
        }

        return optimal_solutions;
    }

} // namespace gapp::detail

#endif // !GAPP_CORE_POPULATION_HPP