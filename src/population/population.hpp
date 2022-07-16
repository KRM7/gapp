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
    FitnessVector minFitness(const FitnessMatrix& pop);

    /* Return the maximum fitness value of the population along each objective axis. */
    FitnessVector maxFitness(const FitnessMatrix& pop);

    /* Return the mean fitness value of the population along each objective axis. */
    FitnessVector fitnessMean(const FitnessMatrix& pop);

    /* Return the standard deviation of the fitness values of the population along each objective axis. */
    FitnessVector fitnessStdDev(const FitnessMatrix& pop);

    /* Return the standard deviation of the fitness values of the population along each objective axis. */
    FitnessVector fitnessStdDev(const FitnessMatrix& pop, const FitnessVector& mean);

    /* Finds the pareto-optimal solutions in a population. */
    template<Gene T>
    Candidates<T> findParetoFront(const Population<T>& pop);

} // namespace genetic_algorithm::detail

namespace genetic_algorithm::detail::_
{
    std::vector<size_t> findParetoFront1D(const FitnessMatrix& fmat);

    std::vector<size_t> findParetoFrontND(const FitnessMatrix& fmat);

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
            _::findParetoFrontND(fitness_matrix); // faster than Kung's algorithm

        return detail::map(optimal_indices, [&pop](size_t idx) { return pop[idx]; });
    }

} // namespace genetic_algorithm::detail

#endif // !GA_POPULATION_HPP