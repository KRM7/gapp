/*
*  MIT License
*
*  Copyright (c) 2021 Krisztián Rugási
*
*  Permission is hereby granted, free of charge, to any person obtaining a copy
*  of this software and associated documentation files (the "Software"), to deal
*  in the Software without restriction, including without limitation the rights
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*  copies of the Software, and to permit persons to whom the Software is
*  furnished to do so, subject to the following conditions:
*
*  The above copyright notice and this permission notice shall be included in all
*  copies or substantial portions of the Software.
*
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
*  SOFTWARE.
*/

#ifndef GA_POPULATION_HPP
#define GA_POPULATION_HPP

#include "candidate.hpp"

#include <vector>

namespace genetic_algorithm
{
    /** The Population type used in all of the genetic algorithms. */
    template<regular_hashable GeneType>
    using Population = std::vector<Candidate<GeneType>>;

    /** A vector of Candidates. */
    template<regular_hashable GeneType>
    using CandidateVec = std::vector<Candidate<GeneType>>;

    /** Return the minimum fitness value of the population along each objective axis. */
    template<regular_hashable GeneType>
    std::vector<double> populationFitnessMin(const Population<GeneType>& pop);

    /** Return the maximum fitness value of the population along each objective axis. */
    template<regular_hashable GeneType>
    std::vector<double> populationFitnessMax(const Population<GeneType>& pop);

    /** Return the mean fitness value of the population along each objective axis. */
    template<regular_hashable GeneType>
    std::vector<double> populationFitnessMean(const Population<GeneType>& pop);

    /** Return the standard deviation of the fitness values of the population along each objective axis. */
    template<regular_hashable GeneType>
    std::vector<double> populationFitnessSD(const Population<GeneType>& pop);

    /** Return all of the pareto-optimal solutions in the population assuming there is only 1 objective function. */
    template<regular_hashable GeneType>
    CandidateVec<GeneType> findParetoFront1D(const Population<GeneType>& pop);

    /** Return all of the pareto-optimal solutions in the population assuming there are more than 1 objective functions. */
    template<regular_hashable GeneType>
    CandidateVec<GeneType> findParetoFrontKung(const Population<GeneType>& pop);

} // namespace genetic_algorithm


/* IMPLEMENTATION */

#include "math.hpp"

#include <algorithm>
#include <numeric>
#include <iterator>
#include <cmath>
#include <cstddef>
#include <cassert>
#include <stdexcept>

namespace genetic_algorithm
{
    template<regular_hashable GeneType>
    std::vector<double> populationFitnessMin(const Population<GeneType>& pop)
    {
        if (pop.empty())
        {
            throw std::invalid_argument("Can't calculate min fitness for an empty population.");
        }
        if (!std::all_of(pop.begin(), pop.end(), [&pop](const Candidate<GeneType>& sol) { return sol.fitness.size() == pop[0].fitness.size(); }))
        {
            throw std::invalid_argument("The fitness vectors of the population must match in size to calculate the minimum.");
        }

        std::vector<double> min_fitness = pop[0].fitness;
        for (size_t i = 1; i < pop.size(); i++)
        {
            for (size_t j = 0; j < min_fitness.size(); j++)
            {
                min_fitness[j] = std::min(min_fitness[j], pop[i].fitness[j]);
            }
        }

        return min_fitness;
    }

    template<regular_hashable GeneType>
    std::vector<double> populationFitnessMax(const Population<GeneType>& pop)
    {
        if (pop.empty())
        {
            throw std::invalid_argument("Can't calculate max fitness for an empty population.");
        }
        if (!std::all_of(pop.begin(), pop.end(), [&pop](const Candidate<GeneType>& sol) { return sol.fitness.size() == pop[0].fitness.size(); }))
        {
            throw std::invalid_argument("The fitness vectors of the population must match in size to calculate the maximum.");
        }

        std::vector<double> max_fitness = pop[0].fitness;
        for (size_t i = 1; i < pop.size(); i++)
        {
            for (size_t j = 0; j < max_fitness.size(); j++)
            {
                max_fitness[j] = std::max(max_fitness[j], pop[i].fitness[j]);
            }
        }

        return max_fitness;
    }

    template<regular_hashable GeneType>
    std::vector<double> populationFitnessMean(const Population<GeneType>& pop)
    {
        if (pop.empty())
        {
            throw std::invalid_argument("Can't calculate mean fitness for an empty population.");
        }
        if (!std::all_of(pop.begin(), pop.end(), [&pop](const Candidate<GeneType>& sol) { return sol.fitness.size() == pop[0].fitness.size(); }))
        {
            throw std::invalid_argument("The fitness vectors of the population must match in size to calculate the mean.");
        }

        std::vector<double> fitness_mean(pop[0].fitness.size(), 0.0);
        for (const auto& sol : pop)
        {
            for (size_t i = 0; i < fitness_mean.size(); i++)
            {
                fitness_mean[i] += sol.fitness[i] / pop.size();
            }
        }

        return fitness_mean;
    }

    template<regular_hashable GeneType>
    std::vector<double> populationFitnessSD(const Population<GeneType>& pop)
    {
        if (pop.empty())
        {
            throw std::invalid_argument("Can't calculate fitness SD for an empty population.");
        }
        if (!std::all_of(pop.begin(), pop.end(), [&pop](const Candidate<GeneType>& sol) { return sol.fitness.size() == pop[0].fitness.size(); }))
        {
            throw std::invalid_argument("The fitness vectors of the population must match in size to calculate the SD.");
        }

        if (pop.size() == 1)
        {
            return std::vector<double>(pop[0].fitness.size(), 0.0);
        }

        auto mean = populationFitnessMean(pop);

        //auto variance = std::vector<long double>(pop[0].fitness.size(), 0.0);
        auto variance = std::vector<double>(pop[0].fitness.size(), 0.0);
        for (const auto& sol : pop)
        {
            for (size_t i = 0; i < variance.size(); i++)
            {
                variance[i] += std::pow(sol.fitness[i] - mean[i], 2) / (pop.size() - 1.0);
            }
        }
        for (auto& elem : variance)
        {
            elem = std::sqrt(elem);
        }

        //return std::vector<double>(variance.begin(), variance.end());
        return variance;
    }

    template<regular_hashable GeneType>
    CandidateVec<GeneType> findParetoFront1D(const Population<GeneType>& pop)
    {
        if (!std::all_of(pop.begin(), pop.end(), [](const Candidate<GeneType>& sol) { return sol.fitness.size() == 1; }))
        {
            throw std::invalid_argument("The size of the fitness vectors of the population must be 1 for this algorithm.");
        }

        if (pop.empty())
        {
            return {};
        }

        CandidateVec<GeneType> optimal_sols;

        /* Even for a single-objective problem, there might be different solutions with the same fitness value. */
        double max_fitness = populationFitnessMax(pop)[0];
        for (const auto& sol : pop)
        {
            //if (sol.fitness[0] == max_fitness) optimal_sols.push_back(sol);

            if (detail::floatIsEqual(max_fitness, sol.fitness[0])) optimal_sols.push_back(sol);
        }

        return optimal_sols;
    }

    template<regular_hashable GeneType>
    CandidateVec<GeneType> findParetoFrontKung(const Population<GeneType>& pop)
    {
        /* See: Kung et al. "On finding the maxima of a set of vectors." Journal of the ACM (JACM) 22.4 (1975): 469-476.*/
        /* Doesn't work for d = 1 (single-objective optimization). */

        using namespace std;
        using iter = vector<size_t>::iterator;

        assert(!pop.empty());
        assert(all_of(pop.begin(), pop.end(), [](const Candidate<GeneType>& sol) { return sol.fitness.size() > 1; }));
        assert(all_of(pop.begin(), pop.end(), [&pop](const Candidate<GeneType>& sol) { return sol.fitness.size() == pop[0].fitness.size(); }));

        size_t dim = pop[0].fitness.size();    /* The number of objectives. */

        /* Find the indices of pareto optimal solutions in the population (Kung's algorithm) assuming fitness maximization. */
        function<vector<size_t>(iter, iter)> pfront = [&pfront, &pop, &dim](iter first, iter last) -> vector<size_t>
        {
            if (distance(first, last) == 1) return { *first };

            vector<size_t> R = pfront(first, first + distance(first, last) / 2);    /* Top half. */
            vector<size_t> S = pfront(first + distance(first, last) / 2, last);     /* Bottom half. */

            /* T = Find all non-dominated elements of the bottom half. */
            vector<size_t> T;
            for (const auto& s : S)
            {
                /* Check if s is dominated by any solution in R. */
                bool is_dominated = false;
                for (const auto& r : R)
                {
                    /* Pareto compare s and r. */
                    /* The first dimension (d = 0) of the fitness vectors doesn't need to be compared since the pop is already sorted. */
                    for (size_t d = 1; d < dim; d++)
                    {
                        if (detail::floatIsLess(pop[r].fitness[d], pop[s].fitness[d]))
                        {
                            is_dominated = false;
                            break;
                        }
                        if (detail::floatIsLess(pop[s].fitness[d], pop[r].fitness[d]))
                        {
                            is_dominated = true;
                        }
                    }
                    if (is_dominated) break;
                }
                if (!is_dominated) T.push_back(s);
            }
            R.insert(R.end(), T.begin(), T.end());

            return R;
        };

        /* Find the indices of the pareto optimal candidates. */
        vector<size_t> indices(pop.size());
        iota(indices.begin(), indices.end(), 0U);

        /* Sort the pop indices into descending order based on first fitness value (needed for Kung's algorithm to work). */
        sort(indices.begin(), indices.end(), [&pop](size_t lidx, size_t ridx) { return pop[lidx].fitness[0] > pop[ridx].fitness[0]; });
        indices = pfront(indices.begin(), indices.end());

        CandidateVec<GeneType> optimal_sols;
        optimal_sols.reserve(indices.size());
        for (const auto& idx : indices)
        {
            optimal_sols.push_back(pop[idx]);
        }

        return optimal_sols;
    }

} // namespace genetic_algorithm

#endif // !GA_POPULATION_HPP