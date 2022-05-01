/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "population.hpp"
#include "../utility/math.hpp"
#include <algorithm>
#include <numeric>
#include <vector>
#include <cmath>
#include <cassert>
#include <cstddef>

namespace genetic_algorithm::detail
{
    FitnessVector toFitnessVector(const FitnessMatrix& fmat)
    {
        assert(std::all_of(fmat.begin(), fmat.end(), [](const FitnessVector& fvec) { return fvec.size() == 1; }));

        return detail::map(fmat, [](const FitnessVector& fvec) { return fvec[0]; });
    }

    double populationFitnessMin(const FitnessVector& pop)
    {
        assert(!pop.empty());

        return *std::min_element(pop.begin(), pop.end());
    }

    FitnessVector populationFitnessMin(const FitnessMatrix& pop)
    {
        assert(!pop.empty());
        assert(std::all_of(pop.begin(), pop.end(), [&pop](const FitnessVector& sol) { return sol.size() == pop[0].size(); }));

        FitnessVector min_fitness = pop[0];
        for (size_t i = 1; i < pop.size(); i++)
        {
            for (size_t j = 0; j < min_fitness.size(); j++)
            {
                min_fitness[j] = std::min(min_fitness[j], pop[i][j]);
            }
        }

        return min_fitness;
    }

    double populationFitnessMax(const FitnessVector& pop)
    {
        assert(!pop.empty());

        return *std::max_element(pop.begin(), pop.end());
    }

    FitnessVector populationFitnessMax(const FitnessMatrix& pop)
    {
        assert(!pop.empty());
        assert(std::all_of(pop.begin(), pop.end(), [&pop](const FitnessVector& sol) { return sol.size() == pop[0].size(); }));

        FitnessVector max_fitness = pop[0];
        for (size_t i = 1; i < pop.size(); i++)
        {
            for (size_t j = 0; j < max_fitness.size(); j++)
            {
                max_fitness[j] = std::max(max_fitness[j], pop[i][j]);
            }
        }

        return max_fitness;
    }

    double populationFitnessMean(const FitnessVector& pop)
    {
        assert(!pop.empty());
        
        return mean(pop);
    }

    FitnessVector populationFitnessMean(const FitnessMatrix& pop)
    {
        assert(!pop.empty());
        assert(std::all_of(pop.begin(), pop.end(), [&pop](const FitnessVector& sol) { return sol.size() == pop[0].size(); }));

        FitnessVector fitness_mean(pop[0].size(), 0.0);
        for (const auto& sol : pop)
        {
            for (size_t i = 0; i < fitness_mean.size(); i++)
            {
                fitness_mean[i] += sol[i] / pop.size();
            }
        }

        return fitness_mean;
    }

    double populationFitnessSD(const FitnessVector& pop)
    {
        assert(!pop.empty());

        return stdDev(pop);
    }

    double populationFitnessSD(const FitnessVector& pop, double mean)
    {
        assert(!pop.empty());

        return stdDev(pop, mean);
    }

    FitnessVector populationFitnessSD(const FitnessMatrix& pop)
    {
        return populationFitnessSD(pop, populationFitnessMean(pop));
    }

    FitnessVector populationFitnessSD(const FitnessMatrix& pop, const FitnessVector& mean)
    {
        assert(!pop.empty());
        assert(std::all_of(pop.begin(), pop.end(), [&pop](const FitnessVector& sol) { return sol.size() == pop[0].size(); }));

        if (pop.size() == 1)
        {
            return FitnessVector(pop[0].size(), 0.0);
        }

        auto variance = FitnessVector(pop[0].size(), 0.0);
        for (const auto& sol : pop)
        {
            for (size_t i = 0; i < variance.size(); i++)
            {
                variance[i] += std::pow(sol[i] - mean[i], 2) / (pop.size() - 1.0);
            }
        }
        for (auto& elem : variance)
        {
            elem = std::sqrt(elem);
        }

        return variance;
    }

} // namespace genetic_algorithm::detail