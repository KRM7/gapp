/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "population.hpp"
#include "../utility/math.hpp"
#include "../utility/algorithm.hpp"
#include <algorithm>
#include <iterator>
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

namespace genetic_algorithm::detail::_
{
    std::vector<size_t> findParetoFront1D(const FitnessMatrix& fmat)
    {
        auto best = std::max_element(fmat.begin(), fmat.end(),
        [](const FitnessVector& lhs, const FitnessVector& rhs) noexcept
        {
            return lhs[0] < rhs[0];
        });

        return detail::find_indices(fmat,
        [best](const FitnessVector& fvec) noexcept
        {
            return detail::floatIsEqual((*best)[0], fvec[0]);
        });
    }

    std::vector<size_t> findParetoFrontND(const FitnessMatrix& fmat)
    {
        auto indices = detail::argsort(fmat.begin(), fmat.end(),
        [](const FitnessVector& lhs, const FitnessVector& rhs) noexcept
        {
            for (size_t i = 0; i < lhs.size(); i++)
            {
                if (lhs[i] != rhs[i]) return lhs[i] > rhs[i];
            }
            return false;
        });

        std::vector<size_t> optimal_indices;

        for (const auto& idx : indices)
        {
            bool dominated = std::any_of(optimal_indices.begin(), optimal_indices.end(),
            [&fmat, idx](size_t optimal_idx) noexcept
            {
                return detail::paretoCompareLess(fmat[idx], fmat[optimal_idx]);
            });
            if (!dominated) optimal_indices.push_back(idx);
        }

        return optimal_indices;
    }

    std::vector<size_t> findParetoFrontKung(const FitnessMatrix& fmat)
    {
        /* See: Kung et al. "On finding the maxima of a set of vectors." Journal of the ACM (JACM) 22.4 (1975): 469-476. */
        /* Doesn't work for d = 1 (single-objective optimization). */

        auto indices = detail::argsort(fmat.begin(), fmat.end(),
        [](const FitnessVector& lhs, const FitnessVector& rhs) noexcept
        {
            for (size_t i = 0; i < lhs.size(); i++)
            {
                if (lhs[i] != rhs[i]) return lhs[i] > rhs[i];
            }
            return false;
        });

        return _::findParetoFrontKungImpl(fmat, indices.begin(), indices.end());
    }

    bool kungCompareLess(const FitnessVector& lhs, const FitnessVector& rhs) noexcept
    {
        bool is_dominated = detail::paretoCompareLess(lhs, rhs, 1);

        bool is_equal = !detail::floatIsEqual(lhs[0], rhs[0]) &&
                        std::equal(lhs.begin() + 1, lhs.end(),
                                   rhs.begin() + 1,
                                   detail::floatIsEqual<double>);

        return is_dominated || is_equal;
    }

    std::vector<size_t> findParetoFrontKungImpl(const FitnessMatrix& fmat, std::vector<size_t>::iterator first, std::vector<size_t>::iterator last)
    {
        if (std::distance(first, last) == 1) return { *first };

        auto middle = first + std::distance(first, last) / 2;

        auto top_half = _::findParetoFrontKungImpl(fmat, first, middle);
        auto bottom_half = _::findParetoFrontKungImpl(fmat, middle, last);

        for (const auto& bad : bottom_half)
        {
            bool is_dominated = false;
            for (const auto& good : top_half)
            {
                if (_::kungCompareLess(fmat[bad], fmat[good]))
                {
                    is_dominated = true;
                    break;
                }
            }
            if (!is_dominated) top_half.push_back(bad);
        }

        return top_half;
    }

} // namespace genetic_algorithm::detail::_