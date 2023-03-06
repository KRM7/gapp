/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "population.hpp"
#include "../utility/math.hpp"
#include "../utility/algorithm.hpp"
#include "../utility/functional.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <numeric>
#include <iterator>
#include <vector>
#include <cmath>
#include <cstddef>

namespace genetic_algorithm::detail
{
    FitnessVector toFitnessVector(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        GA_ASSERT(std::all_of(first, last, detail::is_size(1)));

        std::vector<double> fitness_vector(last - first);
        std::transform(first, last, fitness_vector.begin(), [](const FitnessVector& row) { return row[0]; });

        return fitness_vector;
    }

    FitnessVector minFitness(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        GA_ASSERT(std::distance(first, last) > 0);
        GA_ASSERT(std::all_of(first, last, detail::is_size(first->size())));

        return std::reduce(std::next(first), last, *first, detail::elementwise_min<double>);
    }

    FitnessVector maxFitness(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        GA_ASSERT(std::distance(first, last) > 0);
        GA_ASSERT(std::all_of(first, last, detail::is_size(first->size())));

        return std::reduce(std::next(first), last, *first, detail::elementwise_max<double>);
    }

    FitnessVector fitnessMean(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        GA_ASSERT(std::distance(first, last) > 0);
        GA_ASSERT(std::all_of(first, last, detail::is_size(first->size())));

        const double ninv = 1.0 / (last - first);
        FitnessVector fitness_mean(first->size());

        std::for_each(first, last, [&](const FitnessVector& fvec)
        {
            for (size_t i = 0; i < fvec.size(); i++)
            {
                fitness_mean[i] += fvec[i] * ninv;
            }
        });

        return fitness_mean;
    }

    FitnessVector fitnessVariance(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last, const FitnessVector& fitness_mean)
    {
        GA_ASSERT(std::distance(first, last) > 0);
        GA_ASSERT(std::all_of(first, last, detail::is_size(first->size())));

        const double ninv = 1.0 / (last - first - 1.0);
        FitnessVector fitness_variance(first->size(), 0.0);

        if (std::distance(first, last) == 1) return fitness_variance;

        std::for_each(first, last, [&](const FitnessVector& fitness_vector) noexcept
        {
            for (size_t i = 0; i < fitness_vector.size(); i++)
            {
                fitness_variance[i] += std::pow(fitness_vector[i] - fitness_mean[i], 2) * ninv;
            }
        });

        return fitness_variance;
    }

    FitnessVector fitnessVariance(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        return fitnessVariance(first, last, fitnessMean(first, last));
    }

    FitnessVector fitnessStdDev(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last, const FitnessVector& mean)
    {
        FitnessVector fitness_std_dev = fitnessVariance(first, last, mean);
        for (double& f : fitness_std_dev) { f = std::sqrt(f); }

        return fitness_std_dev;
    }

    FitnessVector fitnessStdDev(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        return fitnessStdDev(first, last, fitnessMean(first, last));
    }

    std::vector<size_t> findParetoFront1D(const FitnessMatrix& fmat)
    {
        GA_ASSERT(std::all_of(fmat.begin(), fmat.end(), detail::is_size(1)));

        const auto best = std::max_element(fmat.begin(), fmat.end(),
        [](const FitnessVector& lhs, const FitnessVector& rhs) noexcept
        {
            return lhs[0] < rhs[0];
        });

        const double max_fitness = (*best)[0];

        return detail::find_indices(fmat, [&](const FitnessVector& fitness_vector) noexcept
        {
            return math::floatIsEqual(max_fitness, fitness_vector[0]);
        });
    }

    std::vector<size_t> findParetoFrontSort(const FitnessMatrix& fmat)
    {
        const auto indices = detail::argsort(fmat.begin(), fmat.end(),
        [](const FitnessVector& lhs, const FitnessVector& rhs) noexcept
        {
            for (size_t i = 0; i < lhs.size(); i++)
            {
                if (lhs[i] != rhs[i]) return lhs[i] > rhs[i];
            }
            return false;
        });

        std::vector<size_t> optimal_indices;

        for (size_t idx : indices)
        {
            bool dominated = std::any_of(optimal_indices.begin(), optimal_indices.end(),
            [&](size_t optimal_idx) noexcept
            {
                return math::paretoCompareLess(fmat[idx], fmat[optimal_idx]);
            });
            if (!dominated) optimal_indices.push_back(idx);
        }

        return optimal_indices;
    }

    std::vector<size_t> findParetoFrontBest(const FitnessMatrix& fmat)
    {
        /* Implementation of the BEST algorithm based on the description in:
         * Godfrey et al. "Algorithms and analyses for maximal vector computation." The VLDB Journal 16, no. 1 (2007): 5-28. */

        if (fmat.empty()) return {};

        auto indices = detail::index_vector(fmat.size());

        std::vector<size_t> optimal_indices;
        optimal_indices.reserve(fmat.size());

        auto first = indices.begin();
        auto last  = indices.end();

        while (first != last)
        {
            auto best = first;
            for (auto it = std::next(first); it < last; ++it)
            {
                auto comp = math::paretoCompare(fmat[*best], fmat[*it]);
                if (comp > 0)
                {
                    std::iter_swap(it--, --last);
                }
                else if (comp < 0)
                {
                    /* Replace and remove best from the index list (can't swap to the back, as that element hasn't been checked yet). */
                    std::iter_swap(best, first++);
                    best = it;
                }
            }
            optimal_indices.push_back(*best);   /* best is definitely optimal */

            /* best was only compared with the elements after it, there could be dominated elements before it still in the index list. */
            for (auto it = first; it < best; ++it)
            {
                if (math::paretoCompareLess(fmat[*it], fmat[*best]))
                {
                    std::iter_swap(it, first++);
                }
            }
            std::iter_swap(best, --last);       /* best shouldn't be selected again. */

            /* None of the remaining indices in [first, last) are dominated by best, but they could be by another element in the list,
               so they can't be added to the optimal list yet. */
        }

        return optimal_indices;
    }

    std::vector<size_t> findParetoFrontKungImpl(const FitnessMatrix& fmat, std::vector<size_t>::const_iterator first, std::vector<size_t>::const_iterator last)
    {
        if (std::distance(first, last) == 1) return { *first };

        auto middle = first + std::distance(first, last) / 2;

        auto top_half    = findParetoFrontKungImpl(fmat, first, middle);
        auto bottom_half = findParetoFrontKungImpl(fmat, middle, last);

        for (const auto& bad : bottom_half)
        {
            bool is_dominated = false;
            for (const auto& good : top_half)
            {
                if (math::paretoCompareLess(fmat[bad], fmat[good]))
                {
                    is_dominated = true;
                    break;
                }
            }
            if (!is_dominated) top_half.push_back(bad);
        }

        return top_half;
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

        return findParetoFrontKungImpl(fmat, indices.cbegin(), indices.cend());
    }

} // namespace genetic_algorithm::detail