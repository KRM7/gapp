/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include "pop_stats.hpp"
#include "../population/candidate.hpp"
#include "../population/population.hpp"
#include "../utility/algorithm.hpp"
#include "../utility/iterators.hpp"
#include "../utility/utility.hpp"
#include <vector>
#include <span>
#include <algorithm>
#include <numeric>
#include <functional>
#include <iterator>
#include <atomic>
#include <cmath>
#include <cstddef>

namespace gapp::detail
{
    FitnessVector minFitness(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        GA_ASSERT(std::distance(first, last) > 0);

        return std::accumulate(std::next(first), last, FitnessVector(*first), detail::elementwise_min<double>);
    }

    FitnessVector maxFitness(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        GA_ASSERT(std::distance(first, last) > 0);

        return std::accumulate(std::next(first), last, FitnessVector(*first), detail::elementwise_max<double>);
    }

    FitnessVector fitnessMean(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        GA_ASSERT(std::distance(first, last) > 0);

        const double ninv = 1.0 / (last - first);
        FitnessVector fitness_mean(first->size());

        std::for_each(first, last, [&](const auto& fvec) noexcept
        {
            for (size_t i = 0; i < fvec.size(); i++)
            {
                fitness_mean[i] += fvec[i] * ninv;
            }
        });

        return fitness_mean;
    }

    FitnessVector fitnessVariance(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last, std::span<const double> fitness_mean)
    {
        GA_ASSERT(std::distance(first, last) > 0);
        GA_ASSERT(first->size() == fitness_mean.size());

        const double ninv = 1.0 / (last - first - 1.0);
        FitnessVector fitness_variance(first->size(), 0.0);

        if (std::distance(first, last) == 1) return fitness_variance;

        std::for_each(first, last, [&](const auto& fitness_vector) noexcept
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

    FitnessVector fitnessStdDev(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last, std::span<const double> mean)
    {
        FitnessVector fitness_std_dev = fitnessVariance(first, last, mean);
        for (double& f : fitness_std_dev) { f = std::sqrt(f); }

        return fitness_std_dev;
    }

    FitnessVector fitnessStdDev(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        return fitnessStdDev(first, last, fitnessMean(first, last));
    }


    /**
    *   HYPERVOLUME
    * 
    * The implementation calculates exact hypervolumes using the WFG algorithm.
    * 
    * See:
    *   While, Lyndon, Lucas Bradstreet, and Luigi Barone. "A fast way of calculating exact hypervolumes."
    *   IEEE Transactions on Evolutionary Computation 16.1 (2011): 86-95.
    */

    struct par { /* tag for parallel execution */ };
    struct seq { /* tag for sequential execution */ };

    static double hypervolume(par, const FitnessMatrix& fmat, std::span<const double> ref_point);
    static double hypervolume(seq, const FitnessMatrix& fmat, std::span<const double> ref_point);

    static inline FitnessMatrix uniqueSortedParetoFront(const FitnessMatrix& fmat)
    {
        const auto optimal_indices = detail::findParetoFrontSort(fmat);

        FitnessMatrix front;
        front.reserve(optimal_indices.size(), fmat.ncols());

        for (size_t optimal_idx : optimal_indices)
        {
            if (front.empty() || fmat[optimal_idx] != front.back()) front.append_row(fmat[optimal_idx]);
        }

        return front;
    }

    static inline FitnessMatrix limitSet(FitnessMatrix fmat, std::span<const double> limit) noexcept
    {
        for (size_t row = 0; row < fmat.nrows(); row++)
        {
            for (size_t col = 0; col < fmat.ncols(); col++)
            {
                fmat[row][col] = std::min(fmat[row][col], limit[col]);
            }
        }
        return fmat;
    }

    static inline double inclusiveHypervolume(std::span<const double> point, std::span<const double> ref_point) noexcept
    {
        GA_ASSERT(point.size() == ref_point.size());
        GA_ASSERT(std::equal(point.begin(), point.end(), ref_point.begin(), ref_point.end(), std::greater_equal{}));
        
        return math::volumeBetween(point, ref_point);
    }

    static inline double exclusiveHypervolume(std::span<const double> point, const FitnessMatrix& rest, std::span<const double> ref_point)
    {
        const double inclusive_hypervolume = inclusiveHypervolume(point, ref_point);
        const double rest_hypervolume = hypervolume(seq{}, limitSet(rest, point), ref_point);

        return inclusive_hypervolume - rest_hypervolume;
    }

    static inline double hypervolume(par, const FitnessMatrix& fmat, std::span<const double> ref_point)
    {
        const FitnessMatrix front = uniqueSortedParetoFront(fmat);

        std::atomic<double> hypervolume = 0.0;
        std::for_each(GA_EXECUTION_UNSEQ, detail::iota_iterator(0_sz), detail::iota_iterator(front.size()), [&](size_t idx)
        {
            const auto point = front[idx];
            const FitnessMatrix rest = { front.begin() + idx + 1, front.end() };

            const double exclusive_hypervolume = exclusiveHypervolume(point, rest, ref_point);

            hypervolume.fetch_add(exclusive_hypervolume, std::memory_order_acq_rel);
        });

        return hypervolume.load(std::memory_order_acquire);
    }

    static inline double hypervolume(seq, const FitnessMatrix& fmat, std::span<const double> ref_point)
    {
        const FitnessMatrix front = uniqueSortedParetoFront(fmat);

        double hypervolume = 0.0;
        for (auto point = front.begin(); point < front.end(); ++point)
        {
            const FitnessMatrix rest{ std::next(point), front.end() };
            hypervolume += exclusiveHypervolume(*point, rest, ref_point);
        }

        return hypervolume;
    }

    double hypervolume(const FitnessMatrix& fmat, std::span<const double> ref_point)
    {
        return hypervolume(par{}, fmat, ref_point);
    }

} // namespace gapp::detail