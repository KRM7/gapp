/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "population.hpp"
#include "../utility/math.hpp"
#include "../utility/algorithm.hpp"
#include "../utility/functional.hpp"
#include <algorithm>
#include <iterator>
#include <vector>
#include <cstddef>

namespace gapp::detail
{
    FitnessMatrix toFitnessMatrix(const PopulationView& pop)
    {
        if (pop.empty()) return {};

        FitnessMatrix fitness_matrix;
        fitness_matrix.reserve(pop.size(), pop[0].fitness.size());

        for (const CandidateInfo& sol : pop)
        {
            fitness_matrix.append_row(sol.fitness);
        }

        return fitness_matrix;
    }

    FitnessVector toFitnessVector(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        FitnessVector fitness_vector(last - first);
        std::transform(first, last, fitness_vector.begin(), detail::element_at(0));

        return fitness_vector;
    }

    FitnessVector toFitnessVector(const PopulationView& pop)
    {
        FitnessVector fitness_vector(pop.size());
        std::transform(pop.begin(), pop.end(), fitness_vector.begin(), [](const CandidateInfo& sol) { return sol.fitness[0]; });

        return fitness_vector;
    }

    small_vector<size_t> findParetoFront(const FitnessMatrix& fmat)
    {
        if (fmat.empty()) return {};

        return (fmat.ncols() == 1) ? findParetoFront1D(fmat) : findParetoFrontBest(fmat);
    }

    small_vector<size_t> findParetoFront1D(const FitnessMatrix& fmat)
    {
        const auto best = std::max_element(fmat.begin(), fmat.end(),
        [](const auto& lhs, const auto& rhs) noexcept
        {
            return lhs[0] < rhs[0];
        });

        const double max_fitness = (*best)[0];

        return detail::find_indices(fmat, [&](const auto& fitness_vector) noexcept
        {
            return math::floatIsEqual(max_fitness, fitness_vector[0]);
        });
    }

    small_vector<size_t> findParetoFrontSort(const FitnessMatrix& fmat)
    {
        const auto indices = detail::argsort(fmat.begin(), fmat.end(), [](const auto& lhs, const auto& rhs) noexcept
        {
            for (size_t i = 0; i < lhs.size(); i++)
            {
                if (lhs[i] != rhs[i]) return lhs[i] > rhs[i];
            }
            return false;
        });

        small_vector<size_t> optimal_indices;

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

    small_vector<size_t> findParetoFrontBest(const FitnessMatrix& fmat)
    {
        // Implementation of the BEST algorithm based on the description in:
        //      Godfrey et al. "Algorithms and analyses for maximal vector computation." The VLDB Journal 16, no. 1 (2007): 5-28.

        auto indices = detail::index_vector(fmat.size());

        small_vector<size_t> optimal_indices;
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
                    // [it] is dominated by [best], so it is removed from the range, but
                    // can't swap to the front of the range here, as that could overwrite the [best] element.
                    // Swapping to the back is fine, but that element hasn't been checked yet, so shouldn't
                    // advance [it] in this case.
                    std::iter_swap(it--, --last);
                }
                else if (comp < 0)
                {
                    // [best] is dominated by [it], so it is removed from the range, but
                    // can't swap to the back of the range here, since that element hasn't been checked yet.
                    // Swapping to the front is fine, since the [best] iterator will have to be updated anyway.
                    std::iter_swap(best, first++);
                    best = it;
                }
            }

            optimal_indices.push_back(*best); // best is definitely optimal at this point

            // best was only compared with the elements after it in the range,
            // there could be dominated elements before it still in the list
            // which have to be eliminated
            for (auto it = first; it < best; ++it)
            {
                if (math::paretoCompareLess(fmat[*it], fmat[*best]))
                {
                    // removing the dominated element by swapping to the front is fine here,
                    // since we know that first != best
                    std::iter_swap(it, first++);
                }
            }
            std::iter_swap(best, --last); // best shouldn't be selected again, remove

            // none of the remaining indices in [first, last) are dominated by best,
            // but they could be dominated by another element in the list,
            // so they can't be added to the optimal list yet, keep going.
        }

        return optimal_indices;
    }

    small_vector<size_t> findParetoFrontKungImpl(const FitnessMatrix& fmat, small_vector<size_t>::const_iterator first, small_vector<size_t>::const_iterator last)
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

    small_vector<size_t> findParetoFrontKung(const FitnessMatrix& fmat)
    {
        /* See: Kung et al. "On finding the maxima of a set of vectors." Journal of the ACM (JACM) 22.4 (1975): 469-476. */
        /* Doesn't work for d = 1 (single-objective optimization). */

        auto indices = detail::argsort(fmat.begin(), fmat.end(),
        [](const auto& lhs, const auto& rhs) noexcept
        {
            for (size_t i = 0; i < lhs.size(); i++)
            {
                if (lhs[i] != rhs[i]) return lhs[i] > rhs[i];
            }
            return false;
        });

        return findParetoFrontKungImpl(fmat, indices.cbegin(), indices.cend());
    }

    FitnessVector findNadirPoint(const FitnessMatrix& fitness_matrix)
    {
        if (fitness_matrix.empty()) return {};

        /* Nadir point estimate = minimum of extreme points along each objective axis. */
        const auto& front_indices = detail::findParetoFront(fitness_matrix);
        FitnessVector nadir_point{ fitness_matrix[front_indices[0]] };

        for (size_t i = 1; i < front_indices.size(); i++)
        {
            detail::elementwise_min(nadir_point, fitness_matrix[front_indices[i]], detail::inplace_t{});
        }

        return nadir_point;
    }

    FitnessVector findFrontNadirPoint(const FitnessMatrix& optimal_points)
    {
        if (optimal_points.empty()) return {};

        /* Nadir point estimate = minimum of extreme points along each objective axis. */
        FitnessVector nadir_point{ optimal_points[0] };
        for (size_t i = 1; i < optimal_points.size(); i++)
        {
            detail::elementwise_min(nadir_point, optimal_points[i], detail::inplace_t{});
        }

        return nadir_point;
    }

} // namespace gapp::detail
