/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "nsga2.hpp"
#include "nd_sort.hpp"
#include "../core/ga_info.hpp"
#include "../core/population.hpp"
#include "../utility/algorithm.hpp"
#include "../utility/functional.hpp"
#include "../utility/math.hpp"
#include "../utility/rng.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <iterator>
#include <vector>
#include <utility>
#include <stdexcept>
#include <cstddef>

namespace gapp::algorithm
{
    using namespace dtl;

    /* Calculate the crowding distances of the solutions in pfronts. */
    static std::vector<double> crowdingDistances(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last, ParetoFronts pfronts)
    {
        GAPP_ASSERT(std::distance(first, last) >= 0);

        std::vector<double> crowding_distances(last - first, 0.0);
        const auto front_bounds = paretoFrontBounds(pfronts);

        for (size_t obj = 0; obj < first->size(); obj++)
        {
            FitnessVector fvec(last - first, 0.0);
            std::transform(first, last, fvec.begin(), detail::element_at(obj));

            for (auto [front_first, front_last] : front_bounds)
            {
                std::sort(front_first, front_last, [&](const FrontInfo& lhs, const FrontInfo& rhs) noexcept
                {
                    return fvec[lhs.idx] < fvec[rhs.idx]; // ascending
                });

                const size_t first_idx = front_first->idx;
                const size_t last_idx = (front_last - 1)->idx;

                const double finterval = std::max(fvec[last_idx] - fvec[first_idx], 1E-6);

                crowding_distances[first_idx] = math::inf<double>;
                crowding_distances[last_idx]  = math::inf<double>;

                for (auto it = front_first + 1; it < front_last - 1; ++it)
                {
                    const size_t next = std::next(it)->idx;
                    const size_t prev = std::prev(it)->idx;
                    crowding_distances[it->idx] += (fvec[next] - fvec[prev]) / finterval;
                }
            }
        }

        return crowding_distances;
    }

    void NSGA2::initializeImpl(const GaInfo& ga)
    {
        GAPP_ASSERT(ga.num_objectives() > 1, "The number of objectives must be greater than 1 for the NSGA-II algorithm.");

        const auto& fmat = ga.fitness_matrix();
        auto pfronts = nonDominatedSort(fmat.begin(), fmat.end());

        ranks_ = paretoRanks(pfronts);
        dists_ = crowdingDistances(fmat.begin(), fmat.end(), std::move(pfronts));
    }

    size_t NSGA2::selectImpl(const GaInfo&, const FitnessMatrix& fmat) const
    {
        GAPP_ASSERT(!fmat.empty() && fmat.size() == ranks_.size());

        const size_t idx1 = rng::randomIdx(fmat);
        const size_t idx2 = rng::randomIdx(fmat);

        // lower ranks and higher crowding distances are better
        const bool first_is_better =
            (ranks_[idx1] != ranks_[idx2]) ?
                ranks_[idx1] < ranks_[idx2] :
                dists_[idx1] > dists_[idx2];

        return first_is_better ? idx1 : idx2;
    }

    std::vector<size_t> NSGA2::nextPopulationImpl(const GaInfo& ga,
                                                  FitnessMatrix::const_iterator parents_first,
                                                  [[maybe_unused]] FitnessMatrix::const_iterator children_first,
                                                  FitnessMatrix::const_iterator children_last)
    {
        const size_t popsize = ga.population_size();

        GAPP_ASSERT(ga.num_objectives() > 1);
        GAPP_ASSERT(size_t(children_first - parents_first) == popsize);
        GAPP_ASSERT(parents_first->size() == ga.num_objectives());

        auto pfronts = nonDominatedSort(parents_first, children_last);
        auto [partial_front_first, partial_front_last] = findPartialFront(pfronts.begin(), pfronts.end(), popsize);

        /* Crowding distances of just the partial front for sorting. */
        dists_ = crowdingDistances(parents_first, children_last, ParetoFronts(partial_front_first, partial_front_last));

        std::sort(partial_front_first, partial_front_last,
        [this](const FrontInfo& lhs, const FrontInfo& rhs) noexcept
        {
            return dists_[lhs.idx] > dists_[rhs.idx]; // descending
        });

        /* Crowding distances of all of the solutions. */
        dists_ = crowdingDistances(parents_first, children_last, ParetoFronts(pfronts.begin(), pfronts.begin() + popsize));
        dists_.resize(popsize);

        /* Add the first popsize elements from pfronts to the next pop. */
        std::vector<size_t> new_pop(popsize);
        for (size_t i = 0; i < popsize; i++)
        {
            new_pop[i] = pfronts[i].idx;
            ranks_[i] = pfronts[i].rank;
        }

        return new_pop;
    }

    std::vector<size_t> NSGA2::optimalSolutionsImpl(const GaInfo&) const
    {
        return detail::find_indices(ranks_, detail::equal_to(0_sz));
    }

} // namespace gapp::algorithm