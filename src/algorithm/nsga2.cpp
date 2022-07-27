/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "nsga2.hpp"
#include "../core/ga_info.hpp"
#include "../population/population.hpp"
#include "../utility/rng.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <iterator>
#include <vector>
#include <limits>
#include <utility>
#include <cstddef>
#include <cassert>

namespace genetic_algorithm::algorithm
{
    using dtl::ParetoFronts;


    std::vector<double> NSGA2::crowdingDistances(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last, ParetoFronts pfronts)
    {
        assert(std::distance(first, last) >= 0);

        std::vector cdistances(size_t(last - first), 0.0);

        using Iter = ParetoFronts::iterator;
        using IterPair = std::pair<Iter, Iter>;

        auto front_bounds = dtl::paretoFrontBounds(pfronts);

        std::for_each(front_bounds.begin(), front_bounds.end(),
        [&cdistances, fmat = first](const IterPair& bounds)
        {
            const auto& [front_first, front_last] = bounds;
            
            /* Calculate distances along each fitness dimension */
            for (size_t dim = 0; dim < fmat->size(); dim++)
            {
                std::sort(front_first, front_last,
                [fmat, dim](const auto& lhs, const auto& rhs) noexcept
                {
                    return fmat[lhs.first][dim] < fmat[rhs.first][dim]; // ascending
                });

                const auto& first_idx = front_first->first;
                const auto& last_idx  = (front_last - 1)->first;

                double finterval = std::max(fmat[last_idx][dim] - fmat[first_idx][dim], 1E-6);

                cdistances[first_idx] = std::numeric_limits<double>::infinity();
                cdistances[last_idx]  = std::numeric_limits<double>::infinity();

                for (auto it = front_first + 1; it < front_last - 1; ++it)
                {
                    size_t this_ = it->first;
                    size_t next = std::next(it)->first;
                    size_t prev = std::prev(it)->first;

                    cdistances[this_] += (fmat[next][dim] - fmat[prev][dim]) / finterval;
                }
            }
        });

        return cdistances;
    }

    void NSGA2::initialize(const GaInfo& ga)
    {
        assert(ga.num_objectives() > 1);
        assert(ga.population_size() != 0);

        auto& fmat = ga.fitness_matrix();
        auto pfronts = dtl::nonDominatedSort(fmat.begin(), fmat.end());

        ranks_ = dtl::paretoRanks(pfronts);
        dists_ = crowdingDistances(fmat.begin(), fmat.end(), std::move(pfronts));
    }

    bool NSGA2::crowdedCompare(size_t lidx, size_t ridx) const noexcept
    {
        return (ranks_[lidx] != ranks_[ridx]) ?
            ranks_[lidx] < ranks_[ridx] :
            dists_[lidx] > dists_[ridx];
    }

    size_t NSGA2::select(const GaInfo&, const FitnessMatrix& pop) const
    {
        assert(!pop.empty() && pop.size() == ranks_.size());

        size_t idx1 = rng::randomIdx(pop);
        size_t idx2 = rng::randomIdx(pop);

        return crowdedCompare(idx1, idx2) ? idx1 : idx2;
    }

    std::vector<size_t> NSGA2::nextPopulation(const GaInfo& ga,
                                              FitnessMatrix::const_iterator parents_first,
                                              FitnessMatrix::const_iterator children_first,
                                              FitnessMatrix::const_iterator children_last)
    {
        const size_t popsize = ga.population_size();

        assert(ga.num_objectives() > 1);
        assert(size_t(children_first - parents_first) == popsize);
        assert(size_t(children_last - parents_first) >= popsize);
        assert(std::all_of(parents_first, children_last, [&ga](const FitnessVector& f) { return f.size() == ga.num_objectives(); }));

        GA_UNUSED(children_first);

        auto pfronts = dtl::nonDominatedSort(parents_first, children_last);
        ranks_ = dtl::paretoRanks(pfronts);
        dists_ = crowdingDistances(parents_first, children_last, pfronts);

        /* Keep track of the details of the candidates added to new_pop to avoid a second sort. */
        std::vector<size_t> new_pop(popsize);
        std::vector<size_t> new_ranks(popsize);
        std::vector<double> new_dists(popsize);

        /* Check if there is going to be a partial front, and find its bounds. */
        auto partial_front_first = pfronts.begin() + popsize;
        auto partial_front_last  = pfronts.begin() + popsize;

        bool has_partial_front = (pfronts.size() > popsize) && (pfronts[popsize - 1].second == pfronts[popsize].second);

        if (has_partial_front)
        {
            size_t partial_front_rank = pfronts[popsize - 1].second;

            partial_front_first = std::find_if(pfronts.begin(), pfronts.end(),
            [partial_front_rank](const auto& sol) noexcept
            { 
                return sol.second == partial_front_rank;
            });

            partial_front_last = dtl::nextFrontBegin(partial_front_first, pfronts.end());
        }

        /* Find the elements that will be added to the next population from the partial front. */
        std::sort(partial_front_first, partial_front_last,
        [this](const auto& lhs, const auto& rhs) noexcept
        {
            size_t lidx = lhs.first;
            size_t ridx = rhs.first;
            /* The ranks will always be the same within a front, so only compare the distances. */
            return dists_[lidx] > dists_[ridx]; // descending
        });

        /* Recalc the crowding distances of the partial front. */
        dtl::ParetoFronts partial_front(partial_front_first, partial_front_last);

        auto changed_dists = crowdingDistances(parents_first, children_last, partial_front);
        for (const auto& [sol_idx, _] : partial_front)
        {
            dists_[sol_idx] = changed_dists[sol_idx];
        }

        /* Add the first popsize elements from pfronts to the next pop. */
        for (size_t i = 0; i < popsize; i++)
        {
            new_pop[i] = pfronts[i].first;
            new_ranks[i] = pfronts[i].second;
            new_dists[i] = dists_[i];
        }

        ranks_ = std::move(new_ranks);
        dists_ = std::move(new_dists);

        return new_pop;
    }

} // namespace genetic_algorithm::algorithm