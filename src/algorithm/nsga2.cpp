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

    size_t NSGA2::select(const GaInfo&, const FitnessMatrix& pop)
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
        assert(size_t(children_first - parents_first) == ga.population_size());
        assert(size_t(children_last - parents_first) >= ga.population_size());
        assert(std::all_of(parents_first, children_last, [&ga](const FitnessVector& f) { return f.size() == ga.num_objectives(); }));

        GA_UNUSED(children_first);

        std::vector<size_t> new_pop;
        new_pop.reserve(ga.population_size());

        auto pfronts = dtl::nonDominatedSort(parents_first, children_last);
        ranks_ = dtl::paretoRanks(pfronts);
        dists_ = crowdingDistances(parents_first, children_last, pfronts);

        /* Keep track of the details of the candidates added to new_pop to avoid a second sort. */
        std::vector<size_t> new_ranks;
        std::vector<double> new_dists;
        new_ranks.reserve(ga.population_size());
        new_dists.reserve(ga.population_size());

        /* Add entire fronts while possible. */
        auto first = pfronts.begin();
        auto last = dtl::nextFrontBegin(first, pfronts.end());
        while (new_pop.size() + std::distance(first, last) <= ga.population_size())
        {
            for (; first != last; first++)
            {
                size_t idx = first->first;
                new_pop.push_back(idx);
                new_ranks.push_back(ranks_[idx]);
                new_dists.push_back(dists_[idx]);
            }
            //first = last; // TODO
            last = dtl::nextFrontBegin(first, pfronts.end());
        }

        /* Add the remaining candidates from the partial front if there is one. */
        if (new_pop.size() != ga.population_size())
        {
            size_t remaining_indices = ga.population_size() - new_pop.size();

            /* Keep track of the candidates added from the partial front,
             * the crowding distances will only need to be updated for these candidates. */
            dtl::ParetoFronts partial_front;
            partial_front.reserve(remaining_indices);

            std::sort(first, last,
            [this](const std::pair<size_t, size_t>& lhs, const std::pair<size_t, size_t>& rhs) noexcept
            {
                size_t lidx = lhs.first;
                size_t ridx = rhs.first;
                return dists_[lidx] > dists_[ridx]; /* The ranks will always be the same */
            });

            for (; new_pop.size() != ga.population_size(); first++)
            {
                size_t idx = first->first;

                new_pop.push_back(idx);
                new_ranks.push_back(ranks_[idx]);
                partial_front.push_back(*first);
            }

            auto changed_dists = crowdingDistances(parents_first, children_last, { partial_front });

            for (size_t i = 0; i < remaining_indices; i++)
            {
                double new_dist = changed_dists[partial_front[i].first];
                new_dists.push_back(new_dist);
            }
        }
        ranks_ = std::move(new_ranks);
        dists_ = std::move(new_dists);

        return new_pop;
    }

} // namespace genetic_algorithm::algorithm