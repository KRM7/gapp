/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "nsga2.hpp"
#include "../core/ga_info.hpp"
#include "../population/population.hpp"
#include "../utility/algorithm.hpp"
#include "../utility/functional.hpp"
#include "../utility/math.hpp"
#include "../utility/rng.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <iterator>
#include <vector>
#include <limits>
#include <utility>
#include <cstddef>
#include <cassert>
#include <stdexcept>

namespace genetic_algorithm::algorithm
{
    using namespace dtl;

    std::vector<double> NSGA2::crowdingDistances(FitnessMatrix::const_iterator fmat, FitnessMatrix::const_iterator last, ParetoFronts pfronts)
    {
        assert(std::distance(fmat, last) >= 0);

        std::vector<double> cdistances(last - fmat);

        using Iter     = ParetoFronts::iterator;
        using IterPair = std::pair<Iter, Iter>;

        auto front_bounds = paretoFrontBounds(pfronts);

        for (size_t dim = 0; dim < fmat->size(); dim++)
        {
            for (const auto& [front_first, front_last] : front_bounds)
            {
                std::sort(front_first, front_last, [&](const FrontInfo& lhs, const FrontInfo& rhs) noexcept
                {
                    return fmat[lhs.idx][dim] < fmat[rhs.idx][dim]; // ascending
                });

                size_t first_idx = front_first->idx;
                size_t last_idx = (front_last - 1)->idx;

                double finterval = std::max(fmat[last_idx][dim] - fmat[first_idx][dim], 1E-6);

                cdistances[first_idx] = math::inf<double>;
                for (auto it = front_first + 1; it < front_last - 1; ++it)
                {
                    size_t next = std::next(it)->idx;
                    size_t prev = std::prev(it)->idx;
                    cdistances[it->idx] += (fmat[next][dim] - fmat[prev][dim]) / finterval;
                }
                cdistances[last_idx] = math::inf<double>;
            }
        }

        return cdistances;
    }

    void NSGA2::initialize(const GaInfo& ga)
    {
        assert(ga.population_size() != 0);

        if (ga.num_objectives() <= 1)
        {
            GA_THROW(std::logic_error, "The number of objectives must be greater than 1 for the NSGA-II algorithm.");
        }

        auto& fmat = ga.fitness_matrix();
        auto pfronts = nonDominatedSort(fmat.begin(), fmat.end());

        ranks_ = paretoRanks(pfronts);
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
        assert(std::all_of(parents_first, children_last, [&ga](const FitnessVector& f) { return f.size() == ga.num_objectives(); }));

        GA_UNUSED(children_first);

        auto pfronts = nonDominatedSort(parents_first, children_last);
        ranks_ = paretoRanks(pfronts);
        dists_ = crowdingDistances(parents_first, children_last, pfronts);

        /* Keep track of the details of the candidates added to new_pop to avoid a second sort. */
        std::vector<size_t> new_pop(popsize);
        std::vector<size_t> new_ranks(popsize);
        std::vector<double> new_dists(popsize);

        auto [partial_front_first, partial_front_last] = findPartialFront(pfronts.begin(), pfronts.end(), popsize);

        std::sort(partial_front_first, partial_front_last,
        [this](const FrontInfo& lhs, const FrontInfo& rhs) noexcept
        {
            /* The ranks will always be the same within a front, so only compare the distances. */
            return dists_[lhs.idx] > dists_[rhs.idx]; // descending
        });

        ParetoFronts partial_front(partial_front_first, pfronts.begin() + popsize);

        /* Recalc the crowding distances of the partial front. */
        auto changed_dists = crowdingDistances(parents_first, children_last, partial_front);
        for (const auto& [sol_idx, _] : partial_front)
        {
            dists_[sol_idx] = changed_dists[sol_idx];
        }

        /* Add the first popsize elements from pfronts to the next pop. */
        for (size_t i = 0; i < popsize; i++)
        {
            new_pop[i]   = pfronts[i].idx;
            new_ranks[i] = pfronts[i].rank;
            new_dists[i] = dists_[i];
        }

        ranks_ = std::move(new_ranks);
        dists_ = std::move(new_dists);

        return new_pop;
    }

    std::optional<std::vector<size_t>> NSGA2::optimalSolutions(const GaInfo&) const
    {
        return detail::find_indices(ranks_, detail::equal_to(0_sz));
    }

} // namespace genetic_algorithm::algorithm