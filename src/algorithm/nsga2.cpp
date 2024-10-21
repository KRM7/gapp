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
    static std::vector<double> crowdingDistances(const FitnessMatrix& fmat, const std::vector<ParetoFrontsRange>& pareto_fronts)
    {
        GAPP_ASSERT(!fmat.empty());
        GAPP_ASSERT(std::none_of(pareto_fronts.begin(), pareto_fronts.end(), detail::is_size(0)));

        std::vector crowding_distances(fmat.size(), 0.0);

        for (size_t obj = 0; obj < fmat.ncols(); obj++)
        {
            const auto fvec = fmat.column(obj);

            for (const auto& front : pareto_fronts)
            {
                std::sort(front.begin(), front.end(), [&](const FrontElement& lhs, const FrontElement& rhs) noexcept
                {
                    return fvec[lhs.idx] < fvec[rhs.idx]; // ascending
                });

                const size_t first = front.front().idx;
                const size_t last  = front.back().idx;

                const double finterval = std::max(fvec[last] - fvec[first], 1E-8);

                crowding_distances[first] = math::inf<double>;
                crowding_distances[last]  = math::inf<double>;

                for (auto it = front.begin() + 1; it < front.end() - 1; ++it)
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

        ParetoFronts pareto_fronts = nonDominatedSort(ga.fitness_matrix());

        ranks_ = pareto_fronts.ranks();
        dists_ = crowdingDistances(ga.fitness_matrix(), pareto_fronts.fronts());
    }

    const CandidateInfo& NSGA2::selectImpl(const GaInfo&, const PopulationView& pop) const
    {
        GAPP_ASSERT(!pop.empty() && pop.size() == ranks_.size());

        const size_t idx1 = rng::randomIndex(pop);
        const size_t idx2 = rng::randomIndex(pop);

        // lower ranks and higher crowding distances are better
        if (ranks_[idx1] < ranks_[idx2]) return pop[idx1];
        if (ranks_[idx1] > ranks_[idx2]) return pop[idx2];
        if (dists_[idx1] > dists_[idx2]) return pop[idx1];
        return pop[idx2];
    }

    CandidatePtrVec NSGA2::nextPopulationImpl(const GaInfo& ga, const PopulationView& pop)
    {
        GAPP_ASSERT(ga.num_objectives() > 1);

        const size_t popsize = ga.population_size();
        const FitnessMatrix fitness_matrix = detail::toFitnessMatrix(pop);

        auto pareto_fronts = nonDominatedSort(fitness_matrix);
        auto partial_front = pareto_fronts.partialFront(popsize);

        if (!partial_front.empty()) { dists_ = crowdingDistances(fitness_matrix, { partial_front }); }
        std::sort(partial_front.begin(), partial_front.end(), [this](const FrontElement& lhs, const FrontElement& rhs)
        {
            return dists_[lhs.idx] > dists_[rhs.idx]; // descending
        });

        pareto_fronts.resize(popsize);
        dists_ = crowdingDistances(fitness_matrix, pareto_fronts.fronts());
        dists_.resize(popsize);

        CandidatePtrVec new_pop(popsize);
        for (size_t i = 0; i < popsize; i++)
        {
            const auto [idx, rank] = pareto_fronts[i];

            new_pop[i] = &pop[idx];
            ranks_[i]  = rank;
        }

        return new_pop;
    }

    small_vector<size_t> NSGA2::optimalSolutionsImpl(const GaInfo&, const PopulationView&) const
    {
        return detail::find_indices(ranks_, detail::equal_to(0_sz));
    }

} // namespace gapp::algorithm
