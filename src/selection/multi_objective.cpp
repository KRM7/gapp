/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "multi_objective.hpp"
#include "selection_dtl.hpp"
#include "../algorithms/ga_info.hpp"
#include "../utility/rng.hpp"
#include <algorithm>
#include <execution>
#include <numeric>
#include <limits>
#include <tuple>
#include <stdexcept>
#include <cassert>

namespace genetic_algorithm::selection::multi_objective
{
    void NSGA2::init(const GaInfo& ga)
    {
        assert(ga.num_objectives() > 1);
        assert(ga.population_size() != 0);

        auto fmat = ga.fitness_matrix();
        auto pfronts = dtl::nonDominatedSort2(fmat);

        ranks_ = dtl::paretoRanks(pfronts);
        dists_ = dtl::crowdingDistances(fmat, std::move(pfronts));
    }

    void NSGA2::prepare(const GaInfo&, const FitnessMatrix&)
    {
        /* Nothing to do, the ranks and distances from the previous nextPopulation call are fine */
    }

    size_t NSGA2::select(const GaInfo&, const FitnessMatrix& pop)
    {
        assert(!pop.empty() && pop.size() == ranks_.size());

        size_t idx1 = rng::randomIdx(pop);
        size_t idx2 = rng::randomIdx(pop);

        return crowdedCompare(idx1, idx2) ? idx1 : idx2;
    }

    std::vector<size_t> NSGA2::nextPopulation(const GaInfo& ga, FitnessMatrix& combined_pop)
    {
        assert(ga.population_size() <= combined_pop.size());
        assert(std::all_of(combined_pop.begin(), combined_pop.end(), [&ga](const FitnessVector& f) { f.size() == ga.num_objectives(); }));

        std::vector<size_t> new_pop;
        new_pop.reserve(ga.population_size());

        auto pfronts = dtl::nonDominatedSort2(combined_pop);
        ranks_ = dtl::paretoRanks(pfronts);
        dists_ = dtl::crowdingDistances(combined_pop, pfronts);

        /* Keep track of the details of the candidates added to new_pop to avoid a second sort. */
        std::vector<size_t> new_ranks;
        std::vector<double> new_dists;
        new_ranks.reserve(ga.population_size());
        new_dists.reserve(ga.population_size());

        /* Add entire fronts while possible. */
        auto first = pfronts.begin();
        auto last  = dtl::nextFrontBegin(first, pfronts.end());
        while (new_pop.size() + std::distance(first, last) <= ga.population_size())
        {
            for (; first != last; first++)
            {
                size_t idx = first->first;
                new_pop.push_back(idx);
                new_ranks.push_back(ranks_[idx]);
                new_dists.push_back(dists_[idx]);
            }
            first = last;
            last  = dtl::nextFrontBegin(first, pfronts.end());
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

            auto changed_dists = dtl::crowdingDistances(combined_pop, { partial_front });

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

    bool NSGA2::crowdedCompare(size_t lidx, size_t ridx) const noexcept
    {
        if (ranks_[lidx] != ranks_[ridx])
        {
            return ranks_[lidx] < ranks_[ridx];
        }
        else
        {
            return dists_[lidx] > dists_[ridx];
        }
    }

    void NSGA3::init(const GaInfo& ga)
    {
        assert(ga.num_objectives() > 1);
        assert(ga.population_size() != 0);

        auto fitness_matrix = ga.fitness_matrix();

        ref_points_ = dtl::generateRefPoints(ga.population_size(), ga.num_objectives());

        ideal_point_ = detail::populationFitnessMax(fitness_matrix);

        extreme_points_ = {};
        updateExtremePoints(extreme_points_, fitness_matrix, ideal_point_);

        nadir_point_ = findNadirPoint(extreme_points_);

        sol_props_ = std::vector<CandidateInfo>(ga.population_size());
        associatePopWithRefs(fitness_matrix);
        ref_niche_counts_ = calcNicheCounts(ga, sol_props_);

        auto pfronts = dtl::nonDominatedSort(fitness_matrix);
        for (size_t i = 0; i < sol_props_.size(); i++)
        {
            sol_props_[i].rank = pfronts.ranks[i];
        }
    }

    void NSGA3::updateIdealPoint(Point& ideal_point, const FitnessMatrix& fmat)
    {
        assert(!fmat.empty());
        assert(ideal_point.size() == fmat[0].size());

        auto fmax = detail::populationFitnessMax(fmat);
        for (size_t i = 0; i < ideal_point.size(); i++)
        {
            ideal_point[i] = std::max(ideal_point[i], fmax[i]);
        }
    }

    std::vector<double> NSGA3::weightVector(size_t dimensions, size_t axis)
    {
        assert(dimensions > axis);

        std::vector weights(dimensions, 1E-6);
        weights[axis] = 1.0;

        return weights;
    }

    void NSGA3::updateExtremePoints(std::vector<Point>& extreme_points, const FitnessMatrix& fmat, const Point& ideal_point)
    {
        assert(!fmat.empty());
        assert(fmat[0].size() == ideal_point.size());

        size_t dim = ideal_point.size();

        std::vector<Point> new_extreme_points;
        new_extreme_points.reserve(extreme_points.size());

        for (size_t i = 0; i < dim; i++)
        {
            auto w = weightVector(dim, i);
            auto ASFi = dtl::ASF(ideal_point, std::move(w));

            auto chebysev_distances_f = detail::map(fmat, ASFi);
            auto fmat_min = std::min_element(chebysev_distances_f.begin(), chebysev_distances_f.end());

            auto chebysev_distances_e = detail::map(extreme_points, ASFi);
            auto ext_min = std::min_element(chebysev_distances_e.begin(), chebysev_distances_e.end());

            if (ext_min != chebysev_distances_e.end() && *ext_min < *fmat_min)
            {
                size_t argmin = size_t(ext_min - chebysev_distances_e.begin());
                new_extreme_points.push_back(extreme_points[argmin]);
            }
            else
            {
                size_t argmin = size_t(fmat_min - chebysev_distances_f.begin());
                new_extreme_points.push_back(fmat[argmin]);
            }
        }
        extreme_points = std::move(new_extreme_points);
    }

    auto NSGA3::findNadirPoint(const std::vector<Point>& extreme_points) -> Point
    {
        Point nadir(extreme_points.size());

        /* Nadir point estimate = minimum of extreme points along each objective axis. */
        for (size_t i = 0; i < extreme_points.size(); i++)
        {
            nadir[i] = extreme_points[0][i];
            for (size_t j = 1; j < extreme_points.size(); j++)
            {
                nadir[i] = std::min(nadir[i], extreme_points[j][i]);
            }
        }

        return nadir;
    }

    void NSGA3::prepare(const GaInfo&, const FitnessMatrix&)
    {
        /* Nothing to do */
    }

    constexpr bool NSGA3::nichedCompare(size_t lidx, size_t ridx) const noexcept
    {
        if (sol_props_[lidx].rank != sol_props_[ridx].rank)
        {
            return sol_props_[lidx].rank < sol_props_[ridx].rank;
        }
        else if (sol_props_[lidx].niche_count != sol_props_[ridx].niche_count)
        {
            return sol_props_[lidx].niche_count < sol_props_[ridx].niche_count;
        }
        else
        {
            return sol_props_[lidx].ref_dist < sol_props_[ridx].ref_dist;
        }
    }

    size_t NSGA3::select(const GaInfo&, const FitnessMatrix& pop)
    {
        assert(!pop.empty());

        size_t idx1 = rng::randomIdx(pop);
        size_t idx2 = rng::randomIdx(pop);

        return nichedCompare(idx1, idx2) ? idx1 : idx2;
    }

    void NSGA3::associatePopWithRefs(const FitnessMatrix& pop)
    {
        assert(!pop.empty());
        assert(std::all_of(pop.begin(), pop.end(), [&pop](const FitnessVector& sol) { return sol.size() == pop[0].size(); }));

        FitnessMatrix fnorms = pop;

        std::for_each(GA_EXECUTION_UNSEQ, fnorms.begin(), fnorms.end(),
        [this](FitnessVector& fnorm)
        {
            for (size_t i = 0; i < fnorm.size(); i++)
            {
                fnorm[i] -= ideal_point_[i];
                fnorm[i] /= std::min(nadir_point_[i] - ideal_point_[i], -1E-6);
            }
        });

        /* Associate each candidate with the closest reference point. */
        sol_props_.resize(pop.size());
        std::transform(GA_EXECUTION_UNSEQ, fnorms.begin(), fnorms.end(), sol_props_.begin(), sol_props_.begin(),
        [this](const std::vector<double>& f, CandidateInfo& info) -> CandidateInfo
        {
            std::tie(info.ref_idx, info.ref_dist) = dtl::findClosestRef(ref_points_, f);

            return info;
        });
    }

    std::vector<size_t> NSGA3::calcNicheCounts(const GaInfo& ga, std::vector<CandidateInfo>& props)
    {
        /* Calculate niche counts for each of the reference points */
        std::vector<size_t> ref_niche_counts(ga.population_size(), 0U);
        for (const auto& info : props)
        {
            ref_niche_counts[info.ref_idx]++;
        }

        /* Assign the niche counts to the candidates too. */
        for (auto& info : props)
        {
            info.niche_count = ref_niche_counts[info.ref_idx];
        }

        return ref_niche_counts;
    }

    std::vector<size_t> NSGA3::nextPopulation(const GaInfo& ga, FitnessMatrix& combined_pop)
    {
        updateIdealPoint(ideal_point_, combined_pop);
        updateExtremePoints(extreme_points_, combined_pop, ideal_point_);
        nadir_point_ = findNadirPoint(extreme_points_);

        std::vector<size_t> new_pop_idxs;
        new_pop_idxs.reserve(ga.population_size());

        sol_props_.resize(combined_pop.size());

        auto pfronts = dtl::nonDominatedSort(combined_pop);
        for (size_t i = 0; i < sol_props_.size(); i++)
        {
            sol_props_[i].rank = pfronts.ranks[i];
        }
        associatePopWithRefs(combined_pop);

        std::vector<CandidateInfo> new_props;
        new_props.reserve(ga.population_size());

        /* Add entire fronts while possible. */
        size_t front_idx = 0;
        for (; new_pop_idxs.size() + pfronts.idxs[front_idx].size() <= ga.population_size(); front_idx++)
        {
            for (const auto& idx : pfronts.idxs[front_idx])
            {
                new_pop_idxs.push_back(idx);
                new_props.push_back(sol_props_[idx]);
            }
        }
        std::vector<size_t> ref_niche_counts = calcNicheCounts(ga, new_props);

        /* Add remaining candidates from the partial front if there is one. */
        std::vector<size_t> partial_front = pfronts.idxs[front_idx];
        while (new_pop_idxs.size() != ga.population_size())
        {
            /* Find the lowest niche count in the partial front. */
            size_t min_count = std::numeric_limits<size_t>::max();
            for (const auto& sol_idx : partial_front)
            {
                size_t ref_idx = sol_props_[sol_idx].ref_idx;
                min_count = std::min(min_count, ref_niche_counts[ref_idx]);
            }

            /* Find the reference points with minimal niche counts, and pick one. */
            std::vector<size_t> refs;
            for (const auto& sol_idx : partial_front)
            {
                size_t ref_idx = sol_props_[sol_idx].ref_idx;
                if (ref_niche_counts[ref_idx] == min_count && std::find(refs.begin(), refs.end(), ref_idx) == refs.end())
                {
                    refs.push_back(ref_idx);
                }
            }
            size_t ref = rng::randomElement(refs);

            /* Find the idx of the closest sol in the partial front associated with this ref point. */
            size_t selected_sol_idx;
            double min_distance = std::numeric_limits<double>::infinity();
            for (const auto& sol_idx : partial_front)
            {
                if (sol_props_[sol_idx].ref_idx == ref && sol_props_[sol_idx].ref_dist < min_distance)
                {
                    min_distance = sol_props_[sol_idx].ref_dist;
                    selected_sol_idx = sol_idx;
                }
            }

            /* Move this candidate to new_pop_idxs and increment the associated niche count. */
            new_pop_idxs.push_back(selected_sol_idx);
            new_props.push_back(sol_props_[selected_sol_idx]);
            partial_front.erase(std::remove(partial_front.begin(), partial_front.end(), selected_sol_idx), partial_front.end());

            ref_niche_counts[ref]++;
        }

        ref_niche_counts_ = calcNicheCounts(ga, new_props);
        sol_props_ = new_props;

        return new_pop_idxs;
    }

} // namespace genetic_algorithm::selection::multi_objective