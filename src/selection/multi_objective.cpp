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

        auto pfronts = dtl::nonDominatedSort(ga.fitness_matrix());
        ranks_ = pfronts.ranks;
        dists_ = dtl::crowdingDistances(ga.fitness_matrix(), pfronts.idxs);
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
        std::vector<size_t> new_pop_idxs;
        new_pop_idxs.reserve(ga.population_size());

        auto pfronts = dtl::nonDominatedSort(combined_pop);
        ranks_ = pfronts.ranks;
        dists_ = dtl::crowdingDistances(combined_pop, pfronts.idxs);

        /* Keep track of the ranks of the candidates added to new_pop_idxs to avoid a second sort. */
        std::vector<size_t> new_ranks;
        new_ranks.reserve(ga.population_size());
        std::vector<double> new_dists;
        new_dists.reserve(ga.population_size());

        /* Add entire fronts while possible. */
        size_t front_idx = 0;
        while (new_pop_idxs.size() + pfronts.idxs[front_idx].size() <= ga.population_size())
        {
            for (const auto& idx : pfronts.idxs[front_idx])
            {
                new_pop_idxs.push_back(idx);
                new_ranks.push_back(ranks_[idx]);
                new_dists.push_back(dists_[idx]);
            }
            front_idx++;
        }

        /* Add the remaining candidates from the partial front if there is one. */
        if (new_pop_idxs.size() != ga.population_size())
        {
            /* Keep track of the candidates added from the partial front,
             * the crowding distances will only need to be updated for these candidates. */
            std::vector<size_t> newly_added_indices(ga.population_size() - new_pop_idxs.size());
            std::iota(newly_added_indices.begin(), newly_added_indices.end(), new_pop_idxs.size());

            std::vector<size_t> partial_front = pfronts.idxs[front_idx];

            std::sort(partial_front.begin(), partial_front.end(),
            [this](size_t lidx, size_t ridx)
            {
                return crowdedCompare(lidx, ridx);
            });

            for (auto it = partial_front.begin(); new_pop_idxs.size() != ga.population_size(); it++)
            {
                new_pop_idxs.push_back(*it);
                new_ranks.push_back(ranks_[*it]);
                new_dists.push_back(dists_[*it]);
            }

            auto changed_dists = dtl::crowdingDistances(combined_pop, { newly_added_indices });
            for (size_t idx : newly_added_indices)
            {
                new_dists[idx] = changed_dists[idx];
            }

        }
        ranks_ = new_ranks;
        dists_ = new_dists;

        return new_pop_idxs;
    }

    constexpr bool NSGA2::crowdedCompare(size_t lidx, size_t ridx) const noexcept
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

        ref_points_ = dtl::generateRefPoints(ga.population_size(), ga.num_objectives());

        ideal_point_ = std::vector(ga.num_objectives(), -std::numeric_limits<double>::max());
        updateIdealPoint(ideal_point_, ga.fitness_matrix());
        extreme_points_ = initExtremePoints(ga.fitness_matrix(), ideal_point_);
        nadir_point_ = nadirPoint(extreme_points_);

        sol_props_ = std::vector<CandidateInfo>(ga.population_size());
        associatePopWithRefs(ga.fitness_matrix());
        ref_niche_counts_ = calcNicheCounts(ga, sol_props_);

        auto pfronts = dtl::nonDominatedSort(ga.fitness_matrix());
        for (size_t i = 0; i < sol_props_.size(); i++)
        {
            sol_props_[i].rank = pfronts.ranks[i];
        }
    }

    void NSGA3::updateIdealPoint(Point& ideal_point, const FitnessMatrix& pop) noexcept
    {
        for (const auto& sol : pop)
        {
            for (size_t i = 0; i < ideal_point.size(); i++)
            {
                ideal_point[i] = std::max(ideal_point[i], sol[i]);
            }
        }
    }

    auto NSGA3::initExtremePoints(const FitnessMatrix& pop, const Point& ideal_point) -> std::vector<Point>
    {
        assert(!pop.empty());

        std::vector<Point> extreme_points(ideal_point.size(), Point(ideal_point.size()));

        /* Identify and update extreme points for each objective axis. */
        for (size_t i = 0; i < extreme_points.size(); i++)
        {
            std::vector<double> weights(extreme_points[i].size(), 1E-6);
            weights[i] = 1.0;

            /* Find the solution with the lowest Chebysev distance to the objective axis. */
            double dmin = std::numeric_limits<double>::infinity();
            std::vector<double> extreme_point = pop[0];
            for (const auto& sol : pop)
            {
                double dist = dtl::ASF(sol, ideal_point, weights);

                if (dist < dmin)
                {
                    dmin = dist;
                    extreme_point = sol;
                }
            }

            extreme_points[i] = extreme_point;
        }

        return extreme_points;
    }

    void NSGA3::updateExtremePoints(std::vector<Point>& extreme_points, const FitnessMatrix& pop, const Point& ideal_point)
    {
        assert(!pop.empty());

        /* Identify and update extreme points for each objective axis. */
        for (size_t i = 0; i < extreme_points.size(); i++)
        {
            std::vector<double> weights(extreme_points[i].size(), 1E-6);
            weights[i] = 1.0;

            /* Find the solution or extreme point with the lowest Chebysev distance to the objective axis. */
            double dmin = std::numeric_limits<double>::infinity();
            std::vector<double> new_extreme_point = extreme_points[i];

            for (const auto& sol : pop)
            {
                double dist = dtl::ASF(sol, ideal_point, weights);

                if (dist < dmin)
                {
                    dmin = dist;
                    new_extreme_point = sol;
                }
            }
            for (const auto& old_extreme_point : extreme_points)
            {
                double dist = dtl::ASF(old_extreme_point, ideal_point, weights);

                if (dist < dmin)
                {
                    dmin = dist;
                    new_extreme_point = old_extreme_point;
                }
            }
            extreme_points[i] = new_extreme_point;
        }
    }

    auto NSGA3::nadirPoint(const std::vector<Point>& extreme_points) -> Point
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
        nadir_point_ = nadirPoint(extreme_points_);

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