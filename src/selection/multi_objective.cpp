/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "multi_objective.hpp"
#include "selection_dtl.hpp"
#include "../core/ga_info.hpp"
#include "../utility/rng.hpp"
#include "../utility/algorithm.hpp"
#include <algorithm>
#include <vector>
#include <set>
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

        auto& fmat = ga.fitness_matrix();
        auto pfronts = dtl::nonDominatedSort(fmat);

        ranks_ = dtl::paretoRanks(pfronts);
        dists_ = dtl::crowdingDistances(fmat, std::move(pfronts));
    }

    size_t NSGA2::select(const GaInfo&, const FitnessMatrix& pop) const
    {
        assert(!pop.empty() && pop.size() == ranks_.size());

        size_t idx1 = rng::randomIdx(pop);
        size_t idx2 = rng::randomIdx(pop);

        return crowdedCompare(idx1, idx2) ? idx1 : idx2;
    }

    std::vector<size_t> NSGA2::nextPopulation(const GaInfo& ga, const FitnessMatrix& fmat)
    {
        assert(ga.population_size() <= fmat.size());
        assert(std::all_of(fmat.begin(), fmat.end(), [&ga](const FitnessVector& f) { return f.size() == ga.num_objectives(); }));

        std::vector<size_t> new_pop;
        new_pop.reserve(ga.population_size());

        auto pfronts = dtl::nonDominatedSort(fmat);
        ranks_ = dtl::paretoRanks(pfronts);
        dists_ = dtl::crowdingDistances(fmat, pfronts);

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
            //first = last;
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

            auto changed_dists = dtl::crowdingDistances(fmat, { partial_front });

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

        auto& fitness_matrix = ga.fitness_matrix();

        auto refs = dtl::generateRefPoints(ga.population_size(), ga.num_objectives());
        ref_points_.reserve(refs.size());
        for (auto& ref : refs)
        {
            ref_points_.emplace_back(std::move(ref));
        }

        ideal_point_ = detail::maxFitness(fitness_matrix.begin(), fitness_matrix.end());
        extreme_points_ = {};

        sol_info_ = std::vector<CandidateInfo>(fitness_matrix.size());
        associatePopWithRefs(sol_info_, fitness_matrix, ref_points_);
        updateNicheCounts(ref_points_, sol_info_);

        auto pfronts = dtl::nonDominatedSort(fitness_matrix);
        for (const auto& elem : pfronts)
        {
            size_t idx = elem.first;
            size_t rank = elem.second;

            sol_info_[idx].rank = rank;
        }
    }

    void NSGA3::updateIdealPoint(Point& ideal_point, const FitnessMatrix& fmat)
    {
        assert(!fmat.empty());
        assert(ideal_point.size() == fmat[0].size());

        auto fmax = detail::maxFitness(fmat.begin(), fmat.end());
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
            auto chebysev_distances_e = detail::map(extreme_points, ASFi);

            auto fmin = std::min_element(chebysev_distances_f.begin(), chebysev_distances_f.end());
            auto emin = std::min_element(chebysev_distances_e.begin(), chebysev_distances_e.end());

            if (emin != chebysev_distances_e.end() && *emin < *fmin)
            {
                size_t argmin = size_t(emin - chebysev_distances_e.begin());
                new_extreme_points.push_back(extreme_points[argmin]);
            }
            else
            {
                size_t argmin = size_t(fmin - chebysev_distances_f.begin());
                new_extreme_points.push_back(fmat[argmin]);
            }
        }
        extreme_points = std::move(new_extreme_points);
    }

    auto NSGA3::findNadirPoint(const std::vector<Point>& extreme_points) -> Point
    {
        assert(!extreme_points.empty());
        assert(std::all_of(extreme_points.begin(), extreme_points.end(), [&extreme_points](const Point& p) { return extreme_points[0].size() == p.size(); }));

        /* Nadir point estimate = minimum of extreme points along each objective axis. */
        Point nadir = extreme_points[0];

        for (size_t i = 1; i < extreme_points.size(); i++)
        {
            for (size_t j = 0; j < nadir.size(); j++)
            {
                nadir[j] = std::min(nadir[j], extreme_points[i][j]);
            }
        }

        return nadir;
    }

    bool NSGA3::nichedCompare(size_t lhs, size_t rhs) const noexcept
    {
        if (sol_info_[lhs].rank != sol_info_[rhs].rank)
        {
            return sol_info_[lhs].rank < sol_info_[rhs].rank;
        }
        else if (nicheCountOf(lhs) != nicheCountOf(rhs))
        {
            return nicheCountOf(lhs) < nicheCountOf(rhs);
        }
        else
        {
            return sol_info_[lhs].ref_dist < sol_info_[rhs].ref_dist;
        }
    }

    size_t NSGA3::select(const GaInfo&, const FitnessMatrix& pop) const
    {
        assert(!pop.empty());

        size_t idx1 = rng::randomIdx(pop);
        size_t idx2 = rng::randomIdx(pop);

        return nichedCompare(idx1, idx2) ? idx1 : idx2;
    }

    auto NSGA3::normalize(const FitnessVector& fvec, const Point& ideal_point, const Point& nadir_point) -> FitnessVector
    {
        FitnessVector fnorm(fvec.size());

        for (size_t i = 0; i < fnorm.size(); i++)
        {
            double f = (ideal_point[i] - fvec[i]) / std::max(ideal_point[i] - nadir_point[i], 1E-6);
            fnorm[i] = f;
        }

        return fnorm;
    }

    void NSGA3::associatePopWithRefs(std::vector<CandidateInfo>& props, const FitnessMatrix& fmat, const std::vector<RefPoint>& refs)
    {
        assert(!fmat.empty());
        assert(std::all_of(fmat.begin(), fmat.end(), [&fmat](const FitnessVector& sol) { return sol.size() == fmat[0].size(); }));
        assert(props.size() == fmat.size());
        // TODO: something is probably wrong in this function, or in something that is called in here
        updateIdealPoint(ideal_point_, fmat);
        updateExtremePoints(extreme_points_, fmat, ideal_point_);
        nadir_point_ = findNadirPoint(extreme_points_);

        auto fnorm = detail::map(fmat,
        [this](const FitnessVector& fvec)
        {
            return normalize(fvec, ideal_point_, nadir_point_);
        });

        /* Associate each candidate with the closest reference point. */
        std::transform(GA_EXECUTION_UNSEQ, fnorm.begin(), fnorm.end(), props.begin(), props.begin(),
        [&refs](const FitnessVector& f, CandidateInfo& info) -> CandidateInfo
        {
            std::tie(info.ref_idx, info.ref_dist) = dtl::findClosestRef(refs, f);

            return info;
        });
    }

    size_t& NSGA3::nicheCountOf(const CandidateInfo& info) noexcept
    {
        return ref_points_[info.ref_idx].niche_count;
    }
    size_t& NSGA3::nicheCountOf(size_t sol_idx) noexcept
    {
        return ref_points_[sol_info_[sol_idx].ref_idx].niche_count;
    }
    const size_t& NSGA3::nicheCountOf(size_t sol_idx) const noexcept
    {
        return ref_points_[sol_info_[sol_idx].ref_idx].niche_count;
    }
    const size_t& NSGA3::nicheCountOf(const CandidateInfo& info) const noexcept
    {
        return ref_points_[info.ref_idx].niche_count;
    }

    void NSGA3::updateNicheCounts(std::vector<RefPoint>& refs, const std::vector<CandidateInfo>& props) noexcept
    {
        assert(!refs.empty());

        for (auto& ref : refs)
        {
            ref.niche_count = 0;
        }
        for (const auto& info : props)
        {
            refs[info.ref_idx].niche_count++;
        }
    }

    std::vector<size_t> NSGA3::nextPopulation(const GaInfo& ga, const FitnessMatrix& fmat)
    {
        assert(fmat.size() >= ga.population_size());

        std::vector<size_t> new_pop;
        new_pop.reserve(ga.population_size());

        auto pfronts = dtl::nonDominatedSort(fmat);

        sol_info_ = std::vector<CandidateInfo>(fmat.size());
        for (size_t i = 0; i < pfronts.size(); i++)
        {
            auto& [idx, rank] = pfronts[i];
            sol_info_[idx].rank = rank; // after this, only the rank is updated, the rest are nonsense (ref idx, ref dist)
        }
        associatePopWithRefs(sol_info_, fmat, ref_points_);

        std::vector<CandidateInfo> new_info;
        new_info.reserve(ga.population_size());

        /* Add entire fronts while possible. */
        auto first = pfronts.begin();
        auto last  = dtl::nextFrontBegin(first, pfronts.end());
        while (new_pop.size() + std::distance(first, last) <= ga.population_size())
        {
            for (; first != last; first++)
            {
                size_t idx = first->first;
                new_pop.push_back(idx);
                new_info.push_back(std::move(sol_info_[idx]));
            }
            last = dtl::nextFrontBegin(first, pfronts.end());
        }
        updateNicheCounts(ref_points_, new_info);

        /* Add remaining candidates from the partial front if there is one. */
        dtl::ParetoFronts partial_front(first, last);
        while (new_pop.size() != ga.population_size())
        {
            /* Find the lowest niche count in the partial front. */
            auto min = std::min_element(partial_front.begin(), partial_front.end(),
            [this](const std::pair<size_t, size_t>& lhs, const std::pair<size_t, size_t>& rhs) noexcept
            {
                size_t lidx = lhs.first;
                size_t ridx = rhs.first;

                return nicheCountOf(lidx) < nicheCountOf(ridx);
            });
            size_t min_cnt = nicheCountOf(min->first);

            /* Find the reference points with minimal niche counts, and pick one. */
            // TODO: find_all_v
            // TODO: this loop takes a long time
            std::vector<size_t> refs;
            for (const auto& elem : partial_front)
            {
                size_t ref_idx = sol_info_[elem.first].ref_idx;
                if (ref_points_[ref_idx].niche_count == min_cnt &&
                    !detail::contains(refs.begin(), refs.end(), ref_idx))
                {
                    refs.push_back(ref_idx);
                }
            }
            size_t ref = rng::randomElement(refs);

            /* Find the idx of the closest sol in the partial front associated with this ref point. */
            size_t selected_idx = partial_front[0].first;
            double dmin = std::numeric_limits<double>::infinity();
            for (const auto& elem : partial_front)
            {
                size_t idx = elem.first;

                if (sol_info_[idx].ref_idx == ref &&
                    sol_info_[idx].ref_dist < dmin)
                {
                    dmin = sol_info_[idx].ref_dist;
                    selected_idx = idx;
                }
            }

            /* Move this candidate to new_pop_idxs and increment the associated niche count. */
            new_pop.push_back(selected_idx);
            new_info.push_back(sol_info_[selected_idx]);

            auto selected = std::find_if(partial_front.begin(), partial_front.end(),
                [selected_idx](const std::pair<size_t, size_t>& elem) { return elem.first == selected_idx; });
            partial_front.erase(selected);

            ref_points_[ref].niche_count++;
        }
        sol_info_ = std::move(new_info);

        return new_pop;
    }

} // namespace genetic_algorithm::selection::multi_objective