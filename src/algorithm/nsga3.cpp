/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "nsga3.hpp"
#include "nd_sort.hpp"
#include "../core/ga_info.hpp"
#include "../population/population.hpp"
#include "../utility/algorithm.hpp"
#include "../utility/functional.hpp"
#include "../utility/math.hpp"
#include "../utility/rng.hpp"
#include "../utility/utility.hpp"
#include <vector>
#include <utility>
#include <algorithm>
#include <functional>
#include <iterator>
#include <cmath>
#include <limits>
#include <cstddef>
#include <cassert>

namespace genetic_algorithm::algorithm
{
    using Point     = NSGA3::Point;
    using RefPoint  = NSGA3::RefPoint;
    using namespace dtl;

    std::vector<Point> NSGA3::generateRefPoints(size_t n, size_t dim)
    {
        assert(n > 0);
        assert(dim > 1);

        /* Generate reference point candidates randomly. */
        size_t ratio = std::max(10_sz, 2 * dim);
        std::vector<Point> candidates(ratio * n - 1);
        std::generate(candidates.begin(), candidates.end(), [&dim]() { return rng::randomSimplexPoint(dim); });

        std::vector<Point> refs;
        refs.reserve(n);
        refs.push_back(rng::randomSimplexPoint(dim));

        std::vector<double> min_distances(candidates.size(), std::numeric_limits<double>::infinity());
        while (refs.size() < n)
        {
            /* Calc the distance of each candidate to the closest ref point. */
            std::transform(GA_EXECUTION_UNSEQ, candidates.begin(), candidates.end(), min_distances.begin(), min_distances.begin(),
            [&refs](const Point& candidate, double dmin) noexcept
            {
                double d = detail::euclideanDistanceSq(candidate, refs.back());
                return std::min(dmin, d);
            });

            /* Add the candidate with highest min_distance to the refs. */
            size_t argmax = detail::argmax(min_distances.begin(), min_distances.end());
            refs.push_back(std::move(candidates[argmax]));

            /* Remove the added candidate and the corresponding min_distance. */
            std::swap(candidates[argmax], candidates.back());
            candidates.pop_back();
            std::swap(min_distances[argmax], min_distances.back());
            min_distances.pop_back();
        }

        return refs;
    }

    std::pair<size_t, double> NSGA3::findClosestRef(const std::vector<RefPoint>& refs, const Point& p)
    {
        assert(!refs.empty());

        auto distances = detail::map(refs, [&p](const RefPoint& ref) { return detail::perpendicularDistanceSq(ref.point, p); });
        auto idx = detail::argmin(distances.begin(), distances.end());

        return { idx, distances[idx] };
    }

    NSGA3::ASF NSGA3::getASF(std::vector<double> ideal_point, std::vector<double> weights) noexcept
    {
        assert(!weights.empty());
        assert(weights.size() == ideal_point.size());

        return [ideal_point = std::move(ideal_point), weights = std::move(weights)]
        (const std::vector<double>& fitness) noexcept
        {
            assert(fitness.size() == weights.size());

            double dmax = (ideal_point[0] - fitness[0]) / weights[0];
            for (size_t i = 0; i < fitness.size(); i++)
            {
                double d = (ideal_point[i] - fitness[i]) / weights[i];
                dmax = std::max(dmax, d);
            }

            return dmax;
        };
    }

    std::vector<double> NSGA3::weightVector(size_t dimensions, size_t axis)
    {
        assert(dimensions > axis);

        std::vector weights(dimensions, 1E-6);
        weights[axis] = 1.0;

        return weights;
    }

    void NSGA3::initialize(const GaInfo& ga)
    {
        assert(ga.num_objectives() > 1);
        assert(ga.population_size() != 0);

        auto& fitness_matrix = ga.fitness_matrix();

        /* Create the reference points. */
        auto refs = generateRefPoints(ga.population_size(), ga.num_objectives());
        ref_points_.reserve(refs.size());
        std::move(refs.begin(), refs.end(), std::back_inserter(ref_points_));

        ideal_point_ = detail::maxFitness(fitness_matrix.begin(), fitness_matrix.end());
        extreme_points_ = {};

        sol_info_.resize(ga.population_size());
        associatePopWithRefs(fitness_matrix.begin(), fitness_matrix.end());
        updateNicheCounts(sol_info_);

        auto pfronts = nonDominatedSort(fitness_matrix.begin(), fitness_matrix.end());
        for (const auto& [idx, rank] : pfronts)
        {
            sol_info_[idx].rank = rank;
        }
    }

    void NSGA3::updateIdealPoint(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        assert(std::distance(first, last) > 0);

        ideal_point_ = detail::max(ideal_point_, detail::maxFitness(first, last));
    }

    void NSGA3::updateExtremePoints(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        assert(std::distance(first, last) > 0);
        assert(first->size() == ideal_point_.size());

        const size_t popsize = size_t(last - first);

        std::vector<Point> new_extreme_points;
        new_extreme_points.reserve(extreme_points_.size());

        for (size_t dim = 0; dim < ideal_point_.size(); dim++)
        {
            auto ASFi = getASF(ideal_point_, weightVector(ideal_point_.size(), dim));

            std::vector<double> chebysev_distances;
            chebysev_distances.reserve(extreme_points_.size() + popsize);

            std::transform(first, last, std::back_inserter(chebysev_distances), ASFi);
            std::transform(extreme_points_.begin(), extreme_points_.end(), std::back_inserter(chebysev_distances), ASFi);

            size_t idx = detail::argmin(chebysev_distances.begin(), chebysev_distances.end());
            
            (idx >= popsize) ?
                new_extreme_points.push_back(extreme_points_[idx - popsize]) :
                new_extreme_points.push_back(first[idx]);
        }
        extreme_points_ = std::move(new_extreme_points);
    }

    auto NSGA3::findNadirPoint(const std::vector<Point>& extreme_points) -> Point
    {
        assert(!extreme_points.empty());

        /* Nadir point estimate = minimum of extreme points along each objective axis. */
        Point nadir_point = extreme_points[0];
        for (size_t i = 1; i < extreme_points.size(); i++)
        {
            nadir_point = detail::min(nadir_point, extreme_points[i]);
        }

        return nadir_point;
    }

    auto NSGA3::normalize(const FitnessVector& fvec) -> FitnessVector
    {
        FitnessVector fnorm(fvec.size());

        for (size_t i = 0; i < fnorm.size(); i++)
        {
            fnorm[i] = (fvec[i] - ideal_point_[i]) / std::max(ideal_point_[i] - nadir_point_[i], 1E-6);
        }

        return fnorm;
    }

    void NSGA3::associatePopWithRefs(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        assert(std::distance(first, last) > 0);
        assert(std::all_of(first, last, [first](const FitnessVector& sol) { return sol.size() == first[0].size(); }));
        assert(sol_info_.size() == size_t(last - first));

        updateIdealPoint(first, last);
        updateExtremePoints(first, last);
        nadir_point_ = findNadirPoint(extreme_points_);

        std::transform(GA_EXECUTION_UNSEQ, first, last, sol_info_.begin(), sol_info_.begin(),
        [this](const FitnessVector& f, CandidateInfo& props)
        {
            std::tie(props.ref_idx, props.ref_dist) = findClosestRef(ref_points_, normalize(f));

            return props;
        });
    }

    void NSGA3::updateNicheCounts(const std::vector<CandidateInfo>& props) noexcept
    {
        assert(!ref_points_.empty());

        std::for_each(ref_points_.begin(), ref_points_.end(), [](RefPoint& ref) { ref.niche_count = 0; });
        std::for_each(props.begin(), props.end(), [this](const CandidateInfo& info) { ref_points_[info.ref_idx].niche_count++; });
    }

    std::vector<size_t> NSGA3::nextPopulation(const GaInfo& ga,
                                              FitnessMatrix::const_iterator parents_first,
                                              FitnessMatrix::const_iterator children_first,
                                              FitnessMatrix::const_iterator children_last)
    {
        assert(ga.num_objectives() > 1);
        assert(size_t(children_first - parents_first) == ga.population_size());
        assert(size_t(children_last - parents_first) >= ga.population_size());
        assert(std::all_of(parents_first, children_last, [&ga](const FitnessVector& f) { return f.size() == ga.num_objectives(); }));

        GA_UNUSED(children_first);

        std::vector<size_t> new_pop;
        new_pop.reserve(ga.population_size());

        auto pfronts = dtl::nonDominatedSort(parents_first, children_last);

        sol_info_ = std::vector<CandidateInfo>(size_t(children_last - parents_first));
        for (const auto& [idx, rank] : pfronts)
        {
            sol_info_[idx].rank = rank;
        } // after this, only the rank is updated, the rest are nonsense (ref idx, ref dist)
        associatePopWithRefs(parents_first, children_last);

        std::vector<CandidateInfo> new_info;
        new_info.reserve(ga.population_size());

        /* Add entire fronts while possible. */
        auto first = pfronts.begin();
        auto last = dtl::nextFrontBegin(first, pfronts.end());
        while (new_pop.size() + std::distance(first, last) <= ga.population_size())
        {
            for (; first != last; first++)
            {
                new_pop.push_back(first->idx);
                new_info.push_back(std::move(sol_info_[first->idx]));
            }
            last = dtl::nextFrontBegin(first, pfronts.end());
        }
        updateNicheCounts(new_info);

        /* Add remaining candidates from the partial front if there is one. */
        dtl::ParetoFronts partial_front(first, last);
        while (new_pop.size() != ga.population_size())
        {
            /* Find the lowest niche count in the partial front. */
            auto min = std::min_element(partial_front.begin(), partial_front.end(),
            [this](const FrontInfo& lhs, const FrontInfo& rhs) noexcept
            {
                return nicheCountOf(lhs.idx) < nicheCountOf(rhs.idx);
            });
            size_t min_cnt = nicheCountOf(min->idx);

            /* Find the reference points with minimal niche counts, and pick one. */
            // TODO: find_all_v
            // TODO: this loop takes a long time
            std::vector<size_t> refs;
            for (const auto& [idx, rank] : partial_front)
            {
                size_t ref_idx = sol_info_[idx].ref_idx;
                if (ref_points_[ref_idx].niche_count == min_cnt &&
                    !detail::contains(refs.begin(), refs.end(), ref_idx))
                {
                    refs.push_back(ref_idx);
                }
            }
            size_t ref = rng::randomElement(refs);

            /* Find the idx of the closest sol in the partial front associated with this ref point. */
            size_t selected_idx = partial_front[0].idx;
            double dmin = std::numeric_limits<double>::infinity();
            for (const auto& [idx, rank] : partial_front)
            {
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
                [selected_idx](const FrontInfo& elem) { return elem.idx == selected_idx; });
            partial_front.erase(selected);

            ref_points_[ref].niche_count++;
        }
        sol_info_ = std::move(new_info);

        return new_pop;
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

    bool NSGA3::nichedCompare(size_t lhs, size_t rhs) const noexcept
    {
        if (sol_info_[lhs].ref_idx == sol_info_[rhs].ref_idx)
        {
            if (sol_info_[lhs].rank != sol_info_[rhs].rank)
            {
                return sol_info_[lhs].rank < sol_info_[rhs].rank;
            }
            else
            {
                return sol_info_[lhs].ref_dist < sol_info_[rhs].ref_dist;
            }
        }

        return rng::randomBool();
    }

    size_t NSGA3::select(const GaInfo&, const FitnessMatrix& pop) const
    {
        assert(!pop.empty());

        size_t idx1 = rng::randomIdx(pop);
        size_t idx2 = rng::randomIdx(pop);

        return nichedCompare(idx1, idx2) ? idx1 : idx2;
    }

} // namespace genetic_algorithm::algorithm