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
#include <iterator>
#include <cmath>
#include <limits>
#include <cstddef>
#include <cassert>

namespace genetic_algorithm::algorithm
{
    using Point     = NSGA3::Point;
    using RefPoint  = NSGA3::RefPoint;

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

        auto distances = detail::map(refs,
        [&p](const RefPoint& ref) noexcept
        {
            return detail::perpendicularDistanceSq(ref.point, p);
        });

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

            double dmax = std::abs(fitness[0] - ideal_point[0]) / weights[0]; // TODO
            for (size_t i = 0; i < fitness.size(); i++)
            {
                double d = std::abs(fitness[i] - ideal_point[i]) / weights[i]; // TODO
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

        auto refs = generateRefPoints(ga.population_size(), ga.num_objectives());
        ref_points_.reserve(refs.size());
        for (auto& ref : refs)
        {
            ref_points_.emplace_back(std::move(ref));
        }

        ideal_point_ = detail::maxFitness(fitness_matrix.begin(), fitness_matrix.end());
        extreme_points_ = {};

        sol_info_ = std::vector<CandidateInfo>(fitness_matrix.size());
        associatePopWithRefs(sol_info_, fitness_matrix.begin(), fitness_matrix.end(), ref_points_);
        updateNicheCounts(ref_points_, sol_info_);

        auto pfronts = dtl::nonDominatedSort(fitness_matrix.begin(), fitness_matrix.end());
        for (const auto& elem : pfronts)
        {
            size_t idx = elem.first;
            size_t rank = elem.second;

            sol_info_[idx].rank = rank;
        }
    }

    void NSGA3::updateIdealPoint(Point& ideal_point, FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        assert(std::distance(first, last) > 0);
        assert(ideal_point.size() == first[0].size());

        auto fmax = detail::maxFitness(first, last);
        for (size_t i = 0; i < ideal_point.size(); i++)
        {
            ideal_point[i] = std::max(ideal_point[i], fmax[i]);
        }
    }

    void NSGA3::updateExtremePoints(std::vector<Point>& extreme_points, const Point& ideal_point, FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        assert(std::distance(first, last) > 0);
        assert(first->size() == ideal_point.size());

        std::vector<Point> new_extreme_points;
        new_extreme_points.reserve(extreme_points.size());

        for (size_t dim = 0; dim < ideal_point.size(); dim++)
        {
            auto w = weightVector(ideal_point.size(), dim);
            auto ASFi = getASF(ideal_point, std::move(w));

            std::vector chebysev_distances_f(size_t(last - first), 0.0);
            std::vector chebysev_distances_e(extreme_points.size(), 0.0);

            std::transform(first, last, chebysev_distances_f.begin(), ASFi);
            std::transform(extreme_points.begin(), extreme_points.end(), chebysev_distances_e.begin(), ASFi);

            auto fmin = std::min_element(chebysev_distances_f.begin(), chebysev_distances_f.end());
            auto emin = std::min_element(chebysev_distances_e.begin(), chebysev_distances_e.end());

            if (!extreme_points.empty() && *emin < *fmin)
            {
                size_t argmin = size_t(emin - chebysev_distances_e.begin());
                new_extreme_points.push_back(extreme_points[argmin]);
            }
            else
            {
                size_t argmin = size_t(fmin - chebysev_distances_f.begin());
                new_extreme_points.push_back(first[argmin]);
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
            for (size_t dim = 0; dim < nadir.size(); dim++)
            {
                nadir[dim] = std::min(nadir[dim], extreme_points[i][dim]);
            }
        }

        return nadir;
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

    void NSGA3::associatePopWithRefs(std::vector<CandidateInfo>& props, FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last, const std::vector<RefPoint>& refs)
    {
        assert(std::distance(first, last) > 0);
        assert(std::all_of(first, last, [first](const FitnessVector& sol) { return sol.size() == first[0].size(); }));
        assert(props.size() == size_t(last - first));
        // TODO: something is probably wrong in this function, or in something that is called in here
        updateIdealPoint(ideal_point_, first, last);
        updateExtremePoints(extreme_points_, ideal_point_, first, last);
        nadir_point_ = findNadirPoint(extreme_points_);

        FitnessMatrix fnorm;
        fnorm.reserve(size_t(last - first));

        std::transform(first, last, std::back_inserter(fnorm),
        [this](const FitnessVector& fvec)
        {
            return normalize(fvec, ideal_point_, nadir_point_);
        });

        /* Associate each candidate with the closest reference point. */
        std::transform(GA_EXECUTION_UNSEQ, fnorm.begin(), fnorm.end(), props.begin(), props.begin(),
        [&refs](const FitnessVector& f, CandidateInfo& info) -> CandidateInfo
        {
            std::tie(info.ref_idx, info.ref_dist) = findClosestRef(refs, f);

            return info;
        });
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

    std::vector<size_t> NSGA3::nextPopulation(const GaInfo& ga,
                                              FitnessMatrix::const_iterator parents_first,
                                              FitnessMatrix::const_iterator children_first,
                                              FitnessMatrix::const_iterator children_last)
    {
        assert(size_t(children_first - parents_first) == ga.population_size());
        assert(size_t(children_last - parents_first) >= ga.population_size());
        assert(std::all_of(parents_first, children_last, [&ga](const FitnessVector& f) { return f.size() == ga.num_objectives(); }));

        std::vector<size_t> new_pop;
        new_pop.reserve(ga.population_size());

        auto pfronts = dtl::nonDominatedSort(parents_first, children_last);

        sol_info_ = std::vector<CandidateInfo>(size_t(children_last - parents_first));
        for (size_t i = 0; i < pfronts.size(); i++)
        {
            auto& [idx, rank] = pfronts[i];
            sol_info_[idx].rank = rank; // after this, only the rank is updated, the rest are nonsense (ref idx, ref dist)
        }
        associatePopWithRefs(sol_info_, parents_first, children_last, ref_points_);

        std::vector<CandidateInfo> new_info;
        new_info.reserve(ga.population_size());

        /* Add entire fronts while possible. */
        auto first = pfronts.begin();
        auto last = dtl::nextFrontBegin(first, pfronts.end());
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
        if (sol_info_[lhs].rank != sol_info_[rhs].rank)
        {
            return sol_info_[lhs].rank < sol_info_[rhs].rank;
        }
        // tiebreak1
        else if (nicheCountOf(lhs) != nicheCountOf(rhs))
        {
            return nicheCountOf(lhs) < nicheCountOf(rhs);
        }
        // tiebreak2
        else
        {
            return sol_info_[lhs].ref_dist < sol_info_[rhs].ref_dist;
        }
    }

    size_t NSGA3::select(const GaInfo&, const FitnessMatrix& pop)
    {
        assert(!pop.empty());

        size_t idx1 = rng::randomIdx(pop);
        size_t idx2 = rng::randomIdx(pop);

        return nichedCompare(idx1, idx2) ? idx1 : idx2;
    }

} // namespace genetic_algorithm::algorithm