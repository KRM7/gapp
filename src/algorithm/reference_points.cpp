/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "reference_points.hpp"
#include "../utility/algorithm.hpp"
#include "../utility/functional.hpp"
#include "../utility/rng.hpp"
#include "../utility/math.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <numeric>
#include <functional>
#include <limits>
#include <cmath>
#include <cassert>

namespace genetic_algorithm::algorithm::dtl
{
    /* Sample a point from a uniform distribution on a unit simplex in dim dimensions. */
    static Point randomSimplexPoint(size_t dim)
    {
        assert(dim > 0);

        Point point(dim);
        std::generate(point.begin(), point.end(), []{ return -std::log(rng::randomReal()); });

        double isum = 1.0 / std::reduce(point.begin(), point.end(), 0.0);
        std::transform(point.begin(), point.end(), point.begin(), [&](double coord) { return isum * coord; });

        return point;
    }

    /* Generate n random reference points in dim dimensions. */
    static std::vector<Point> generateRandomRefpoints(size_t dim, size_t n)
    {
        assert(dim > 0);

        std::vector<Point> reference_points;
        reference_points.reserve(n);

        while (n--) reference_points.push_back(randomSimplexPoint(dim));

        return reference_points;
    }

    /* Generate reference points by picking n points from a set of randomly generated points. */
    static std::vector<Point> generateRandomRefpointsPick(size_t dim, size_t n)
    {
        assert(dim > 0);
        assert(n > 0);

        const size_t npoints = n * std::max(10_sz, 2 * dim);
        std::vector<Point> candidate_points = generateRandomRefpoints(dim, npoints);

        std::vector<Point> points;
        points.reserve(n);
        points.push_back(candidate_points.back());
        candidate_points.pop_back();

        auto min_distances = detail::map(candidate_points, [&](const Point& p) { return detail::euclideanDistanceSq(p, points.back()); });

        while (points.size() < n)
        {
            size_t idx = detail::argmax(min_distances.begin(), min_distances.begin(), min_distances.end());
            points.push_back(std::move(candidate_points[idx]));

            /* Remove the added candidate and the corresponding min_distance. */
            std::swap(candidate_points[idx], candidate_points.back());
            candidate_points.pop_back();
            std::swap(min_distances[idx], min_distances.back());
            min_distances.pop_back();

            /* Calc the distance of each candidate to the closest ref point. */
            std::transform(GA_EXECUTION_UNSEQ, candidate_points.begin(), candidate_points.end(), min_distances.begin(), min_distances.begin(),
            [&last_point = points.back()](const Point& candidate, double current_min) noexcept
            {
                double dist = detail::euclideanDistanceSq(candidate, last_point);
                return std::min(current_min, dist);
            });
        }

        return points;
    }


    ReferencePoints generateReferencePoints(size_t dim, size_t n)
    {
        std::vector<Point> points = generateRandomRefpointsPick(dim, n);

        ReferencePoints refs;
        refs.reserve(points.size());
        std::move(points.begin(), points.end(), std::back_inserter(refs));

        return refs;
    }

    std::pair<size_t, double> findClosestRef(const std::vector<ReferencePoint>& refs, const Point& p)
    {
        assert(!refs.empty());
        assert(std::all_of(refs.begin(), refs.end(), [&](const ReferencePoint& ref) { return ref.point.size() == p.size(); }));

        auto distances = detail::map(refs, [&](const ReferencePoint& ref) { return detail::perpendicularDistanceSq(p, ref.point); });
        auto closest_idx = detail::argmin(distances.begin(), distances.begin(), distances.end());

        return { closest_idx, distances[closest_idx] };
    }

} // namespace genetic_algorithm::algorithm::dtl