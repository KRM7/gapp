/* Copyright (c) 2022 Kriszti�n Rug�si. Subject to the MIT License. */

#include "reference_lines.hpp"
#include "../utility/algorithm.hpp"
#include "../utility/functional.hpp"
#include "../utility/rng.hpp"
#include "../utility/qrng.hpp"
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

        const double sum = std::reduce(point.begin(), point.end(), 0.0);
        std::transform(point.begin(), point.end(), point.begin(), detail::divide_by(sum));

        return point;
    }

    /* Generate n random points on a unit simplex in dim dimensions. */
    static std::vector<Point> randomSimplexPoints(size_t dim, size_t n)
    {
        assert(dim > 0);

        std::vector<Point> reference_points;
        reference_points.reserve(n);

        while (n--) reference_points.push_back(randomSimplexPoint(dim));

        return reference_points;
    }

    /* Generate n quasirandom points on a unit simplex in dim dimensions. */
    static std::vector<Point> quasirandomSimplexPoints(size_t dim, size_t n)
    {
        assert(dim > 0);

        detail::QuasiRandom qrng(dim - 1);

        std::vector<Point> points;

        for (size_t p = 0; p < n; p++)
        {
            Point point(dim + 1, 0.0);
            Point qrand = qrng();
            std::copy(qrand.begin(), qrand.end(), point.begin() + 1);
            point.back() = 1.0;

            for (size_t i = dim - 1; i > 1; i--)
            {
                for (size_t j = 1; j < i; j++)
                {
                    if (point[i] < point[j])
                    {
                        for (size_t k = j; k < i + 1; k++)
                        {
                            point[k] = point[j - 1] + point[i + 1] - point[k];
                        }
                    }
                }
            }

            for (size_t i = 0; i < dim; i++)
            {
                point[i] = point[i + 1] - point[i];
            }
            point.pop_back();

            points.push_back(point);
        }

        return points;
    }

    /* Generate reference points by picking n points from a set of randomly generated points. */
    static std::vector<Point> generateRandomRefpointsPick(size_t dim, size_t n)
    {
        assert(dim > 0);
        assert(n > 0);

        const size_t npoints = n * std::max(10_sz, 2 * dim);
        std::vector<Point> candidate_points = randomSimplexPoints(dim, npoints);

        std::vector<Point> points;
        points.reserve(n);
        points.push_back(candidate_points.back());
        candidate_points.pop_back();

        auto min_distances = detail::map(candidate_points, [&](const Point& p) { return math::euclideanDistanceSq(p, points.back()); });

        while (points.size() < n)
        {
            size_t idx = detail::argmax(min_distances.begin(), min_distances.end());
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
                const double dist = math::euclideanDistanceSq(candidate, last_point);
                return std::min(current_min, dist);
            });
        }

        return points;
    }


    std::vector<Point> generateReferencePoints(size_t dim, size_t n)
    {
        //std::vector<Point> points = generateRandomRefpointsPick(dim, n);
        std::vector<Point> points = quasirandomSimplexPoints(dim, n);
        std::transform(points.begin(), points.end(), points.begin(), [](Point p) noexcept { return math::normalizeVector(std::move(p)); });

        return points;
    }

} // namespace genetic_algorithm::algorithm::dtl