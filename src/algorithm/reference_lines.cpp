/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "reference_lines.hpp"
#include "../utility/algorithm.hpp"
#include "../utility/functional.hpp"
#include "../utility/qrng.hpp"
#include "../utility/math.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <execution>
#include <numeric>
#include <functional>
#include <iterator>
#include <utility>
#include <cmath>

namespace genetic_algorithm::algorithm::reflines
{
    /*
    * The unit-hypercube -> unit-simplex transformations used for the quasirandom points are based on:
    * 
    *   Pillards, Tim, and Ronald Cools. "Transforming low-discrepancy sequences from a cube to a simplex."
    *   Journal of computational and applied mathematics 174, no. 1 (2005): 29-42.
    */

    /* Transform a point from the n-dimensional unit hypercube to the n-dimensional unit simplex. */
    static inline Point simplexMappingLog(Point&& point)
    {
        GA_ASSERT(std::all_of(point.begin(), point.end(), detail::between(0.0, 1.0)));

        std::transform(point.begin(), point.end(), point.begin(), [](double p) { return -std::log(std::max(p, math::small<double>)); });
        const double sum = std::reduce(point.begin(), point.end(), 0.0);
        std::transform(point.begin(), point.end(), point.begin(), detail::divide_by(sum));

        return point;
    }

    /* Transform a point from the n-dimensional unit hypercube to the (n+1)-dimensional unit simplex. */
    static inline Point simplexMappingSort(Point&& point)
    {
        GA_ASSERT(std::all_of(point.begin(), point.end(), detail::between(0.0, 1.0)));

        point.push_back(1.0);
        std::sort(point.begin(), point.end());
        std::adjacent_difference(point.begin(), point.end(), point.begin());

        return point;
    }

    /* Transform a point from the n-dimensional unit hypercube to the (n+1)-dimensional unit simplex. */
    static inline Point simplexMappingRoot(Point&& point)
    {
        GA_ASSERT(std::all_of(point.begin(), point.end(), detail::between(0.0, 1.0)));

        if (point.empty()) return { 1.0 };

        point.back() = std::pow(point.back(), 1.0 / point.size());

        for (auto it = std::next(point.rbegin()); it != point.rend(); ++it)
        {
            *it = *std::prev(it) * std::pow(*it, 1.0 / std::distance(it, point.rend()));
        }
        point.push_back(1.0);

        std::adjacent_difference(point.begin(), point.end(), point.begin());

        return point;
    }

    /* Transform a point from the n-dimensional unit hypercube to the (n+1)-dimensional unit simplex. */
    static inline Point simplexMappingMirror(Point&& point)
    {
        GA_ASSERT(std::all_of(point.begin(), point.end(), detail::between(0.0, 1.0)));

        point.push_back(1.0);
        if (point.size() == 1) return point;

        for (auto last = std::prev(point.end(), 2); ; )
        {
            bool has_lower = false;

            for (auto first = point.begin(); first != last; ++first)
            {
                if (*first > *last)
                {
                    has_lower = true;
                    const double high = *std::next(last);
                    const double low  = (first != point.begin()) ? *std::prev(first) : 0.0;

                    *first = low + high - *first;
                    *last  = low + high - *last;
                }
            }
            if (!has_lower && last == point.begin()) break;
            if (!has_lower) --last;
        }
        std::adjacent_difference(point.begin(), point.end(), point.begin());

        return point;
    }


    template<auto Mapping>
    struct SimplexMappingTraits
    {
        static constexpr size_t input_dim(size_t output_dim) noexcept { return output_dim - 1; }
    };

    template<>
    struct SimplexMappingTraits<simplexMappingLog>
    {
        static constexpr size_t input_dim(size_t output_dim) noexcept { return output_dim; }
    };


    template<auto Mapping>
    static std::vector<Point> quasirandomSimplexPoints(size_t dim, size_t num_points)
    {
        std::vector<Point> points(num_points);
        size_t input_dim = SimplexMappingTraits<Mapping>::input_dim(dim);

        if (dim == 0) return points;

        rng::QuasiRandom qrng{ input_dim };

        for (size_t i = 0; i < num_points; i++)
        {
            Point hypercubePoint = qrng();
            Point simplexPoint = Mapping(std::move(hypercubePoint));

            points[i] = std::move(simplexPoint);
        }

        return points;
    }


    std::vector<Point> quasirandomSimplexPointsSort(size_t dim, size_t num_points)
    {
        return quasirandomSimplexPoints<simplexMappingSort>(dim, num_points);
    }

    std::vector<Point> quasirandomSimplexPointsRoot(size_t dim, size_t num_points)
    {
        return quasirandomSimplexPoints<simplexMappingRoot>(dim, num_points);
    }

    std::vector<Point> quasirandomSimplexPointsMirror(size_t dim, size_t num_points)
    {
        return quasirandomSimplexPoints<simplexMappingMirror>(dim, num_points);
    }

    std::vector<Point> quasirandomSimplexPointsLog(size_t dim, size_t num_points)
    {
        return quasirandomSimplexPoints<simplexMappingLog>(dim, num_points);
    }


    std::vector<Point> pickSparseSubset(size_t dim, size_t num_points, RefLineGenerator generator, Positive<size_t> k)
    {
        if (num_points == 0) return {};

        std::vector<Point> candidate_points = generator(dim, k * num_points);

        std::vector<Point> points;
        points.reserve(num_points);
        points.push_back(candidate_points.back());
        candidate_points.pop_back();

        auto min_distances = detail::map(candidate_points, [&](const Point& p) noexcept { return math::euclideanDistanceSq(p, points.back()); });

        while (points.size() < num_points)
        {
            const size_t idx = detail::argmax(min_distances.begin(), min_distances.end());
            points.push_back(std::move(candidate_points[idx]));

            /* Remove the added candidate and the corresponding min_distance. */
            std::swap(candidate_points[idx], candidate_points.back());
            candidate_points.pop_back();
            std::swap(min_distances[idx], min_distances.back());
            min_distances.pop_back();

            /* Calc the distance of each candidate to the closest ref point. */
            std::transform(GA_EXECUTION_UNSEQ, candidate_points.begin(), candidate_points.end(), min_distances.begin(), min_distances.begin(),
            [&](const Point& candidate, double current_min) noexcept
            {
                const double dist = math::euclideanDistanceSq(candidate, points.back());
                return std::min(current_min, dist);
            });
        }

        return points;
    }

} // namespace genetic_algorithm::algorithm::reflines