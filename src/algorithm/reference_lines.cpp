/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "reference_lines.hpp"
#include "../utility/algorithm.hpp"
#include "../utility/functional.hpp"
#include "../utility/matrix.hpp"
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

namespace gapp::algorithm::reflines
{
    using SimplexMapping = void(*)(FitnessVector&);

    /*
    * The unit-hypercube -> unit-simplex transformations used for the quasirandom points are based on:
    * 
    *   Pillards, Tim, and Ronald Cools. "Transforming low-discrepancy sequences from a cube to a simplex."
    *   Journal of computational and applied mathematics 174, no. 1 (2005): 29-42.
    */

    /* Transform a point from the n-dimensional unit hypercube to the n-dimensional unit simplex. */
    static inline void simplexMappingLog(FitnessVector& point)
    {
        GAPP_ASSERT(std::all_of(point.begin(), point.end(), detail::between(0.0, 1.0)));

        std::transform(point.begin(), point.end(), point.begin(), [](double p) { return -std::log(std::max(p, math::small<double>)); });
        const double sum = std::reduce(point.begin(), point.end(), 0.0);
        std::transform(point.begin(), point.end(), point.begin(), detail::divide_by(sum));
    }

    /* Transform a point from the n-dimensional unit hypercube to the (n+1)-dimensional unit simplex. */
    static inline void simplexMappingSort(FitnessVector& point)
    {
        GAPP_ASSERT(std::all_of(point.begin(), point.end(), detail::between(0.0, 1.0)));

        point.push_back(1.0);
        std::sort(point.begin(), point.end());
        std::adjacent_difference(point.begin(), point.end(), point.begin());
    }

    /* Transform a point from the n-dimensional unit hypercube to the (n+1)-dimensional unit simplex. */
    static inline void simplexMappingRoot(FitnessVector& point)
    {
        GAPP_ASSERT(std::all_of(point.begin(), point.end(), detail::between(0.0, 1.0)));

        if (point.empty())
        {
            point = { 1.0 };
            return;
        }

        point.back() = std::pow(point.back(), 1.0 / point.size());

        for (auto it = std::next(point.rbegin()); it != point.rend(); ++it)
        {
            *it = *std::prev(it) * std::pow(*it, 1.0 / std::distance(it, point.rend()));
        }
        point.push_back(1.0);

        std::adjacent_difference(point.begin(), point.end(), point.begin());
    }

    /* Transform a point from the n-dimensional unit hypercube to the (n+1)-dimensional unit simplex. */
    static inline void simplexMappingMirror(FitnessVector& point)
    {
        GAPP_ASSERT(std::all_of(point.begin(), point.end(), detail::between(0.0, 1.0)));

        point.push_back(1.0);
        if (point.size() == 1) return;

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
    }


    template<SimplexMapping>
    struct SimplexMappingTraits
    {
        static constexpr size_t input_dim(size_t output_dim) noexcept { return output_dim - 1; }
    };

    template<>
    struct SimplexMappingTraits<simplexMappingLog>
    {
        static constexpr size_t input_dim(size_t output_dim) noexcept { return output_dim; }
    };


    template<SimplexMapping F>
    static FitnessMatrix quasirandomSimplexPoints(size_t dim, size_t num_points)
    {
        FitnessMatrix points(num_points, dim);
        size_t input_dim = SimplexMappingTraits<F>::input_dim(dim);

        if (dim == 0) return points;

        rng::QuasiRandom qrng{ input_dim };

        for (size_t i = 0; i < num_points; i++)
        {
            auto point = qrng();
            std::invoke(F, point);
            points[i] = std::move(point);
        }

        return points;
    }


    FitnessMatrix quasirandomSimplexPointsSort(size_t dim, size_t num_points)
    {
        return quasirandomSimplexPoints<simplexMappingSort>(dim, num_points);
    }

    FitnessMatrix quasirandomSimplexPointsRoot(size_t dim, size_t num_points)
    {
        return quasirandomSimplexPoints<simplexMappingRoot>(dim, num_points);
    }

    FitnessMatrix quasirandomSimplexPointsMirror(size_t dim, size_t num_points)
    {
        return quasirandomSimplexPoints<simplexMappingMirror>(dim, num_points);
    }

    FitnessMatrix quasirandomSimplexPointsLog(size_t dim, size_t num_points)
    {
        return quasirandomSimplexPoints<simplexMappingLog>(dim, num_points);
    }


    FitnessMatrix pickSparseSubset(size_t dim, size_t num_points, RefLineGenerator generator, Positive<size_t> k)
    {
        if (dim * num_points == 0) return {};

        FitnessMatrix candidate_points = generator(dim, k * num_points);

        FitnessMatrix points;
        points.reserve(num_points, dim);
        points.append_row(candidate_points.back());
        candidate_points.pop_back();

        std::vector<double> min_distances(candidate_points.size());
        std::transform(candidate_points.begin(), candidate_points.end(), min_distances.begin(), std::bind_front(math::euclideanDistanceSq, points.back()));
        
        while (points.size() < num_points)
        {
            const size_t idx = detail::argmax(min_distances.begin(), min_distances.end());
            points.append_row(candidate_points[idx]);

            /* Remove the added candidate and the corresponding min_distance. */
            using std::swap;
            swap(candidate_points[idx], candidate_points.back());
            candidate_points.pop_back();
            swap(min_distances[idx], min_distances.back());
            min_distances.pop_back();

            /* Calc the distance of each candidate to the closest ref point. */
            for (size_t i = 0; i < candidate_points.size(); i++)
            {
                double dist = math::euclideanDistanceSq(candidate_points[i], points.back());
                min_distances[i] = std::min(min_distances[i], dist);
            }
        }

        return points;
    }

} // namespace gapp::algorithm::reflines