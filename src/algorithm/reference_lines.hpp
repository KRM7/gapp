/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_REFERENCE_LINES_HPP
#define GA_ALGORITHM_REFERENCE_LINES_HPP

#include "../utility/math.hpp"
#include "../utility/bounded_value.hpp"
#include <vector>
#include <cstddef>

/** Methods for generating reference lines for the NSGA3 algorithm. */
namespace genetic_algorithm::algorithm::reflines
{
    using math::Point;

    using RefLineGenerator = std::vector<Point>(*)(size_t, size_t);

    /**
    * Generates a set of points on the unit simplex by mapping a set of quasirandom
    * points generated in a unit hypercube onto the unit simplex.
    *
    * @param dim The dimension of the generated points.
    * @param num_points The number of points to generate.
    * @returns The generated simplex points.
    */
    std::vector<Point> quasirandomSimplexPointsMirror(size_t dim, size_t num_points);

    /**
    * Generates a set of points on the unit simplex by mapping a set of quasirandom
    * points generated in a unit hypercube onto the unit simplex.
    *
    * @param dim The dimension of the generated points.
    * @param num_points The number of points to generate.
    * @returns The generated simplex points.
    */
    std::vector<Point> quasirandomSimplexPointsSort(size_t dim, size_t num_points);

    /**
    * Generates a set of points on the unit simplex by mapping a set of quasirandom
    * points generated in a unit hypercube onto the unit simplex.
    * 
    * @param dim The dimension of the generated points.
    * @param num_points The number of points to generate.
    * @returns The generated simplex points.
    */
    std::vector<Point> quasirandomSimplexPointsRoot(size_t dim, size_t num_points);

    /**
    * Generates a set of points on the unit simplex by mapping a set of quasirandom
    * points generated in a unit hypercube onto the unit simplex.
    *
    * @param dim The dimension of the generated points.
    * @param num_points The number of points to generate.
    * @returns The generated simplex points.
    */
    std::vector<Point> quasirandomSimplexPointsLog(size_t dim, size_t num_points);


    /**
    * Generate a set of reference points by picking a subset of the points generated
    * by another simplex point generator function.
    * 
    * @param dim The dimension of the generated points.
    * @param num_points The number of points to generate.
    * @param generator The simplex point generator function used for creating the initial set of points.
    * @param k The multiple of num_points to use for the size of the initial point set generated using the given generator function.
    * @returns The generated simplex points.
    */
    std::vector<Point> pickSparseSubset(size_t dim, size_t num_points, RefLineGenerator generator, Positive<size_t> k = 10);

} // namespace genetic_algorithm::algorithm::reflines

#endif // !GA_ALGORITHM_REFERENCE_LINES_HPP