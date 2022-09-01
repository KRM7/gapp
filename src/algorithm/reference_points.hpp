/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_REFERENCE_POINTS_HPP
#define GA_ALGORITHM_REFERENCE_POINTS_HPP

#include <vector>
#include <utility>
#include <cstddef>

namespace genetic_algorithm::algorithm::dtl
{
    using Point = std::vector<double>;

    struct ReferencePoint
    {
        const Point point;
        size_t niche_count;

        ReferencePoint(Point p) :
            point(std::move(p)),
            niche_count(0) {}

        friend bool operator==(const ReferencePoint& lhs, const ReferencePoint& rhs) = default;
    };

    using ReferencePoints = std::vector<ReferencePoint>;


    /* Generate n reference points in dim dimensions. */
    ReferencePoints generateReferencePoints(size_t dim, size_t n);

    /* Find the index and distance of the closest reference line to the point p. */
    std::pair<size_t, double> findClosestRef(const ReferencePoints& refs, const Point& p);

} // namespace genetic_algorithm::algorithm::dtl

#endif // !GA_ALGORITHM_REFERENCE_POINTS_HPP