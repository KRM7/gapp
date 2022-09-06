/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_REFERENCE_LINES_HPP
#define GA_ALGORITHM_REFERENCE_LINES_HPP

#include <vector>
#include <utility>
#include <cstddef>
#include "../utility/math.hpp"

namespace genetic_algorithm::algorithm::dtl
{
    using math::Point;

    /* A reference direction for the NSGA-III algorithm. */
    struct RefLine
    {
        const Point direction;
        size_t niche_count;

        RefLine(Point p) noexcept :
            direction(std::move(p)),
            niche_count(0)
        {}
    };

    /* Generate n reference points in dim dimensions. */
    std::vector<RefLine> generateReferencePoints(size_t dim, size_t n);

} // namespace genetic_algorithm::algorithm::dtl

#endif // !GA_ALGORITHM_REFERENCE_LINES_HPP