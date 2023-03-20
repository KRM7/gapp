/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/catch_approx.hpp>
#include "algorithm/reference_lines.hpp"
#include "utility/functional.hpp"
#include <algorithm>
#include <numeric>
#include <vector>
#include <cstddef>

using Catch::Approx;
using namespace genetic_algorithm::detail;
using namespace genetic_algorithm::math;
using namespace genetic_algorithm::algorithm::reflines;


TEMPLATE_TEST_CASE_SIG("reference_lines", "[pareto_front]", ((auto F), F), quasirandomSimplexPointsMirror, quasirandomSimplexPointsSort, quasirandomSimplexPointsRoot, quasirandomSimplexPointsLog)
{
    const size_t num_points = GENERATE(0, 1, 10);
    const size_t dim = GENERATE(0, 1, 2, 3, 100);

    INFO("Dimensions: " << dim);

    std::vector<Point> points = F(dim, num_points);

    REQUIRE(points.size() == num_points);
    REQUIRE(std::all_of(points.begin(), points.end(), is_size(dim)));

    if (dim == 0) return;

    for (const Point& point : points)
    {
        CAPTURE(point);
        REQUIRE(std::accumulate(point.begin(), point.end(), 0.0) == Approx(1.0).margin(1E-4));
        REQUIRE(std::all_of(point.begin(), point.end(), greater_eq_than(0.0)));
    }
}

TEMPLATE_TEST_CASE_SIG("reference_lines_subset", "[pareto_front]", ((auto F), F), quasirandomSimplexPointsMirror, quasirandomSimplexPointsSort, quasirandomSimplexPointsRoot, quasirandomSimplexPointsLog)
{
    const size_t num_points = GENERATE(0, 1, 10);
    const size_t dim = GENERATE(0, 1, 2, 3, 100);

    INFO("Dimensions: " << dim);

    std::vector<Point> points = pickSparseSubset(dim, num_points, F);

    REQUIRE(points.size() == num_points);
    REQUIRE(std::all_of(points.begin(), points.end(), is_size(dim)));

    if (dim == 0) return;

    for (const Point& point : points)
    {
        CAPTURE(point);
        REQUIRE(std::accumulate(point.begin(), point.end(), 0.0) == Approx(1.0).margin(1E-4));
        REQUIRE(std::all_of(point.begin(), point.end(), greater_eq_than(0.0)));
    }
}