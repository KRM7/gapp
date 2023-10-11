/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include "algorithm/reference_lines.hpp"
#include "core/population.hpp"
#include "metrics/pop_stats.hpp"
#include "utility/math.hpp"

using namespace gapp;
using namespace gapp::algorithm::reflines;


template<auto F>
static double hypervolumeOf(size_t dim, size_t count)
{
    auto points = F(dim, count);

    FitnessMatrix fmat;
    fmat.reserve(count, dim);
    for (const auto& point : points) { fmat.append_row(point); }

    return detail::hypervolume(fmat, math::Point(dim, 0.0));
}


TEST_CASE("ref_points_count", "[benchmark]")
{
    constexpr size_t dim = 3;
    const size_t count = GENERATE(40, 200, 1500);

    WARN("Number of generated points: " << count);

    BENCHMARK("quasirandom_mirror")
    {
        return quasirandomSimplexPointsMirror(dim, count);
    };

    BENCHMARK("quasirandom_sort")
    {
        return quasirandomSimplexPointsSort(dim, count);
    };

    BENCHMARK("quasirandom_root")
    {
        return quasirandomSimplexPointsRoot(dim, count);
    };

    BENCHMARK("quasirandom_log")
    {
        return quasirandomSimplexPointsLog(dim, count);
    };

    WARN("Mirror hypervolume: " << hypervolumeOf<quasirandomSimplexPointsMirror>(dim, count));
    WARN("Sort hypervolume: " << hypervolumeOf<quasirandomSimplexPointsSort>(dim, count));
    WARN("Root hypervolume: " << hypervolumeOf<quasirandomSimplexPointsRoot>(dim, count));
    WARN("Log hypervolume: " << hypervolumeOf<quasirandomSimplexPointsLog>(dim, count));
}

TEST_CASE("ref_points_dimensions", "[benchmark]")
{
    constexpr size_t count = 100;
    const size_t dimensions = GENERATE(3, 15, 100);

    WARN("Number of dimensions: " << dimensions);

    BENCHMARK("quasirandom_mirror")
    {
        return quasirandomSimplexPointsMirror(dimensions, count);
    };

    BENCHMARK("quasirandom_sort")
    {
        return quasirandomSimplexPointsSort(dimensions, count);
    };

    BENCHMARK("quasirandom_root")
    {
        return quasirandomSimplexPointsRoot(dimensions, count);
    };

    BENCHMARK("quasirandom_log")
    {
        return quasirandomSimplexPointsLog(dimensions, count);
    };

    if (dimensions < 20) // takes too long otherwise
    {
        WARN("Mirror hypervolume: " << hypervolumeOf<quasirandomSimplexPointsMirror>(dimensions, 30));
        WARN("Sort hypervolume: " << hypervolumeOf<quasirandomSimplexPointsSort>(dimensions, 30));
        WARN("Root hypervolume: " << hypervolumeOf<quasirandomSimplexPointsRoot>(dimensions, 30));
        WARN("Log hypervolume: " << hypervolumeOf<quasirandomSimplexPointsLog>(dimensions, 30));
    }
}
