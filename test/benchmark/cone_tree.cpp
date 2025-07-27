/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include "utility/cone_tree.hpp"
#include "utility/rng.hpp"
#include "utility/algorithm.hpp"

using namespace gapp;
using namespace gapp::detail;
using namespace Catch;


static ConeTree::Point randomPoint(size_t dim)
{
    ConeTree::Point point(dim);
    for (double& elem : point) elem = rng::randomReal();
    return point;
}

static std::vector<ConeTree::Point> randomPoints(size_t n, size_t dim)
{
    std::vector<ConeTree::Point> points(n);
    for (auto& point : points) point = randomPoint(dim);
    return points;
}

static auto linearFind(const std::vector<ConeTree::Point>& lines, const ConeTree::Point& point)
{
    auto idistance = [&](const auto& line) { return std::inner_product(point.begin(), point.end(), line.begin(), 0.0); };

    return detail::max_element(lines.begin(), lines.end(), idistance);
}


TEST_CASE("cone_tree_ctor", "[cone_tree]")
{
    constexpr size_t ndim = 3;

    BENCHMARK_ADVANCED("small")(Benchmark::Chronometer meter)
    {
        auto points = randomPoints(100, ndim);
        meter.measure([&] { return ConeTree{ points }; });
    };

    BENCHMARK_ADVANCED("medium")(Benchmark::Chronometer meter)
    {
        auto points = randomPoints(1000, ndim);
        meter.measure([&] { return ConeTree{ points }; });
    };

    BENCHMARK_ADVANCED("large")(Benchmark::Chronometer meter)
    {
        auto points = randomPoints(10000, ndim);
        meter.measure([&] { return ConeTree{ points }; });
    };
}

TEST_CASE("cone_tree_lookup_size", "[cone_tree]")
{
    constexpr size_t ndim = 3;
    const size_t size = GENERATE(100, 1000, 10000);

    WARN("Number of points: " << size);

    auto points = randomPoints(size, ndim);
    auto tree = ConeTree{ points };

    BENCHMARK("cone_tree") { return tree.findBestMatch(randomPoint(ndim)); };
    BENCHMARK("linsearch") { return linearFind(points, randomPoint(ndim)); };
}

TEST_CASE("cone_tree_lookup_dim", "[cone_tree]")
{
    const size_t ndim = GENERATE(3, 15, 100);
    constexpr size_t size = 10000;

    WARN("Number of dimensions: " << ndim);

    auto points = randomPoints(size, ndim);
    auto tree = ConeTree{ points };

    BENCHMARK("cone_tree") { return tree.findBestMatch(randomPoint(ndim)); };
    BENCHMARK("linsearch") { return linearFind(points, randomPoint(ndim)); };
}
