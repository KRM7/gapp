#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include "utility/algorithm.hpp"
#include "utility/rng.hpp"
#include <algorithm>

using namespace gapp;

std::vector<double> randomVector(std::size_t size)
{
    std::vector vec(size, 0.0);
    for (double& elem : vec) elem = rng::randomReal();
    std::sort(vec.begin(), vec.end());
    return vec;
}


TEST_CASE("binary_search", "[benchmark]")
{
    const auto vlen = GENERATE(100, 1000, 10000);
    WARN("Number of elements: " << vlen);

    auto v = randomVector(vlen);

    BENCHMARK("std::find_if") { return std::find_if(v.begin(), v.end(), detail::greater_eq_than(rng::randomReal())); };
    BENCHMARK("std::lower_bound") { return std::lower_bound(v.begin(), v.end(), rng::randomReal()); };
    BENCHMARK("detail::lower_bound") { return detail::lower_bound(v.begin(), v.end(), rng::randomReal()); };
}
