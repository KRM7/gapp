/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include "utility/rng.hpp"

using namespace gapp::rng;

TEST_CASE("sample_unique", "[benchmark]")
{
    BENCHMARK("small range - select few") { return sampleUnique(0, 100, 2); };
    BENCHMARK("small range - select half") { return sampleUnique(0, 100, 50); };
    BENCHMARK("small range - select many") { return sampleUnique(0, 100, 95); };

    BENCHMARK("large range - select few") { return sampleUnique(0, 50000, 3); };
    BENCHMARK("large range - select half") { return sampleUnique(0, 50000, 25000); };
    BENCHMARK("large range - select many") { return sampleUnique(0, 50000, 49750); };
}
