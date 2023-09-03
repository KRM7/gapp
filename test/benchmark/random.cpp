/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include "utility/rng.hpp"

using namespace gapp::rng;

TEST_CASE("prng", "[benchmark]")
{
    BENCHMARK("randomBool") { return randomBool(); };
    BENCHMARK("randomInt")  { return randomInt(0, 100); };
    BENCHMARK("randomReal") { return randomReal(0.0, 1.0); };
    BENCHMARK("randomNorm") { return randomNormal(0.0, 10.0); };
}
