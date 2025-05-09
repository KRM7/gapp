/* Copyright (c) 2024 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <thread>
#include "gapp.hpp"

using namespace gapp;

TEST_CASE("reproducibility_single_thread", "[reproducibility]")
{
    RCGA ga{ 10 };
    problems::Sphere f{ 3 };

    execution_threads(1);

    rng::prng.seed(0x9e3779b97f4a7c14);
    const auto solutions1 = ga.solve(f, f.bounds(), 5);

    rng::prng.seed(0x9e3779b97f4a7c14);
    const auto solutions2 = ga.solve(f, f.bounds(), 5);

    REQUIRE(solutions1 == solutions2);

    execution_threads(std::thread::hardware_concurrency());
}

TEST_CASE("reproducibility_multi_thread", "[reproducibility]")
{
    const size_t nthreads = GENERATE(2, 7, 16, 27, 128);

    RCGA ga{ 10 };
    problems::Sphere f{ 3 };

    execution_threads(nthreads);

    rng::prng.seed(0x9e3779b97f4a7c15);
    const auto solutions1 = ga.solve(f, f.bounds(), 5);

    rng::prng.seed(0x9e3779b97f4a7c15);
    const auto solutions2 = ga.solve(f, f.bounds(), 5);

    REQUIRE(solutions1 == solutions2);
}
