/* Copyright (c) 2025 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include "utility/rng.hpp"
#include "utility/distribution.hpp"
#include "utility/math.hpp"
#include <limits>
#include <cstdint>

using namespace gapp;


TEST_CASE("uniform_bool_distribution", "[distribution]")
{
    rng::prng.seed(std::random_device{}());

    rng::uniform_bool_distribution dist;

    REQUIRE(dist.min() == 0);
    REQUIRE(dist.max() == 1);

    REQUIRE(dist == rng::uniform_bool_distribution{});

    std::ignore = dist(rng::prng);

    REQUIRE(dist != rng::uniform_bool_distribution{});
}

TEST_CASE("uniform_int_distribution", "[distribution]")
{
    rng::prng.seed(std::random_device{}());

    rng::uniform_int_distribution<int64_t> dist1(0, 0);

    REQUIRE(dist1.min() == 0);
    REQUIRE(dist1.max() == 0);

    for (size_t i = 0; i < 1000; i++)
    {
        REQUIRE(dist1(rng::prng) == 0);
    }


    rng::uniform_int_distribution<uint64_t> dist2(0, 0);

    REQUIRE(dist2.min() == 0);
    REQUIRE(dist2.max() == 0);

    for (size_t i = 0; i < 1000; i++)
    {
        REQUIRE(dist2(rng::prng) == 0);
    }


    rng::uniform_int_distribution<int64_t> dist3(-100, 100);

    REQUIRE(dist3.min() == -100);
    REQUIRE(dist3.max() == 100);

    for (size_t i = 0; i < 1000; i++)
    {
        const int64_t n = dist3(rng::prng);

        REQUIRE(n >= -100);
        REQUIRE(n <= 100);
    }


    rng::uniform_int_distribution<uint64_t> dist4(0, 100);

    REQUIRE(dist4.min() == 0);
    REQUIRE(dist4.max() == 100);

    for (size_t i = 0; i < 1000; i++)
    {
        const uint64_t n = dist4(rng::prng);

        REQUIRE(n <= 100);
    }


    rng::uniform_int_distribution<int64_t> dist5(std::numeric_limits<int64_t>::min(), 0);

    REQUIRE(dist5.min() == std::numeric_limits<int64_t>::min());
    REQUIRE(dist5.max() == 0);

    for (size_t i = 0; i < 1000; i++)
    {
        const int64_t n = dist5(rng::prng);

        REQUIRE(n <= 0);
    }


    rng::uniform_int_distribution<int64_t> dist6(std::numeric_limits<int64_t>::min(), std::numeric_limits<int64_t>::max());

    REQUIRE(dist6.min() == std::numeric_limits<int64_t>::min());
    REQUIRE(dist6.max() == std::numeric_limits<int64_t>::max());

    REQUIRE_NOTHROW(dist6(rng::prng));


    rng::uniform_int_distribution<uint64_t> dist7(std::numeric_limits<uint64_t>::min(), std::numeric_limits<uint64_t>::max());

    REQUIRE(dist7.min() == std::numeric_limits<uint64_t>::min());
    REQUIRE(dist7.max() == std::numeric_limits<uint64_t>::max());

    REQUIRE_NOTHROW(dist7(rng::prng));


    rng::uniform_int_distribution<int8_t> dist8(1, 6);

    REQUIRE(dist8.min() == 1);
    REQUIRE(dist8.max() == 6);

    for (size_t i = 0; i < 1000; i++)
    {
        const int64_t n = dist8(rng::prng);

        REQUIRE(n >= 1);
        REQUIRE(n <= 6);
    }
}

TEST_CASE("generate_canonical", "[distribution]")
{
    rng::prng.seed(std::random_device{}());

    for (size_t i = 0; i < 1000; i++)
    {
        const float f = rng::generate_canonical<float>(rng::prng);

        REQUIRE(f >= 0.0);
        REQUIRE(f < 1.0);
    }

    for (size_t i = 0; i < 1000; i++)
    {
        const double f = rng::generate_canonical<double>(rng::prng);

        REQUIRE(f >= 0.0);
        REQUIRE(f < 1.0);
    }
}

TEST_CASE("uniform_real_distribution", "[distribution]")
{
    rng::prng.seed(std::random_device{}());

    rng::uniform_real_distribution<double> dist1;

    REQUIRE(dist1.min() == 0.0);
    REQUIRE(dist1.max() == 1.0);

    for (size_t i = 0; i < 1000; i++)
    {
        const double f = dist1(rng::prng);

        REQUIRE(f >= 0.0);
        REQUIRE(f < 1.0);
    }

    rng::uniform_real_distribution<double> dist2(-100.0, 100.0);

    REQUIRE(dist2.min() == -100.0);
    REQUIRE(dist2.max() == 100.0);

    for (size_t i = 0; i < 1000; i++)
    {
        const double f = dist2(rng::prng);

        REQUIRE(f >= -100.0);
        REQUIRE(f < 100.0);
    }

    REQUIRE(dist1 != dist2);
}

TEST_CASE("exponential_distribution", "[distribution]")
{
    rng::prng.seed(std::random_device{}());

    rng::exponential_distribution<double> dist;

    REQUIRE(dist.min() == 0.0);
    REQUIRE(dist.max() == math::inf<double>);

    for (size_t i = 0; i < 1000; i++)
    {
        REQUIRE(dist(rng::prng) >= 0.0);
    }
}

TEST_CASE("normal_distribution", "[distribution]")
{
    rng::prng.seed(std::random_device{}());

    rng::normal_distribution<double> dist;

    REQUIRE(dist.mean() == 0.0);
    REQUIRE(dist.stddev() == 1.0);

    REQUIRE(dist.min() == -math::inf<double>);
    REQUIRE(dist.max() == math::inf<double>);

    REQUIRE_NOTHROW(dist(rng::prng));
}

TEST_CASE("small_poisson_distribution", "[distribution]")
{
    rng::prng.seed(std::random_device{}());

    rng::small_poisson_distribution<int64_t> dist(2.0);

    REQUIRE(dist.min() == 0);
    REQUIRE(dist.max() == std::numeric_limits<int64_t>::max());

    for (size_t i = 0; i < 1000; i++)
    {
        REQUIRE(dist(rng::prng) >= 0);
    }
}

TEST_CASE("symmetric_binomial_distribution", "[distribution]")
{
    rng::prng.seed(std::random_device{}());

    rng::symmetric_binomial_distribution<int64_t> dist1(100);

    REQUIRE(dist1.min() == 0);
    REQUIRE(dist1.max() == 100);

    for (size_t i = 0; i < 1000; i++)
    {
        REQUIRE(dist1(rng::prng) >= 0);
    }


    rng::symmetric_binomial_distribution<int64_t> dist2(1000);

    REQUIRE(dist2.min() == 0);
    REQUIRE(dist2.max() == 1000);

    for (size_t i = 0; i < 1000; i++)
    {
        REQUIRE(dist2(rng::prng) >= 0);
    }

    REQUIRE(dist1 != dist2);
}

TEST_CASE("binomial_distribution", "[distribution]")
{
    rng::prng.seed(std::random_device{}());

    rng::binomial_distribution<int64_t> dist1(100, 0.02);

    REQUIRE(dist1.min() == 0);
    REQUIRE(dist1.max() == 100);

    for (size_t i = 0; i < 1000; i++)
    {
        const int64_t n = dist1(rng::prng);

        REQUIRE(n >= 0);
        REQUIRE(n <= 100);
    }


    rng::binomial_distribution<int64_t> dist2(100, 0.35);

    REQUIRE(dist1.min() == 0);
    REQUIRE(dist1.max() == 100);

    for (size_t i = 0; i < 1000; i++)
    {
        const int64_t n = dist2(rng::prng);

        REQUIRE(n >= 0);
        REQUIRE(n <= 100);
    }

    REQUIRE(dist1 != dist2);
}
