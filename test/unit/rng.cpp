/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/catch_approx.hpp>
#include "utility/rng.hpp"
#include "utility/functional.hpp"
#include <algorithm>
#include <limits>
#include <cstdint>

using namespace gapp;
using namespace gapp::rng;


TEST_CASE("random_int", "[rng]")
{
    REQUIRE(randomInt(-1, -1) == -1);

    for (size_t i = 0; i < 10'000; i++)
    {
        int n = randomInt(-1, 1);
        REQUIRE((-1 <= n && n <= 1));
    }
}

TEST_CASE("random_real", "[rng]")
{
    REQUIRE(randomReal(1.0, 1.0) == Catch::Approx(1.0));

    for (size_t i = 0; i < 10'000; i++)
    {
        double n1 = randomReal();
        REQUIRE((0.0 <= n1 && n1 <= 1.0));

        double n2 = randomReal(-1.0, 1.0);
        REQUIRE((-1.0 <= n2 && n2 <= 1.0));
    }
}

TEST_CASE("random_binomial", "[rng]")
{
    REQUIRE(randomBinomial(10u, 0.0) == 0);
    REQUIRE(randomBinomial(10u, 1.0) == 10);

    REQUIRE(randomBinomial(0u, 0.0) == 0);
    REQUIRE(randomBinomial(0u, 1.0) == 0);

    for (size_t i = 0; i < 10'000; i++)
    {
        REQUIRE(randomBinomial(1u, 0.5) <= 1);
        REQUIRE(randomBinomial(5u, 0.8) <= 5);
    }
}

TEST_CASE("random_index", "[rng]")
{
    REQUIRE(randomIndex(std::vector{ true }) == 0);
    REQUIRE(randomIndex(std::vector{ 0.12, 0.32 }) <= 1);
}

TEST_CASE("random_element", "[rng]")
{
    const std::vector single{ 2 };
    const std::vector vector{ 0, 3 };

    REQUIRE(randomElement(single.begin(), single.end()) == single.begin());

    for (size_t i = 0; i < 100; i++)
    {
        REQUIRE(randomElement(vector.begin(), vector.end()) != vector.end());
    }
}

TEST_CASE("random_element_cdf", "[rng]")
{
    const std::vector single{ 2 };

    REQUIRE(randomElement(single, std::vector{ 1.0 }) == 2);

    const std::vector vec1{ 0, 3, 9 };
    const std::vector cdf1{ 0.3, 0.4, 1.0 };

    for (size_t i = 0; i < 100; i++)
    {
        const int elt = randomElement(vec1, cdf1);
        REQUIRE((elt == 0 || elt == 3 || elt == 9));
    }

    const std::vector vec2{ 0, 3, 9 };
    const std::vector cdf2{ 0.0, 0.0, 1.0 };

    for (size_t i = 0; i < 100; i++)
    {
        const int elt = randomElement(vec2, cdf2);
        REQUIRE(elt == 9);
    }
}

TEST_CASE("sample_unique", "[rng]")
{
    auto [lbound, ubound] = GENERATE(std::pair{ -60, -10 }, std::pair{ -20, 30 }, std::pair{ 10, 60 });
    constexpr size_t count = 25;

    for (size_t i = 0; i < 100; i++)
    {
        auto nums = sampleUnique(lbound, ubound, count);

        REQUIRE(nums.size() == count);

        REQUIRE(std::all_of(nums.begin(), nums.end(), detail::between(lbound, ubound - 1)));

        REQUIRE((
            std::sort(nums.begin(), nums.end()),
            std::unique(nums.begin(), nums.end()) == nums.end()
        ));
    }
}

TEMPLATE_TEST_CASE("sample_unique_bounds_8", "[rng]", std::int8_t, std::uint8_t)
{
    using IntegralType = TestType;
    constexpr IntegralType low  = std::numeric_limits<IntegralType>::min();
    constexpr IntegralType high = std::numeric_limits<IntegralType>::max();
    const size_t count = GENERATE(2, 130, 250);

    auto nums = sampleUnique<IntegralType>(low, high, count);

    REQUIRE(nums.size() == count);
    REQUIRE(std::all_of(nums.begin(), nums.end(), detail::between<IntegralType>(low, high - 1)));
    REQUIRE((
        std::sort(nums.begin(), nums.end()),
        std::unique(nums.begin(), nums.end()) == nums.end()
    ));
}

TEMPLATE_TEST_CASE("sample_unique_bounds_64", "[rng]", std::int64_t, std::uint64_t)
{
    using IntegralType = TestType;
    constexpr IntegralType low = std::numeric_limits<IntegralType>::min();
    constexpr IntegralType high = std::numeric_limits<IntegralType>::max();
    constexpr size_t count = 3;

    auto nums = sampleUnique<IntegralType>(low, high, count);

    REQUIRE(nums.size() == count);
    REQUIRE(std::all_of(nums.begin(), nums.end(), detail::between<IntegralType>(low, high - 1)));
    REQUIRE((
        std::sort(nums.begin(), nums.end()),
        std::unique(nums.begin(), nums.end()) == nums.end()
    ));
}

TEST_CASE("sample_cdf", "[rng]")
{
    const std::vector cdf1 = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };

    for (size_t i = 0; i < 100; i++)
    {
        const size_t idx = rng::sampleCdf(cdf1);
        REQUIRE(idx < cdf1.size());
    }

    const std::vector cdf2 = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0 };

    for (size_t i = 0; i < 100; i++)
    {
        const size_t idx = rng::sampleCdf(cdf2);
        REQUIRE(idx > 4);
        REQUIRE(idx < cdf1.size());
    }
}
