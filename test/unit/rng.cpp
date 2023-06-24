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

    for (size_t i = 0; i < 100; i++)
    {
        int n = randomInt(-1, 1);
        REQUIRE((-1 <= n && n <= 1));
    }
}

TEST_CASE("random_real", "[rng]")
{
    REQUIRE(randomReal(1.0, 1.0) == Catch::Approx(1.0));

    for (size_t i = 0; i < 100; i++)
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

    for (size_t i = 0; i < 100; i++)
    {
        REQUIRE(randomBinomial(1u, 0.5) <= 1);
        REQUIRE(randomBinomial(5u, 0.8) <= 5);
    }
}

TEST_CASE("random_idx", "[rng]")
{
    REQUIRE(randomIdx(std::vector{ true }) == 0);
    REQUIRE(randomIdx(std::vector{ 0.12, 0.32 }) <= 1);
}

TEST_CASE("random_element", "[rng]")
{
    std::vector<int> single{ 2 };
    std::vector<int> vector{ 0, 3 };

    REQUIRE(randomElement(single.begin(), single.end()) == single.begin());

    for (size_t i = 0; i < 100; i++)
    {
        REQUIRE(randomElement(vector.begin(), vector.end()) != vector.end());
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
