/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include "utility/bounded_value.hpp"

using namespace gapp;

TEST_CASE("probability", "[probability]")
{
    SECTION("constructor")
    {
        [[maybe_unused]] auto p1 = Probability(0.123);
        [[maybe_unused]] auto p2 = Probability(0.0);
        [[maybe_unused]] auto p3 = 1.0_p;
        SUCCEED();
    }

    SECTION("constexpr")
    {
        constexpr Probability p = 0.2;

        STATIC_REQUIRE(p == 0.2);
    }

    SECTION("assignment")
    {
        Probability p = 0.3;
        REQUIRE(p == 0.3);

        p = 0.1;
        REQUIRE(p == 0.1);
    }

    SECTION("comparisons")
    {
        Probability p1 = 0.2;
        Probability p2 = 0.7;

        REQUIRE(p1 != p2);
        REQUIRE(p1 < p2);
        REQUIRE(p2 > p1);

        REQUIRE(p1 < 0.5);
        REQUIRE(p2 >= 0.5);
    }
}