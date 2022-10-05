#include <catch2/catch_test_macros.hpp>
#include "utility/probability.hpp"

using namespace genetic_algorithm;

TEST_CASE("probability", "[probability]")
{
    SECTION("good values")
    {
        REQUIRE_NOTHROW(Probability(0.123));
        REQUIRE_NOTHROW(Probability(0.0));
        REQUIRE_NOTHROW(Probability(1.0));
    }

    SECTION("bad values")
    {
        REQUIRE_THROWS(Probability(-0.3));
        REQUIRE_THROWS(Probability(1.3));
        REQUIRE_THROWS(Probability(100.0/0.0));
    }

    SECTION("compile time constant")
    {
        constexpr Probability p = 0.2;

        STATIC_REQUIRE(p == 0.2);
        STATIC_REQUIRE(*p == 0.2);
    }

    SECTION("assignment")
    {
        Probability p = 0.3;

        REQUIRE_NOTHROW(p = 0.1);
        REQUIRE_THROWS(p = -3.0);

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