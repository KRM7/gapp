#include <catch2/catch_test_macros.hpp>
#include "utility/probability.hpp"

using namespace genetic_algorithm;

TEST_CASE("probability", "[probability]")
{
    REQUIRE_NOTHROW(Probability(0.123));
    REQUIRE_NOTHROW(Probability(0.0));
    REQUIRE_NOTHROW(Probability(1.0));

    REQUIRE_THROWS(Probability(-0.3));
    REQUIRE_THROWS(Probability(1.3));
    REQUIRE_THROWS(Probability(100.0/0.0));

    constexpr Probability p1 = 0.2;

    STATIC_REQUIRE(p1 == 0.2);
    STATIC_REQUIRE(*p1 == 0.2);


    Probability p2 = 0.3;

    REQUIRE_NOTHROW(p2 = 0.1);
    REQUIRE_THROWS(p2 = -3.0);

    REQUIRE(p2 == 0.1);
}