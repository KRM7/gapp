#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "utility/math.hpp"
#include <vector>
#include <limits>
#include <numbers>

using namespace genetic_algorithm::math;
using namespace Catch;

TEST_CASE("fp_compare")
{
    constexpr double INF   = std::numeric_limits<double>::infinity();
    constexpr double BIG   = std::numeric_limits<double>::max();
    constexpr double SMALL = std::numeric_limits<double>::denorm_min();
    constexpr double NaN   = std::numeric_limits<double>::quiet_NaN();

    SECTION("equality")
    {
        REQUIRE(floatIsEqual(0.0, 0.0));
        REQUIRE(floatIsEqual(1702.17, 1702.17));
        REQUIRE(floatIsEqual(SMALL, SMALL));
        REQUIRE(floatIsEqual(BIG, BIG));
        REQUIRE(floatIsEqual(INF, INF));
        REQUIRE(!floatIsEqual(INF, -INF));
        REQUIRE(!floatIsEqual(NaN, NaN));

        REQUIRE(!floatIsEqual(0.0, INF));
        REQUIRE(!floatIsEqual(0.0, BIG));
        REQUIRE(!floatIsEqual(0.0, NaN));

        REQUIRE(!floatIsEqual(SMALL, BIG));
        REQUIRE(!floatIsEqual(SMALL, INF));
        REQUIRE(!floatIsEqual(SMALL, NaN));

        REQUIRE(!floatIsEqual(BIG, INF));
        REQUIRE(!floatIsEqual(BIG, NaN));

        REQUIRE(!floatIsEqual(INF, NaN));
    }

    SECTION("less")
    {
        REQUIRE(!floatIsLess(0.0, 0.0));
        REQUIRE(!floatIsLess(SMALL, SMALL));
        REQUIRE(!floatIsLess(BIG, BIG));
        REQUIRE(!floatIsLess(INF, INF));
        REQUIRE(!floatIsLess(NaN, NaN));

        REQUIRE(floatIsLess(0.0, INF));
        REQUIRE(!floatIsLess(INF, 0.0));

        REQUIRE(floatIsLess(0.0, BIG));
        REQUIRE(!floatIsLess(BIG, 0.0));

        REQUIRE(!floatIsLess(0.0, NaN));
        REQUIRE(!floatIsLess(NaN, 0.0));

        REQUIRE(floatIsLess(SMALL, BIG));
        REQUIRE(!floatIsLess(BIG, SMALL));

        REQUIRE(floatIsLess(SMALL, INF));
        REQUIRE(!floatIsLess(INF, SMALL));

        REQUIRE(!floatIsLess(SMALL, NaN));
        REQUIRE(!floatIsLess(NaN, SMALL));

        REQUIRE(floatIsLess(BIG, INF));
        REQUIRE(!floatIsLess(INF, BIG));

        REQUIRE(!floatIsLess(BIG, NaN));
        REQUIRE(!floatIsLess(NaN, BIG));

        REQUIRE(!floatIsLess(INF, NaN));
        REQUIRE(!floatIsLess(NaN, INF));
    }

    SECTION("approx equal")
    {
        LocalTolerances _(10, 1E-12);

        REQUIRE(floatIsEqual(0.0, 1E-13));
        REQUIRE(!floatIsEqual(0.0, 1E-11));

        REQUIRE(floatIsEqual(1.28E+32, 1.28E+32 + 1E+15));
        REQUIRE(!floatIsEqual(1.28E+32, 1.29E+32));
    }
}

TEST_CASE("pareto_compare_less", "[math]")
{
    const std::vector vec = { 3.0, 2.0, 1.0 };

    SECTION("self")
    {
        REQUIRE(!paretoCompareLess(vec, vec));
    }

    SECTION("bigger")
    {
        const std::vector other = { 3.0, 3.0, 1.0 };

        REQUIRE(paretoCompareLess(vec, other));
        REQUIRE(!paretoCompareLess(other, vec));
    }

    SECTION("equal")
    {
        const std::vector other = { 4.0, 5.1, 0.0 };

        REQUIRE(!paretoCompareLess(vec, other));
        REQUIRE(!paretoCompareLess(other, vec));
    }

    SECTION("1D")
    {
        REQUIRE(!paretoCompareLess({ 1.0 }, { 1.0 }));

        REQUIRE(paretoCompareLess({ 1.0 }, { 2.3 }));
        REQUIRE(!paretoCompareLess({ 2.3 }, { 1.0 }));
    }

    SECTION("empty")
    {
        REQUIRE(!paretoCompareLess({}, {}));
    }
}

TEST_CASE("pareto_compare_three_way", "[math]")
{
    const std::vector vec = { 3.0, 2.0, 1.0 };

    SECTION("self")
    {
        REQUIRE(paretoCompare(vec, vec) == 0);
    }

    SECTION("bigger")
    {
        const std::vector other = { 3.0, 3.0, 1.0 };

        REQUIRE(paretoCompare(vec, other) < 0);
        REQUIRE(paretoCompare(other, vec) > 0);
    }
    
    SECTION("equal")
    {
        const std::vector other = { 4.0, 5.1, 0.0 };

        REQUIRE(paretoCompare(vec, other) == 0);
        REQUIRE(paretoCompare(other, vec) == 0);
    }

    SECTION("1D")
    {
        REQUIRE(paretoCompare({ 1.0 }, { 1.0 }) == 0);

        REQUIRE(paretoCompare({ 1.0 }, { 2.3 }) < 0);
        REQUIRE(paretoCompare({ 2.3 }, { 1.0 }) > 0);
    }

    SECTION("empty")
    {
        REQUIRE(paretoCompare({}, {}) == 0);
    }
}

TEST_CASE("euclidean_norm", "[math]")
{
    REQUIRE(euclideanNorm({}) == 0.0);
    REQUIRE(euclideanNorm({ 1.0, 4.5, 3.2, 0.3 }) == Approx(5.62).margin(0.01));
}

TEST_CASE("normalize", "[math]")
{
    const std::vector vec            = { 1.0,  4.5,  3.2,  0.3 };
    const std::vector normalized_vec = { 0.18, 0.80, 0.57, 0.05 }; // not exact

    REQUIRE_THAT(normalizeVector(vec), Matchers::Approx(normalized_vec).margin(0.01));  
    REQUIRE_THAT(normalizeVector(normalized_vec), Matchers::Approx(normalized_vec).margin(0.01));

    REQUIRE(normalizeVector({}) == std::vector<double>{});
}

TEST_CASE("euclidean_distance", "[math]")
{
    REQUIRE(euclideanDistanceSq({}, {}) == Approx(0.0));

    REQUIRE(euclideanDistanceSq({ 3.0 }, { 1.0 }) == Approx(4.0).margin(0.01));
    REQUIRE(euclideanDistanceSq({ 1.0, 0.0 }, { 2.0, 1.0 }) == Approx(2.0).margin(0.01));
    REQUIRE(euclideanDistanceSq({ 1.0, 2.0, 0.0 }, { 3.0, 0.0, 1.0 }) == Approx(9.0).margin(0.01));
}

TEST_CASE("perpendicular_distance", "[math]")
{
    REQUIRE(perpendicularDistanceSq({}, {}) == Approx(0.0).margin(1E-8));

    REQUIRE(perpendicularDistanceSq({ 3.1 }, { 0.95 }) == Approx(0.0).margin(1E-8));
    REQUIRE(perpendicularDistanceSq({ 0.4, 0.9 }, { 2.5, 1.0 }) == Approx(3.53).margin(0.01));
}

TEST_CASE("integral_sin_pow", "[math]")
{
    SECTION("n = 0")
    {
        double x = GENERATE(0.0, 0.13, 0.5, 1.85, 3.72, -1.2);

        REQUIRE(integralSinPow(0, x) == Approx(x).margin(1E-6));
    }

    SECTION("n = 1")
    {
        REQUIRE(integralSinPow(1, 0.0) == Approx(-1.0).margin(0.01));
        REQUIRE(integralSinPow(1, 0.75) == Approx(-0.73).margin(0.01));
        REQUIRE(integralSinPow(1, 3.14) == Approx(1.0).margin(0.01));
        REQUIRE(integralSinPow(1, 4.6) == Approx(0.11).margin(0.01));
    }

    SECTION("n > 1")
    {
        REQUIRE(integralSinPow(2, 0.0) == Approx(0.0).margin(0.01));
        REQUIRE(integralSinPow(2, 0.4) == Approx(0.02).margin(0.01));

        REQUIRE(integralSinPow(3, 0.0) == Approx(-0.66).margin(0.01));
        REQUIRE(integralSinPow(3, 1.2) == Approx(-0.35).margin(0.01));

        REQUIRE(integralSinPow(5, 0.0) == Approx(-0.53).margin(0.01));
        REQUIRE(integralSinPow(5, 7.6) == Approx(-0.24).margin(0.01));

        REQUIRE(integralSinPow(12, 0.0) == Approx(0.0).margin(0.01));
        REQUIRE(integralSinPow(12, 3.1) == Approx(0.71).margin(0.01));
    }
}