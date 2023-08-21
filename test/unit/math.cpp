/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "utility/math.hpp"
#include <vector>
#include <limits>
#include <numbers>
#include <utility>

using namespace gapp::math;
using namespace Catch;

TEST_CASE("fp_compare", "[math]")
{
    constexpr double INF   = std::numeric_limits<double>::infinity();
    constexpr double BIG   = std::numeric_limits<double>::max();
    constexpr double SMALL = std::numeric_limits<double>::denorm_min();
    constexpr double NaN   = std::numeric_limits<double>::quiet_NaN();

    SECTION("is_equal")
    {
        auto [abs, rel] = GENERATE(std::pair{ 0.0, 0.0 }, std::pair{ 1E-12, 10 * eps<double> });

        ScopedTolerances _(abs, rel);
        INFO("Relative tolerance eps: " << rel << ", absolute tolerance: " << abs);

        REQUIRE(floatIsEqual(0.0, 0.0));
        REQUIRE(floatIsEqual(0.0, -0.0));
        REQUIRE(floatIsEqual(-0.0, 0.0));
        REQUIRE(floatIsEqual(-0.0, -0.0));

        REQUIRE(floatIsEqual(1702.17, 1702.17));

        REQUIRE(floatIsEqual(SMALL, SMALL));
        REQUIRE(floatIsEqual(BIG, BIG));

        REQUIRE(floatIsEqual(INF, INF));
        REQUIRE(floatIsEqual(-INF, -INF));
        REQUIRE(!floatIsEqual(INF, -INF));
        REQUIRE(!floatIsEqual(-INF, INF));

        REQUIRE(!floatIsEqual(NaN, NaN));

        REQUIRE(!floatIsEqual(0.0, INF));
        REQUIRE(!floatIsEqual(0.0, BIG));
        REQUIRE(!floatIsEqual(0.0, NaN));
        REQUIRE(!floatIsEqual(INF, 0.0));
        REQUIRE(!floatIsEqual(BIG, 0.0));
        REQUIRE(!floatIsEqual(NaN, 0.0));

        REQUIRE(!floatIsEqual(SMALL, BIG));
        REQUIRE(!floatIsEqual(SMALL, INF));
        REQUIRE(!floatIsEqual(SMALL, NaN));
        REQUIRE(!floatIsEqual(BIG, SMALL));
        REQUIRE(!floatIsEqual(INF, BIG));
        REQUIRE(!floatIsEqual(NaN, SMALL));

        REQUIRE(!floatIsEqual(BIG, INF));
        REQUIRE(!floatIsEqual(BIG, NaN));
        REQUIRE(!floatIsEqual(INF, BIG));
        REQUIRE(!floatIsEqual(NaN, BIG));

        REQUIRE(!floatIsEqual(INF, NaN));
        REQUIRE(!floatIsEqual(NaN, INF));
    }

    SECTION("approx is_equal")
    {
        ScopedTolerances _(1E-12, 10 * eps<double>);

        REQUIRE(floatIsEqual(0.0, 1E-13));
        REQUIRE(!floatIsEqual(0.0, 1E-11));

        REQUIRE(floatIsEqual(1.28E+32, 1.28E+32 + 1E+15));
        REQUIRE(!floatIsEqual(1.28E+32, 1.29E+32));
    }

    SECTION("is_less")
    {
        auto [abs, rel] = GENERATE(std::pair{ 0.0, 0.0 }, std::pair{ 1E-12, 10 * eps<double> });

        ScopedTolerances _(abs, rel);
        INFO("Relative tolerance eps: " << rel << ", absolute tolerance: " << abs);

        REQUIRE(!floatIsLess(0.0, 0.0));
        REQUIRE(!floatIsLess(0.0, -0.0));
        REQUIRE(!floatIsLess(-0.0, 0.0));
        REQUIRE(!floatIsLess(-0.0, -0.0));

        REQUIRE(!floatIsLess(4.0, 4.0));
        REQUIRE(floatIsLess(0.0, 4.0));
        REQUIRE(!floatIsLess(4.0, 0.0));

        REQUIRE(!floatIsLess(SMALL, SMALL));
        REQUIRE(!floatIsLess(BIG, BIG));
        REQUIRE(!floatIsLess(INF, INF));
        REQUIRE(!floatIsLess(NaN, NaN));

        REQUIRE(floatIsLess(-INF, INF));
        REQUIRE(!floatIsLess(INF, -INF));
        REQUIRE(!floatIsLess(INF, INF));
        REQUIRE(!floatIsLess(-INF, -INF));

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

    SECTION("approx is_less")
    {
        ScopedTolerances _(1E-12, 10 * eps<double>);

        REQUIRE(!floatIsLess(0.0, 1E-13));
        REQUIRE(floatIsLess(0.0, 1E-11));

        REQUIRE(!floatIsLess(1.28E+32, 1.28E+32 + 1E+15));
        REQUIRE(floatIsLess(1.28E+32, 1.29E+32));
    }

    SECTION("three-way comparison")
    {
        auto [abs, rel] = GENERATE(std::pair{ 0.0, 0.0 }, std::pair{ 1E-12, 10 * eps<double> });

        ScopedTolerances _(abs, rel);
        INFO("Relative tolerance eps: " << rel << ", absolute tolerance: " << abs);

        REQUIRE((floatCompare(0.0, 0.0) == 0));
        REQUIRE((floatCompare(0.0, -0.0) == 0));
        REQUIRE((floatCompare(-0.0, 0.0) == 0));
        REQUIRE((floatCompare(-0.0, -0.0) == 0));

        REQUIRE((floatCompare(4.0, 4.0) == 0));
        REQUIRE((floatCompare(0.0, 4.0) < 0));
        REQUIRE((floatCompare(4.0, 0.0) > 0));

        REQUIRE((floatCompare(SMALL, SMALL) == 0));
        REQUIRE((floatCompare(BIG, BIG) == 0));
        REQUIRE((floatCompare(INF, INF) == 0));

        REQUIRE((floatCompare(-INF, INF) < 0));
        REQUIRE((floatCompare(INF, -INF) > 0));
        REQUIRE((floatCompare(INF, INF) == 0));
        REQUIRE((floatCompare(-INF, -INF) == 0));

        REQUIRE((floatCompare(0.0, INF) < 0));
        REQUIRE((floatCompare(INF, 0.0) > 0));

        REQUIRE((floatCompare(0.0, BIG) < 0));
        REQUIRE((floatCompare(BIG, 0.0) > 0));

        REQUIRE((floatCompare(SMALL, BIG) < 0));
        REQUIRE((floatCompare(BIG, SMALL) > 0));

        REQUIRE((floatCompare(SMALL, INF) < 0));
        REQUIRE((floatCompare(INF, SMALL) > 0));

        REQUIRE((floatCompare(BIG, INF) < 0));
        REQUIRE((floatCompare(INF, BIG) > 0));
    }

    SECTION("is_less_not_greater")
    {
        auto [abs, rel] = GENERATE(std::pair{ 0.0, 0.0 }, std::pair{ 1E-12, 10 * eps<double> });

        ScopedTolerances _(abs, rel);
        INFO("Relative tolerance eps: " << rel << ", absolute tolerance: " << abs);

        REQUIRE(!floatIsLessAssumeNotGreater(0.0, 0.0));
        REQUIRE(!floatIsLessAssumeNotGreater(0.0, -0.0));
        REQUIRE(!floatIsLessAssumeNotGreater(-0.0, 0.0));
        REQUIRE(!floatIsLessAssumeNotGreater(-0.0, -0.0));

        REQUIRE(!floatIsLessAssumeNotGreater(4.0, 4.0));
        REQUIRE(floatIsLessAssumeNotGreater(0.0, 4.0));

        REQUIRE(!floatIsLessAssumeNotGreater(SMALL, SMALL));
        REQUIRE(!floatIsLessAssumeNotGreater(BIG, BIG));
        REQUIRE(!floatIsLessAssumeNotGreater(INF, INF));
        REQUIRE(!floatIsLessAssumeNotGreater(NaN, NaN));

        REQUIRE(floatIsLessAssumeNotGreater(-INF, INF));
        REQUIRE(!floatIsLessAssumeNotGreater(INF, INF));
        REQUIRE(!floatIsLessAssumeNotGreater(-INF, -INF));

        REQUIRE(!floatIsLessAssumeNotGreater(0.0, NaN));
        REQUIRE(!floatIsLessAssumeNotGreater(NaN, 0.0));

        REQUIRE(floatIsLessAssumeNotGreater(0.0, BIG));
        REQUIRE(floatIsLessAssumeNotGreater(0.0, INF));
        REQUIRE(floatIsLessAssumeNotGreater(SMALL, BIG));
        REQUIRE(floatIsLessAssumeNotGreater(SMALL, INF));
        REQUIRE(floatIsLessAssumeNotGreater(BIG, INF));

        REQUIRE(!floatIsLessAssumeNotGreater(SMALL, NaN));
        REQUIRE(!floatIsLessAssumeNotGreater(NaN, SMALL));

        REQUIRE(!floatIsLessAssumeNotGreater(BIG, NaN));
        REQUIRE(!floatIsLessAssumeNotGreater(NaN, BIG));

        REQUIRE(!floatIsLessAssumeNotGreater(INF, NaN));
        REQUIRE(!floatIsLessAssumeNotGreater(NaN, INF));
    }
}

TEST_CASE("pareto_compare_less", "[math]")
{
    auto [abs, rel] = GENERATE(std::pair{ 0.0, 0.0 }, std::pair{ 1E-12, 10 * eps<double> });

    ScopedTolerances _(abs, rel);
    INFO("Relative tolerance eps: " << rel << ", absolute tolerance: " << abs);

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
        REQUIRE(!paretoCompareLess(std::vector{ 1.0 }, std::vector{ 1.0 }));

        REQUIRE(paretoCompareLess(std::vector{ 1.0 }, std::vector{ 2.3 }));
        REQUIRE(!paretoCompareLess(std::vector{ 2.3 }, std::vector{ 1.0 }));
    }

    SECTION("empty")
    {
        REQUIRE(!paretoCompareLess({}, {}));
    }
}

TEST_CASE("pareto_compare_three_way", "[math]")
{
    auto [abs, rel] = GENERATE(std::pair{ 0.0, 0.0 }, std::pair{ 1E-12, 10 * eps<double> });

    ScopedTolerances _(abs, rel);
    INFO("Relative tolerance eps: " << rel << ", absolute tolerance: " << abs);

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
        REQUIRE(paretoCompare(std::vector{ 1.0 }, std::vector{ 1.0 }) == 0);

        REQUIRE(paretoCompare(std::vector{ 1.0 }, std::vector{ 2.3 }) < 0);
        REQUIRE(paretoCompare(std::vector{ 2.3 }, std::vector{ 1.0 }) > 0);
    }

    SECTION("empty")
    {
        REQUIRE(paretoCompare({}, {}) == 0);
    }
}

TEST_CASE("euclidean_norm", "[math]")
{
    REQUIRE(euclideanNorm({}) == 0.0);
    REQUIRE(euclideanNorm(std::vector{ 1.0, 4.5, 3.2, 0.3 }) == Approx(5.62).margin(0.01));
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

    REQUIRE(euclideanDistanceSq(std::vector{ 3.0 }, std::vector{ 1.0 }) == Approx(4.0).margin(0.01));
    REQUIRE(euclideanDistanceSq(std::vector{ 1.0, 0.0 }, std::vector{ 2.0, 1.0 }) == Approx(2.0).margin(0.01));
    REQUIRE(euclideanDistanceSq(std::vector{ 1.0, 2.0, 0.0 }, std::vector{ 3.0, 0.0, 1.0 }) == Approx(9.0).margin(0.01));
}

TEST_CASE("perpendicular_distance", "[math]")
{
    REQUIRE(perpendicularDistanceSq({}, {}) == Approx(0.0).margin(1E-8));

    REQUIRE(perpendicularDistanceSq(std::vector{ 3.1 }, std::vector{ 0.95 }) == Approx(0.0).margin(1E-8));
    REQUIRE(perpendicularDistanceSq(std::vector{ 0.4, 0.9 }, std::vector{ 2.5, 1.0 }) == Approx(3.53).margin(0.01));
}

TEST_CASE("volume_between", "[math]")
{
    REQUIRE(volumeBetween(std::vector{ 0.0 }, std::vector{ 1.0 }) == Approx(1.0).margin(1E-8));
    REQUIRE(volumeBetween(std::vector{ 1.0 }, std::vector{ 0.0 }) == Approx(1.0).margin(1E-8));

    REQUIRE(volumeBetween(std::vector{ -1.0, -1.0 }, std::vector{ 1.0, 1.0 }) == Approx(4.0).margin(1E-8));
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