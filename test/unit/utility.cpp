/* Copyright (c) 2024 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include "utility/utility.hpp"
#include <type_traits>
#include <limits>
#include <cstdint>

using namespace gapp;
using namespace gapp::detail;


TEMPLATE_TEST_CASE("range_length_signed", "[utility]", std::int8_t, std::int16_t, std::int32_t, std::int64_t)
{
    using IntType  = TestType;
    using UIntType = std::make_unsigned_t<IntType>;

    constexpr IntType small = std::numeric_limits<IntType>::min();
    constexpr IntType large = std::numeric_limits<IntType>::max();

    STATIC_REQUIRE(range_length(IntType{ 0 }, IntType{ 0 }) == 0);
    STATIC_REQUIRE(range_length(IntType{ -1 }, IntType{ -1 }) == 0);
    STATIC_REQUIRE(range_length(IntType{ 1 }, IntType{ 1 }) == 0);

    STATIC_REQUIRE(range_length(IntType{ 0 }, IntType{ 1 }) == 1);
    STATIC_REQUIRE(range_length(IntType{ -1 }, IntType{ 0 }) == 1);

    STATIC_REQUIRE(range_length(IntType{ -1 }, IntType{ 1 }) == 2);

    STATIC_REQUIRE(range_length(small, small) == 0);
    STATIC_REQUIRE(range_length(large, large) == 0);

    STATIC_REQUIRE(range_length(small, IntType{ 0 }) == (UIntType(large) + 1));
    STATIC_REQUIRE(range_length(IntType{ 0 }, large) == UIntType(large));
    STATIC_REQUIRE(range_length(small, large) == std::numeric_limits<UIntType>::max());
}

TEMPLATE_TEST_CASE("range_length_unsigned", "[utility]", std::uint8_t, std::uint16_t, std::uint32_t, std::uint64_t)
{
    using IntType = TestType;

    constexpr IntType small = std::numeric_limits<IntType>::min();
    constexpr IntType large = std::numeric_limits<IntType>::max();

    STATIC_REQUIRE(range_length(IntType{ 0 }, IntType{ 0 }) == 0);
    STATIC_REQUIRE(range_length(IntType{ 1 }, IntType{ 1 }) == 0);

    STATIC_REQUIRE(range_length(IntType{ 0 }, IntType{ 1 }) == 1);

    STATIC_REQUIRE(range_length(small, small) == 0);
    STATIC_REQUIRE(range_length(large, large) == 0);

    STATIC_REQUIRE(range_length(small, large) == large);
}

TEST_CASE("next_mod", "[algorithm]")
{
    REQUIRE(detail::next_mod(0, 3) == 1);
    REQUIRE(detail::next_mod(1, 3) == 2);
    REQUIRE(detail::next_mod(2, 3) == 0);
}

TEST_CASE("prev_mod", "[algorithm]")
{
    REQUIRE(detail::prev_mod(0, 3) == 2);
    REQUIRE(detail::prev_mod(1, 3) == 0);
    REQUIRE(detail::prev_mod(2, 3) == 1);
}

TEST_CASE("increment_mod", "[algorithm]")
{
    int n = 0;

    detail::increment_mod(n, 3);
    REQUIRE(n == 1);

    detail::increment_mod(n, 3);
    REQUIRE(n == 2);

    detail::increment_mod(n, 3);
    REQUIRE(n == 0);
}

TEST_CASE("decrement_mod", "[algorithm]")
{
    int n = 0;

    detail::decrement_mod(n, 3);
    REQUIRE(n == 2);

    detail::decrement_mod(n, 3);
    REQUIRE(n == 1);

    detail::decrement_mod(n, 3);
    REQUIRE(n == 0);
}
