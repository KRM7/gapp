/* Copyright (c) 2025 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include "utility/bit.hpp"
#include <utility>
#include <cstdint>

using namespace gapp::detail;

TEST_CASE("float_bits", "[bit]")
{
    STATIC_REQUIRE(exponent_bits<float> == 8);
    STATIC_REQUIRE(exponent_bits<double> == 11);

    STATIC_REQUIRE(mantissa_bits<float> == 23);
    STATIC_REQUIRE(mantissa_bits<double> == 52);

    STATIC_REQUIRE(implicit_mantissa_bits<float> == 24);
    STATIC_REQUIRE(implicit_mantissa_bits<double> == 53);
}

TEST_CASE("set_sign_bit", "[bit]")
{
    STATIC_REQUIRE(set_sign_bit(3.1f, 0) == 3.1f);
    STATIC_REQUIRE(set_sign_bit(3.1f, 1) == -3.1f);

    STATIC_REQUIRE(set_sign_bit(3.1, 0) == 3.1);
    STATIC_REQUIRE(set_sign_bit(3.1, 1) == -3.1);

    STATIC_REQUIRE(set_sign_bit(-1.2f, 0) == 1.2f);
    STATIC_REQUIRE(set_sign_bit(-1.2f, 1) == -1.2f);

    STATIC_REQUIRE(set_sign_bit(-1.2, 0) == 1.2);
    STATIC_REQUIRE(set_sign_bit(-1.2, 1) == -1.2);
}

TEST_CASE("is_nth_bit_set", "[bit]")
{
    STATIC_REQUIRE(is_nth_bit_set(1, 0));
    STATIC_REQUIRE(!is_nth_bit_set(0, 0));

    STATIC_REQUIRE(is_nth_bit_set(0b0101, 2));
    STATIC_REQUIRE(!is_nth_bit_set(0b0101, 3));
}

TEST_CASE("msb", "[bit]")
{
    STATIC_REQUIRE(msb(0) == 0);
    STATIC_REQUIRE(msb(-1) == 1);
}

TEST_CASE("lsb", "[bit]")
{
    STATIC_REQUIRE(lsb(0) == 0);
    STATIC_REQUIRE(lsb(-1) == 1);
}

TEST_CASE("mask_right_n", "[bit]")
{
    STATIC_REQUIRE(mask_right_n<int64_t>(0) == 0);
    STATIC_REQUIRE(mask_right_n<int64_t>(64) == -1);
    STATIC_REQUIRE(mask_right_n<int64_t>(4) == 0b1111);
}

TEST_CASE("mask_left_n", "[bit]")
{
    STATIC_REQUIRE(mask_left_n<int64_t>(0) == 0);
    STATIC_REQUIRE(mask_left_n<int64_t>(64) == -1);
    STATIC_REQUIRE(mask_left_n<int64_t>(4) == int64_t(0xF000'0000'0000'0000));
}

TEST_CASE("block_of", "[bit]")
{
    STATIC_REQUIRE(block_of<int64_t>(0) == 0);
    STATIC_REQUIRE(block_of<int64_t>(1) == -1);
}

TEST_CASE("extract_bits", "[bit]")
{
    STATIC_REQUIRE(extract_bits<63, 64>(int64_t(-1)) == 1);
    STATIC_REQUIRE(extract_bits<63, 64>(int64_t(0)) == 0);

    STATIC_REQUIRE(extract_bits<60, 64>(int64_t(-1)) == 0b1111);
    STATIC_REQUIRE(extract_bits<13, 41>(int64_t(0)) == 0);

    STATIC_REQUIRE(extract_bits<1, 7>(int64_t(0b1111'1010)) == 0b0011'1101);
}
