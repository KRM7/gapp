/* Copyright (c) 2024 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include "utility/dynamic_bitset.hpp"
#include <utility>

using gapp::detail::dynamic_bitset;


TEST_CASE("constructor_default", "[dynamic_bitset]")
{
    dynamic_bitset bitset;

    REQUIRE(bitset.empty());
    REQUIRE(bitset.size() == 0);
}

TEST_CASE("constructor_size", "[dynamic_bitset]")
{
    dynamic_bitset bitset1(0);

    REQUIRE(bitset1.empty());
    REQUIRE(bitset1.size() == 0);

    dynamic_bitset bitset2(3);

    REQUIRE(!bitset2.empty());
    REQUIRE(bitset2.size() == 3);
}

TEST_CASE("constructor_size_value", "[dynamic_bitset]")
{
    dynamic_bitset bitset(100, true);

    REQUIRE(bitset[0] == true);
    REQUIRE(bitset.size() == 100);
}

TEST_CASE("index_operator", "[dynamic_bitset]")
{
    dynamic_bitset bitset(100, false);

    bitset[80] = true;

    REQUIRE(bitset[79] == false);
    REQUIRE(bitset[80] == true);
    REQUIRE(bitset[81] == false);

    REQUIRE(std::as_const(bitset)[80]);

    bitset[80].flip();

    REQUIRE(bitset[79] == false);
    REQUIRE(bitset[80] == false);
    REQUIRE(bitset[81] == false);
}

TEST_CASE("clear", "[dynamic_bitset]")
{
    dynamic_bitset bitset(100, true);
    bitset.clear();

    REQUIRE(bitset.empty());
}

TEST_CASE("resize", "[dynamic_bitset]")
{
    dynamic_bitset bitset(100, true);

    bitset.resize(99);

    REQUIRE(bitset.size() == 99);
    REQUIRE(bitset[98] == true);

    bitset.resize(100, false);

    REQUIRE(bitset.size() == 100);
    REQUIRE(bitset[98] == true);
    REQUIRE(bitset[99] == false);

    bitset.resize(101, true);

    REQUIRE(bitset.size() == 101);
    REQUIRE(bitset[98] == true);
    REQUIRE(bitset[99] == false);
    REQUIRE(bitset[100] == true);
}

TEST_CASE("find_first_true", "[dynamic_bitset]")
{
    dynamic_bitset bitset(100, false);
    REQUIRE(bitset.find_first(true) == 100);

    bitset[99] = true;
    REQUIRE(bitset.find_first(true) == 99);

    bitset[64] = true;
    REQUIRE(bitset.find_first(true) == 64);

    bitset[0] = true;
    REQUIRE(bitset.find_first(true) == 0);
}

TEST_CASE("find_first_false", "[dynamic_bitset]")
{
    dynamic_bitset bitset(100, true);
    REQUIRE(bitset.find_first(false) == 100);

    bitset[99] = false;
    REQUIRE(bitset.find_first(false) == 99);

    bitset[64] = false;
    REQUIRE(bitset.find_first(false) == 64);

    bitset[0] = false;
    REQUIRE(bitset.find_first(false) == 0);
}

TEST_CASE("popcount", "[dynamic_bitset]")
{
    dynamic_bitset bitset(100, false);
    REQUIRE(bitset.popcount() == 0);

    bitset.resize(101, true);
    REQUIRE(bitset.popcount() == 1);

    bitset.resize(100, true);
    REQUIRE(bitset.popcount() == 0);

    bitset[0] = true;
    REQUIRE(bitset.popcount() == 1);

    bitset[99] = true;
    REQUIRE(bitset.popcount() == 2);
}

TEST_CASE("fill", "[dynamic_bitset]")
{
    dynamic_bitset bitset(100, true);

    bitset.fill(false);
    REQUIRE(bitset.popcount() == 0);

    bitset.fill(true);
    REQUIRE(bitset.popcount() == 100);
}

TEST_CASE("any_set", "[dynamic_bitset]")
{
    dynamic_bitset bitset(100, false);
    REQUIRE(!bitset.any_set());

    bitset[99] = true;
    REQUIRE(bitset.any_set());
}

TEST_CASE("all_set", "[dynamic_bitset]")
{
    dynamic_bitset bitset(100, false);
    REQUIRE(!bitset.all_set());

    bitset[99] = true;
    REQUIRE(!bitset.all_set());

    bitset.fill(true);
    REQUIRE(bitset.all_set());
}

TEST_CASE("none_set", "[dynamic_bitset]")
{
    dynamic_bitset bitset(100, false);
    REQUIRE(bitset.none_set());

    bitset[99] = true;
    REQUIRE(!bitset.none_set());
}

TEST_CASE("operator~", "[dynamic_bitset]")
{
    dynamic_bitset bitset1(100, false);
    REQUIRE(bitset1.popcount() == 0);

    dynamic_bitset bitset2 = ~bitset1;
    REQUIRE(bitset2.popcount() == 100);
}
