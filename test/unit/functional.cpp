/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include "utility/functional.hpp"
#include <algorithm>
#include <iterator>
#include <type_traits>


template<typename T>
constexpr T square(const T& n) noexcept { return n * n; }

template<typename T>
constexpr T increment(const T& n) noexcept { return n + 1; }

using namespace gapp::detail;


TEST_CASE("map", "[functional]")
{
    const std::vector nums = { 0.0, 1.2, 5.0, 2.5 };

    const auto res = map(nums, [](double n) { return n + 1.5; });
    REQUIRE(res == std::vector{ 1.5, 2.7, 6.5, 4.0 });

    REQUIRE(map(std::vector<int>{}, std::identity{}).empty());
}

TEST_CASE("flatten", "[functional]")
{
    std::vector<std::pair<int, int>> num_pairs = { { 0, 1 }, { 1, 3 }, { 5, 2 }};

    REQUIRE(flatten(num_pairs) == std::vector{ 0, 1, 1, 3, 5, 2 });
    REQUIRE(flatten(std::move(num_pairs)) == std::vector{ 0, 1, 1, 3, 5, 2 });
}

TEST_CASE("arithmetic_funcs", "[functional]")
{
    const std::vector nums = { 1, 2, 4, 2, 9 };
    std::vector<int> out;

    SECTION("multiply")
    {
        std::transform(nums.begin(), nums.end(), std::back_inserter(out), multiply_by(2));
        REQUIRE(out == std::vector{ 2, 4, 8, 4, 18 });
    }

    SECTION("divide")
    {
        std::transform(nums.begin(), nums.end(), std::back_inserter(out), divide_by(2));
        REQUIRE(out == std::vector{ 0, 1, 2, 1, 4 });
    }

    SECTION("add")
    {
        std::transform(nums.begin(), nums.end(), std::back_inserter(out), add(3));
        REQUIRE(out == std::vector{ 4, 5, 7, 5, 12 });
    }

    SECTION("subtract")
    {
        std::transform(nums.begin(), nums.end(), std::back_inserter(out), subtract(1));
        REQUIRE(out == std::vector{ 0, 1, 3, 1, 8 });
    }

    SECTION("multiply add")
    {
        std::transform(nums.begin(), nums.end(), std::back_inserter(out), multiply_add(2, 1));
        REQUIRE(out == std::vector{ 3, 5, 9, 5, 19 });
    }
}

TEST_CASE("comparison_funcs", "[functional]")
{
    const std::vector nums = { 1, 1, 3, 2, 4, 6, 9 };

    SECTION("equal_to")
    {
        auto it = std::find_if(nums.begin(), nums.end(), equal_to(9));
        REQUIRE(it == nums.end() - 1);
    }

    SECTION("not_equal_to")
    {
        auto it = std::find_if(nums.begin(), nums.end(), not_equal_to(1));
        REQUIRE(it == nums.begin() + 2);
    }

    SECTION("greater_than")
    {
        auto it = std::find_if(nums.begin(), nums.end(), greater_than(3.2));
        REQUIRE(it == nums.begin() + 4);
    }

    SECTION("greater_eq_than")
    {
        auto it = std::find_if(nums.begin(), nums.end(), greater_eq_than(6.0));
        REQUIRE(it == nums.begin() + 5);
    }

    SECTION("less_than")
    {
        auto it = std::find_if(nums.begin(), nums.end(), less_than(0));
        REQUIRE(it == nums.end());
    }

    SECTION("less_eq_than")
    {
        auto it = std::find_if(nums.begin(), nums.end(), less_eq_than(-1.0));
        REQUIRE(it == nums.end());
    }

    SECTION("between")
    {
        size_t n = std::count_if(nums.begin(), nums.end(), between(3, 6));
        REQUIRE(n == 3);
    }
}

TEST_CASE("is_size", "[functional]")
{
    std::vector<double> empty;

    REQUIRE(is_size(0)(empty));
    REQUIRE(!is_size(1)(empty));
}

TEST_CASE("element_at", "[functional]")
{
    std::vector<double> vec{ 4.0, 2.0, 3.0 };

    REQUIRE(element_at(0)(vec) == 4.0);
    REQUIRE(element_at(2)(vec) == 3.0);
}