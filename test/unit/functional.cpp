#include <catch2/catch_test_macros.hpp>
#include "utility/functional.hpp"
#include <algorithm>
#include <iterator>
#include <type_traits>


template<typename T>
constexpr T square(const T& n) noexcept { return n * n; }

template<typename T>
constexpr T increment(const T& n) noexcept { return n + 1; }

using namespace genetic_algorithm::detail;

TEST_CASE("lforward", "[functional]")
{
    int n = 2;

    STATIC_REQUIRE(std::is_same_v<int&&, decltype(lforward<int>(n))>);
    STATIC_REQUIRE(std::is_same_v<int&&, decltype(lforward<int>(1))>);
    STATIC_REQUIRE(std::is_same_v<int&&, decltype(lforward<int>(std::move(n)))>);
    STATIC_REQUIRE(std::is_same_v<std::reference_wrapper<int>, decltype(lforward<int&>(n))>);
}

TEST_CASE("compose", "[functional]")
{
    SECTION("functions return by value")
    {
        auto f = compose(square<double>, increment<double>, increment<double>);
        using F = decltype(f);

        STATIC_REQUIRE(std::is_nothrow_invocable_v<F, double>);
        STATIC_REQUIRE(!std::is_invocable_v<F, void*>);
        STATIC_REQUIRE(std::is_same_v<std::invoke_result_t<F, double>, double>);
        STATIC_REQUIRE(std::is_same_v<std::invoke_result_t<F, double&>, double>);

        REQUIRE(f(2.0) == 6.0);
        REQUIRE(f(5) == 27.0);
    }

    SECTION("functions return references")
    {
        auto f = compose(std::identity{}, std::identity{}, std::identity{});
        using F = decltype(f);

        STATIC_REQUIRE(std::is_nothrow_invocable_v<F, double>);
        STATIC_REQUIRE(std::is_invocable_v<F, void*>);
        STATIC_REQUIRE(std::is_same_v<std::invoke_result_t<F, double>, double>);
        STATIC_REQUIRE(std::is_same_v<std::invoke_result_t<F, double&>, double&>);

        REQUIRE(f(1.5) == 1.5);
        REQUIRE(f(nullptr) == nullptr);

        auto in = 1;
        auto& out = f(in);
        REQUIRE(&in == &out);
    }

    SECTION("functions return by value and ref")
    {
        auto f = compose(square<double>, std::identity{}, square<double>, std::identity{});
        using F = decltype(f);

        STATIC_REQUIRE(std::is_nothrow_invocable_v<F, double>);
        STATIC_REQUIRE(!std::is_invocable_v<F, void*>);
        STATIC_REQUIRE(std::is_same_v<std::invoke_result_t<F, double>, double>);
        STATIC_REQUIRE(std::is_same_v<std::invoke_result_t<F, double&>, double>);

        REQUIRE(f(2.0) == 5.0);
    }
}

TEST_CASE("map", "[functional]")
{
    const std::vector nums = { 0.0, 1.2, 5.0, 2.5 };

    const auto res = map(nums, [](double n) { return n + 1.5; });
    REQUIRE(res == std::vector{ 1.5, 2.7, 6.5, 4.0 });

    REQUIRE(map(std::vector<int>{}, std::identity{}) == std::vector<int>{});
}

TEST_CASE("flatten", "[functional]")
{
    std::vector<std::pair<int, int>> num_pairs = { { 0, 1 }, { 1, 3 }, { 5, 2 }};

    SECTION("lvalue")
    {
        REQUIRE(flatten(num_pairs) == std::vector{ 0, 1, 1, 3, 5, 2 });
    }
    SECTION("rvalue")
    {
        REQUIRE(flatten(std::move(num_pairs)) == std::vector{ 0, 1, 1, 3, 5, 2 });
    }
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
}