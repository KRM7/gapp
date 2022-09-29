#include <catch2/catch_test_macros.hpp>
#include "utility/algorithm.hpp"

using namespace genetic_algorithm;

TEST_CASE("index_vector", "[algorithm]")
{
    auto indices = detail::index_vector(3);

    REQUIRE(indices.size() == 3);
    REQUIRE(indices[0] == 0);
    REQUIRE(indices[2] == 2);

    indices = detail::index_vector(4, 1);

    REQUIRE(indices.size() == 4);
    REQUIRE(indices[0] == 1);
    REQUIRE(indices[3] == 4);
}

TEST_CASE("argsort", "[algorithm]")
{
    const std::vector nums = { 4.0, 0.0, 2.0, 1.0 };

    SECTION("simple ascending")
    {
        const auto indices = detail::argsort(nums.begin(), nums.end());

        REQUIRE(indices ==  std::vector<size_t>{ 1, 3, 2, 0 });
    }

    SECTION("descending with reverse_iterators")
    {
        const auto indices = detail::argsort(nums.rbegin(), nums.rend());

        REQUIRE(indices ==  std::vector<size_t>{ 0, 2, 3, 1 });
    }

    SECTION("descending with comparison function")
    {
        const auto indices = detail::argsort(nums.begin(), nums.end(), std::greater<>{});

        REQUIRE(indices ==  std::vector<size_t>{ 0, 2, 3, 1 });
    }

    SECTION("empty range")
    {
        const auto indices = detail::argsort(nums.begin(), nums.begin());

        REQUIRE(indices.empty());
    }

    SECTION("invalid range")
    {
        decltype(nums) another_range;

        REQUIRE_THROWS(detail::argsort(nums.end(), nums.begin()));
        REQUIRE_THROWS(detail::argsort(nums.begin(), another_range.begin()));
    }
}

TEST_CASE("partial_argsort", "[algorithm]")
{
    
}