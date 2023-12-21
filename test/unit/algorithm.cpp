/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include <vector>
#include <tuple>
#include <cstddef>
#include "utility/algorithm.hpp"
#include "utility/rng.hpp"
#include "utility/utility.hpp"


static constexpr auto is_odd = [](int n) { return n % 2; };
static constexpr auto is_big = [](int n) { return n > 10; };

static constexpr auto always_true = [](auto) { return true; };
static constexpr auto always_false = [](auto) { return false; };


using namespace gapp;

TEST_CASE("increment_mod", "[algorithm]")
{
    int n = 0;

    detail::increment_mod(n, 2);
    REQUIRE(n == 1);

    detail::increment_mod(n, 2);
    REQUIRE(n == 0);
}

TEST_CASE("index_vector", "[algorithm]")
{
    REQUIRE(detail::index_vector(3) == std::vector<size_t>{ 0, 1, 2 });
    REQUIRE(detail::index_vector(4, 2) == std::vector<size_t>{ 2, 3, 4, 5 });
}

TEST_CASE("argsort", "[algorithm]")
{
    const std::vector nums = { 4.0, 0.0, 2.0, 1.0 };

    SECTION("iterators")
    {
        const auto indices = detail::argsort(nums.begin(), nums.end());
        REQUIRE(indices == std::vector<size_t>{ 1, 3, 2, 0 });
    }

    SECTION("reverse iterators")
    {
        const auto indices = detail::argsort(nums.rbegin(), nums.rend());
        REQUIRE(indices == std::vector<size_t>{ 0, 2, 3, 1 });
    }

    SECTION("compare function")
    {
        const auto indices = detail::argsort(nums.begin(), nums.end(), std::greater<>{});
        REQUIRE(indices == std::vector<size_t>{ 0, 2, 3, 1 });
    }

    SECTION("empty range")
    {
        const auto indices = detail::argsort(nums.begin(), nums.begin());
        REQUIRE(indices.empty());
    }
}

TEST_CASE("partial_argsort", "[algorithm]")
{
    const std::vector nums = { 4.0, 0.0, 2.0, 1.0, 5.0 };

    SECTION("iterators")
    {
        const auto indices = detail::partial_argsort(nums.begin(), nums.begin() + 2, nums.end());

        REQUIRE(indices.size() == nums.size());
        REQUIRE(indices[0] == 1);
        REQUIRE(indices[1] == 3);
    }

    SECTION("reverse iterators")
    {
        const auto indices = detail::partial_argsort(nums.rbegin(), nums.rbegin() + 2, nums.rend());

        REQUIRE(indices.size() == nums.size());
        REQUIRE(indices[0] == 4);
        REQUIRE(indices[1] == 0);
    }

    SECTION("compare function")
    {
        const auto indices = detail::partial_argsort(nums.begin(), nums.begin() + 2, nums.end(), std::greater<>{});

        REQUIRE(indices.size() == nums.size());
        REQUIRE(indices[0] == 4);
        REQUIRE(indices[1] == 0);
    }

    SECTION("empty range")
    {
        auto indices = detail::partial_argsort(nums.begin(), nums.begin(), nums.begin());
        REQUIRE(indices.empty());
    }
}

TEST_CASE("max_element", "[algorithm]")
{
    const std::vector nums = { 4.0, 0.0, 2.0, 5.0, 1.0 };

    REQUIRE(*detail::max_element(nums.begin(), nums.end()) == 5.0);
    REQUIRE(*detail::max_element(nums.rbegin(), nums.rend()) == 5.0);

    REQUIRE(detail::max_element(nums.begin(), nums.begin()) == nums.begin());

    REQUIRE(*detail::max_element(nums.begin(), nums.end(), std::negate{}) == 0.0);
}

TEST_CASE("min_element", "[algorithm]")
{
    const std::vector nums = { 4.0, 0.0, 2.0, 5.0, 1.0 };

    REQUIRE(*detail::min_element(nums.begin(), nums.end()) == 0.0);
    REQUIRE(*detail::min_element(nums.rbegin(), nums.rend()) == 0.0);

    REQUIRE(detail::min_element(nums.begin(), nums.begin()) == nums.begin());

    REQUIRE(*detail::min_element(nums.begin(), nums.end(), std::negate{}) == 5.0);
}

TEST_CASE("argmax", "[algorithm]")
{
    const std::vector nums = { 4.0, 0.0, 2.0, 5.0, 1.0 };

    REQUIRE(detail::argmax(nums.begin(), nums.end()) == 3);
    REQUIRE(detail::argmax(nums.rbegin(), nums.rend()) == 3);

    REQUIRE(detail::argmax(nums.begin(), nums.end(), std::negate{}) == 1);

    REQUIRE(detail::argmax(nums.begin(), nums.begin() + 3) == 0);
    REQUIRE(detail::argmax(nums.begin() + 1, nums.end()) == 2);
}

TEST_CASE("argmin", "[algorithm]")
{
    const std::vector nums = { 4.0, 0.0, 2.0, 5.0, 1.0 };

    REQUIRE(detail::argmin(nums.begin(), nums.end()) == 1);
    REQUIRE(detail::argmin(nums.rbegin(), nums.rend()) == 1);

    REQUIRE(detail::argmin(nums.begin(), nums.end(), std::negate{}) == 3);

    REQUIRE(detail::argmin(nums.begin() + 2, nums.end()) == 2);
}

TEST_CASE("max", "[algorithm]")
{
    REQUIRE(detail::max(1, 2) == 2);
    REQUIRE(detail::max(0, 6, -10) == 6);

    REQUIRE(detail::max<int>(-1, 3u) == 3);
}

TEST_CASE("min", "[algorithm]")
{
    REQUIRE(detail::min(1, 2) == 1);
    REQUIRE(detail::min(0, 6, -10) == -10);

    REQUIRE(detail::min<int>(-1, 3u) == -1);
}

TEST_CASE("partial_shuffle", "[algorithm]")
{
    std::vector nums = { 4.0, 0.0, 2.0, 5.0, 1.0 };

    SECTION("Shuffle empty subrange")
    {
        detail::partial_shuffle(nums.begin(), nums.begin(), nums.end(), rng::prng);
        REQUIRE_THAT(nums, Catch::Matchers::Equals(std::vector{ 4.0, 0.0, 2.0, 5.0, 1.0 }));
    }

    SECTION("Shuffle subrange")
    {
        detail::partial_shuffle(nums.begin(), nums.end() - 2, nums.end(), rng::prng);
        REQUIRE_THAT(nums, Catch::Matchers::UnorderedEquals(std::vector{ 4.0, 0.0, 2.0, 5.0, 1.0 }));
    }

    SECTION("Shuffle entire range")
    {
        detail::partial_shuffle(nums.begin(), nums.end(), nums.end(), rng::prng);
        REQUIRE_THAT(nums, Catch::Matchers::UnorderedEquals(std::vector{ 4.0, 0.0, 2.0, 5.0, 1.0 }));
    }
}

TEST_CASE("contains", "[algorithm]")
{
    const std::vector nums = { 4.0, 0.0, 2.0, 5.0, 1.0 };

    REQUIRE(detail::contains(nums.begin(), nums.end(), 0.0));
    REQUIRE(detail::contains(nums.begin(), nums.end(), 1.0));

    REQUIRE(!detail::contains(nums.begin(), nums.end() - 1, 1.0));
    REQUIRE(!detail::contains(nums.begin(), nums.end(), 0.001));
}

TEST_CASE("find_all", "[algorithm]")
{
    const std::vector nums = { 4, 0, 2, 5, 1 };

    const auto odd_nums = detail::find_all(nums.begin(), nums.end(), is_odd);
    REQUIRE(odd_nums == std::vector{ 5, 1 });

    const auto big_nums = detail::find_all(nums.begin(), nums.end(), is_big);
    REQUIRE(big_nums.empty());

    REQUIRE(detail::find_all(nums.begin(), nums.end(), always_true).size() == nums.size());
    REQUIRE(detail::find_all(nums.begin(), nums.end(), always_false).empty());
}

TEST_CASE("find_indices", "[algorithm]")
{
    const std::vector nums = { 4, 0, 2, 5, 1 };

    const auto odd_num_idxs = detail::find_indices(nums, is_odd);
    REQUIRE(odd_num_idxs == std::vector{ 3_sz, 4_sz });

    const auto big_num_idxs = detail::find_indices(nums, is_big);
    REQUIRE(big_num_idxs.empty());

    const auto all = detail::find_indices(nums, always_true);
    REQUIRE(all == std::vector{ 0_sz, 1_sz, 2_sz, 3_sz, 4_sz });

    const auto none = detail::find_indices(nums, always_false);
    REQUIRE(none.empty());
}

TEST_CASE("index_of", "[algorithm]")
{
    const std::vector nums = { 4, 0, 2, 5, 1 };

    REQUIRE(detail::index_of(nums, 4) == 0_sz);
    REQUIRE(detail::index_of(nums, 2) == 2_sz);
    REQUIRE(detail::index_of(nums, 1) == 4_sz);
    REQUIRE(!detail::index_of(nums, 7).has_value());
}

TEST_CASE("find_index", "[algorithm]")
{
    const std::vector nums = { 4, 0, 2, 5, 1 };

    const auto first_idx = detail::find_index(nums, always_true);
    REQUIRE(first_idx == 0_sz);

    const auto none = detail::find_index(nums, always_false);
    REQUIRE(!none.has_value());

    const auto first_odd_idx = detail::find_index(nums, is_odd);
    REQUIRE(first_odd_idx == 3_sz);

    const auto six = detail::find_index(nums, [](int i) { return i == 6; });
    REQUIRE(!six.has_value());
}

TEST_CASE("elementwise_min", "[algorithm]")
{
    std::vector nums1 = { 4, 0, 2, 5, 1 };
    std::vector nums2 = { 2, 3, 1, 6, 0 };

    detail::elementwise_min(nums1, nums2, detail::inplace_t{});
    REQUIRE(nums1 == std::vector{ 2, 0, 1, 5, 0 });
}

TEST_CASE("elementwise_max", "[algorithm]")
{
    std::vector nums1 = { 4, 0, 2, 5, 1 };
    std::vector nums2 = { 2, 3, 1, 6, 0 };

    detail::elementwise_max(nums1, nums2, detail::inplace_t{});
    REQUIRE(nums1 == std::vector{ 4, 3, 2, 6, 1 });
}

TEST_CASE("erase_first_stable", "[algorithm]")
{
    std::vector nums = { 4, 0, 2, 5, 1, 3, 1 };

    REQUIRE(detail::erase_first_stable(nums, 0));
    REQUIRE(nums == std::vector{ 4, 2, 5, 1, 3, 1 });

    REQUIRE(detail::erase_first_stable(nums, 1));
    REQUIRE(nums == std::vector{ 4, 2, 5, 3, 1 });

    REQUIRE(!detail::erase_first_stable(nums, 7));
    REQUIRE(nums == std::vector{ 4, 2, 5, 3, 1 });

    std::vector<int> empty_vec;
    REQUIRE(!detail::erase_first_stable(empty_vec, 3));
    REQUIRE(empty_vec.empty());
}

TEST_CASE("select", "[algorithm]")
{
    const std::vector nums = { 4, 0, 2, 5, 1, 3, 1 };

    auto selected = detail::select(nums, { 0, 1, 4 });
    REQUIRE(selected == std::vector{ 4, 0, 1 });

    selected = detail::select(selected, { 2 });
    REQUIRE(selected == std::vector{ 1 });

    selected = detail::select(nums, { });
    REQUIRE(selected.empty());

    selected = detail::select(std::vector{ 1, 3, 5 }, { 0, 1 });
    REQUIRE(selected == std::vector{ 1, 3 });
}

TEST_CASE("erase_duplicates", "[algorithm]")
{
    std::vector nums = { 1, 0, 1, 5, 1, 3, 1 };

    detail::erase_duplicates(nums);
    REQUIRE_THAT(nums, Catch::Matchers::UnorderedEquals(std::vector{ 0, 1, 3, 5 }));
}