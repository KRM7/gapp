#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include <vector>
#include <tuple>
#include <cstddef>
#include "utility/algorithm.hpp"
#include "utility/rng.hpp"

using namespace genetic_algorithm;

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

    SECTION("bad range")
    {
        //REQUIRE_THROWS(detail::argsort(nums.end(), nums.begin()));
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

    SECTION("bad range")
    {
        //REQUIRE_THROWS(detail::partial_argsort(nums.begin(), nums.end(), nums.begin()));
        //REQUIRE_THROWS(detail::partial_argsort(nums.end(), nums.begin(), nums.end()));
    }
}

TEST_CASE("argmax", "[algorithm]")
{
    const std::vector nums = { 4.0, 0.0, 2.0, 5.0, 1.0 };

    REQUIRE(detail::argmax(nums.begin(), nums.end()) == 3);
    REQUIRE(detail::argmax(nums.rbegin(), nums.rend()) == 3);

    REQUIRE(detail::argmax(nums.begin(), nums.end(), std::greater<>{}) == 1);

    REQUIRE(detail::argmax(nums.begin(), nums.begin() + 3) == 0);
    REQUIRE(detail::argmax(nums.begin() + 1, nums.end()) == 2);

    //REQUIRE_THROWS(detail::argmax(nums.begin(), nums.begin()));
    //REQUIRE_THROWS(detail::argmax(nums.end(), nums.begin()));
}

TEST_CASE("argmin", "[algorithm]")
{
    const std::vector nums = { 4.0, 0.0, 2.0, 5.0, 1.0 };

    REQUIRE(detail::argmin(nums.begin(), nums.end()) == 1);
    REQUIRE(detail::argmin(nums.rbegin(), nums.rend()) == 1);

    REQUIRE(detail::argmin(nums.begin(), nums.end(), std::greater<>{}) == 3);

    REQUIRE(detail::argmin(nums.begin() + 2, nums.end()) == 2);

    //REQUIRE_THROWS(detail::argmin(nums.begin(), nums.begin()));
    //REQUIRE_THROWS(detail::argmin(nums.end(), nums.begin()));
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

    SECTION("bad range")
    {
        //REQUIRE_THROWS(detail::partial_shuffle(nums.end(), nums.begin(), nums.end(), rng::prng));
        //REQUIRE_THROWS(detail::partial_shuffle(nums.begin(), nums.end(), nums.begin(), rng::prng));
    }
}

TEST_CASE("contains", "[algorithm]")
{
    const std::vector nums = { 4.0, 0.0, 2.0, 5.0, 1.0 };

    REQUIRE(detail::contains(nums.begin(), nums.end(), 0.0));
    REQUIRE(detail::contains(nums.begin(), nums.end(), 1.0));

    REQUIRE(!detail::contains(nums.begin(), nums.end() - 1, 1.0));
    REQUIRE(!detail::contains(nums.begin(), nums.end(), 0.001));

    //REQUIRE_THROWS(detail::contains(nums.end(), nums.begin(), 0.0));
}

TEST_CASE("find_all", "[algorithm]")
{
    const std::vector nums = { 4, 0, 2, 5, 1 };

    const auto odd = detail::find_all(nums.begin(), nums.end(), [](int n) { return n % 2; });
    REQUIRE(odd == std::vector{ nums.end() - 2, nums.end() - 1 });

    const auto big = detail::find_all(nums.begin(), nums.end(), [](int n) { return n > 10; });
    REQUIRE(big.empty());

    REQUIRE(detail::find_all(nums.begin(), nums.begin(), [](int) { return true; }).empty());
    REQUIRE(detail::find_all(nums.begin(), nums.end(), [](int) { return false; }).empty());

    //REQUIRE_THROWS(detail::find_all(nums.end(), nums.begin(), std::logical_not<>{}));
}

TEST_CASE("find_all_v", "[algorithm]")
{
    const std::vector nums = { 4, 0, 2, 5, 1 };

    const auto odd = detail::find_all_v(nums.begin(), nums.end(), [](int n) { return n % 2; });
    REQUIRE(odd == std::vector{ 5, 1 });

    const auto big = detail::find_all_v(nums.begin(), nums.end(), [](int n) { return n > 10; });
    REQUIRE(big.empty());

    REQUIRE(detail::find_all_v(nums.begin(), nums.begin(), [](int) { return true; }).empty());
    REQUIRE(detail::find_all_v(nums.begin(), nums.end(), [](int) { return false; }).empty());

    //REQUIRE_THROWS(detail::find_all_v(nums.end(), nums.begin(), std::logical_not<>{}));
}

TEST_CASE("find_indices", "[algorithm]")
{
    const std::vector nums = { 4, 0, 2, 5, 1 };

    const auto odd = detail::find_indices(nums, [](int n) { return n % 2; });
    REQUIRE(odd == std::vector<size_t>{ 3, 4 });

    const auto big = detail::find_indices(nums, [](int n) { return n > 10; });
    REQUIRE(big.empty());

    const auto all = detail::find_indices(nums, [](int) { return true; });
    REQUIRE(all == std::vector<size_t>{ 0, 1, 2, 3, 4 });

    const auto none = detail::find_indices(nums, [](int) { return false; });
    REQUIRE(none.empty());
}

TEST_CASE("index_of", "[algorithm]")
{
    const std::vector nums = { 4, 0, 2, 5, 1 };

    REQUIRE(detail::index_of(nums, 4) == 0);
    REQUIRE(detail::index_of(nums, 2) == 2);
    REQUIRE(detail::index_of(nums, 1) == 4);
    REQUIRE(!detail::index_of(nums, 7).has_value());
}

TEST_CASE("find_index", "[algorithm]")
{
    const std::vector nums = { 4, 0, 2, 5, 1 };

    auto first = detail::find_index(nums, [](int) { return true; });
    REQUIRE(first == 0);

    auto none = detail::find_index(nums, [](int) { return false; });
    REQUIRE(!none.has_value());

    auto odd = detail::find_index(nums, [](int i) { return i % 2; });
    REQUIRE(odd == 3);

    auto six = detail::find_index(nums, [](int i) { return i == 6; });
    REQUIRE(!six.has_value());
}

TEST_CASE("elementwise_min", "[algorithm]")
{
    const std::vector nums1 = { 4, 0, 2, 5, 1 };
    const std::vector nums2 = { 2, 3, 1, 6, 0 };

    auto min = detail::elementwise_min(nums1, nums2);
    REQUIRE(min == std::vector{ 2, 0, 1, 5, 0 });

    //REQUIRE_THROWS(detail::elementwise_min(nums1, std::vector{ 1, 2, 3}));
}

TEST_CASE("elementwise_max", "[algorithm]")
{
    const std::vector nums1 = { 4, 0, 2, 5, 1 };
    const std::vector nums2 = { 2, 3, 1, 6, 0 };

    auto max = detail::elementwise_max(nums1, nums2);
    REQUIRE(max == std::vector{ 4, 3, 2, 6, 1 });

    //REQUIRE_THROWS(detail::elementwise_max(nums1, std::vector{ 1, 2, 3}));
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

    selected = detail::select(nums, {});
    REQUIRE(selected.empty());

    selected = detail::select(std::vector{ 1, 3, 5 }, { 0, 1 });
    REQUIRE(selected == std::vector{ 1, 3 });

    selected = detail::select<int>({}, {});
    REQUIRE(selected.empty());

    //REQUIRE_THROWS(detail::select({}, { 1 }));
    //REQUIRE_THROWS(detail::select(nums, { 71 }));
}

TEST_CASE("erase_duplicates", "[algorithm]")
{
    std::vector nums = { 1, 0, 1, 5, 1, 3, 1 };

    detail::erase_duplicates(nums);
    REQUIRE_THAT(nums, Catch::Matchers::UnorderedEquals(std::vector{ 0, 1, 3, 5 }));
}

TEST_CASE("transform_reduce", "[algorithm]")
{
    constexpr std::tuple vals = { 0, 3.14, 'a', "abcd", size_t{ 2 } };
    
    constexpr int num_arithmetic = detail::transform_reduce(vals, 0, std::identity{},
    [](int acc, const auto& val)
    {
        if constexpr (std::is_arithmetic_v<std::remove_reference_t<decltype(val)>>)
        {
            return acc + 1;
        }
        else
        {
            return acc;
        }
    });

    STATIC_REQUIRE(num_arithmetic == 4);

    const int sum = detail::transform_reduce(std::tuple{}, 1, std::identity{}, std::plus<>{});
    REQUIRE(sum == 1);
}