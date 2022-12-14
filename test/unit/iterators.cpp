/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include "utility/iterators.hpp"
#include <vector>
#include <algorithm>
#include <type_traits>

using namespace genetic_algorithm::detail;

TEST_CASE("stable_iterator", "[iterators]")
{
    using iterator = stable_iterator<std::vector<int>>;
    using const_iterator = const_stable_iterator<std::vector<int>>;

    SECTION("member typedefs")
    {
        STATIC_REQUIRE(std::is_same_v<iterator::value_type, int>);
        STATIC_REQUIRE(std::is_same_v<iterator::reference, int&>);
        STATIC_REQUIRE(std::is_same_v<iterator::pointer, int*>);

        STATIC_REQUIRE(std::is_same_v<const_iterator::value_type, int>);
        STATIC_REQUIRE(std::is_same_v<const_iterator::reference, const int&>);
        STATIC_REQUIRE(std::is_same_v<const_iterator::pointer, const int*>);
    }

    SECTION("constructors")
    {
        STATIC_REQUIRE(std::is_default_constructible_v<iterator>);
        STATIC_REQUIRE(std::is_copy_constructible_v<iterator>);
        STATIC_REQUIRE(std::is_copy_assignable_v<iterator>);
        STATIC_REQUIRE(std::is_nothrow_move_constructible_v<iterator>);
        STATIC_REQUIRE(std::is_nothrow_move_assignable_v<iterator>);

        STATIC_REQUIRE(std::is_default_constructible_v<const_iterator>);
        STATIC_REQUIRE(std::is_copy_constructible_v<const_iterator>);
        STATIC_REQUIRE(std::is_copy_assignable_v<const_iterator>);
        STATIC_REQUIRE(std::is_nothrow_move_constructible_v<const_iterator>);
        STATIC_REQUIRE(std::is_nothrow_move_assignable_v<const_iterator>);
    }

    std::vector nums = { 1, 3, 4, 2, 8 };

    iterator first(nums, 0);
    iterator last(nums, nums.size());
    iterator it;

    const_iterator cfirst(nums, 0);
    const_iterator clast(nums, nums.size());
    const_iterator cit;

    SECTION("factory functions")
    {
        const std::vector<int>& cnums = nums;

        REQUIRE(first == stable_begin(nums));
        REQUIRE(last == stable_end(nums));

        REQUIRE(cfirst == stable_begin(cnums));
        REQUIRE(clast == stable_end(cnums));

        REQUIRE(cfirst == stable_cbegin(nums));
        REQUIRE(clast == stable_cend(nums));
    }

    SECTION("const conversion")
    {
        const_iterator const_copy = first;

        REQUIRE(const_copy == first);
        REQUIRE(*const_copy == *first);
    }

    SECTION("dereference")
    {
        REQUIRE(*first == 1);
        REQUIRE(*cfirst == 1);

        //REQUIRE_THROWS(*last);
        //REQUIRE_THROWS(*clast);

        //REQUIRE_THROWS(*it);
        //REQUIRE_THROWS(*cit);
    }

    SECTION("assignment")
    {
        *first = 7;
        REQUIRE(nums[0] == 7);
    }

    SECTION("comparisons")
    {
        REQUIRE(cfirst == first);
        REQUIRE(it == cit);

        REQUIRE(first < clast);
        REQUIRE(first <= last);
        REQUIRE(last > first);
        REQUIRE(clast >= cfirst);
        REQUIRE(first != last);

        //REQUIRE_THROWS(first == it);
        //REQUIRE_THROWS(cfirst < it);
    }

    SECTION("advance")
    {
        REQUIRE(*++first == 3);
        REQUIRE(*++first == 4);
        REQUIRE(*--first == 3);
        REQUIRE(*first-- == 3);
        REQUIRE(*first == 1);

        //REQUIRE_THROWS(--first);
        //REQUIRE_THROWS(cfirst--);
        //REQUIRE_THROWS(++last);
        //REQUIRE_THROWS(clast++);

        //REQUIRE_THROWS(++it);
        //REQUIRE_THROWS(--it);
    }

    SECTION("arithmetic")
    {
        REQUIRE(*(first + 2) == 4);
        REQUIRE(last - nums.size() == first);
        REQUIRE(*(first += 3) == 2);
        REQUIRE(*first == 2);

        REQUIRE(size_t(clast - cfirst) == nums.size());
        REQUIRE(size_t(last - cfirst) == nums.size());

        //REQUIRE_THROWS(it + 3);
        //REQUIRE_THROWS(cit + 0);

        //REQUIRE_THROWS(first + 21);
        //REQUIRE_THROWS(first - 8);
        //REQUIRE_THROWS(clast += 174);
    }

    SECTION("algorithms")
    {
        std::sort(first, last);
        REQUIRE(nums == std::vector{ 1, 2, 3, 4, 8 });

        auto res = std::find(cfirst, clast, 8);
        REQUIRE(res == --last);
    }
}

TEST_CASE("iota_iterator", "[iterators]")
{
    iota_iterator first(1);
    iota_iterator last(5);

    SECTION("dereference")
    {
        REQUIRE(*first == 1);
        REQUIRE(*last == 5);

        REQUIRE(*iota_iterator<int>{} == 0);
    }

    SECTION("comparison")
    {
        REQUIRE(first != last);
        REQUIRE(first < last);
        REQUIRE(last >= first);
    }

    SECTION("advance")
    {
        REQUIRE(*++first == 2);
        REQUIRE(*first++ == 2);
        REQUIRE(*first == 3);
    }

    SECTION("arithmetic")
    {
        REQUIRE((first + 4) == last);
        REQUIRE(*(first += 2) == 3);
        REQUIRE(*(last - 2) == *first);
    }

    SECTION("algorithms")
    {
        auto r1 = std::find(first, last, 3);
        auto r2 = std::find(first, last, 7);

        REQUIRE(*r1 == 3);
        REQUIRE(*r2 == 5);
    }
}