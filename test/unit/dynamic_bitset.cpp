/* Copyright (c) 2024 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include "utility/dynamic_bitset.hpp"
#include <algorithm>
#include <utility>

using gapp::dynamic_bitset;


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

    bitset.resize(256, true);

    REQUIRE(bitset[255] == true);

    dynamic_bitset empty;
    empty.resize(31, true);

    REQUIRE(empty[30] == true);
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

    REQUIRE(dynamic_bitset{}.find_first(true) == 0);
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

    REQUIRE(dynamic_bitset{}.find_first(false) == 0);
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

    bitset.resize(64);
    bitset.fill(true);

    REQUIRE(bitset.popcount() == 64);

    dynamic_bitset empty;
    REQUIRE(empty.popcount() == 0);
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

    REQUIRE(!dynamic_bitset{}.any_set());
}

TEST_CASE("all_set", "[dynamic_bitset]")
{
    dynamic_bitset bitset(100, false);
    REQUIRE(!bitset.all_set());

    bitset[99] = true;
    REQUIRE(!bitset.all_set());

    bitset.fill(true);
    REQUIRE(bitset.all_set());

    REQUIRE(dynamic_bitset{}.all_set());
}

TEST_CASE("none_set", "[dynamic_bitset]")
{
    dynamic_bitset bitset(100, false);
    REQUIRE(bitset.none_set());

    bitset[99] = true;
    REQUIRE(!bitset.none_set());

    REQUIRE(dynamic_bitset{}.none_set());
}

TEST_CASE("operator~", "[dynamic_bitset]")
{
    dynamic_bitset bitset1(100, false);
    REQUIRE(bitset1.popcount() == 0);

    dynamic_bitset bitset2 = ~bitset1;
    REQUIRE(bitset2.popcount() == 100);
}

TEST_CASE("flip", "[dynamic_bitset]")
{
    dynamic_bitset bitset(100, false);
    bitset[17] = true;

    REQUIRE(bitset.popcount() == 1);

    bitset.flip();
    REQUIRE(bitset.popcount() == 99);
}

TEST_CASE("operator&", "[dynamic_bitset]")
{
    const dynamic_bitset all_ones(100, true);
    const dynamic_bitset all_zeros(100, false);

    REQUIRE((all_ones & all_zeros) == all_zeros);
    REQUIRE((all_zeros & all_ones) == all_zeros);
    REQUIRE((all_ones & all_ones) == all_ones);
    REQUIRE((all_zeros & all_zeros) == all_zeros);
}

TEST_CASE("operator&=", "[dynamic_bitset]")
{
    dynamic_bitset all_ones(100, true);
    dynamic_bitset all_zeros(100, false);

    all_ones &= all_ones;
    REQUIRE(all_ones.popcount() == 100);

    all_zeros &= all_zeros;
    REQUIRE(all_zeros.popcount() == 0);

    all_zeros &= all_ones;
    REQUIRE(all_zeros.popcount() == 0);

    all_ones &= all_zeros;
    REQUIRE(all_ones.popcount() == 0);
}

TEST_CASE("operator|", "[dynamic_bitset]")
{
    const dynamic_bitset all_ones(100, true);
    const dynamic_bitset all_zeros(100, false);

    REQUIRE((all_ones | all_zeros) == all_ones);
    REQUIRE((all_zeros | all_ones) == all_ones);
    REQUIRE((all_ones | all_ones) == all_ones);
    REQUIRE((all_zeros | all_zeros) == all_zeros);
}

TEST_CASE("operator|=", "[dynamic_bitset]")
{
    dynamic_bitset all_ones(100, true);
    dynamic_bitset all_zeros(100, false);

    all_ones |= all_ones;
    REQUIRE(all_ones.popcount() == 100);

    all_zeros |= all_zeros;
    REQUIRE(all_zeros.popcount() == 0);

    all_ones |= all_zeros;
    REQUIRE(all_ones.popcount() == 100);

    all_zeros |= all_ones;
    REQUIRE(all_zeros.popcount() == 100);
}

TEST_CASE("operator^", "[dynamic_bitset]")
{
    const dynamic_bitset all_ones(100, true);
    const dynamic_bitset all_zeros(100, false);

    REQUIRE((all_ones ^ all_zeros) == all_ones);
    REQUIRE((all_zeros ^ all_ones) == all_ones);
    REQUIRE((all_ones ^ all_ones) == all_zeros);
    REQUIRE((all_zeros ^ all_zeros) == all_zeros);
}

TEST_CASE("operator^=", "[dynamic_bitset]")
{
    dynamic_bitset all_ones(100, true);
    dynamic_bitset all_zeros(100, false);

    all_zeros ^= all_zeros;
    REQUIRE(all_zeros.popcount() == 0);

    all_ones ^= all_zeros;
    REQUIRE(all_ones.popcount() == 100);

    all_ones ^= all_ones;
    REQUIRE(all_ones.popcount() == 0);

    all_ones.fill(true);

    all_zeros ^= all_ones;
    REQUIRE(all_zeros.popcount() == 100);
}

TEST_CASE("operator-", "[dynamic_bitset]")
{
    const dynamic_bitset all_ones(100, true);
    const dynamic_bitset all_zeros(100, false);

    REQUIRE((all_ones - all_zeros) == all_ones);
    REQUIRE((all_zeros - all_ones) == all_zeros);
    REQUIRE((all_ones - all_ones) == all_zeros);
    REQUIRE((all_zeros - all_zeros) == all_zeros);
}

TEST_CASE("operator-=", "[dynamic_bitset]")
{
    dynamic_bitset all_ones(100, true);
    dynamic_bitset all_zeros(100, false);

    all_zeros -= all_zeros;
    REQUIRE(all_zeros.popcount() == 0);

    all_ones -= all_zeros;
    REQUIRE(all_ones.popcount() == 100);

    all_zeros -= all_ones;
    REQUIRE(all_zeros.popcount() == 0);

    all_ones -= all_ones;
    REQUIRE(all_ones.popcount() == 0);
}

TEST_CASE("swap", "[dynamic_bitset]")
{
    dynamic_bitset all_ones(100, true);
    dynamic_bitset all_zeros(50, false);

    all_ones.swap(all_zeros);

    REQUIRE(all_ones.size() == 50);
    REQUIRE(all_ones.popcount() == 0);

    REQUIRE(all_zeros.size() == 100);
    REQUIRE(all_zeros.popcount() == 100);

    all_zeros.swap(all_ones);

    REQUIRE(all_ones.size() == 100);
    REQUIRE(all_ones.popcount() == 100);

    REQUIRE(all_zeros.size() == 50);
    REQUIRE(all_zeros.popcount() == 0);

    using std::swap;
    swap(all_ones, all_zeros);

    REQUIRE(all_ones.size() == 50);
    REQUIRE(all_ones.popcount() == 0);

    REQUIRE(all_zeros.size() == 100);
    REQUIRE(all_zeros.popcount() == 100);

    swap(all_ones, all_zeros);

    REQUIRE(all_ones.size() == 100);
    REQUIRE(all_ones.popcount() == 100);

    REQUIRE(all_zeros.size() == 50);
    REQUIRE(all_zeros.popcount() == 0);
}

TEST_CASE("hash", "[dynamic_bitset]")
{
    const dynamic_bitset bitset(100, true);

    std::hash<dynamic_bitset>{}(bitset);
    SUCCEED();
}

TEST_CASE("reference_copy", "[dynamic_bitset]")
{
    dynamic_bitset bitset(100, false);

    const auto ref = bitset[10];
    const auto ref_copy = ref;

    ref_copy = true;

    REQUIRE(ref == true);
    REQUIRE(bitset[10] == true);
}

TEST_CASE("reference_operator&=", "[dynamic_bitset]")
{
    dynamic_bitset ones(100, true);
    dynamic_bitset zeros(100, false);

    REQUIRE(ones[10] & ones[2]);
    REQUIRE(ones[10] & true);

    REQUIRE((ones[10] &= ones[9]));
    REQUIRE(!(ones[10] &= false));
}

TEST_CASE("reference_operator|=", "[dynamic_bitset]")
{
    dynamic_bitset ones(100, true);
    dynamic_bitset zeros(100, false);

    REQUIRE(ones[10] | zeros[2]);
    REQUIRE(ones[10] | false);

    REQUIRE(!(zeros[10] |= false));
    REQUIRE((zeros[10] |= ones[9]));
}

TEST_CASE("reference_operator^=", "[dynamic_bitset]")
{
    dynamic_bitset ones(100, true);
    dynamic_bitset zeros(100, false);

    REQUIRE(ones[10] ^ zeros[2]);
    REQUIRE(ones[10] ^ false);

    REQUIRE(!(zeros[10] ^= false));
    REQUIRE((zeros[10] ^= ones[9]));
}

TEST_CASE("reference_operator-=", "[dynamic_bitset]")
{
    dynamic_bitset ones(100, true);
    dynamic_bitset zeros(100, false);

    REQUIRE(ones[10] - zeros[2]);
    REQUIRE(ones[10] - false);

    REQUIRE(!(zeros[10] -= true));
    REQUIRE(!(zeros[10] -= false));
    REQUIRE(!(ones[10] -= ones[9]));
}

TEST_CASE("reference_swap", "[dynamic_bitset]")
{
    dynamic_bitset ones(100, true);
    dynamic_bitset zeros(100, false);

    ones[10].swap(zeros[0]);

    REQUIRE(!ones[10]);
    REQUIRE(zeros[0]);

    ones[10].swap(ones[9]);

    REQUIRE(ones[10]);
    REQUIRE(!ones[9]);

    using std::swap;
    swap(ones[9], ones[10]);

    REQUIRE(ones[9]);
    REQUIRE(!ones[10]);

    swap(ones[32], ones[33]);

    REQUIRE(ones[32]);
    REQUIRE(ones[33]);
}

TEST_CASE("push_back", "[dynamic_bitset]")
{
    dynamic_bitset bitset;

    REQUIRE(bitset.empty());

    bitset.push_back(false);
    bitset.push_back(true);
    bitset.push_back(false);

    REQUIRE(bitset.size() == 3);

    REQUIRE(bitset[0] == false);
    REQUIRE(bitset[1] == true);
    REQUIRE(bitset[2] == false);

    bitset.resize(128);
    bitset.push_back(true);

    REQUIRE(bitset[128] == true);
}

TEST_CASE("equality_comparison", "[dynamic_bitset]")
{
    {
        dynamic_bitset bitset1;

        REQUIRE(bitset1 == bitset1);
    }

    {
        dynamic_bitset bitset1(100, true);
        dynamic_bitset bitset2(100, true);

        REQUIRE(bitset1 == bitset2);
    }

    {
        dynamic_bitset bitset1(100, true);
        dynamic_bitset bitset2(100, false);

        REQUIRE(bitset1 != bitset2);
    }

    {
        dynamic_bitset bitset1(100, true);
        dynamic_bitset bitset2(100, true);

        bitset2[99] = false;

        REQUIRE(bitset1 != bitset2);
    }

    {
        dynamic_bitset bitset1(100, true);
        dynamic_bitset bitset2(101, true);

        REQUIRE(bitset1 != bitset2);
    }
}

TEST_CASE("three_way_comparison", "[dynamic_bitset]")
{
    {
        dynamic_bitset bitset1;

        REQUIRE((bitset1 <=> bitset1) == std::strong_ordering::equal);
    }

    {
        dynamic_bitset bitset1(100, true);
        dynamic_bitset bitset2(100, true);

        REQUIRE((bitset1 <=> bitset2) == std::strong_ordering::equal);
        REQUIRE((bitset2 <=> bitset1) == std::strong_ordering::equal);
    }

    {
        dynamic_bitset bitset1(100, true);
        dynamic_bitset bitset2(100, false);

        REQUIRE((bitset1 <=> bitset2) == std::strong_ordering::greater);
        REQUIRE((bitset2 <=> bitset1) == std::strong_ordering::less);
    }

    {
        dynamic_bitset bitset1(100, true);
        dynamic_bitset bitset2(256, false);

        REQUIRE((bitset1 <=> bitset2) == std::strong_ordering::greater);
        REQUIRE((bitset2 <=> bitset1) == std::strong_ordering::less);
    }

    {
        dynamic_bitset bitset1(100, true);
        dynamic_bitset bitset2(100, true);
        bitset2[98] = false;

        REQUIRE((bitset1 <=> bitset2) == std::strong_ordering::greater);
        REQUIRE((bitset2 <=> bitset1) == std::strong_ordering::less);
    }

    {
        dynamic_bitset bitset1(100, true);
        dynamic_bitset bitset2(101, true);
        bitset2[100] = false;

        REQUIRE((bitset1 <=> bitset2) == std::strong_ordering::less);
        REQUIRE((bitset2 <=> bitset1) == std::strong_ordering::greater);
    }

    {
        dynamic_bitset bitset1(100, true);
        dynamic_bitset bitset2(101, true);

        REQUIRE((bitset1 <=> bitset2) == std::strong_ordering::less);
        REQUIRE((bitset2 <=> bitset1) == std::strong_ordering::greater);
    }
}

TEST_CASE("iterators", "[dynamic_bitset]")
{
    dynamic_bitset bitset(100, false);

    std::transform(bitset.begin(), bitset.end(), bitset.begin(),
    [cnt = 0](dynamic_bitset::reference) mutable
    {
        return bool(cnt++ % 2);
    });

    REQUIRE(std::count(bitset.crbegin(), bitset.crend(), true) == 50);
    REQUIRE(bitset.popcount() == 50);
}
