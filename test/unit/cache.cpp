/* Copyright (c) 2024 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "utility/cache.hpp"
#include <utility>
#include <cstddef>

using gapp::detail::fifo_cache;


TEST_CASE("constructor", "[fifo_cache]")
{
    fifo_cache<int, int> cache1;
    REQUIRE(cache1.capacity() == 0);

    fifo_cache<int, int> cache2(4);
    REQUIRE(cache2.capacity() == 4);

    fifo_cache<int, int> cache3 = cache2;
    REQUIRE(cache3.capacity() == 4);

    fifo_cache<int, int> cache4 = std::move(cache2);
    REQUIRE(cache4.capacity() == 4);
    REQUIRE(cache2.empty());
}

TEST_CASE("copy_complex", "[fifo_cache]")
{
    fifo_cache<int, int> cache1(4);

    cache1.insert(1, 2);
    cache1.insert(3, 6);
    cache1.insert(2, 4);
    cache1.insert(4, 8);

    fifo_cache<int, int> cache2 = cache1;

    REQUIRE(cache2.size() == 4);
    REQUIRE(cache2.capacity() == 4);

    cache2.insert(5, 10);

    REQUIRE(cache2.get(1) == nullptr);
    REQUIRE(*cache2.get(5) == 10);

    cache2.insert(6, 12);

    REQUIRE(cache2.get(3) == nullptr);
    REQUIRE(*cache2.get(6) == 12);
}

TEST_CASE("size/capacity", "[fifo_cache]")
{
    fifo_cache<int, int> cache(5);

    REQUIRE(cache.size() == 0);
    REQUIRE(cache.capacity() == 5);

    cache.insert(1, 2);

    REQUIRE(cache.size() == 1);
    REQUIRE(cache.capacity() == 5);

    cache.insert(2, 4);
    cache.insert(3, 6);
    cache.insert(4, 8);
    cache.insert(5, 10);

    REQUIRE(cache.size() == 5);
    REQUIRE(cache.capacity() == 5);

    cache.insert(6, 12);

    REQUIRE(cache.size() == 5);
    REQUIRE(cache.capacity() == 5);
}

TEST_CASE("full/empty", "[fifo_cache]")
{
    fifo_cache<int, int> cache(4);

    REQUIRE(cache.empty());
    REQUIRE(!cache.full());

    cache.insert(1, 2);

    REQUIRE(!cache.empty());
    REQUIRE(!cache.full());

    cache.insert(2, 4);
    cache.insert(3, 6);
    cache.insert(4, 8);

    REQUIRE(!cache.empty());
    REQUIRE(cache.full());

    cache.insert(5, 10);

    REQUIRE(!cache.empty());
    REQUIRE(cache.full());
}

TEST_CASE("insert/get", "[fifo_cache]")
{
    fifo_cache<int, int> cache(4);

    REQUIRE(cache.get(1) == nullptr);

    cache.insert(1, 2);

    REQUIRE(*cache.get(1) == 2);
    REQUIRE(cache.get(2) == nullptr);

    cache.insert(1, -1);

    REQUIRE(*cache.get(1) == -1);

    cache.insert(2, 4);
    cache.insert(3, 6);
    cache.insert(4, 8);

    REQUIRE(*cache.get(1) == -1);
    REQUIRE(*cache.get(3) == 6);
    REQUIRE(*cache.get(4) == 8);

    cache.insert(5, 10);

    REQUIRE(*cache.get(5) == 10);
    REQUIRE(cache.get(1) == nullptr);

    cache.insert(6, 12);

    REQUIRE(*cache.get(6) == 12);
    REQUIRE(cache.get(2) == nullptr);


    fifo_cache<int, int> empty;

    REQUIRE(empty.empty());
    REQUIRE(empty.capacity() == 0);

    empty.insert(1, 1);

    REQUIRE(empty.empty());
}

TEST_CASE("try_insert", "[fifo_cache]")
{
    fifo_cache<int, int> cache(4);

    REQUIRE(cache.get(1) == nullptr);

    cache.try_insert(1, 2);

    REQUIRE(*cache.get(1) == 2);
    REQUIRE(cache.get(2) == nullptr);

    cache.try_insert(1, -1);

    REQUIRE(*cache.get(1) == 2);

    cache.try_insert(2, 4);
    cache.try_insert(3, 6);
    cache.try_insert(4, 8);

    REQUIRE(*cache.get(1) == 2);
    REQUIRE(*cache.get(3) == 6);
    REQUIRE(*cache.get(4) == 8);

    cache.try_insert(5, 10);

    REQUIRE(*cache.get(5) == 10);
    REQUIRE(cache.get(1) == nullptr);

    cache.try_insert(6, 12);

    REQUIRE(*cache.get(6) == 12);
    REQUIRE(cache.get(2) == nullptr);


    fifo_cache<int, int> empty;

    REQUIRE(empty.empty());
    REQUIRE(empty.capacity() == 0);

    empty.try_insert(1, 1);

    REQUIRE(empty.empty());
}

TEST_CASE("insert_range", "[fifo_cache]")
{
    const std::vector keys{ 1, 2, 3, 4 };

    fifo_cache<int, int> cache1(4);
    cache1.insert(keys.begin(), keys.end(), [](int n) { return n * 2; });

    REQUIRE(cache1.size() == 4);
    REQUIRE(*cache1.get(1) == 2);
    REQUIRE(*cache1.get(3) == 6);

    fifo_cache<int, int> cache2(2);
    cache2.insert(keys.begin(), keys.end(), [](int n) { return n * 2; });

    REQUIRE(cache2.size() == 2);
    REQUIRE(*cache2.get(3) == 6);
    REQUIRE(*cache2.get(4) == 8);
}

TEST_CASE("contains", "[fifo_cache]")
{
    fifo_cache<int, int> cache(4);

    REQUIRE(!cache.contains(3));
    REQUIRE(!cache.contains(2));

    cache.insert(3, 2);

    REQUIRE(cache.contains(3));
    REQUIRE(!cache.contains(2));
}

TEST_CASE("clear", "[fifo_cache]")
{
    fifo_cache<int, int> cache(3);

    cache.insert(1, 2);
    cache.insert(2, 4);

    REQUIRE(cache.size() == 2);
    REQUIRE(cache.capacity() == 3);

    cache.clear();

    REQUIRE(cache.empty());
    REQUIRE(cache.capacity() == 3);
}

TEST_CASE("reset", "[fifo_cache]")
{
    const size_t new_capacity = GENERATE(2, 3, 5);

    fifo_cache<int, int> cache(3);

    cache.insert(1, 2);
    cache.insert(2, 4);

    REQUIRE(cache.size() == 2);
    REQUIRE(cache.capacity() == 3);

    cache.reset(new_capacity);

    REQUIRE(cache.empty());
    REQUIRE(cache.capacity() == new_capacity);
}

TEST_CASE("swap", "[fifo_cache]")
{
    using std::swap;

    fifo_cache<int, int> cache1(4);
    fifo_cache<int, int> cache2(5);

    REQUIRE(cache1 == cache2);

    swap(cache1, cache2);

    REQUIRE(cache1 == cache2);

    cache1.insert(1, 2);
    cache1.insert(2, 4);

    cache2.insert(1, 3);

    REQUIRE(cache1 != cache2);

    REQUIRE(cache1.size() == 2);
    REQUIRE(cache1.capacity() == 5);

    REQUIRE(cache2.size() == 1);
    REQUIRE(cache2.capacity() == 4);

    swap(cache1, cache2);

    REQUIRE(cache1.size() == 1);
    REQUIRE(cache1.capacity() == 4);

    REQUIRE(cache2.size() == 2);
    REQUIRE(cache2.capacity() == 5);

    REQUIRE(*cache1.get(1) == 3);
}

TEST_CASE("comparison", "[fifo_cache]")
{
    fifo_cache<int, int> cache1(3);
    fifo_cache<int, int> cache2(4);

    REQUIRE(cache1 == cache2);

    cache1.insert(1, 2);
    cache2.insert(1, 2);

    REQUIRE(cache1 == cache2);

    cache1.insert(2, 4);
    cache2.insert(3, 6);

    REQUIRE(cache1 != cache2);

    cache1.insert(3, 6);
    cache2.insert(2, 4);

    REQUIRE(cache1 != cache2);
}
