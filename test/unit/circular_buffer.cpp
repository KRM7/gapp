/* Copyright (c) 2024 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "utility/circular_buffer.hpp"
#include <algorithm>
#include <utility>
#include <cstddef>

using gapp::detail::circular_buffer;


TEST_CASE("constructor", "[circular_buffer]")
{
    const size_t capacity = GENERATE(1, 2, 3, 5, 100);

    circular_buffer<int> buffer(capacity);

    REQUIRE(buffer.capacity() == capacity);
    REQUIRE(buffer.size() == 0);

    circular_buffer<int> buffer_copy = buffer;

    REQUIRE(buffer.capacity() == capacity);
    REQUIRE(buffer.size() == 0);

    REQUIRE(buffer_copy.capacity() == capacity);
    REQUIRE(buffer_copy.size() == 0);

    circular_buffer<int> buffer_moved = std::move(buffer);

    REQUIRE(buffer.capacity() == 0);
    REQUIRE(buffer.size() == 0);

    REQUIRE(buffer_moved.capacity() == capacity);
    REQUIRE(buffer_moved.size() == 0);
}

TEST_CASE("emplace_back", "[circular_buffer]")
{
    circular_buffer<int> buffer(4);

    REQUIRE(buffer.capacity() == 4);

    buffer.emplace_back(1);
    REQUIRE(buffer.size() == 1);

    buffer.emplace_back(2);
    buffer.emplace_back(3);
    REQUIRE(buffer.size() == 3);

    buffer.emplace_back(4);
    REQUIRE(buffer.size() == 4);

    buffer.emplace_back(5);
    REQUIRE(buffer.size() == 4);

    buffer.emplace_back(6);
    REQUIRE(buffer.size() == 4);

    REQUIRE(buffer.capacity() == 4);
}

TEST_CASE("empty/full", "[circular_buffer]")
{
    circular_buffer<int> buffer(4);

    REQUIRE(buffer.empty());
    REQUIRE(!buffer.full());

    buffer.emplace_back(1);
    REQUIRE(!buffer.empty());
    REQUIRE(!buffer.full());

    buffer.emplace_back(2);
    buffer.emplace_back(3);
    buffer.emplace_back(4);
    REQUIRE(!buffer.empty());
    REQUIRE(buffer.full());

    buffer.emplace_back(5);
    REQUIRE(!buffer.empty());
    REQUIRE(buffer.full());
}

TEST_CASE("front/back", "[circular_buffer]")
{
    circular_buffer<int> buffer(4);

    buffer.emplace_back(1);
    REQUIRE(buffer.back() == 1);
    REQUIRE(buffer.front() == 1);

    buffer.emplace_back(2);
    REQUIRE(buffer.back() == 2);
    REQUIRE(buffer.front() == 1);

    buffer.emplace_back(3);
    buffer.emplace_back(4);
    REQUIRE(buffer.back() == 4);
    REQUIRE(buffer.front() == 1);

    buffer.emplace_back(5);
    REQUIRE(buffer.back() == 5);
    REQUIRE(buffer.front() == 2);

    buffer.emplace_back(6);
    REQUIRE(buffer.back() == 6);
    REQUIRE(buffer.front() == 3);
}

TEST_CASE("element_access", "[circular_buffer]")
{
    circular_buffer<int> buffer(4);

    buffer.emplace_back(1);
    buffer.emplace_back(2);
    buffer.emplace_back(3);
    buffer.emplace_back(4);

    REQUIRE(buffer[0] == 1);
    REQUIRE(buffer[1] == 2);
    REQUIRE(buffer[2] == 3);
    REQUIRE(buffer[3] == 4);

    buffer.emplace_back(5);
    buffer.emplace_back(6);

    REQUIRE(buffer[0] == 3);
    REQUIRE(buffer[1] == 4);
    REQUIRE(buffer[2] == 5);
    REQUIRE(buffer[3] == 6);

    buffer.emplace_back(7);
    buffer.emplace_back(8);
    buffer.emplace_back(9);

    REQUIRE(buffer[0] == 6);
    REQUIRE(buffer[1] == 7);
    REQUIRE(buffer[2] == 8);
    REQUIRE(buffer[3] == 9);
}

TEST_CASE("emplace_front", "[circular_buffer]")
{
    circular_buffer<int> buffer(4);

    buffer.emplace_front(1);
    REQUIRE(buffer.size() == 1);
    REQUIRE(buffer.front() == 1);

    buffer.emplace_front(2);
    REQUIRE(buffer.size() == 2);
    REQUIRE(buffer.front() == 2);

    buffer.emplace_front(3);
    buffer.emplace_front(4);
    REQUIRE(buffer.size() == 4);
    REQUIRE(buffer.front() == 4);
    REQUIRE(buffer.back() == 1);

    buffer.emplace_front(5);
    REQUIRE(buffer.size() == 4);
    REQUIRE(buffer.front() == 5);
    REQUIRE(buffer.back() == 2);
}

TEST_CASE("pop_front/back", "[circular_buffer]")
{
    circular_buffer<int> buffer(4);

    buffer.emplace_back(1);
    buffer.emplace_back(2);

    SECTION("pop_front")
    {
        buffer.pop_front();
        REQUIRE(buffer.size() == 1);
        REQUIRE(buffer.front() == 2);
        REQUIRE(buffer.back() == 2);

        buffer.pop_front();
        REQUIRE(buffer.empty());
    }

    SECTION("pop_back")
    {
        buffer.pop_back();
        REQUIRE(buffer.size() == 1);
        REQUIRE(buffer.front() == 1);
        REQUIRE(buffer.back() == 1);

        buffer.pop_back();
        REQUIRE(buffer.empty());
    }

    SECTION("pop_both")
    {
        buffer.pop_front();
        REQUIRE(buffer.size() == 1);
        REQUIRE(buffer.front() == 2);
        REQUIRE(buffer.back() == 2);

        buffer.pop_back();
        REQUIRE(buffer.empty());
    }
}

TEST_CASE("push_pop", "[circular_buffer]")
{
    circular_buffer<int> buffer(4);

    buffer.push_back(1);
    buffer.push_back(2);
    buffer.pop_front();

    REQUIRE(buffer.size() == 1);
    REQUIRE(buffer.front() == 2);
    REQUIRE(buffer.back() == 2);

    buffer.push_back(3);

    REQUIRE(buffer.size() == 2);
    REQUIRE(buffer.back() == 3);

    buffer.push_back(4);
    buffer.pop_front();

    REQUIRE(buffer.size() == 2);
    REQUIRE(buffer.back() == 4);

    buffer.push_back(5);
    buffer.push_back(6);

    REQUIRE(buffer.size() == 4);
    REQUIRE(buffer.front() == 3);
    REQUIRE(buffer.back() == 6);
}

TEST_CASE("set_capacity", "[circular_buffer]")
{
    circular_buffer<int> buffer(4, { 1, 2, 3 });

    REQUIRE(buffer.size() == 3);
    REQUIRE(buffer.capacity() == 4);

    buffer.set_capacity(4);

    REQUIRE(buffer.size() == 3);
    REQUIRE(buffer.capacity() == 4);

    REQUIRE(buffer.front() == 1);
    REQUIRE(buffer.back() == 3);

    buffer.set_capacity(10);

    REQUIRE(buffer.size() == 3);
    REQUIRE(buffer.capacity() == 10);

    REQUIRE(buffer.front() == 1);
    REQUIRE(buffer.back() == 3);

    buffer.set_capacity(1);

    REQUIRE(buffer.size() == 1);
    REQUIRE(buffer.capacity() == 1);

    REQUIRE(buffer.front() == 1);
}

TEST_CASE("clear", "[circular_buffer]")
{
    circular_buffer<int> buffer(4);

    buffer.emplace_back(1);
    buffer.emplace_back(2);
    buffer.emplace_back(3);

    buffer.clear();

    REQUIRE(buffer.empty());
    REQUIRE(buffer.capacity() == 4);

    buffer.emplace_front(1);

    REQUIRE(buffer.size() == 1);
}

TEST_CASE("reset", "[circular_buffer]")
{
    const size_t new_capacity = GENERATE(1, 3, 5);

    circular_buffer<int> buffer(4, { 1, 2, 3 });

    REQUIRE(buffer.size() == 3);
    REQUIRE(buffer.capacity() == 4);

    buffer.reset(new_capacity);

    REQUIRE(buffer.empty());
    REQUIRE(buffer.capacity() == new_capacity);
}

TEST_CASE("comparisons", "[circular_buffer]")
{
    circular_buffer<int> buffer1(4);

    REQUIRE(buffer1 == circular_buffer<int>{});

    buffer1.push_back(1);
    buffer1.push_back(2);
    buffer1.push_back(3);
    buffer1.push_back(4);

    circular_buffer<int> buffer2(4);

    buffer2.push_back(0);
    buffer2.push_back(1);
    buffer2.push_back(2);
    buffer2.push_back(3);

    REQUIRE(buffer1 != buffer2);

    buffer2.push_back(4);

    REQUIRE(buffer1 == buffer2);

    REQUIRE(buffer1 != circular_buffer<int>{});
}

TEST_CASE("iterators", "[circular_buffer]")
{
    circular_buffer<int> buffer(4);

    buffer.push_back(0);
    buffer.push_back(2);
    buffer.push_back(1);
    buffer.push_back(3);
    buffer.push_back(0);

    REQUIRE(!std::is_sorted(buffer.cbegin(), buffer.cend()));

    std::sort(buffer.begin(), buffer.end());

    REQUIRE(std::is_sorted(buffer.rbegin(), buffer.rend(), std::greater{}));
}

TEST_CASE("swap", "[circular_buffer]")
{
    circular_buffer<int> buffer1(4, { 1, 2, 3 });
    circular_buffer<int> buffer2(5, { 4, 5 });

    REQUIRE(buffer1 != buffer2);

    REQUIRE(buffer1.size() == 3);
    REQUIRE(buffer1.capacity() == 4);

    REQUIRE(buffer2.size() == 2);
    REQUIRE(buffer2.capacity() == 5);

    using std::swap;
    swap(buffer1, buffer2);

    REQUIRE(buffer1.size() == 2);
    REQUIRE(buffer1.capacity() == 5);

    REQUIRE(buffer2.size() == 3);
    REQUIRE(buffer2.capacity() == 4);

    REQUIRE(buffer1.front() == 4);
    REQUIRE(buffer2.front() == 1);
}
