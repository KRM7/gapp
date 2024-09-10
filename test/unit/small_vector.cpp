/* Copyright (c) 2024 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "utility/small_vector.hpp"
#include <algorithm>
#include <vector>
#include <iterator>
#include <memory>
#include <type_traits>
#include <string>
#include <sstream>
#include <utility>
#include <cstddef>

using namespace gapp;

inline constexpr size_t EMPTY = 0;
inline constexpr size_t SMALL_SIZE = 4;
inline constexpr size_t LARGE_SIZE = 100;

using TrivialType = int;

struct NonDefaultConstructibleType
{
    NonDefaultConstructibleType(int i) : i_(i) {}

    friend bool operator==(const NonDefaultConstructibleType&, const NonDefaultConstructibleType&) = default;

    int i_ = 0;
};

struct MoveOnlyType
{
    MoveOnlyType() = default;
    MoveOnlyType(int i) : i_(i) {}

    MoveOnlyType(const MoveOnlyType&) = delete;
    MoveOnlyType(MoveOnlyType&& o) : i_(o.i_) {}

    MoveOnlyType& operator=(const MoveOnlyType&) = delete;
    MoveOnlyType& operator=(MoveOnlyType&& o) { i_ = o.i_; return *this; }

    ~MoveOnlyType() noexcept {}

    friend bool operator==(const MoveOnlyType&, const MoveOnlyType&) = default;

    int i_ = 0;
};

struct ImmovableType
{
    ImmovableType() = default;
    ImmovableType(int i) : i_(i) {}

    ImmovableType(const ImmovableType&) = delete;
    ImmovableType(ImmovableType&&) = delete;

    ImmovableType& operator=(const ImmovableType&) = delete;
    ImmovableType& operator=(ImmovableType&&) = delete;

    friend bool operator==(const ImmovableType&, const ImmovableType&) = default;

    int i_ = 0;
};

struct NonTrivialType
{
    NonTrivialType() {}
    NonTrivialType(int i) : i_(i) {}
    NonTrivialType(const NonTrivialType& o) : i_(o.i_) {}
    NonTrivialType(NonTrivialType&& o) : i_(o.i_) {}
    NonTrivialType& operator=(const NonTrivialType& o) { i_ = o.i_; return *this; }
    NonTrivialType& operator=(NonTrivialType&& o) { i_ = o.i_; return *this; }
    ~NonTrivialType() noexcept {}

    friend bool operator==(const NonTrivialType&, const NonTrivialType&) = default;

    int i_ = 0;
};

template<typename T>
constexpr auto equal_to(T rhs)
{
    return [rhs = std::move(rhs)](const auto& lhs) { return lhs == rhs; };
}


    //-----------------------------------//
    //           OBJECT LAYOUT           //
    //-----------------------------------//


TEST_CASE("small_vector_size", "[object_layout][!mayfail]") // fails under clang-cl because of no support for no_unique_address
{
    STATIC_REQUIRE(std::is_standard_layout_v<small_vector<int>>);

    CHECK(sizeof(small_vector<int>) == detail::cache_line_size);
}


    //-----------------------------------//
    //            CONSTRUCTORS           //
    //-----------------------------------//


TEMPLATE_TEST_CASE("small_vector()", "[constructor]", TrivialType, MoveOnlyType, ImmovableType, NonTrivialType)
{
    small_vector<TestType> vec;

    REQUIRE(vec.empty());
    REQUIRE(vec.begin() == vec.end());
    REQUIRE(vec.cbegin() == vec.end());

    REQUIRE(vec.size() == 0);
    REQUIRE(vec.capacity() > 0);
}

TEMPLATE_TEST_CASE("small_vector(Alloc)", "[constructor]", TrivialType, MoveOnlyType, ImmovableType, NonTrivialType)
{
    small_vector<TestType> vec(std::allocator<TestType>{});

    REQUIRE(vec.empty());
    REQUIRE(vec.begin() == vec.end());
    REQUIRE(vec.cbegin() == vec.end());

    REQUIRE(vec.size() == 0);
    REQUIRE(vec.capacity() > 0);
}

TEMPLATE_TEST_CASE("small_vector(size_t)", "[constructor]", TrivialType, MoveOnlyType, ImmovableType, NonTrivialType)
{
    const size_t size = GENERATE(EMPTY, SMALL_SIZE, LARGE_SIZE);

    small_vector<TestType> vec(size);

    REQUIRE(vec.size() == size);
    REQUIRE(vec.capacity() >= size);
}

TEMPLATE_TEST_CASE("small_vector(size_t, Alloc)", "[constructor]", TrivialType, MoveOnlyType, ImmovableType, NonTrivialType)
{
    const size_t size = GENERATE(EMPTY, SMALL_SIZE, LARGE_SIZE);

    small_vector<TestType> vec(size, std::allocator<TestType>{});

    REQUIRE(vec.size() == size);
    REQUIRE(vec.capacity() >= size);
}

TEMPLATE_TEST_CASE("small_vector(size_t, const T&)", "[constructor]", TrivialType, NonTrivialType, NonDefaultConstructibleType)
{
    const size_t size = GENERATE(EMPTY, SMALL_SIZE, LARGE_SIZE);

    small_vector vec(size, TestType{ 0 });

    REQUIRE(vec.size() == size);
    REQUIRE(vec.capacity() >= size);
    REQUIRE(std::all_of(vec.begin(), vec.end(), equal_to<TestType>(0)));
}

TEMPLATE_TEST_CASE("small_vector(size_t, const T&, Alloc)", "[constructor]", TrivialType, NonTrivialType, NonDefaultConstructibleType)
{
    const size_t size = GENERATE(EMPTY, SMALL_SIZE, LARGE_SIZE);

    small_vector vec(size, TestType{ 1 }, std::allocator<TestType>{});

    REQUIRE(vec.size() == size);
    REQUIRE(vec.capacity() >= size);
    REQUIRE(std::all_of(vec.begin(), vec.end(), equal_to<TestType>(1)));
}

TEMPLATE_TEST_CASE("small_vector(Iter, Iter)", "[constructor]", TrivialType, MoveOnlyType, ImmovableType, NonTrivialType, NonDefaultConstructibleType)
{
    const size_t size = GENERATE(EMPTY, SMALL_SIZE, LARGE_SIZE);

    const std::vector<int> source(size, 2);

    small_vector<TestType> vec(source.begin(), source.end());

    REQUIRE(vec.size() == source.size());
    REQUIRE(vec.capacity() >= source.size());
}

TEMPLATE_TEST_CASE("small_vector(FwdIter, FwdIter, Alloc)", "[constructor]", TrivialType, MoveOnlyType, ImmovableType, NonTrivialType, NonDefaultConstructibleType)
{
    const size_t size = GENERATE(EMPTY, SMALL_SIZE, LARGE_SIZE);

    const std::vector<int> source(size, 3);

    small_vector<TestType> vec(source.begin(), source.end(), std::allocator<TestType>{});

    REQUIRE(vec.size() == source.size());
    REQUIRE(vec.capacity() >= source.size());
}

TEST_CASE("small_vector(InputIter, InputIter)", "[constructor]")
{
    const size_t size = GENERATE(EMPTY, SMALL_SIZE, LARGE_SIZE);

    std::istringstream source(std::string(size, 'c'));

    small_vector<char> vec{ std::istream_iterator<char>(source), std::istream_iterator<char>() };

    REQUIRE(vec.size() == size);
    REQUIRE(vec.capacity() >= size);
    REQUIRE(std::all_of(vec.begin(), vec.end(), equal_to('c')));
}

TEST_CASE("small_vector(nullptr, nullptr)", "[constructor]")
{
    const small_vector<int> vec(static_cast<int*>(nullptr), static_cast<int*>(nullptr));

    REQUIRE(vec.empty());
    REQUIRE(vec.capacity() != 0);
}

TEMPLATE_TEST_CASE("small_vector(initializer_list)", "[constructor]", TrivialType, NonTrivialType, NonDefaultConstructibleType)
{
    small_vector vec{ TestType{ 1 }, TestType{ 4 }, TestType{ 2 } };

    REQUIRE(vec.is_small());
    REQUIRE(vec.size() == 3);
    REQUIRE(vec.capacity() >= 3);
    REQUIRE(vec.back() == 2);
}

TEMPLATE_TEST_CASE("small_vector(initializer_list, Alloc)", "[constructor]", TrivialType, NonTrivialType, NonDefaultConstructibleType)
{
    small_vector vec({ TestType{ 1 }, TestType{ 4 }, TestType{ 2 } }, std::allocator<TestType>{});

    REQUIRE(vec.is_small());
    REQUIRE(vec.size() == 3);
    REQUIRE(vec.capacity() >= 3);
    REQUIRE(vec.back() == 2);
}

TEMPLATE_TEST_CASE("small_vector(const small_vector&)", "[constructor]", TrivialType, NonTrivialType, NonDefaultConstructibleType)
{
    const size_t size = GENERATE(EMPTY, SMALL_SIZE, LARGE_SIZE);

    const small_vector source(size, TestType{ 26 });

    small_vector vec(source);

    REQUIRE(vec.size() == source.size());
    REQUIRE(vec == source);
    REQUIRE(vec.capacity() != 0);
}

TEMPLATE_TEST_CASE("small_vector(const small_vector&, Alloc)", "[constructor]", TrivialType, NonTrivialType, NonDefaultConstructibleType)
{
    const size_t size = GENERATE(EMPTY, SMALL_SIZE, LARGE_SIZE);

    const small_vector source(size, TestType{ 26 });

    small_vector vec(source, std::allocator<TestType>{});

    REQUIRE(vec.size() == source.size());
    REQUIRE(vec == source);
    REQUIRE(vec.capacity() != 0);
}

TEMPLATE_TEST_CASE("small_vector(small_vector&&)", "[constructor]", TrivialType, NonTrivialType, NonDefaultConstructibleType)
{
    const size_t size = GENERATE(EMPTY, SMALL_SIZE, LARGE_SIZE);

    small_vector source(size, TestType{ 26 });
    small_vector source_copy(source);

    small_vector vec(std::move(source));

    REQUIRE(vec.size() == source_copy.size());
    REQUIRE(vec == source_copy);

    REQUIRE(source.empty());
    REQUIRE(source.capacity() != 0);

    source.push_back(TestType{ 11 });
    REQUIRE(source.size() == 1);
}

TEMPLATE_TEST_CASE("small_vector<MoveOnlyType>(small_vector&&)", "[constructor]", MoveOnlyType)
{
    const size_t size = GENERATE(EMPTY, SMALL_SIZE, LARGE_SIZE);

    small_vector<TestType> source(size);
    small_vector<TestType> source_copy(size);

    small_vector vec(std::move(source));

    REQUIRE(vec.size() == source_copy.size());
    REQUIRE(vec == source_copy);

    REQUIRE(source.empty());
    REQUIRE(source.capacity() != 0);

    source.push_back(TestType{ 11 });
    REQUIRE(source.size() == 1);
}

    //-----------------------------------//
    //             ASSIGNMENT            //
    //-----------------------------------//

TEMPLATE_TEST_CASE("assign(size_t count, const T& val)", "[assignment]", TrivialType, NonTrivialType, NonDefaultConstructibleType)
{
    const size_t count = GENERATE(EMPTY, SMALL_SIZE - 1, LARGE_SIZE + 1);
    const size_t dest_size = GENERATE(EMPTY, SMALL_SIZE, LARGE_SIZE);

    small_vector dest(dest_size, TestType{ 2 });

    dest.assign(count, TestType{ 3 });

    REQUIRE(dest.size() == count);
    REQUIRE(std::all_of(dest.begin(), dest.end(), equal_to<TestType>(3)));
}

TEMPLATE_TEST_CASE("assign(FwdIter, FwdIter)", "[assignment]", TrivialType, NonTrivialType, NonDefaultConstructibleType)
{
    const size_t src_size = GENERATE(EMPTY, SMALL_SIZE - 1, LARGE_SIZE + 1);
    const size_t dest_size = GENERATE(EMPTY, SMALL_SIZE, LARGE_SIZE);

    const small_vector source(src_size, TestType{ 4 });
    small_vector dest(dest_size, TestType{ 3 });

    dest.assign(source.begin(), source.end());

    REQUIRE(dest.size() == source.size());
    REQUIRE(dest == source);


    small_vector dest2 = dest;
    dest2.reserve(2 * source.size());
    dest2.assign(source.begin(), source.end());

    REQUIRE(dest.size() == source.size());
    REQUIRE(dest == source);
}

TEST_CASE("assign(InputIter, InputIter)", "[assignment]")
{
    const size_t src_size = GENERATE(EMPTY, SMALL_SIZE - 1, LARGE_SIZE + 1);
    const size_t dst_size = GENERATE(EMPTY, SMALL_SIZE, LARGE_SIZE);

    std::istringstream source(std::string(src_size, 'c'));
    small_vector<char> dest(dst_size);

    dest.assign(std::istream_iterator<char>(source), std::istream_iterator<char>());

    REQUIRE(dest.size() == src_size);
    REQUIRE(dest.capacity() >= src_size);
    REQUIRE(std::all_of(dest.begin(), dest.end(), equal_to('c')));
}

TEST_CASE("assign(nullptr, nullptr)", "[assignment]")
{
    small_vector<int> vec(0, 2);
    vec.assign(static_cast<int*>(nullptr), static_cast<int*>(nullptr));

    REQUIRE(vec.empty());
}

TEMPLATE_TEST_CASE("operator=(const small_vector&)", "[assignment]", TrivialType, NonTrivialType, NonDefaultConstructibleType)
{
    const size_t src_size = GENERATE(EMPTY, SMALL_SIZE - 1, LARGE_SIZE + 1);
    const size_t dest_size = GENERATE(EMPTY, SMALL_SIZE, LARGE_SIZE);

    const small_vector source(src_size, TestType{ 4 });
    small_vector dest(dest_size, TestType{ 3 });

    dest = source;

    REQUIRE(dest.size() == source.size());
    REQUIRE(dest == source);


    small_vector dest2 = dest;
    dest2.reserve(2 * source.size());
    dest2 = source;

    REQUIRE(dest.size() == source.size());
    REQUIRE(dest == source);
}

TEMPLATE_TEST_CASE("operator=(small_vector&&)", "[assignment]", TrivialType, NonTrivialType, NonDefaultConstructibleType)
{
    const size_t src_size = GENERATE(EMPTY, SMALL_SIZE - 1, LARGE_SIZE + 1);
    const size_t dest_size = GENERATE(EMPTY, SMALL_SIZE, LARGE_SIZE);

    small_vector source(src_size, TestType{ 4 });
    const small_vector src_copy(source);
    small_vector dest(dest_size, TestType{ 3 });

    dest = std::move(source);

    REQUIRE(dest.size() == src_copy.size());
    REQUIRE(dest == src_copy);

    REQUIRE(source.capacity() != 0);

    source = dest;
    REQUIRE(source == dest);
}

TEMPLATE_TEST_CASE("operator=(initializer_list)", "[assignment]", TrivialType, NonTrivialType, NonDefaultConstructibleType)
{
    const size_t dest_size = GENERATE(SMALL_SIZE, LARGE_SIZE);

    const small_vector source{ TestType{ 1 }, TestType{ 2 }, TestType{ 4 } };
    small_vector dest(dest_size, TestType{ 0 });

    dest = { TestType{ 1 }, TestType{ 2 }, TestType{ 4 } };

    REQUIRE(dest.size() == source.size());
    REQUIRE(dest == source);
}

    //-----------------------------------//
    //             ITERATORS             //
    //-----------------------------------//

TEMPLATE_TEST_CASE("forward_iteration", "[iterators]", TrivialType, NonTrivialType)
{
    const size_t size = GENERATE(SMALL_SIZE, LARGE_SIZE);
    const small_vector vec(size, TestType{ 1 });

    REQUIRE(std::all_of(vec.begin(), vec.end(), equal_to(TestType{ 1 })));
}

TEMPLATE_TEST_CASE("reverse_iteration", "[iterators]", TrivialType, NonTrivialType)
{
    const size_t size = GENERATE(SMALL_SIZE, LARGE_SIZE);
    const small_vector vec(size, TestType{ 1 });

    REQUIRE(std::all_of(vec.rbegin(), vec.rend(), equal_to(TestType{ 1 })));
}

    //-----------------------------------//
    //           ELEMENT ACCESS          //
    //-----------------------------------//

TEMPLATE_TEST_CASE("operator[]", "[element_access]", TrivialType, NonTrivialType, MoveOnlyType, ImmovableType)
{
    const size_t size = GENERATE(SMALL_SIZE, LARGE_SIZE);
    const small_vector<TestType> vec(size);

    REQUIRE(vec[0] == vec[1]);
    REQUIRE(vec[1] == TestType{});
}

TEMPLATE_TEST_CASE("at", "[element_access]", TrivialType, NonTrivialType, MoveOnlyType, ImmovableType)
{
    const size_t size = GENERATE(SMALL_SIZE, LARGE_SIZE);
    const small_vector<TestType> vec(size);

    REQUIRE(vec.at(0) == vec.at(1));
    REQUIRE(vec.at(1) == TestType{});
    REQUIRE_THROWS(vec.at(size));
}

TEMPLATE_TEST_CASE("front/back", "[element_access]", TrivialType, NonTrivialType)
{
    const size_t size = GENERATE(SMALL_SIZE, LARGE_SIZE);
    small_vector vec(size, TestType{ 2 });

    vec.back() = TestType{ 3 };
    vec.front() = TestType{ 0 };

    REQUIRE(vec.front() == TestType{ 0 });
    REQUIRE(vec.back() == TestType{ 3 });
}

    //-----------------------------------//
    //              CAPACITY             //
    //-----------------------------------//


TEMPLATE_TEST_CASE("is_small", "[capacity]", TrivialType, NonTrivialType)
{
    const size_t size = GENERATE(SMALL_SIZE, LARGE_SIZE);
    small_vector vec(size, TestType{ 2 });

    REQUIRE(vec.is_small() == (size == SMALL_SIZE));
}

TEMPLATE_TEST_CASE("reserve", "[capacity]", TrivialType, NonTrivialType)
{
    const size_t size = GENERATE(SMALL_SIZE, LARGE_SIZE);
    small_vector vec(size, TestType{ 2 });

    vec.reserve(2 * LARGE_SIZE);
    REQUIRE(vec.capacity() >= 2 * LARGE_SIZE);

    vec.reserve(SMALL_SIZE);
    REQUIRE(vec.capacity() >= 2 * LARGE_SIZE);
}

TEMPLATE_TEST_CASE("shrink_to_fit", "[capacity]", TrivialType, NonTrivialType)
{
    small_vector vec(LARGE_SIZE, TestType{ 2 });

    vec.reserve(2 * LARGE_SIZE);
    REQUIRE(vec.capacity() >= 2 * LARGE_SIZE);

    vec.shrink_to_fit();
    REQUIRE(vec.capacity() == LARGE_SIZE);
}

    //-----------------------------------//
    //             MODIFIERS             //
    //-----------------------------------//

TEMPLATE_TEST_CASE("clear", "[modifiers]", TrivialType, NonTrivialType)
{
    const size_t size = GENERATE(SMALL_SIZE, LARGE_SIZE);
    small_vector vec(size, TestType{ 1 });

    vec.clear();

    REQUIRE(vec.empty());
}

TEMPLATE_TEST_CASE("swap", "[modifiers]", TrivialType, NonTrivialType)
{
    const size_t left_size = GENERATE(SMALL_SIZE, LARGE_SIZE);
    const size_t right_size = GENERATE(SMALL_SIZE, LARGE_SIZE);

    small_vector left(left_size, TestType{ 1 });
    const small_vector old_left(left);
    small_vector right(right_size, TestType{ 2 });
    const small_vector old_right(right);

    using std::swap;
    swap(left, right);

    REQUIRE(right == old_left);
    REQUIRE(left == old_right);

    swap(left, right);

    REQUIRE(left == old_left);
    REQUIRE(right == old_right);
}

TEMPLATE_TEST_CASE("push_back(const T&)", "[modifiers]", TrivialType, NonTrivialType, NonDefaultConstructibleType)
{
    const size_t size = GENERATE(SMALL_SIZE, LARGE_SIZE);

    small_vector vec(size, TestType{ 0 });
    const TestType elem(1);

    for (size_t i = 0; i < LARGE_SIZE + 1; i++) vec.push_back(elem);

    REQUIRE(vec.size() == size + LARGE_SIZE + 1);
    REQUIRE(vec.front() == TestType{ 0 });
    REQUIRE(vec.back() == TestType{ 1 });
}

TEMPLATE_TEST_CASE("push_back(T&&)", "[modifiers]", TrivialType, NonTrivialType, MoveOnlyType)
{
    const size_t size = GENERATE(SMALL_SIZE, LARGE_SIZE);

    small_vector<TestType> vec(size);

    for (size_t i = 0; i < LARGE_SIZE + 1; i++) vec.push_back(TestType{ 1 });

    REQUIRE(vec.size() == size + LARGE_SIZE + 1);
    REQUIRE(vec.front() == TestType{});
    REQUIRE(vec.back() == TestType{ 1 });
}

TEMPLATE_TEST_CASE("emplace_back(...)", "[modifiers]", TrivialType, NonTrivialType, MoveOnlyType)
{
    const size_t size = GENERATE(SMALL_SIZE, LARGE_SIZE);

    small_vector<TestType> vec(size);

    const TestType& val = vec.emplace_back();
    REQUIRE(val == TestType{});

    for (size_t i = 0; i < LARGE_SIZE; i++) vec.emplace_back(1);

    REQUIRE(vec.size() == size + LARGE_SIZE + 1);
    REQUIRE(vec.front() == TestType{});
    REQUIRE(vec.back() == TestType{ 1 });
}

TEMPLATE_TEST_CASE("pop_back()", "[modifiers]", TrivialType, NonTrivialType, MoveOnlyType, ImmovableType)
{
    const size_t size = GENERATE(SMALL_SIZE, LARGE_SIZE);
    small_vector<TestType> vec(size);

    vec.pop_back();
    REQUIRE(vec.size() == (size - 1));

    vec.pop_back();
    REQUIRE(vec.size() == (size - 2));
}

TEMPLATE_TEST_CASE("resize(n)", "[modifiers]", TrivialType, NonTrivialType, MoveOnlyType)
{
    const size_t size = GENERATE(SMALL_SIZE, LARGE_SIZE);
    small_vector<TestType> vec(size);

    vec.resize(2 * LARGE_SIZE);
    REQUIRE(vec.size() == 2 * LARGE_SIZE);
    REQUIRE(vec.back() == TestType{});

    vec.resize(0);
    REQUIRE(vec.empty());
}

TEMPLATE_TEST_CASE("resize(n, value)", "[modifiers]", TrivialType, NonTrivialType, NonDefaultConstructibleType)
{
    const size_t size = GENERATE(SMALL_SIZE, LARGE_SIZE);
    small_vector vec(size, TestType{ 1 });

    vec.resize(2 * LARGE_SIZE, TestType{ 2 });
    REQUIRE(vec.size() == 2 * LARGE_SIZE);
    REQUIRE(vec.back() == TestType{ 2 });

    vec.resize(0, TestType{ 3 });
    REQUIRE(vec.empty());
}

TEMPLATE_TEST_CASE("erase(pos)", "[modifiers]", TrivialType, NonTrivialType, MoveOnlyType)
{
    const size_t size = GENERATE(SMALL_SIZE, LARGE_SIZE);

    small_vector<TestType> vec(size);
    vec.front() = TestType{ 2 };

    auto it1 = vec.erase(vec.begin());

    REQUIRE(it1 == vec.begin());
    REQUIRE(vec.size() == size - 1);
    REQUIRE(vec.front() == TestType{});

    auto it2 = vec.erase(vec.end() - 1);

    REQUIRE(it2 == vec.end());
    REQUIRE(vec.size() == size - 2);
}

TEMPLATE_TEST_CASE("erase(first, last)", "[modifiers]", TrivialType, NonTrivialType)
{
    small_vector vec{ TestType{ 0 }, TestType{ 1 }, TestType{ 2 }, TestType{ 3 }, TestType{ 4 } };

    SECTION("front")
    {
        auto it = vec.erase(vec.begin(), vec.begin() + 1);

        REQUIRE(vec == small_vector{ TestType{ 1 }, TestType{ 2 }, TestType{ 3 }, TestType{ 4 } });
        REQUIRE(*it == TestType{ 1 });
    }
    SECTION("back")
    {
        auto it = vec.erase(vec.end() - 1, vec.end());

        REQUIRE(vec == small_vector{ TestType{ 0 }, TestType{ 1 }, TestType{ 2 }, TestType{ 3 } });
        REQUIRE(it == vec.end());
    }
    SECTION("nothing")
    {
        auto it = vec.erase(vec.begin(), vec.begin());

        REQUIRE(vec == small_vector{ TestType{ 0 }, TestType{ 1 }, TestType{ 2 }, TestType{ 3 }, TestType{ 4 } });
        REQUIRE(it == vec.begin());
    }
    SECTION("everything")
    {
        auto it = vec.erase(vec.begin(), vec.end());

        REQUIRE(vec.empty());
        REQUIRE(it == vec.end());

    }
    SECTION("middle")
    {
        auto it = vec.erase(vec.begin() + 1, vec.begin() + 3);

        REQUIRE(vec == small_vector{ TestType{ 0 }, TestType{ 3 }, TestType{ 4 } });
        REQUIRE(*it == TestType{ 3 });
    }
}

TEMPLATE_TEST_CASE("insert(pos, const T&)", "[modifiers]", TrivialType, NonTrivialType)
{
    small_vector vec{ TestType{ 0 }, TestType{ 1 }, TestType{ 2 }, TestType{ 3 } };
    const TestType value{ 21 };

    SECTION("front")
    {
        auto it = vec.insert(vec.begin(), value);

        REQUIRE(vec == small_vector{ TestType{ 21 }, TestType{ 0 }, TestType{ 1 }, TestType{ 2 }, TestType{ 3 } });
        REQUIRE(*it == value);
    }
    SECTION("back")
    {
        auto it = vec.insert(vec.end(), value);

        REQUIRE(vec == small_vector{ TestType{ 0 }, TestType{ 1 }, TestType{ 2 }, TestType{ 3 }, TestType{ 21 } });
        REQUIRE(*it == value);
    }
    SECTION("middle")
    {
        auto it = vec.insert(vec.begin() + 2, value);

        REQUIRE(vec == small_vector{ TestType{ 0 }, TestType{ 1 }, TestType{ 21 }, TestType{ 2 }, TestType{ 3 } });
        REQUIRE(*it == value);
    }
    SECTION("many")
    {
        for (size_t i = 0; i < 100; i++) vec.insert(vec.begin(), value);

        REQUIRE(vec.size() == 104);
        REQUIRE(vec.front() == TestType{ 21 });
        REQUIRE(vec.back() == TestType{ 3 });
    }
}

TEMPLATE_TEST_CASE("insert(pos, T&&)", "[modifiers]", TrivialType, NonTrivialType)
{
    small_vector vec{ TestType{ 0 }, TestType{ 1 }, TestType{ 2 }, TestType{ 3 } };

    SECTION("front")
    {
        auto it = vec.insert(vec.begin(), TestType{ 21 });

        REQUIRE(vec == small_vector{ TestType{ 21 }, TestType{ 0 }, TestType{ 1 }, TestType{ 2 }, TestType{ 3 } });
        REQUIRE(*it == TestType{ 21 });

    }
    SECTION("back")
    {
        auto it = vec.insert(vec.end(), TestType{ 21 });

        REQUIRE(vec == small_vector{ TestType{ 0 }, TestType{ 1 }, TestType{ 2 }, TestType{ 3 }, TestType{ 21 } });
        REQUIRE(*it == TestType{ 21 });
    }
    SECTION("middle")
    {
        auto it = vec.insert(vec.begin() + 2, TestType{ 21 });

        REQUIRE(vec == small_vector{ TestType{ 0 }, TestType{ 1 }, TestType{ 21 }, TestType{ 2 }, TestType{ 3 } });
        REQUIRE(*it == TestType{ 21 });
    }
    SECTION("many")
    {
        for (size_t i = 0; i < 100; i++) vec.insert(vec.begin(), TestType{ 21 });

        REQUIRE(vec.size() == 104);
        REQUIRE(vec.front() == TestType{ 21 });
        REQUIRE(vec.back() == TestType{ 3 });
    }
}

TEMPLATE_TEST_CASE("insert(pos, count, const T&)", "[modifiers]", TrivialType, NonTrivialType)
{
    small_vector dest{ TestType{ 0 }, TestType{ 1 } };

    SECTION("front")
    {
        auto it = dest.insert(dest.begin(), 3, TestType{ 4 });

        REQUIRE(dest == small_vector{ TestType{ 4 }, TestType{ 4 }, TestType{ 4 }, TestType{ 0 }, TestType{ 1 } });
        REQUIRE(it == dest.begin());
        REQUIRE(*it == TestType{ 4 });
    }
    SECTION("back")
    {
        auto it = dest.insert(dest.end(), 3, TestType{ 4 });

        REQUIRE(dest == small_vector{ TestType{ 0 }, TestType{ 1 }, TestType{ 4 }, TestType{ 4 }, TestType{ 4 } });
        REQUIRE(it == (dest.begin() + 2));
        REQUIRE(*it == TestType{ 4 });
    }
    SECTION("middle")
    {
        auto it = dest.insert(dest.begin() + 1, 3, TestType{ 4 });

        REQUIRE(dest == small_vector{ TestType{ 0 }, TestType{ 4 }, TestType{ 4 }, TestType{ 4 }, TestType{ 1 } });
        REQUIRE(it == (dest.begin() + 1));
        REQUIRE(*it == TestType{ 4 });
    }
    SECTION("nothing")
    {
        auto it = dest.insert(dest.end(), 0, TestType{ 4 });

        REQUIRE(dest == small_vector{ TestType{ 0 }, TestType{ 1 } });
        REQUIRE(it == dest.end());
    }
    SECTION("many")
    {
        for (size_t i = 0; i < 100; i++) dest.insert(dest.begin(), 3, TestType{ 4 });

        REQUIRE(dest.size() == 302);
        REQUIRE(dest.front() == TestType{ 4 });
        REQUIRE(dest.back() == TestType{ 1 });
    }
}

TEMPLATE_TEST_CASE("insert(pos, FwdIter, FwdIter)", "[modifiers]", TrivialType, NonTrivialType)
{
    small_vector dest{ TestType{ 0 }, TestType{ 1 } };
    const small_vector src{ TestType{ 2 }, TestType{ 3 }, TestType{ 4 } };

    SECTION("front")
    {
        auto it = dest.insert(dest.begin(), src.begin(), src.end());

        REQUIRE(dest == small_vector{ TestType{ 2 }, TestType{ 3 }, TestType{ 4 }, TestType{ 0 }, TestType{ 1 } });
        REQUIRE(*it == TestType{ 2 });
    }
    SECTION("back")
    {
        auto it = dest.insert(dest.end(), src.begin(), src.end());

        REQUIRE(dest == small_vector{ TestType{ 0 }, TestType{ 1 }, TestType{ 2 }, TestType{ 3 }, TestType{ 4 } });
        REQUIRE(*it == TestType{ 2 });
    }
    SECTION("middle")
    {
        auto it = dest.insert(dest.begin() + 1, src.begin(), src.end());

        REQUIRE(dest == small_vector{ TestType{ 0 }, TestType{ 2 }, TestType{ 3 }, TestType{ 4 }, TestType{ 1 } });
        REQUIRE(*it == TestType{ 2 });
    }
    SECTION("empty")
    {
        auto it = dest.insert(dest.end(), src.begin(), src.begin());

        REQUIRE(dest == small_vector{ TestType{ 0 }, TestType{ 1 } });
        REQUIRE(it == dest.end());
    }
    SECTION("many")
    {
        for (size_t i = 0; i < 100; i++) dest.insert(dest.begin(), src.begin(), src.end());

        REQUIRE(dest.size() == 100 * src.size() + 2);
        REQUIRE(dest.front() == TestType{ 2 });
        REQUIRE(dest.back() == TestType{ 1 });
    }
}

TEST_CASE("insert(pos, InputIter, InputIter)", "[modifiers]")
{
    small_vector<char> dest(2, 'a');
    std::istringstream src(std::string(3, 'c'));

    SECTION("front")
    {
        auto it = dest.insert(dest.begin(), std::istream_iterator<char>(src), std::istream_iterator<char>());

        REQUIRE(dest == small_vector{ 'c', 'c', 'c', 'a', 'a' });
        REQUIRE(it == dest.begin());
    }
    SECTION("back")
    {
        auto it = dest.insert(dest.end(), std::istream_iterator<char>(src), std::istream_iterator<char>());

        REQUIRE(dest == small_vector{ 'a', 'a', 'c', 'c', 'c' });
        REQUIRE(it == (dest.begin() + 2));
        REQUIRE(*it == 'c');
    }
    SECTION("middle")
    {
        auto it = dest.insert(dest.begin() + 1, std::istream_iterator<char>(src), std::istream_iterator<char>());

        REQUIRE(dest == small_vector{ 'a', 'c', 'c', 'c', 'a' });
        REQUIRE(it == (dest.begin() + 1));
        REQUIRE(*it == 'c');
    }
    SECTION("empty")
    {
        auto it = dest.insert(dest.end(), std::istream_iterator<char>(), std::istream_iterator<char>());

        REQUIRE(dest == small_vector{ 'a', 'a' });
        REQUIRE(it == dest.end());
    }
    SECTION("many")
    {
        for (size_t i = 0; i < 100; i++)
        {
            std::istringstream input(std::string(3, 'c'));
            dest.insert(dest.begin(), std::istream_iterator<char>(input), std::istream_iterator<char>());
        }

        REQUIRE(dest.size() == 302);
        REQUIRE(dest.front() == 'c');
        REQUIRE(dest.back() == 'a');
    }
}

TEMPLATE_TEST_CASE("emplace(pos, Args&&...)", "[modifiers]", TrivialType, NonTrivialType, MoveOnlyType)
{
    small_vector<TestType> vec(2);

    auto it1 = vec.emplace(vec.begin(), 1);

    REQUIRE(*it1 == 1);
    REQUIRE(vec.size() == 3);
    REQUIRE(vec.front() == TestType{ 1 });

    auto it2 = vec.emplace(vec.end(), 2);

    REQUIRE(*it2 == 2);
    REQUIRE(vec.size() == 4);
    REQUIRE(vec.back() == TestType{ 2 });

    auto it3 = vec.emplace(vec.begin() + 1);

    REQUIRE(*it3 == TestType{});
    REQUIRE(vec.size() == 5);

    for (size_t i = 0; i < 100; i++) vec.emplace(vec.begin());
    REQUIRE(vec.size() == 105);
    REQUIRE(vec.front() == TestType{});
    REQUIRE(vec.back() == TestType{ 2 });
}

    //-----------------------------------//
    //       ALLOCATOR PROPAGATION       //
    //-----------------------------------//

template<typename T>
struct DummyAllocator
{
    using value_type = T;

    T* allocate(std::size_t size) const
    {
        if (size == 0) return nullptr;

        const std::size_t bytes_needed = size * sizeof(value_type);
        void* const storage = ::operator new[](bytes_needed);

        return static_cast<T*>(storage);
    }

    void deallocate(T* storage, std::size_t) const noexcept
    {
        ::operator delete[](storage);
    }

    using propagate_on_container_copy_assignment = std::true_type;
    using propagate_on_container_move_assignment = std::true_type;
    using propagate_on_container_swap = std::true_type;

    using is_always_equal = std::false_type;
    friend bool operator==(const DummyAllocator&, const DummyAllocator&) { return false; }
};


template<typename T>
using small_vector2 = small_vector<T, 8, DummyAllocator<T>>;


TEMPLATE_TEST_CASE("propagate_on_copy_assignment", "[allocators]", TrivialType, NonTrivialType)
{
    const size_t src_size = GENERATE(SMALL_SIZE, LARGE_SIZE + 1);
    const size_t dest_size = GENERATE(SMALL_SIZE, LARGE_SIZE);

    const small_vector2<TestType> source(src_size, 4);
    small_vector2<TestType> dest(dest_size, 3);

    dest = source;

    REQUIRE(dest.size() == source.size());
    REQUIRE(dest == source);
}

TEMPLATE_TEST_CASE("propagate_on_move_assignment", "[allocators]", TrivialType, NonTrivialType)
{
    const size_t src_size = GENERATE(SMALL_SIZE, LARGE_SIZE + 1);
    const size_t dest_size = GENERATE(SMALL_SIZE, LARGE_SIZE);

    small_vector2<TestType> source(src_size, 4);
    const small_vector2<TestType> src_copy(source);
    small_vector2<TestType> dest(dest_size, 3);

    dest = std::move(source);

    REQUIRE(dest.size() == src_copy.size());
    REQUIRE(dest == src_copy);

    REQUIRE(source.capacity() != 0);

    source = dest;
    REQUIRE(source == dest);
}

TEMPLATE_TEST_CASE("propagate_on_swap", "[allocators]", TrivialType, NonTrivialType)
{
    const size_t left_size = GENERATE(SMALL_SIZE, LARGE_SIZE);
    const size_t right_size = GENERATE(SMALL_SIZE, LARGE_SIZE);

    small_vector2<TestType> left(left_size, 1);
    const small_vector2<TestType> old_left(left);
    small_vector2<TestType> right(right_size, 2);
    const small_vector2<TestType> old_right(right);

    using std::swap;
    swap(left, right);

    REQUIRE(right == old_left);
    REQUIRE(left == old_right);

    swap(left, right);

    REQUIRE(left == old_left);
    REQUIRE(right == old_right);
}
