/* Copyright (c) 2025 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include "utility/hash.hpp"
#include <vector>
#include <tuple>
#include <string>

using namespace gapp;

template<typename T>
struct A {};

struct B : A<int> {};

TEST_CASE("hashable", "[hash]")
{
    STATIC_REQUIRE(detail::hashable<int>);
    STATIC_REQUIRE(detail::hashable<std::string>);
    
    STATIC_REQUIRE(!detail::hashable<A<double>>);
    STATIC_REQUIRE(!detail::hashable<std::vector<int>>);
}

TEST_CASE("hash", "[hash]")
{
    std::ignore = detail::hash(12);
    std::ignore = detail::hash_combine(1, 2, 3);
    std::ignore = detail::hash_range(std::vector{ 1, 2, 3 });

    SUCCEED();
}
