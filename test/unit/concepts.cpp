/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include "utility/concepts.hpp"
#include <vector>
#include <string>

using namespace gapp;

template<typename T>
struct A {};

struct B : A<int> {};

TEST_CASE("hashable", "[concepts]")
{
    STATIC_REQUIRE(detail::hashable<int>);
    STATIC_REQUIRE(detail::hashable<std::string>);
    
    STATIC_REQUIRE(!detail::hashable<A<double>>);
    STATIC_REQUIRE(!detail::hashable<std::vector<int>>);
}