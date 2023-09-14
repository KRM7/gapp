﻿/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include "utility/concepts.hpp"
#include <vector>
#include <array>
#include <set>
#include <string>

using namespace gapp;

template<typename T>
struct A {};

struct B : A<int> {};

TEST_CASE("hashable", "[concepts]")
{
    STATIC_REQUIRE(detail::Hashable<int>);
    STATIC_REQUIRE(detail::Hashable<std::string>);
    
    STATIC_REQUIRE(!detail::Hashable<A<double>>);
    STATIC_REQUIRE(!detail::Hashable<std::vector<int>>);
}

TEST_CASE("container", "[concepts]")
{
    STATIC_REQUIRE(detail::Container<std::vector<int>>);
    STATIC_REQUIRE(detail::Container<std::array<int, 3>>);
    STATIC_REQUIRE(detail::Container<std::set<std::string>>);

    STATIC_REQUIRE(detail::Container<std::vector<bool>>);

    STATIC_REQUIRE(!detail::Container<int*>);
    STATIC_REQUIRE(!detail::Container<void>);
}

TEST_CASE("indexable_container", "[concepts]")
{
    STATIC_REQUIRE(detail::IndexableContainer<std::vector<int>>);
    STATIC_REQUIRE(detail::IndexableContainer<std::array<int, 3>>);
    STATIC_REQUIRE(detail::IndexableContainer<std::string>);

    STATIC_REQUIRE(detail::IndexableContainer<std::vector<bool>>);

    STATIC_REQUIRE(!detail::IndexableContainer<int*>);
    STATIC_REQUIRE(!detail::IndexableContainer<void>);
}