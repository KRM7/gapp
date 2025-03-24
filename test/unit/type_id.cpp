/* Copyright (c) 2025 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include "utility/type_id.hpp"

using namespace gapp;

TEST_CASE("type_id", "[type_id]")
{
    REQUIRE(detail::type_id<int>() == detail::type_id<int>());
    REQUIRE(detail::type_id<int>() != detail::type_id<const int>());
}
