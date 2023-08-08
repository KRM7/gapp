/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include "utility/scope_exit.hpp"

using namespace gapp::detail;


TEST_CASE("scope_exit", "[utility]")
{
    int n = 0;

    {
        ScopeExit on_exit{ [&] { n = 2; } };
        REQUIRE(n == 0);
    }
    REQUIRE(n == 2);

    {
        ScopeExit on_exit{ [&] { n = 3; } };
        on_exit.release();
    }   
    REQUIRE(n == 2);
}

TEST_CASE("restore_on_exit", "[utility]")
{
    int n = 3;
    double f = 2.5;

    {
        RestoreOnExit on_exit(n, f);

        n = 10;
        f = 0.2;

        REQUIRE(n == 10);
        REQUIRE(f == 0.2);
    }

    REQUIRE(n == 3);
    REQUIRE(f == 2.5);
}
