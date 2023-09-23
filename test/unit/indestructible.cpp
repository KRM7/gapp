/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include "utility/indestructible.hpp"
#include "utility/utility.hpp"
#include <stdexcept>

using namespace gapp::detail;

struct Widget
{
    Widget() : n(3) {}
    Widget(int val) : n(val) {}
    ~Widget() noexcept(false) { GAPP_THROW(std::logic_error, ""); }

    int n = 0;
};

TEST_CASE("indestructible", "[indestructible]")
{
    Indestructible<Widget> widget1;
    REQUIRE(widget1->n == 3);

    Indestructible<Widget> widget2(4);
    REQUIRE(widget2->n == 4);
}
