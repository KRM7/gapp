/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include "utility/cone_tree.hpp"
#include <type_traits>

using namespace gapp::detail;

TEST_CASE("cone_tree constructors", "[cone_tree]")
{
    STATIC_REQUIRE(std::is_nothrow_default_constructible_v<ConeTree>);
    STATIC_REQUIRE(std::is_copy_constructible_v<ConeTree>);
    STATIC_REQUIRE(std::is_copy_assignable_v<ConeTree>);
    STATIC_REQUIRE(std::is_nothrow_move_constructible_v<ConeTree>);
    STATIC_REQUIRE(std::is_nothrow_move_assignable_v<ConeTree>);

    std::vector<ConeTree::Point> points = { { 0.0, 1.0 }, { 1.0, 0.0 }, { 1.0, 1.0 }, { 0.0, 0.0 } };
    ConeTree tree(points);

    REQUIRE(tree.size() == points.size());
}

TEST_CASE("cone_tree lookup", "[cone_tree]")
{
    std::vector<ConeTree::Point> points = {
        { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.2 }, { 0.0, 0.2, 0.0 }, { 0.0, 0.2, 0.2 },
        { 0.0, 0.2, 0.4 }, { 0.0, 0.2, 0.6 }, { 0.0, 0.2, 0.8 }, { 0.0, 0.4, 0.2 },
        { 0.0, 0.4, 0.4 }, { 0.0, 0.4, 0.6 }, { 0.0, 0.4, 0.8 }, { 0.0, 0.6, 0.2 },
        { 0.0, 0.6, 0.4 }, { 0.0, 0.6, 0.6 }, { 0.0, 0.6, 0.8 }, { 0.0, 0.8, 0.2 },
        { 0.0, 0.8, 0.4 }, { 0.0, 0.8, 0.6 }, { 0.0, 0.8, 0.8 }, { 0.2, 0.0, 0.0 },
        { 0.2, 0.0, 0.2 }, { 0.2, 0.0, 0.4 }, { 0.2, 0.0, 0.6 }, { 0.2, 0.0, 0.8 },
        { 0.2, 0.2, 0.0 }, { 0.2, 0.2, 0.2 }, { 0.2, 0.2, 0.4 }, { 0.2, 0.2, 0.6 },
        { 0.2, 0.2, 0.8 }, { 0.2, 0.4, 0.0 }, { 0.2, 0.4, 0.2 }, { 0.2, 0.4, 0.4 },
        { 0.2, 0.4, 0.6 }, { 0.2, 0.4, 0.8 }, { 0.2, 0.6, 0.0 }, { 0.2, 0.6, 0.2 },
        { 0.2, 0.6, 0.4 }, { 0.2, 0.6, 0.6 }, { 0.2, 0.6, 0.8 }, { 0.2, 0.8, 0.0 },
        { 0.2, 0.8, 0.2 }, { 0.2, 0.8, 0.4 }, { 0.2, 0.8, 0.6 }, { 0.2, 0.8, 0.8 },
        { 0.4, 0.0, 0.2 }, { 0.4, 0.0, 0.4 }, { 0.4, 0.0, 0.6 }, { 0.4, 0.0, 0.8 },
        { 0.4, 0.2, 0.0 }, { 0.4, 0.2, 0.2 }, { 0.4, 0.2, 0.4 }, { 0.4, 0.2, 0.6 },
        { 0.4, 0.2, 0.8 }, { 0.4, 0.4, 0.2 }, { 0.4, 0.4, 0.6 }, { 0.4, 0.6, 0.0 },
        { 0.4, 0.6, 0.2 }, { 0.4, 0.6, 0.4 }, { 0.4, 0.6, 0.6 }, { 0.4, 0.6, 0.8 },
        { 0.4, 0.8, 0.0 }, { 0.4, 0.8, 0.2 }, { 0.4, 0.8, 0.6 }, { 0.6, 0.0, 0.2 },
        { 0.6, 0.0, 0.4 }, { 0.6, 0.0, 0.8 }, { 0.6, 0.2, 0.0 }, { 0.6, 0.2, 0.2 },
        { 0.6, 0.2, 0.4 }, { 0.6, 0.2, 0.6 }, { 0.6, 0.2, 0.8 }, { 0.6, 0.4, 0.0 },
        { 0.6, 0.4, 0.2 }, { 0.6, 0.4, 0.4 }, { 0.6, 0.4, 0.6 }, { 0.6, 0.4, 0.8 },
        { 0.6, 0.6, 0.2 }, { 0.6, 0.6, 0.4 }, { 0.6, 0.6, 0.8 }, { 0.6, 0.8, 0.0 },
        { 0.6, 0.8, 0.2 }, { 0.6, 0.8, 0.4 }, { 0.6, 0.8, 0.6 }, { 0.6, 0.8, 0.8 },
        { 0.8, 0.0, 0.2 }, { 0.8, 0.0, 0.6 }, { 0.8, 0.2, 0.0 }, { 0.8, 0.2, 0.2 },
        { 0.8, 0.2, 0.4 }, { 0.8, 0.2, 0.6 }, { 0.8, 0.2, 0.8 }, { 0.8, 0.4, 0.2 },
        { 0.8, 0.4, 0.6 }, { 0.8, 0.6, 0.0 }, { 0.8, 0.6, 0.2 }, { 0.8, 0.6, 0.4 },
        { 0.8, 0.6, 0.6 }, { 0.8, 0.6, 0.8 }, { 0.8, 0.8, 0.2 }, { 0.8, 0.8, 0.6 },
    };

    ConeTree tree(points);

    auto best = tree.findBestMatch(ConeTree::Point{ 1.0, 1.0, 0.1 });
    REQUIRE(*best.elem == ConeTree::Point{ 0.8, 0.8, 0.6 });

    best = tree.findBestMatch(ConeTree::Point{ 0.1, 0.5, 0.8 });
    REQUIRE(*best.elem == ConeTree::Point{ 0.6, 0.8, 0.8 });
}

TEST_CASE("empty_cone_tree", "[cone_tree]")
{
    ConeTree tree;

    auto best = tree.findBestMatch(ConeTree::Point{ 1.0, 1.0 });

    REQUIRE(best.elem == tree.end());
}
