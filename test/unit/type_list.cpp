/* Copyright (c) 2025 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include "utility/type_list.hpp"
#include "utility/type_id.hpp"
#include <type_traits>

using namespace gapp::detail;

using empty_type_list = type_list<>;
using test_type_list = type_list<void, int, double>;

TEST_CASE("tuple_to_list", "[type_list]")
{
    STATIC_REQUIRE(std::is_same_v<tuple_to_list_t<std::tuple<>>, type_list<>>);
    STATIC_REQUIRE(std::is_same_v<tuple_to_list_t<std::tuple<void, int>>, type_list<void, int>>);
}

TEST_CASE("args_to_list", "[type_list]")
{
    STATIC_REQUIRE(std::is_same_v<args_to_list_t<std::tuple<>>, type_list<>>);
    STATIC_REQUIRE(std::is_same_v<args_to_list_t<std::tuple<void, int>>, type_list<void, int>>);
}

TEST_CASE("to_tuple", "[type_list]")
{
    STATIC_REQUIRE(std::is_same_v<type_list<>::to_tuple, std::tuple<>>);
    STATIC_REQUIRE(std::is_same_v<type_list<void>::to_tuple, std::tuple<void>>);
    STATIC_REQUIRE(std::is_same_v<type_list<int, float, long>::to_tuple, std::tuple<int, float, long>>);
}

TEST_CASE("type_list_size", "[type_list]")
{
    STATIC_REQUIRE(empty_type_list::size == 0);
    STATIC_REQUIRE(test_type_list::size == 3);
}

TEST_CASE("type_list_contains", "[type_list]")
{
    STATIC_REQUIRE(!empty_type_list::contains<void>);
    STATIC_REQUIRE(test_type_list::contains<void>);

    STATIC_REQUIRE(test_type_list::contains<double>);

    STATIC_REQUIRE(!test_type_list::contains<const double>);
    STATIC_REQUIRE(!test_type_list::contains<float>);
}

TEST_CASE("type_list_index_of", "[type_list]")
{
    STATIC_REQUIRE(test_type_list::index_of<void> == 0);
    STATIC_REQUIRE(test_type_list::index_of<int> == 1);
    STATIC_REQUIRE(test_type_list::index_of<double> == 2);
}

TEST_CASE("type_list_to_tuple", "[type_list]")
{
    STATIC_REQUIRE(std::is_same_v<std::tuple<>, empty_type_list::to_tuple>);
    STATIC_REQUIRE(std::is_same_v<std::tuple<void, int, double>, test_type_list::to_tuple>);
}

TEST_CASE("filter_type_list", "[type_list]")
{
    STATIC_REQUIRE(empty_type_list::filter_types_t<std::is_void>::size == 0);
    STATIC_REQUIRE(test_type_list::filter_types_t<std::is_void>::size == 1);

    STATIC_REQUIRE(test_type_list::filter_types_t<std::is_arithmetic>::size == 2);
}

TEST_CASE("type_list_apply", "[type_list]")
{
    STATIC_REQUIRE(empty_type_list::apply([]<typename... Ts>() { return sizeof...(Ts); }) == 0);
    STATIC_REQUIRE(test_type_list::apply([]<typename... Ts>() { return sizeof...(Ts); }) == 3);
}

TEST_CASE("type_list_for_each", "[type_list]")
{
    STATIC_REQUIRE((empty_type_list::for_each([]<typename T>(size_t) {}), true));

    empty_type_list::for_each([]<typename T>(size_t) { FAIL(); });
    test_type_list::for_each([]<typename T>(size_t i) { if (i > 2) FAIL(); });
}

TEST_CASE("type_list_find_index", "[type_list]")
{
    STATIC_REQUIRE(!empty_type_list::find_index([]<typename T>() { return true; }));
    STATIC_REQUIRE(!empty_type_list::find_index([]<typename T>() { return false; }));

    STATIC_REQUIRE(test_type_list::find_index([]<typename T>() { return true; }) == 0);
    STATIC_REQUIRE(!test_type_list::find_index([]<typename T>() { return false; }));

    STATIC_REQUIRE(test_type_list::find_index([]<typename T>() { return std::is_same_v<T, double>; }) == 2);
}

TEST_CASE("type_list_find_type_id", "[type_list]")
{
    REQUIRE(!empty_type_list::index_of_id(type_id<void>()));

    REQUIRE(!test_type_list::index_of_id(type_id<const int>()));
    REQUIRE(test_type_list::index_of_id(type_id<int>()) == 1);
}
