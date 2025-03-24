/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include "utility/type_traits.hpp"
#include <vector>
#include <list>

using namespace gapp::detail;

template<typename T>
class Vec : public std::vector<T> {};

TEST_CASE("concat_tup", "[type_traits]")
{
    STATIC_REQUIRE(std::is_same_v<concat_tup_t<std::tuple<int, double>, std::tuple<>>, std::tuple<int, double>>);
    STATIC_REQUIRE(std::is_same_v<concat_tup_t<std::tuple<>, std::tuple<int, double>>, std::tuple<int, double>>);

    STATIC_REQUIRE(std::is_same_v<concat_tup_t<std::tuple<int, double>, std::tuple<float>>, std::tuple<int, double, float>>);

    STATIC_REQUIRE(std::is_same_v<concat_tup_t<std::tuple<int, double>, float>, std::tuple<int, double, float>>);
    STATIC_REQUIRE(std::is_same_v<concat_tup_t<float, std::tuple<int, double>>, std::tuple<float, int, double>>);
}

TEST_CASE("filter_types", "[type_traits]")
{
    STATIC_REQUIRE(std::is_same_v<filter_types_t<std::is_floating_point>, std::tuple<>>);
    STATIC_REQUIRE(std::is_same_v<filter_types_t<std::is_floating_point, int, void, float, double, long>, std::tuple<float, double>>);
}

TEST_CASE("map_types", "[type_traits]")
{
    STATIC_REQUIRE(std::is_same_v<map_types_t<std::tuple, std::tuple<>>, std::tuple<>>);
    STATIC_REQUIRE(std::is_same_v<map_types_t<std::tuple, std::tuple<int, double, long>>, std::tuple<std::tuple<int>, std::tuple<double>, std::tuple<long>>>);
}

TEST_CASE("is_same_template", "[type_traits]")
{
    STATIC_REQUIRE(is_same_template_v<std::vector, std::vector>);
    STATIC_REQUIRE(!is_same_template_v<std::vector, std::list>);
    STATIC_REQUIRE(!is_same_template_v<std::vector, Vec>);
}

TEST_CASE("is_one_of_templates", "[type_traits]")
{
    STATIC_REQUIRE(is_one_of_templates_v<std::vector, std::vector, Vec>);
    STATIC_REQUIRE(is_one_of_templates_v<std::vector, std::vector>);

    STATIC_REQUIRE(!is_one_of_templates_v<Vec, std::vector>);
    STATIC_REQUIRE(!is_one_of_templates_v<std::vector>);
}

TEST_CASE("unique_types", "[type_traits]")
{
    STATIC_REQUIRE(unique_types_v<>);
    STATIC_REQUIRE(unique_types_v<void>);
    STATIC_REQUIRE(unique_types_v<void, int, double>);
    STATIC_REQUIRE(unique_types_v<const int, int>);

    STATIC_REQUIRE(!unique_types_v<int, int>);
    STATIC_REQUIRE(!unique_types_v<void, int, float, double, void>);
}

TEST_CASE("nth_type", "[type_traits]")
{
    STATIC_REQUIRE(std::is_same_v<nth_type_t<1, int, void>, void>);
    STATIC_REQUIRE(std::is_same_v<nth_type_t<3, int, void, double, float>, float>);
    STATIC_REQUIRE(std::is_same_v<nth_type_t<0, int, void>, int>);
}

TEST_CASE("index_of_type", "[type_traits]")
{
    STATIC_REQUIRE(index_of_type_v<int, int> == 0);
    STATIC_REQUIRE(index_of_type_v<int, int, int, int> == 0);
    STATIC_REQUIRE(index_of_type_v<double, int, double, int> == 1);
    STATIC_REQUIRE(index_of_type_v<void, int, float, long, void> == 3);
}

TEST_CASE("is_derived_from_spec_of", "[type_traits]")
{
    STATIC_REQUIRE(is_derived_from_spec_of_v<Vec<int>, std::vector>);
    STATIC_REQUIRE(is_derived_from_spec_of_v<std::vector<int>, std::vector>);

    STATIC_REQUIRE(!is_derived_from_spec_of_v<void, std::vector>);
    STATIC_REQUIRE(!is_derived_from_spec_of_v<int*, std::vector>);
}

TEST_CASE("is_specialization_of", "[type_traits]")
{
    STATIC_REQUIRE(is_specialization_of_v<std::vector<int>, std::vector>);
    STATIC_REQUIRE(!is_specialization_of_v<void, std::vector>);
}

TEST_CASE("is_reverse_iterator", "[type_traits]")
{
    using Iter = std::vector<double>::iterator;
    using RevIter = std::vector<double>::reverse_iterator;

    STATIC_REQUIRE(!is_reverse_iterator_v<Iter>);
    STATIC_REQUIRE(is_reverse_iterator_v<RevIter>);
}

TEST_CASE("dereference", "[type_traits]")
{
    using Iter = std::vector<double>::const_iterator;

    STATIC_REQUIRE(std::is_same_v<dereference_t<int*>, int&>);
    STATIC_REQUIRE(std::is_same_v<dereference_t<Iter>, const double&>);
}

TEST_CASE("remove_rvalue_ref", "[type_traits]")
{
    STATIC_REQUIRE(std::is_same_v<remove_rvalue_ref_t<int>,   int>);
    STATIC_REQUIRE(std::is_same_v<remove_rvalue_ref_t<int&>,  int&>);
    STATIC_REQUIRE(std::is_same_v<remove_rvalue_ref_t<int&&>, int>);
    STATIC_REQUIRE(std::is_same_v<remove_rvalue_ref_t<int*>,  int*>);
}

TEST_CASE("copy_const_t", "[type_traits]")
{
    STATIC_REQUIRE(std::is_same_v<copy_const_t<const int, double>, const double>);
    STATIC_REQUIRE(std::is_same_v<copy_const_t<int, double>, double>);
    STATIC_REQUIRE(std::is_same_v<copy_const_t<int, const double>, const double>);
    STATIC_REQUIRE(std::is_same_v<copy_const_t<const int, const double>, const double>);
}

TEST_CASE("copy_volatile_t", "[type_traits]")
{
    STATIC_REQUIRE(std::is_same_v<copy_volatile_t<volatile int, double>, volatile double>);
    STATIC_REQUIRE(std::is_same_v<copy_volatile_t<int, double>, double>);
    STATIC_REQUIRE(std::is_same_v<copy_volatile_t<int, volatile double>, volatile double>);
    STATIC_REQUIRE(std::is_same_v<copy_volatile_t<volatile int, volatile double>, volatile double>);
}

TEST_CASE("copy_cv_t", "[type_traits]")
{
    STATIC_REQUIRE(std::is_same_v<copy_cv_t<const int, volatile double>, const volatile double>);
    STATIC_REQUIRE(std::is_same_v<copy_cv_t<const volatile int, double>, const volatile double>);
    STATIC_REQUIRE(std::is_same_v<copy_cv_t<int, double>, double>);
    STATIC_REQUIRE(std::is_same_v<copy_cv_t<int, const double>, const double>);
    STATIC_REQUIRE(std::is_same_v<copy_cv_t<volatile int, const double>, const volatile double>);
    STATIC_REQUIRE(std::is_same_v<copy_cv_t<const int, const double>, const double>);
    STATIC_REQUIRE(std::is_same_v<copy_cv_t<const volatile int, const volatile double>, const volatile double>);
}
