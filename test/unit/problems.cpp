﻿/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_all.hpp>
#include "problems/problems.hpp"
#include "utility/math.hpp"
#include "utility/rng.hpp"
#include <algorithm>
#include <vector>
#include <cstddef>


using namespace gapp;
using namespace gapp::problems;
using namespace gapp::rng;
using namespace gapp::math;
using namespace Catch::Matchers;

Candidate<RealGene> randomSolution(const BoundsVector<RealGene>& bounds)
{
    Chromosome<RealGene> point;
    point.reserve(bounds.size());

    for (const auto& bound : bounds) { point.push_back( randomReal(bound.lower(), bound.upper()) ); }

    return Candidate<RealGene>{ std::move(point) };
}

Candidate<BinaryGene> randomSolution(const BoundsVector<BinaryGene>& bounds)
{
    Chromosome<BinaryGene> point(bounds.size());
    std::generate(point.begin(), point.end(), randomBool);

    return Candidate<BinaryGene>{ std::move(point) };
}


TEMPLATE_TEST_CASE("single_objective_problems", "[problems]", Sphere, Rastrigin, Rosenbrock, Schwefel, Griewank, Ackley, Levy)
{
    const size_t var_count = GENERATE(1, 10, 100, 1000);

    INFO("Number of variables: " + std::to_string(var_count));

    TestType func(var_count);

    REQUIRE_THAT( func(func.optimum()).std_vec(), Approx(func.optimal_value().std_vec()).margin(1E-6) );

    const auto random_sol = randomSolution(func.bounds());

    REQUIRE( func.bounds().size() == var_count );
    REQUIRE( !paretoCompareLess(func.optimal_value(), func(random_sol)) );
}

TEST_CASE("kursawe", "[problems]")
{
    const size_t var_count = GENERATE(2, 10, 100, 1000);

    INFO("Number of variables: " + std::to_string(var_count));

    Kursawe func(var_count);

    REQUIRE_THAT( func(func.optimum()).std_vec(), Approx(func.optimal_value().std_vec()).margin(1E-6));

    REQUIRE( !paretoCompareLess(func.ideal_point(), func.nadir_point()) );
    REQUIRE( !paretoCompareLess(func.optimal_value(), func.nadir_point()) );
    REQUIRE( !paretoCompareLess(func.ideal_point(), func.optimal_value()) );

    REQUIRE( func.bounds().size() == func.num_vars());

    const auto random_sol = randomSolution(func.bounds());

    REQUIRE( !paretoCompareLess(func.optimal_value(), func(random_sol)) );
    REQUIRE( !paretoCompareLess(func.ideal_point(), func(random_sol)) );
}

TEMPLATE_TEST_CASE("zdt_suite", "[problems]", ZDT1, ZDT2, ZDT3, ZDT4, ZDT5, ZDT6)
{
    const size_t var_count = GENERATE(2, 3, 10, 100, 1000);

    INFO("Number of variables: " + std::to_string(var_count));

    TestType func(var_count);

    REQUIRE_THAT(func(func.optimum()).std_vec(), Approx(func.optimal_value().std_vec()).margin(1E-6));

    REQUIRE( !paretoCompareLess(func.ideal_point(), func.nadir_point()) );
    REQUIRE( !paretoCompareLess(func.optimal_value(), func.nadir_point()) );
    REQUIRE( !paretoCompareLess(func.ideal_point(), func.optimal_value()) );

    REQUIRE( func.bounds().size() == func.num_vars() );

    const auto random_sol = randomSolution(func.bounds());

    REQUIRE( !paretoCompareLess(func.optimal_value(), func(random_sol)) );
    REQUIRE( !paretoCompareLess(func.ideal_point(), func(random_sol)) );
}

TEMPLATE_TEST_CASE("dtlz_suite", "[problems]", DTLZ1, DTLZ2, DTLZ3, DTLZ4, DTLZ5, DTLZ6, DTLZ7)
{
    const size_t num_obj = GENERATE(2, 3, 10, 100, 1000);

    INFO("Number of objectives: " + std::to_string(num_obj));

    TestType func(num_obj);

    REQUIRE_THAT( func(func.optimum()).std_vec(), Approx(func.optimal_value().std_vec()).margin(1E-6) );

    REQUIRE( !paretoCompareLess(func.ideal_point(), func.nadir_point()) );
    REQUIRE( !paretoCompareLess(func.optimal_value(), func.nadir_point()) );
    REQUIRE( !paretoCompareLess(func.ideal_point(), func.optimal_value()) );

    REQUIRE( func.bounds().size() == func.num_vars() );

    const auto random_sol = randomSolution(func.bounds());

    REQUIRE( !paretoCompareLess(func.optimal_value(), func(random_sol)) );
    REQUIRE( !paretoCompareLess(func.ideal_point(), func(random_sol)) );
}
