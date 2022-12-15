/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_all.hpp>
#include "problems/problems.hpp"
#include "utility/math.hpp"
#include "utility/rng.hpp"
#include <algorithm>
#include <vector>
#include <cstddef>


using namespace genetic_algorithm::problems;
using namespace genetic_algorithm::rng;
using namespace genetic_algorithm::math;
using namespace Catch::Matchers;

std::vector<double> randomPoint(const auto& bounds)
{
    std::vector<double> point;
    point.reserve(bounds.size());

    for (const auto& [lbound, ubound] : bounds) { point.push_back(randomReal(lbound, ubound)); }

    return point;
}


TEMPLATE_TEST_CASE("single_objective_problems", "[problems]", Sphere, Rastrigin, Rosenbrock, Schwefel, Griewank, Ackley, Levy)
{
    REQUIRE_THROWS(TestType(0));

    const size_t var_count = GENERATE(1, 10, 100, 1000);

    INFO("Number of variables: " + std::to_string(var_count));

    TestType func(var_count);

    REQUIRE_THAT( func(func.optimum()), Approx(func.optimal_value()).margin(1E-6) );

    const auto random_point = randomPoint(func.bounds());

    REQUIRE( func.bounds().size() == var_count );
    REQUIRE( !paretoCompareLess(func.optimal_value(), func(random_point)) );
}

// TODO kursawe

TEMPLATE_TEST_CASE("zdt_suite", "[problems]", ZDT1, ZDT2, ZDT3, ZDT4, ZDT6)
{
    // TODO ideal point, nadir point?
    // TODO optimum, optimal point members

    REQUIRE_THROWS(TestType(0));

    const size_t var_count = GENERATE(2, 10, 100, 1000); // TODO fix 1 var

    INFO("Number of variables: " + std::to_string(var_count));

    TestType func(var_count);

    const auto optimum = std::vector(func.num_vars(), 0.0);
    const auto random_point = randomPoint(func.bounds());

    REQUIRE( func.bounds().size() == var_count );
    REQUIRE( !paretoCompareLess(func(optimum), func(random_point)) );
}

// TODO zdt5

TEMPLATE_TEST_CASE("dtlz_suite", "[problems]", DTLZ1, DTLZ2, DTLZ3, DTLZ4, DTLZ5, DTLZ6)
{
    REQUIRE_THROWS(TestType(0));
    REQUIRE_THROWS(TestType(1));

    const size_t num_obj = GENERATE(2, 10, 100, 1000);

    INFO("Number of objectives: " + std::to_string(num_obj));

    TestType func(num_obj);

    std::vector optimum(func.num_vars(), 0.5);
    std::fill(optimum.begin(), optimum.begin() + num_obj, 0.0);

    const auto random_point = randomPoint(func.bounds());

    REQUIRE( !paretoCompareLess(func(optimum), func(random_point)) );
}

// TODO dtlz7