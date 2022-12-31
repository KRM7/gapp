/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_all.hpp>
#include "problems/problems.hpp"
#include "utility/math.hpp"
#include "utility/rng.hpp"
#include <algorithm>
#include <vector>
#include <cstddef>


using namespace genetic_algorithm;
using namespace genetic_algorithm::problems;
using namespace genetic_algorithm::rng;
using namespace genetic_algorithm::math;
using namespace Catch::Matchers;

std::vector<double> randomPoint(const typename BenchmarkFunction<RealGene>::BoundsVec& bounds)
{
    std::vector<double> point;
    point.reserve(bounds.size());

    for (const auto& [lbound, ubound] : bounds) { point.push_back(randomReal(lbound, ubound)); }

    return point;
}

std::vector<char> randomPoint(const typename BenchmarkFunction<BinaryGene>::BoundsVec& bounds)
{
    std::vector<char> point(bounds.size());
    std::generate(point.begin(), point.end(), randomBool);

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

TEST_CASE("kursawe", "[problems]")
{
    REQUIRE_THROWS(Kursawe(0));
    REQUIRE_THROWS(Kursawe(1));

    const size_t var_count = GENERATE(2, 10, 100, 1000);

    INFO("Number of variables: " + std::to_string(var_count));

    Kursawe func(var_count);

    REQUIRE_THAT( func(func.optimum()), Approx(func.optimal_value()).margin(1E-6) );

    REQUIRE( !paretoCompareLess(func.ideal_point(), func.nadir_point()) );
    REQUIRE( !paretoCompareLess(func.optimal_value(), func.nadir_point()) );
    REQUIRE( !paretoCompareLess(func.ideal_point(), func.optimal_value()) );

    REQUIRE( func.bounds().size() == func.num_vars());

    const auto random_point = randomPoint(func.bounds());

    REQUIRE( !paretoCompareLess(func.optimal_value(), func(random_point)) );
    REQUIRE( !paretoCompareLess(func.ideal_point(), func(random_point)) );
}

TEMPLATE_TEST_CASE("zdt_suite", "[problems]", ZDT1, ZDT2, ZDT3, ZDT4, ZDT5, ZDT6)
{
    REQUIRE_THROWS(TestType(0));
    REQUIRE_THROWS(TestType(1));

    const size_t var_count = GENERATE(2, 3, 10, 100, 1000);

    INFO("Number of variables: " + std::to_string(var_count));

    TestType func(var_count);

    REQUIRE_THAT(func(func.optimum()), Approx(func.optimal_value()).margin(1E-6));

    REQUIRE( !paretoCompareLess(func.ideal_point(), func.nadir_point()) );
    REQUIRE( !paretoCompareLess(func.optimal_value(), func.nadir_point()) );
    REQUIRE( !paretoCompareLess(func.ideal_point(), func.optimal_value()) );

    REQUIRE( func.bounds().size() == func.num_vars() );

    const auto random_point = randomPoint(func.bounds());

    REQUIRE( !paretoCompareLess(func.optimal_value(), func(random_point)) );
    REQUIRE( !paretoCompareLess(func.ideal_point(), func(random_point)) );
}

TEMPLATE_TEST_CASE("dtlz_suite", "[problems]", DTLZ1, DTLZ2, DTLZ3, DTLZ4, DTLZ5, DTLZ6, DTLZ7)
{
    REQUIRE_THROWS(TestType(0));
    REQUIRE_THROWS(TestType(1));

    const size_t num_obj = GENERATE(2, 3, 10, 100, 1000);

    INFO("Number of objectives: " + std::to_string(num_obj));

    TestType func(num_obj);

    REQUIRE_THAT( func(func.optimum()), Approx(func.optimal_value()).margin(1E-6) );

    REQUIRE( !paretoCompareLess(func.ideal_point(), func.nadir_point()) );
    REQUIRE( !paretoCompareLess(func.optimal_value(), func.nadir_point()) );
    REQUIRE( !paretoCompareLess(func.ideal_point(), func.optimal_value()) );

    REQUIRE( func.bounds().size() == func.num_vars() );

    const auto random_point = randomPoint(func.bounds());

    REQUIRE( !paretoCompareLess(func.optimal_value(), func(random_point)) );
    REQUIRE( !paretoCompareLess(func.ideal_point(), func(random_point)) );
}
