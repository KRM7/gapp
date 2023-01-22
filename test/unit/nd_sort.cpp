/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#include "algorithm/nd_sort.hpp"
#include "utility/math.hpp"
#include "utility/utility.hpp"
#include <vector>
#include <algorithm>
#include <cstddef>
#include <limits>

using namespace genetic_algorithm;
using namespace genetic_algorithm::algorithm::dtl;
using namespace genetic_algorithm::detail;
using namespace Catch;

static constexpr double inf = std::numeric_limits<double>::infinity();

TEMPLATE_TEST_CASE_SIG("nd_sort", "[pareto_front]", ((auto F), F), fastNonDominatedSort, dominanceDegreeSort)
{
    const FitnessMatrix fmat = {
        {  3.0, 3.0  },  // p2 - 0
        {  4.0, 4.0  },  // p1 - 1
        {  5.0, 6.0  },  // p0 - 2
        {  7.0, 2.0  },  // p0 - 3
        {  4.0, 2.0  },  // p2 - 4
        {  1.0, 4.0  },  // p2 - 5
        {  1.0, 2.0  },  // p4 - 6
        {  7.0, 0.0  },  // p1 - 7
        {  2.0, 2.0  },  // p3 - 8
        {  1.0, 6.0  },  // p1 - 9
        {  6.0, 4.0  },  // p0 - 10
        {  3.0, 1.0  },  // p3 - 11
        {  3.0, 7.0  },  // p0 - 12
        {  1.0, 1.0  },  // p5 - 13
        {  2.0, 5.0  },  // p1 - 14
        {  6.0, 1.0  },  // p1 - 15
        { -1.0, 0.0  },  // p6 - 16
        { -2.0, inf  },  // p0 - 17
        {  2.9, 0.9  },  // p4 - 18 ~
    };

    ParetoFronts pareto_fronts = F(fmat.begin(), fmat.end());

    ParetoFronts expected_fronts = {
        { 0,  2 },
        { 1,  1 },
        { 2,  0 },
        { 3,  0 },
        { 4,  2 },
        { 5,  2 },
        { 6,  4 },
        { 7,  1 },
        { 8,  3 },
        { 9,  1 },
        { 10, 0 },
        { 11, 3 },
        { 12, 0 },
        { 13, 5 },
        { 14, 1 },
        { 15, 1 },
        { 16, 6 },
        { 17, 0 },
        { 18, 4 }
    };

    REQUIRE_THAT(pareto_fronts, Matchers::UnorderedEquals(expected_fronts));

    REQUIRE(pareto_fronts.front().rank == 0);
    REQUIRE(pareto_fronts.back().rank  == 6);

    REQUIRE(std::is_sorted(pareto_fronts.begin(), pareto_fronts.end(), [](auto lhs, auto rhs) { return lhs.rank < rhs.rank; }));

    REQUIRE(std::adjacent_find(pareto_fronts.begin(), pareto_fronts.end(), [](auto lhs, auto rhs) { return (rhs.rank - lhs.rank) > 1; }) == pareto_fronts.end());


    math::ScopedTolerances _(0, 0.11);

    ParetoFronts pareto_fronts_approx = F(fmat.begin(), fmat.end());
    expected_fronts[18].rank = 3;

    REQUIRE_THAT(pareto_fronts_approx, Matchers::UnorderedEquals(expected_fronts));
}