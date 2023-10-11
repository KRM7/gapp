/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#include "core/population.hpp"
#include "utility/utility.hpp"
#include <vector>

using namespace gapp;
using namespace gapp::detail;
using namespace gapp::math;
using namespace Catch;

TEST_CASE("find_pareto_front_1D", "[pareto_front]")
{
    FitnessMatrix fmat = { { 0.0 }, { -1.2 }, { 3.5 }, { 3.41 }, { 2.3 }, { 3.499 }, { -112.0 }, { 3.5 }, { 2.7 }, { 0.0 }, { 3.5 } };

    SECTION("single optimum")
    {
        fmat[1] = { 12.0 };

        auto optimal_indices = findParetoFront1D(fmat);

        REQUIRE(optimal_indices == std::vector{ 1_sz });
    }

    SECTION("multiple optimum")
    {
        ScopedTolerances _(0.0, 0.0);

        auto optimal_indices = findParetoFront1D(fmat);

        REQUIRE(optimal_indices == std::vector{ 2_sz, 7_sz, 10_sz });
    }

    SECTION("multiple optimum approx")
    {
        ScopedTolerances _(0.1, 0.0);

        auto optimal_indices = findParetoFront1D(fmat);

        REQUIRE(optimal_indices == std::vector{ 2_sz, 3_sz, 5_sz, 7_sz, 10_sz });
    }
}

TEMPLATE_TEST_CASE_SIG("find_pareto_front_nd", "[pareto_front]", ((auto F), F), findParetoFrontSort, findParetoFrontBest, findParetoFrontKung)
{
    FitnessMatrix fmat = {
        { 0.0,   0.0  },
        { 1.0,   2.0  },
        { 1.0,   3.0  },
        { 1.0,   4.0  },
        { 1.0,   4.97 }, // ~
        { 1.0,   5.0  }, //
        { 2.0,   3.0  },
        { 3.0,   1.0  },
        { 3.0,   2.0  },
        { 3.0,   3.0  }, //
        { 4.0,   1.0  },
        { 4.0,   2.0  }, //
        { 5.0,  -1.0  }, // ~
        { 5.01, -0.99 }, //
        { 5.01, -0.99 }
    };

    SECTION("single optimum")
    {
        fmat[3] = { 31.0, 7.0 };

        auto optimal_indices = F(fmat);

        REQUIRE(optimal_indices == std::vector{ 3_sz });
    }

    SECTION("multiple optimum")
    {
        ScopedTolerances _(0.0, 0.0);

        auto optimal_indices = F(fmat);

        REQUIRE_THAT(optimal_indices, Matchers::UnorderedEquals(std::vector<size_t>{ 5, 9, 11, 13, 14 }));
    }

    SECTION("multiple optimum approx")
    {
        ScopedTolerances _(0.1, 0.0);

        auto optimal_indices = F(fmat);

        REQUIRE_THAT(optimal_indices, Matchers::UnorderedEquals(std::vector<size_t>{ 4, 5, 9, 11, 12, 13, 14 }));
    }
}