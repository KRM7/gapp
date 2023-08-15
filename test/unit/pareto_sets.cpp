/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#include "population/population.hpp"
#include "utility/math.hpp"
#include <algorithm>

using namespace gapp;
using namespace gapp::math;
using namespace gapp::detail;
using namespace Catch;


static Population<int> fmatToPopulation(const FitnessMatrix& fmat)
{
    Population<int> pop(fmat.size(), Candidate<int>(0));

    for (size_t i = 0; i < fmat.size(); i++) { pop[i].fitness = FitnessVector(fmat[i]); }

    return pop;
}

template<typename T>
static constexpr bool fcomp(const Candidate<T>& lhs, const Candidate<T>& rhs) { return lhs.fitness == rhs.fitness; }


TEST_CASE("merge_pareto_sets", "[pareto_front]")
{
    ScopedTolerances _(0.0, 0.0);

    FitnessMatrix front1 = {
        { 10.0, -1.0 }, //
        {  8.0, 1.0  }, //
        {  7.0, 3.0  },
        {  6.0, 4.0  },
        {  5.0, 5.0  },
        {  4.0, 6.0  },
        {  3.0, 8.0  }, //
        {  2.0, 11.0 }, //
        {  1.0, 12.0 }, //
        { -1.0, 10.0 },
    };

    FitnessMatrix front2 = {
        {  9.0, -1.0 },
        {  8.0, 1.0  }, //
        {  7.0, 4.0  }, //
        {  6.0, 6.0  }, //
        {  5.0, 7.0  }, //
        {  3.0, 8.0  }, //
        {  2.0, 9.0  },
        {  1.0, 10.0 },
        { -1.0, 14.0 }, //
    };

    FitnessMatrix pareto_set = {
        { 10.0, -1.0 },
        {  8.0, 1.0  },
        {  8.0, 1.0  },
        {  7.0, 4.0  },
        {  6.0, 6.0  },
        {  5.0, 7.0  },
        {  3.0, 8.0  },
        {  3.0, 8.0  },
        {  2.0, 11.0 },
        {  1.0, 12.0 },
        { -1.0, 14.0 },
    };

    Population<int> pop1 = fmatToPopulation(front1);
    Population<int> pop2 = fmatToPopulation(front2);
    Population<int> expected_pop = fmatToPopulation(pareto_set);

    Population<int> optimal_pop = mergeParetoSets(pop1, pop2);

    REQUIRE(optimal_pop.size() == expected_pop.size());

    REQUIRE(std::is_permutation(optimal_pop.begin(), optimal_pop.end(), expected_pop.begin(), fcomp<int>));
}