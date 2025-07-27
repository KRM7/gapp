/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#include "algorithm/soga_replacement.hpp"
#include "encoding/encoding.hpp"
#include "utility/math.hpp"
#include "test_utils.hpp"

using namespace gapp;
using namespace gapp::replacement;

static constexpr size_t POPSIZE = 10;
static const BinaryGA context = []{ BinaryGA GA(POPSIZE); GA.solve(DummyFitnessFunction<BinaryGene>(10), 1); return GA; }();

static const Population<BinaryGene> pop = []
{
    Population<BinaryGene> p(20);
    // parents
    p[0].fitness = { math::inf<double> };
    p[1].fitness = { math::large<double> };
    p[2].fitness = { 0.0 };
    p[3].fitness = { -math::inf<double> };
    p[4].fitness = { 1.0 };
    p[5].fitness = { math::inf<double> };
    p[6].fitness = { 0.0 };
    p[7].fitness = { -math::large<double> };
    p[8].fitness = { -math::inf<double> };
    p[9].fitness = { math::small<double> };
    // children
    p[10].fitness = { 0.0 };
    p[11].fitness = { math::large<double> };
    p[12].fitness = { math::small<double> };
    p[13].fitness = { -1.0 };
    p[14].fitness = { -math::inf<double> };
    p[15].fitness = { 500.0 };
    p[16].fitness = { math::large<double> };
    p[17].fitness = { math::inf<double> };
    p[18].fitness = { 0.0 };
    p[19].fitness = { -math::large<double> };

    return p;
}();

TEST_CASE("replacement_best", "[replacement][single-objective]")
{
    math::ScopedTolerances _(0.0, 0.0);

    std::unique_ptr<Replacement> replacement = std::make_unique<KeepBest>();
    const auto indices = replacement->nextPopulationImpl(context, pop);

    const std::vector<size_t> expected = { 0, 1, 4, 5, 9, 11, 12, 15, 16, 17 };

    REQUIRE_THAT(indices.std_vec(), Catch::Matchers::UnorderedEquals(expected));
}

TEST_CASE("replacement_children", "[replacement][single-objective]")
{
    std::unique_ptr<Replacement> replacement = std::make_unique<KeepChildren>();
    const auto indices = replacement->nextPopulationImpl(context, pop);

    const std::vector<size_t> expected = { 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 };

    REQUIRE_THAT(indices.std_vec(), Catch::Matchers::UnorderedEquals(expected));
}

TEST_CASE("replacement_elitism", "[replacement][single-objective]")
{
    std::unique_ptr<Replacement> replacement = std::make_unique<Elitism>(2);
    const auto indices = replacement->nextPopulationImpl(context, pop);

    const std::vector<size_t> expected = { 0, 5, 10, 11, 12, 13, 14, 15, 16, 17 };

    REQUIRE_THAT(indices.std_vec(), Catch::Matchers::UnorderedEquals(expected));
}
