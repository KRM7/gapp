/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#include "algorithm/soga_replacement.hpp"
#include "encoding/encoding.hpp"
#include "utility/math.hpp"
#include "test_utils.hpp"

using namespace genetic_algorithm;
using namespace genetic_algorithm::replacement;


const size_t POPSIZE = 10;
const BinaryGA context{ DummyFitnessFunction<BinaryGene>(10), POPSIZE };

const FitnessMatrix fmat = {
    // parents
    { math::inf<double> },
    { math::large<double> },
    { 0.0 },
    { -math::inf<double> },
    { 1.0 },
    { math::inf<double> },
    { 0.0 },
    { -math::large<double> },
    { -math::inf<double> },
    { math::small<double> },
    // children
    { 0.0 },
    { math::large<double> },
    { math::small<double> },
    { -1.0 },
    { -math::inf<double> },
    { 500.0 },
    { math::large<double> },
    { math::inf<double> },
    { 0.0 },
    { -math::large<double> }
};


TEST_CASE("replacement_best", "[replacement][single-objective]")
{
    math::ScopedTolerances _(0, 0.0);

    std::unique_ptr<Replacement> replacement = std::make_unique<KeepBest>();
    const auto indices = replacement->nextPopulationImpl(context, fmat.begin(), fmat.begin() + POPSIZE, fmat.end());

    const std::vector<size_t> expected = { 0, 1, 4, 5, 9, 11, 12, 15, 16, 17 };

    REQUIRE_THAT(indices, Catch::Matchers::UnorderedEquals(expected));
}

TEST_CASE("replacement_children", "[replacement][single-objective]")
{
    std::unique_ptr<Replacement> replacement = std::make_unique<KeepChildren>();
    const auto indices = replacement->nextPopulationImpl(context, fmat.begin(), fmat.begin() + POPSIZE, fmat.end());

    const std::vector<size_t> expected = { 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 };

    REQUIRE_THAT(indices, Catch::Matchers::UnorderedEquals(expected));
}

TEST_CASE("replacement_elitism", "[replacement][single-objective]")
{
    std::unique_ptr<Replacement> replacement = std::make_unique<Elitism>(2);
    const auto indices = replacement->nextPopulationImpl(context, fmat.begin(), fmat.begin() + POPSIZE, fmat.end());

    const std::vector<size_t> expected = { 0, 5, 10, 11, 12, 13, 14, 15, 16, 17 };

    REQUIRE_THAT(indices, Catch::Matchers::UnorderedEquals(expected));
}
