/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include "algorithm/soga_selection.hpp"
#include "encoding/encoding.hpp"
#include "utility/math.hpp"
#include "test_utils.hpp"

using namespace genetic_algorithm;
using namespace genetic_algorithm::selection;

const BinaryGA context{ DummyFitnessFunction<BinaryGene>(10), 10 };


TEST_CASE("roulette_selection", "[selection][single-objective]")
{
    FitnessMatrix fmat = {
        { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }
    };

    std::unique_ptr<Selection> selection = std::make_unique<Roulette>();
    selection->initializeImpl(context);

    selection->prepareSelectionsImpl(context, fmat);
    REQUIRE(selection->selectImpl(context, fmat) < fmat.size());

    fmat[3] = { math::small<double> };
    selection->prepareSelectionsImpl(context, fmat);
    REQUIRE(selection->selectImpl(context, fmat) == 3);

    fmat[3] = { math::large<double> };
    fmat[4] = { math::large<double> };
    selection->prepareSelectionsImpl(context, fmat);
    size_t idx = selection->selectImpl(context, fmat);
    REQUIRE((idx == 3 || idx == 4));

    fmat[0] = { -math::large<double> };
    selection->prepareSelectionsImpl(context, fmat);
    REQUIRE(selection->selectImpl(context, fmat) < fmat.size());

    fmat[3] = { -math::large<double> };
    fmat[4] = { -math::large<double> };
    selection->prepareSelectionsImpl(context, fmat);
    REQUIRE(selection->selectImpl(context, fmat) < fmat.size());
}

TEST_CASE("tournament_selection", "[selection][single-objective]")
{
    FitnessMatrix fmat = {
        { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }
    };

    std::unique_ptr<Selection> selection = std::make_unique<Tournament>();
    selection->initializeImpl(context);

    selection->prepareSelectionsImpl(context, fmat);
    REQUIRE(selection->selectImpl(context, fmat) < fmat.size());

    fmat[0] = { -math::inf<double> };
    selection->prepareSelectionsImpl(context, fmat);
    REQUIRE(selection->selectImpl(context, fmat) != 0);

    fmat[1] = { math::large<double> };
    selection->prepareSelectionsImpl(context, fmat);
    REQUIRE(selection->selectImpl(context, fmat) < fmat.size());

    fmat[4] = { math::inf<double> };
    selection->prepareSelectionsImpl(context, fmat);
    REQUIRE(selection->selectImpl(context, fmat) < fmat.size());
}

TEST_CASE("rank_selection", "[selection][single-objective]")
{
    FitnessMatrix fmat = {
        { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }
    };

    std::unique_ptr<Selection> selection = std::make_unique<Rank>();
    selection->initializeImpl(context);

    selection->prepareSelectionsImpl(context, fmat);
    REQUIRE(selection->selectImpl(context, fmat) < fmat.size());

    fmat[0] = { -math::inf<double> };
    selection->prepareSelectionsImpl(context, fmat);
    REQUIRE(selection->selectImpl(context, fmat) != 0);

    fmat[1] = { math::large<double> };
    selection->prepareSelectionsImpl(context, fmat);
    REQUIRE(selection->selectImpl(context, fmat) < fmat.size());

    fmat[4] = { math::inf<double> };
    selection->prepareSelectionsImpl(context, fmat);
    REQUIRE(selection->selectImpl(context, fmat) < fmat.size());
}

TEST_CASE("sigma_selection" "[selection][single-objective]")
{
    FitnessMatrix fmat = {
        { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }
    };

    std::unique_ptr<Selection> selection = std::make_unique<Sigma>();
    selection->initializeImpl(context);

    selection->prepareSelectionsImpl(context, fmat);
    REQUIRE(selection->selectImpl(context, fmat) < fmat.size());

    fmat[3] = { math::small<double> };
    selection->prepareSelectionsImpl(context, fmat);
    REQUIRE(selection->selectImpl(context, fmat) < fmat.size());

    fmat[3] = { math::large<double> };
    fmat[4] = { math::large<double> };
    selection->prepareSelectionsImpl(context, fmat);
    REQUIRE(selection->selectImpl(context, fmat) < fmat.size());

    fmat[0] = { -math::large<double> };
    selection->prepareSelectionsImpl(context, fmat);
    REQUIRE(selection->selectImpl(context, fmat) < fmat.size());

    fmat[3] = { -math::large<double> };
    fmat[4] = { -math::large<double> };
    selection->prepareSelectionsImpl(context, fmat);
    REQUIRE(selection->selectImpl(context, fmat) < fmat.size());
}

TEST_CASE("boltzmann_selection" "[selection][single-objective]")
{
    FitnessMatrix fmat = {
        { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }
    };

    std::unique_ptr<Selection> selection = std::make_unique<Boltzmann>();
    selection->initializeImpl(context);

    selection->prepareSelectionsImpl(context, fmat);
    REQUIRE(selection->selectImpl(context, fmat) < fmat.size());

    fmat[3] = { math::small<double> };
    selection->prepareSelectionsImpl(context, fmat);
    REQUIRE(selection->selectImpl(context, fmat) < fmat.size());

    fmat[3] = { math::large<double> };
    fmat[4] = { math::large<double> };
    selection->prepareSelectionsImpl(context, fmat);
    REQUIRE(selection->selectImpl(context, fmat) < fmat.size());

    fmat[0] = { -math::large<double> };
    selection->prepareSelectionsImpl(context, fmat);
    REQUIRE(selection->selectImpl(context, fmat) < fmat.size());

    fmat[3] = { -math::large<double> };
    fmat[4] = { -math::large<double> };
    selection->prepareSelectionsImpl(context, fmat);
    REQUIRE(selection->selectImpl(context, fmat) < fmat.size());
}