/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include "algorithm/soga_selection.hpp"
#include "encoding/encoding.hpp"
#include "utility/math.hpp"
#include "test_utils.hpp"

using namespace gapp;
using namespace gapp::selection;

static constexpr size_t POPSIZE = 10;
static const BinaryGA context = []{ BinaryGA GA(POPSIZE); GA.solve(DummyFitnessFunction<BinaryGene>(10), 1); return GA; }();
static const Candidate<BinaryGene> sol = []{ Candidate<BinaryGene> c; c.fitness = { 0.0 }; return c; }();

TEST_CASE("roulette_selection", "[selection][single-objective]")
{
    Population<BinaryGene> pop(10, sol);

    std::unique_ptr<Selection> selection = std::make_unique<Roulette>();
    selection->initializeImpl(context);

    selection->prepareSelectionsImpl(context, pop);
    REQUIRE(selection->selectImpl(context, pop) < pop.size());

    pop[3].fitness = { math::small<double> };
    selection->prepareSelectionsImpl(context, pop);
    REQUIRE(selection->selectImpl(context, pop) == 3);

    pop[3].fitness = { math::large<double> };
    pop[4].fitness = { math::large<double> };
    selection->prepareSelectionsImpl(context, pop);
    const size_t idx = selection->selectImpl(context, pop);
    REQUIRE((idx == 3 || idx == 4));

    pop[0].fitness = { -math::large<double> };
    selection->prepareSelectionsImpl(context, pop);
    REQUIRE(selection->selectImpl(context, pop) < pop.size());

    pop[3].fitness = { -math::large<double> };
    pop[4].fitness = { -math::large<double> };
    selection->prepareSelectionsImpl(context, pop);
    REQUIRE(selection->selectImpl(context, pop) < pop.size());
}

TEST_CASE("tournament_selection", "[selection][single-objective]")
{
    Population<BinaryGene> pop(10, sol);

    std::unique_ptr<Selection> selection = std::make_unique<Tournament>();
    selection->initializeImpl(context);

    selection->prepareSelectionsImpl(context, pop);
    REQUIRE(selection->selectImpl(context, pop) < pop.size());

    pop[0].fitness = { -math::inf<double> };
    selection->prepareSelectionsImpl(context, pop);
    REQUIRE(selection->selectImpl(context, pop) < pop.size());

    pop[1].fitness = { math::large<double> };
    selection->prepareSelectionsImpl(context, pop);
    REQUIRE(selection->selectImpl(context, pop) < pop.size());

    pop[4].fitness = { math::inf<double> };
    selection->prepareSelectionsImpl(context, pop);
    REQUIRE(selection->selectImpl(context, pop) < pop.size());
}

TEST_CASE("rank_selection", "[selection][single-objective]")
{
    Population<BinaryGene> pop(10, sol);

    std::unique_ptr<Selection> selection = std::make_unique<Rank>(0.0);
    selection->initializeImpl(context);

    selection->prepareSelectionsImpl(context, pop);
    REQUIRE(selection->selectImpl(context, pop) < pop.size());

    pop[0].fitness = { -math::inf<double> };
    selection->prepareSelectionsImpl(context, pop);
    REQUIRE(selection->selectImpl(context, pop) != 0);

    pop[1].fitness = { math::large<double> };
    selection->prepareSelectionsImpl(context, pop);
    REQUIRE(selection->selectImpl(context, pop) < pop.size());

    pop[4].fitness = { math::inf<double> };
    selection->prepareSelectionsImpl(context, pop);
    REQUIRE(selection->selectImpl(context, pop) < pop.size());
}

TEST_CASE("sigma_selection" "[selection][single-objective]")
{
    Population<BinaryGene> pop(10, sol);

    std::unique_ptr<Selection> selection = std::make_unique<Sigma>();
    selection->initializeImpl(context);

    selection->prepareSelectionsImpl(context, pop);
    REQUIRE(selection->selectImpl(context, pop) < pop.size());

    pop[3].fitness = { math::small<double> };
    selection->prepareSelectionsImpl(context, pop);
    REQUIRE(selection->selectImpl(context, pop) < pop.size());

    pop[3].fitness = { math::large<double> };
    pop[4].fitness = { math::large<double> };
    selection->prepareSelectionsImpl(context, pop);
    REQUIRE(selection->selectImpl(context, pop) < pop.size());

    pop[0].fitness = { -math::large<double> };
    selection->prepareSelectionsImpl(context, pop);
    REQUIRE(selection->selectImpl(context, pop) < pop.size());

    pop[3].fitness = { -math::large<double> };
    pop[4].fitness = { -math::large<double> };
    selection->prepareSelectionsImpl(context, pop);
    REQUIRE(selection->selectImpl(context, pop) < pop.size());
}

TEST_CASE("boltzmann_selection" "[selection][single-objective]")
{
    Population<BinaryGene> pop(10, sol);

    std::unique_ptr<Selection> selection = std::make_unique<Boltzmann>();
    selection->initializeImpl(context);

    selection->prepareSelectionsImpl(context, pop);
    REQUIRE(selection->selectImpl(context, pop) < pop.size());

    pop[3].fitness = { math::small<double> };
    selection->prepareSelectionsImpl(context, pop);
    REQUIRE(selection->selectImpl(context, pop) < pop.size());

    pop[3].fitness = { math::large<double> };
    pop[4].fitness = { math::large<double> };
    selection->prepareSelectionsImpl(context, pop);
    REQUIRE(selection->selectImpl(context, pop) < pop.size());

    pop[0].fitness = { -math::large<double> };
    selection->prepareSelectionsImpl(context, pop);
    REQUIRE(selection->selectImpl(context, pop) < pop.size());

    pop[3].fitness = { -math::large<double> };
    pop[4].fitness = { -math::large<double> };
    selection->prepareSelectionsImpl(context, pop);
    REQUIRE(selection->selectImpl(context, pop) < pop.size());
}
