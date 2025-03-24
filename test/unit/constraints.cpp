/* Copyright (c) 2024 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include "gapp.hpp"
#include "test_utils.hpp"

using namespace gapp;

constexpr size_t population_size = 10;
constexpr size_t generation_count = 5;
constexpr size_t chromosome_length = 10;

TEST_CASE("unconstrained_problem", "[constraints]")
{
    RCGA ga{ population_size };

    ga.repair_function([&](const GaInfo&, const Candidate<RealGene>& sol)
    {
        REQUIRE(sol.constraint_violation.empty());
        REQUIRE(!sol.has_constraint_violation());
        return false;
    });

    const auto solutions = ga.solve(DummyFitnessFunction<RealGene>{ chromosome_length }, Bounds{ -1.0, 1.0 }, generation_count);
    REQUIRE(solutions[0].constraint_violation.empty());

    REQUIRE(ga.num_constraints() == 0);
}

TEST_CASE("constrained_problem", "[constraints]")
{
    RCGA ga{ population_size };

    ga.constraints_function([](const GaInfo&, const Candidate<RealGene>&)
    {
        return CVVector{ 1.0, 0.0 };
    });

    ga.repair_function([&](const GaInfo&, Candidate<RealGene>& sol)
    {
        REQUIRE(sol.constraint_violation == CVVector{ 1.0, 0.0 });
        REQUIRE(sol.has_constraint_violation());
        return false;
    });

    const auto solutions = ga.solve(DummyFitnessFunction<RealGene>{ chromosome_length }, Bounds{ -1.0, 1.0 }, generation_count);
    REQUIRE(solutions[0].constraint_violation.size() == 2);

    REQUIRE(ga.num_constraints() == 2);
}
