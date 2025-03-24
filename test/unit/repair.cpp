/* Copyright (c) 2024 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <atomic>
#include "gapp.hpp"
#include "test_utils.hpp"

using namespace gapp;

TEST_CASE("repair_function", "[repair_function]")
{
    constexpr size_t population_size = 10;
    constexpr size_t generation_count = 5;
    constexpr size_t chromosome_length = 10;

    RCGA ga{ population_size };

    std::atomic<size_t> repair_count = 0;

    ga.repair_function([&](const GaInfo&, Candidate<RealGene>& sol)
    {
        for (RealGene& gene : sol) gene = std::max(0.0, gene);
        repair_count++;

        return true;
    });

    const auto solutions = ga.solve(DummyFitnessFunction<RealGene>{ chromosome_length }, Bounds{ -1.0, 1.0 }, generation_count);

    for (const auto& sol : solutions)
    {
        REQUIRE(std::all_of(sol.begin(), sol.end(), detail::greater_eq_than(0.0)));
    }

    REQUIRE(repair_count == population_size * generation_count);
}
