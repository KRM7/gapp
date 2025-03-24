/* Copyright (c) 2025 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include "core/fitness_function.hpp"
#include <ranges>

using namespace gapp;

TEST_CASE("simple_fitness_function", "[fitness_function]")
{
    detail::FitnessLambda<RealGene> f(3, [](const auto&)
    {
        return FitnessVector{ 0.0, 0.0 };
    });

    REQUIRE(f.chrom_lens()[0] == 3);
    REQUIRE(f.chrom_len<RealGene>() == 3);

    REQUIRE(!f.is_dynamic());

    REQUIRE(std::ranges::equal(f(Candidate<RealGene>{}), FitnessVector{ 0.0, 0.0 }));
}

TEST_CASE("mixed_fitness_function", "[fitness_function]")
{
    detail::FitnessLambda<MixedGene<RealGene, BinaryGene>> f({ 3, 2 }, [](const auto&)
    {
        return FitnessVector{ 0.0, 0.0 };
    });

    REQUIRE(std::ranges::equal(f.chrom_lens(), std::array{ 3, 2 }));

    REQUIRE(f.chrom_len<RealGene>() == 3);
    REQUIRE(f.chrom_len<BinaryGene>() == 2);

    REQUIRE(!f.is_dynamic());

    REQUIRE(std::ranges::equal(f(Candidate<MixedGene<RealGene, BinaryGene>>{}), FitnessVector{ 0.0, 0.0 }));
}
