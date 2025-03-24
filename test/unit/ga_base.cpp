/* Copyright (c) 2025 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include "gapp.hpp"

using namespace gapp;

detail::FitnessLambda<RealGene> rf(3, [](const auto&) { return FitnessVector{ 0.0 }; });
detail::FitnessLambda<MixedGene<RealGene, BinaryGene>> mf({ 3, 4 }, [](const auto&) { return FitnessVector{ 0.0 }; });

TEST_CASE("ga_info_simple", "[ga_base]")
{
    RCGA ga;
    GaInfo& ga_info = ga;

    REQUIRE(ga_info.fitness_function() == nullptr);

    REQUIRE(ga_info.chrom_lens().empty());
    REQUIRE(ga_info.chrom_len<RealGene>() == 0);

    REQUIRE(ga_info.num_objectives() == 0);
    REQUIRE(ga_info.num_constraints() == 0);

    REQUIRE(ga_info.population_view().empty());
    REQUIRE(ga_info.fitness_matrix().empty());

    REQUIRE(ga_info.num_fitness_evals() == 0);
    REQUIRE(ga_info.generation_cntr() == 0);

    REQUIRE(ga_info.crossover_method<RealGene>());
    REQUIRE(ga_info.mutation_method<RealGene>());

    REQUIRE(!ga_info.crossover_method<BinaryGene>());
    REQUIRE(!ga_info.mutation_method<BinaryGene>());

    ga_info.crossover_method<RealGene>()->crossover_rate(0.3);
    REQUIRE(ga_info.crossover_method<RealGene>()->crossover_rate() == 0.3);

    ga_info.mutation_method<RealGene>()->mutation_rate(0.3);
    REQUIRE(ga_info.mutation_method<RealGene>()->mutation_rate() == 0.3);


    ga.solve(rf, Bounds<RealGene>{ 0.0, 1.0 }, 2);

    REQUIRE(ga_info.fitness_function());
    
    REQUIRE(ga_info.chrom_lens().size() == 1);
    REQUIRE(ga_info.chrom_len<RealGene>() == 3);

    REQUIRE(ga_info.num_objectives() == 1);
    REQUIRE(ga_info.num_constraints() == 0);

    REQUIRE(ga_info.population_view().size() == ga_info.population_size());
    REQUIRE(ga_info.fitness_matrix().size() == ga_info.population_size());

    REQUIRE(ga_info.num_fitness_evals() > 0);
    REQUIRE(ga_info.generation_cntr() == 1);
}

TEST_CASE("ga_info_mixed", "[ga_base]")
{
    MixedGA<RealGene, BinaryGene> ga;
    GaInfo& ga_info = ga;

    REQUIRE(ga_info.chrom_lens().empty());
    REQUIRE(ga_info.chrom_len<RealGene>() == 0);
    REQUIRE(ga_info.chrom_len<BinaryGene>() == 0);

    REQUIRE(ga_info.crossover_method<MixedGene<RealGene, BinaryGene>>());
    REQUIRE(ga_info.crossover_method<RealGene>());
    REQUIRE(ga_info.crossover_method<BinaryGene>());
    REQUIRE(!ga_info.crossover_method<IntegerGene>());

    REQUIRE(ga_info.mutation_method<MixedGene<RealGene, BinaryGene>>());
    REQUIRE(ga_info.mutation_method<RealGene>());
    REQUIRE(ga_info.mutation_method<BinaryGene>());
    REQUIRE(!ga_info.mutation_method<IntegerGene>());

    ga_info.mutation_method<BinaryGene>()->mutation_rate(0.1);
    REQUIRE(ga_info.mutation_method<BinaryGene>()->mutation_rate() == 0.1);

    ga.solve(mf, Bounds<RealGene>{ 0.0, 1.0 }, 2);

    REQUIRE(ga_info.chrom_lens().size() == 2);
    REQUIRE(ga_info.chrom_len<RealGene>() == 3);
    REQUIRE(ga_info.chrom_len<BinaryGene>() == 4);
}
