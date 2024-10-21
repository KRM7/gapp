/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "core/candidate.hpp"
#include "mutation/mutation.hpp"
#include "encoding/encoding.hpp"
#include "utility/functional.hpp"
#include "test_utils.hpp"
#include <algorithm>
#include <numeric>

using namespace gapp;
using namespace gapp::mutation;

TEST_CASE("mutation_fitness_eval", "[crossover]")
{
    BinaryGA context;
    context.solve(DummyFitnessFunction<BinaryGene>(10), 1);

    binary::Flip mutation{ 0.0 };

    Candidate<BinaryGene> candidate{ { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } };
    candidate.fitness = { 0.0 };

    SECTION("evaluated, unchanged")
    {
        mutation.mutation_rate(0.0);

        mutation(context, candidate);

        REQUIRE(candidate.is_evaluated());
        REQUIRE(candidate.fitness == FitnessVector{ 0.0 });
    }

    SECTION("evaluated, changed")
    {
        mutation.mutation_rate(1.0);

        mutation(context, candidate);

        REQUIRE(!candidate.is_evaluated());
    }

    SECTION("unevaluated")
    {
        candidate.fitness.clear();
        mutation.mutation_rate(0.1);

        mutation(context, candidate);

        REQUIRE(!candidate.is_evaluated());
    }
}

TEMPLATE_TEST_CASE("binary_mutation", "[mutation]", binary::Flip)
{
    using Mutation = TestType;

    BinaryGA context;
    context.solve(DummyFitnessFunction<BinaryGene>(10), 1);

    Candidate<BinaryGene> candidate{ { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } };
    candidate.fitness = { 0.0 };

    Candidate<BinaryGene> old_candidate = candidate;

    SECTION("mutation probability = 0.0")
    {
        constexpr Mutation mutation{ 0.0 };

        mutation(context, candidate);

        REQUIRE(candidate.is_evaluated());
        REQUIRE(candidate.chromosome == old_candidate.chromosome);
        REQUIRE(candidate.fitness == old_candidate.fitness);
    }

    SECTION("mutation probability = 1.0")
    {
        constexpr Mutation mutation{ 1.0 };

        mutation(context, candidate);

        REQUIRE(!candidate.is_evaluated());
        REQUIRE(candidate.chromosome != old_candidate.chromosome);
    }

    REQUIRE(
        candidate.chromosome.size() == old_candidate.chromosome.size()
    );

    REQUIRE(
        std::all_of(candidate.chromosome.begin(), candidate.chromosome.end(), detail::between(0, 1))
    );
}

TEMPLATE_TEST_CASE("real_mutation", "[mutation]", real::Boundary, real::Gauss, real::NonUniform, real::Polynomial, real::Uniform)
{
    using Mutation = TestType;
    const Bounds bounds = { -1.0, 1.0 };

    RCGA context;
    context.solve(DummyFitnessFunction<RealGene>(10), bounds, 1);

    Candidate<RealGene> candidate{ { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } };
    candidate.fitness = { 0.0 };

    Candidate<RealGene> old_candidate = candidate;

    SECTION("mutation probability = 0.0")
    {
        constexpr Mutation mutation{ 0.0 };

        mutation(context, candidate);

        REQUIRE(candidate.is_evaluated());
        REQUIRE(candidate.chromosome == old_candidate.chromosome);
        REQUIRE(candidate.fitness == old_candidate.fitness);
    }

    SECTION("mutation probability = 1.0")
    {
        constexpr Mutation mutation{ 1.0 };

        mutation(context, candidate);

        REQUIRE(( candidate.is_evaluated() == (candidate.chromosome == old_candidate.chromosome)));
    }

    REQUIRE(
        candidate.chromosome.size() == old_candidate.chromosome.size()
    );

    REQUIRE(
        std::all_of(candidate.chromosome.begin(), candidate.chromosome.end(), detail::between(bounds.lower(), bounds.upper()))
    );
}

TEMPLATE_TEST_CASE("perm_mutation", "[mutation]", perm::Inversion, perm::Shift, perm::Shuffle, perm::Swap2, perm::Swap3)
{
    using Mutation = TestType;

    PermutationGA context;
    context.solve(DummyFitnessFunction<PermutationGene>(10), 1);

    Candidate<PermutationGene> candidate{ { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 } };
    candidate.fitness = { 0.0 };

    Candidate<PermutationGene> old_candidate = candidate;

    SECTION("mutation probability = 0.0")
    {
        constexpr Mutation mutation{ 0.0 };

        mutation(context, candidate);

        REQUIRE(candidate.is_evaluated());
        REQUIRE(candidate.chromosome == old_candidate.chromosome);
        REQUIRE(candidate.fitness == old_candidate.fitness);
    }

    SECTION("mutation probability = 1.0")
    {
        constexpr Mutation mutation{ 1.0 };

        mutation(context, candidate);

        if constexpr (!std::is_same_v<Mutation, perm::Shuffle>) // shuffle could return the same sequence
        {
            REQUIRE(!candidate.is_evaluated());
            REQUIRE(candidate.chromosome != old_candidate.chromosome);
        }
    }

    REQUIRE(
        candidate.chromosome.size() == old_candidate.chromosome.size()
    );

    REQUIRE(
        std::all_of(candidate.chromosome.begin(), candidate.chromosome.end(), detail::between(0, 9))
    );

    REQUIRE((
        std::sort(candidate.chromosome.begin(), candidate.chromosome.end()),
        std::unique(candidate.chromosome.begin(), candidate.chromosome.end()) == candidate.chromosome.end()
    ));
}

TEMPLATE_TEST_CASE("integer_mutation", "[mutation]", integer::Uniform)
{
    using Mutation = TestType;
    const Bounds<IntegerGene> bounds = { 0, 3 };

    IntegerGA context;
    context.solve(DummyFitnessFunction<IntegerGene>(10), bounds, 1);

    Candidate<IntegerGene> candidate{ { 0, 1, 2, 3, 3, 1, 0, 1, 0, 2 } };
    candidate.fitness = { 0.0 };

    Candidate<IntegerGene> old_candidate = candidate;

    SECTION("mutation probability = 0.0")
    {
        constexpr Mutation mutation{ 0.0 };

        mutation(context, candidate);

        REQUIRE(candidate.is_evaluated());
        REQUIRE(candidate.chromosome == old_candidate.chromosome);
        REQUIRE(candidate.fitness == old_candidate.fitness);
    }

    SECTION("mutation probability = 1.0")
    {
        constexpr Mutation mutation{ 1.0 };

        mutation(context, candidate);

        REQUIRE(!candidate.is_evaluated());
        REQUIRE(candidate.chromosome != old_candidate.chromosome);
    }

    REQUIRE(
        candidate.chromosome.size() == old_candidate.chromosome.size()
    );

    REQUIRE(
        std::all_of(candidate.chromosome.begin(), candidate.chromosome.end(), detail::between(bounds.lower(), bounds.upper()))
    );
}
