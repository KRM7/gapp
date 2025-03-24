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

TEMPLATE_TEST_CASE("binary_mutation", "[mutation]", binary::Flip)
{
    using Mutation = TestType;

    BinaryGA context;
    context.solve(DummyFitnessFunction<BinaryGene>(10), 1);

    Candidate<BinaryGene> candidate{ { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } };
    candidate.fitness = { 0.0 };

    const Candidate<BinaryGene> old_candidate = candidate;

    SECTION("mutation probability = 0.0")
    {
        Mutation mutation{ 0.0 };

        mutation(context, candidate);

        REQUIRE(candidate.chromosome == old_candidate.chromosome);
    }

    SECTION("mutation probability = 1.0")
    {
        Mutation mutation{ 1.0 };

        mutation(context, candidate);

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

    constexpr size_t chrom_len = 10;
    const BoundsVector<RealGene> bounds(chrom_len, { -1.0, 1.0 });

    RCGA context;
    context.solve(DummyFitnessFunction<RealGene>(chrom_len), bounds, 1);

    Candidate<RealGene> candidate{ { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, bounds };
    candidate.fitness = { 0.0 };

    const Candidate<RealGene> old_candidate = candidate;

    SECTION("mutation probability = 0.0")
    {
        Mutation mutation{ 0.0_p };

        mutation(context, candidate);

        REQUIRE(candidate.chromosome == old_candidate.chromosome);
    }

    SECTION("mutation probability = 1.0")
    {
        Mutation mutation{ 1.0_p };

        mutation(context, candidate);
    }

    REQUIRE(
        candidate.chromosome.size() == old_candidate.chromosome.size()
    );

    REQUIRE(
        std::all_of(candidate.chromosome.begin(), candidate.chromosome.end(), detail::between(-1.0, 1.0))
    );
}

TEMPLATE_TEST_CASE("perm_mutation", "[mutation]", perm::Inversion, perm::Shift, perm::Shuffle, perm::Swap2, perm::Swap3)
{
    using Mutation = TestType;

    PermutationGA context;
    context.solve(DummyFitnessFunction<PermutationGene>(10), 1);

    Candidate<PermutationGene> candidate{ { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 } };
    candidate.fitness = { 0.0 };

    const Candidate<PermutationGene> old_candidate = candidate;

    SECTION("mutation probability = 0.0")
    {
        Mutation mutation{ 0.0 };

        mutation(context, candidate);

        REQUIRE(candidate.chromosome == old_candidate.chromosome);
    }

    SECTION("mutation probability = 1.0")
    {
        Mutation mutation{ 1.0 };

        mutation(context, candidate);

        // Shuffle mutation could result in the same sequence as the input
        if constexpr (!std::is_same_v<Mutation, perm::Shuffle>)
        {
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

    constexpr size_t chrom_len = 10;
    const BoundsVector<IntegerGene> bounds(chrom_len, { 0, 3 });

    IntegerGA context;
    context.solve(DummyFitnessFunction<IntegerGene>(chrom_len), bounds, 1);

    Candidate<IntegerGene> candidate{ { 0, 1, 2, 3, 3, 1, 0, 1, 0, 2 }, bounds };
    candidate.fitness = { 0.0 };

    const Candidate<IntegerGene> old_candidate = candidate;

    SECTION("mutation probability = 0.0")
    {
        Mutation mutation{ 0.0 };

        mutation(context, candidate);

        REQUIRE(candidate.chromosome == old_candidate.chromosome);
    }

    SECTION("mutation probability = 1.0")
    {
        Mutation mutation{ 1.0 };

        mutation(context, candidate);

        REQUIRE(candidate.chromosome != old_candidate.chromosome);
    }

    REQUIRE(
        candidate.chromosome.size() == old_candidate.chromosome.size()
    );

    REQUIRE(
        std::all_of(candidate.chromosome.begin(), candidate.chromosome.end(), detail::between(0, 3))
    );
}

TEST_CASE("mixed_mutation", "[mutation]")
{
    Mixed mutation{ binary::Flip{ 0.0 }, real::Boundary{ 0.0 } };

    REQUIRE(mutation.mutation_rates() == std::array<Probability, 2>{ 0.0, 0.0 });

    mutation.mutation_rates({ 0.1, 0.2 });
    REQUIRE(mutation.mutation_rates() == std::array<Probability, 2>{ 0.1, 0.2 });

    mutation.mutation_rates(1.0);
    REQUIRE(mutation.mutation_rates() == std::array<Probability, 2>{ 1.0, 1.0 });

    mutation.mutation_rate<RealGene>(0.5);
    REQUIRE(mutation.mutation_rate<RealGene>() == 0.5);

    REQUIRE(mutation.allow_variable_chrom_length<BinaryGene>());
    REQUIRE(!mutation.allow_variable_chrom_length<RealGene>());

    REQUIRE(mutation.component<BinaryGene>().mutation_rate() == 1.0);
    REQUIRE(mutation.component<RealGene>().mutation_rate() == 0.5);


    const std::array chrom_lens = { 3_sz, 4_sz };
    const BoundsVector<RealGene> bounds(chrom_lens[1], { 0.0,1.0 });

    MixedGA<BinaryGene, RealGene> context;
    context.solve(DummyFitnessFunction<MixedGene<BinaryGene, RealGene>>(chrom_lens), bounds, 1);

    Candidate<MixedGene<BinaryGene, RealGene>> candidate{ { Chromosome<BinaryGene>{ 0, 0, 0 }, Chromosome<RealGene>{ 0.0, 0.0, 0.0, 0.0 } }, bounds };

    mutation(context, candidate);

    REQUIRE(candidate.chrom_len<BinaryGene>() == chrom_lens[0]);
    REQUIRE(candidate.chrom_len<RealGene>() == chrom_lens[1]);

    REQUIRE(std::all_of(candidate.chrom<BinaryGene>().begin(), candidate.chrom<BinaryGene>().end(), detail::between(0, 1)));
    REQUIRE(std::all_of(candidate.chrom<RealGene>().begin(), candidate.chrom<RealGene>().end(), detail::between(0.0, 1.0)));

    REQUIRE(detail::equal(candidate.bounds<RealGene>(), BoundsView<RealGene>(bounds)));
}

TEST_CASE("mutation_bounds", "[mutation]")
{
    const BoundsVector<RealGene> bounds(10, { 0.0, 1.0 });

    RCGA context;
    context.solve(DummyFitnessFunction<RealGene>(10), bounds, 1);

    real::Boundary mutation;
    Candidate<RealGene> candidate{ { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, bounds };

    mutation(context, candidate);

    REQUIRE(detail::equal(candidate.bounds<RealGene>(), BoundsView<RealGene>(bounds)));
}

TEST_CASE("mutation_lambda", "[mutation]")
{
    mutation::Lambda<BinaryGene> mutation([](const GaInfo&, const auto&, auto&)
    {
        return;
    });

    mutation.mutation_rate(0.1);
    REQUIRE(mutation.mutation_rate() == 0.1);
}

TEST_CASE("mutation_callable", "[mutation]")
{
    RCGA ga;
    ga.mutation_method([](const GaInfo&, const Candidate<RealGene>&, Chromosome<RealGene>&) {});

    REQUIRE(ga.mutation_method().use_default_mutation_rate());
}
