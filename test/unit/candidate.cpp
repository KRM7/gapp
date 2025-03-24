/* Copyright (c) 2025 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include "core/candidate.hpp"
#include <ranges>
#include <unordered_map>

using namespace gapp;

TEST_CASE("bounds", "[candidate]")
{
    Bounds<RealGene> bounds{ 0.0, 1.0 };

    REQUIRE(bounds.lower() == 0.0);
    REQUIRE(bounds.upper() == 1.0);

    REQUIRE(bounds == bounds);
}

TEST_CASE("simple_encoding", "[candidate]")
{
    STATIC_REQUIRE(Candidate<BinaryGene>::NumChroms == 1);
    STATIC_REQUIRE(Candidate<RealGene>::NumChroms == 1);

    const BoundsVector<RealGene> bounds(4, { 0.0, 1.0 });

    Candidate<BinaryGene> unbounded_candidate({ 1, 1, 1, 1, 1 });
    Candidate<RealGene> bounded_candidate({ 0.0, 0.0, 0.0, 0.0 }, bounds);

    REQUIRE(unbounded_candidate.chrom_len() == 5);
    REQUIRE(bounded_candidate.chrom_len() == 4);

    REQUIRE(unbounded_candidate.chrom_len<BinaryGene>() == 5);
    REQUIRE(bounded_candidate.chrom_len<RealGene>() == 4);

    REQUIRE(unbounded_candidate.chrom() == unbounded_candidate.chromosome);
    REQUIRE(bounded_candidate.chrom() == bounded_candidate.chromosome);

    REQUIRE(unbounded_candidate.chrom<BinaryGene>() == Chromosome<BinaryGene>{ 1, 1, 1, 1, 1 });
    REQUIRE(bounded_candidate.chrom<RealGene>() == Chromosome<RealGene>{ 0.0, 0.0, 0.0, 0.0 });

    REQUIRE(std::ranges::equal(unbounded_candidate, unbounded_candidate.chromosome));
    REQUIRE(std::ranges::equal(bounded_candidate, bounded_candidate.chromosome));


    REQUIRE(std::ranges::equal(bounded_candidate.bounds(), bounds));
    REQUIRE(std::ranges::equal(bounded_candidate.bounds<RealGene>(), bounded_candidate.gene_bounds));


    REQUIRE(!unbounded_candidate.is_evaluated());
    REQUIRE(!bounded_candidate.is_evaluated());

    REQUIRE(unbounded_candidate.num_objectives() == 0);
    REQUIRE(bounded_candidate.num_objectives() == 0);

    unbounded_candidate.fitness = { 0.0, 0.0 };
    bounded_candidate.fitness = { 1.0 };

    REQUIRE(unbounded_candidate.is_evaluated());
    REQUIRE(bounded_candidate.is_evaluated());

    REQUIRE(unbounded_candidate.num_objectives() == 2);
    REQUIRE(bounded_candidate.num_objectives() == 1);


    REQUIRE(unbounded_candidate.num_constraints() == 0);
    REQUIRE(bounded_candidate.num_constraints() == 0);

    REQUIRE(!unbounded_candidate.has_constraint_violation());
    REQUIRE(!bounded_candidate.has_constraint_violation());

    unbounded_candidate.constraint_violation = { 1.0, 1.0 };
    bounded_candidate.constraint_violation = { 0.0 };

    REQUIRE(unbounded_candidate.num_constraints() == 2);
    REQUIRE(bounded_candidate.num_constraints() == 1);

    REQUIRE(unbounded_candidate.has_constraint_violation());
    REQUIRE(!bounded_candidate.has_constraint_violation());


    REQUIRE(unbounded_candidate == unbounded_candidate);
    REQUIRE(bounded_candidate == bounded_candidate);
}

TEST_CASE("mixed_encoding", "[candidate]")
{
    STATIC_REQUIRE(Candidate<MixedGene<BinaryGene, PermutationGene>>::NumChroms == 2);
    STATIC_REQUIRE(Candidate<MixedGene<BinaryGene, PermutationGene, IntegerGene>>::NumChroms == 3);

    const BoundsVector<RealGene> real_bounds(4, { 0.0, 1.0 });
    const BoundsVector<IntegerGene> int_bounds(3, { 1, 4 });

    const Chromosome<BinaryGene> bin_chrom{ 0, 1, 0, 1, 0 };
    const Chromosome<RealGene> real_chrom{ 0.0, 1.0, 1.0, 0.0 };
    const Chromosome<IntegerGene> int_chrom{ 1, 2, 1 };
    const Chromosome<PermutationGene> perm_chrom{ 0, 1, 2 };

    Candidate<MixedGene<BinaryGene, PermutationGene>> unbounded_candidate{ std::tuple{ bin_chrom, perm_chrom } };
    Candidate<MixedGene<RealGene, IntegerGene>> bounded_candidate{ { real_chrom, int_chrom }, { real_bounds, int_bounds } };
    Candidate<MixedGene<BinaryGene, RealGene>> partially_bounded_candidate{ { bin_chrom, real_chrom }, real_bounds };

    REQUIRE(unbounded_candidate.chrom_len<BinaryGene>() == bin_chrom.size());
    REQUIRE(unbounded_candidate.chrom_len<PermutationGene>() == perm_chrom.size());

    REQUIRE(bounded_candidate.chrom_len<RealGene>() == real_chrom.size());
    REQUIRE(bounded_candidate.chrom_len<IntegerGene>() == int_chrom.size());

    REQUIRE(partially_bounded_candidate.chrom_len<BinaryGene>() == bin_chrom.size());
    REQUIRE(partially_bounded_candidate.chrom_len<RealGene>() == real_chrom.size());

    REQUIRE(unbounded_candidate.chrom<BinaryGene>() == bin_chrom);
    REQUIRE(unbounded_candidate.chrom<PermutationGene>() == perm_chrom);

    REQUIRE(bounded_candidate.chrom<RealGene>() == real_chrom);
    REQUIRE(bounded_candidate.chrom<IntegerGene>() == int_chrom);

    REQUIRE(partially_bounded_candidate.chrom<BinaryGene>() == bin_chrom);
    REQUIRE(partially_bounded_candidate.chrom<RealGene>() == real_chrom);


    REQUIRE(std::ranges::equal(bounded_candidate.bounds<RealGene>(), real_bounds));
    REQUIRE(std::ranges::equal(bounded_candidate.bounds<IntegerGene>(), int_bounds));
    REQUIRE(std::ranges::equal(partially_bounded_candidate.bounds<RealGene>(), real_bounds));


    REQUIRE(!partially_bounded_candidate.is_evaluated());
    REQUIRE(partially_bounded_candidate.num_objectives() == 0);

    REQUIRE(partially_bounded_candidate.num_constraints() == 0);
    REQUIRE(!partially_bounded_candidate.has_constraint_violation());

    REQUIRE(static_cast<Candidate<BinaryGene>&>(unbounded_candidate) == static_cast<Candidate<BinaryGene>&>(partially_bounded_candidate));
    REQUIRE(static_cast<Candidate<RealGene>&>(bounded_candidate) == static_cast<Candidate<RealGene>&>(partially_bounded_candidate));


    REQUIRE(!static_cast<CandidateInfo&>(bounded_candidate).is_evaluated());
    REQUIRE(static_cast<CandidateInfo&>(bounded_candidate).fitness.empty());


    REQUIRE_NOTHROW(dynamic_cast<decltype(bounded_candidate)&>(static_cast<Candidate<RealGene>&>(bounded_candidate)));
    REQUIRE_NOTHROW(dynamic_cast<decltype(bounded_candidate)&>(static_cast<CandidateInfo&>(bounded_candidate)));
}

TEST_CASE("mixed_candidate_move", "[candidate]")
{
    Candidate<MixedGene<BinaryGene, PermutationGene>> c1;

    c1.fitness = FitnessVector(10, 1.0);
    c1.constraint_violation = CVVector(20, 0.0);

    Candidate<MixedGene<BinaryGene, PermutationGene>> c2 = std::move(c1);

    REQUIRE(c2.num_objectives() == 10);
    REQUIRE(c2.num_constraints() == 20);

    c1 = std::move(c2);

    REQUIRE(c1.num_objectives() == 10);
    REQUIRE(c1.num_constraints() == 20);
}

TEST_CASE("candidate_hash", "[candidate]")
{
    std::unordered_map<Candidate<RealGene>, int> map1;
    std::unordered_map<Candidate<MixedGene<BinaryGene, RealGene>>, int> map2;

    SUCCEED();
}
