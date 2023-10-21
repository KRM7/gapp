/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include "core/candidate.hpp"
#include "encoding/encoding.hpp"
#include "crossover/crossover.hpp"
#include "crossover/crossover_impl.hpp"
#include "test_utils.hpp"

using namespace gapp;
using namespace gapp::crossover;
using namespace gapp::crossover::dtl;

TEST_CASE("single_point_crossover", "[crossover]")
{
    Candidate<int> parent1{ { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } };
    Candidate<int> parent2{ { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 } };
                                            //

    auto [child1, child2] = singlePointCrossoverImpl(parent1, parent2, 5);

    REQUIRE(child1.chromosome == Chromosome<int>{ { 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0 } });
    REQUIRE(child2.chromosome == Chromosome<int>{ { 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1 } });
}

TEST_CASE("two_point_crossover", "[crossover]")
{
    Candidate<int> parent1{ { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } };
    Candidate<int> parent2{ { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 } };
                                       //                //

    auto [child1, child2] = twoPointCrossoverImpl(parent1, parent2, { 9, 3 });

    REQUIRE(child1.chromosome == Chromosome<int>{ { 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0 } });
    REQUIRE(child2.chromosome == Chromosome<int>{ { 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1 } });
}

TEST_CASE("npoint_crossover", "[crossover]")
{
    Candidate<char> parent1{ { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } };
    Candidate<char> parent2{ { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 } };
                                 //    //          //          //

    auto [child1, child2] = nPointCrossoverImpl(parent1, parent2, { 1, 3, 7, 11 });

    REQUIRE(child1.chromosome == Chromosome<char>{ { 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1 } });
    REQUIRE(child2.chromosome == Chromosome<char>{ { 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0 } });
}


TEMPLATE_TEST_CASE("order1_crossover", "[crossover]", int, unsigned)
{
    Candidate<TestType> parent1{ { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 } };
    Candidate<TestType> parent2{ { 4, 5, 0, 6, 1, 2, 8, 3, 9, 7 } };

    auto child1 = order1CrossoverImpl(parent1, parent2, 4, 8);
    auto child2 = order1CrossoverImpl(parent2, parent1, 4, 8);

    REQUIRE(child1.chromosome == Chromosome<TestType>{ { 1, 2, 8, 3, 4, 5, 6, 7, 9, 0 } });
    REQUIRE(child2.chromosome == Chromosome<TestType>{ { 4, 5, 6, 7, 1, 2, 8, 3, 9, 0 } });
}

TEMPLATE_TEST_CASE("order2_crossover", "[crossover]", int, unsigned)
{
    Candidate<TestType> parent1{ { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 } };
    Candidate<TestType> parent2{ { 4, 5, 0, 6, 1, 2, 8, 3, 9, 7 } };

    auto child1 = order2CrossoverImpl(parent1, parent2, 4, 8);
    auto child2 = order2CrossoverImpl(parent2, parent1, 4, 8);

    REQUIRE(child1.chromosome == Chromosome<TestType>{ { 0, 1, 2, 8, 4, 5, 6, 7, 3, 9 } });
    REQUIRE(child2.chromosome == Chromosome<TestType>{ { 0, 4, 5, 6, 1, 2, 8, 3, 7, 9 } });
}

TEMPLATE_TEST_CASE("position_crossover", "[crossover]", int, unsigned)
{
    Candidate<TestType> parent1{ { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 } };
    Candidate<TestType> parent2{ { 4, 5, 0, 6, 1, 2, 8, 3, 9, 7 } };

    auto child1 = positionCrossoverImpl(parent1, parent2, { { 0, 3, 4, 7 } });
    auto child2 = positionCrossoverImpl(parent2, parent1, { { 0, 3, 4, 7 } });

    REQUIRE(child1.chromosome == Chromosome<TestType>{ { 0, 5, 6, 3, 4, 1, 2, 7, 8, 9 } });
    REQUIRE(child2.chromosome == Chromosome<TestType>{ { 4, 0, 2, 6, 1, 5, 7, 3, 8, 9 } });
}

TEST_CASE("cycle_crossover", "[crossover]")
{
    Candidate<int> parent1{ { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 } };
    Candidate<int> parent2{ { 4, 5, 0, 6, 1, 2, 8, 3, 9, 7 } };

    // cycle0 : 0 - 4 - 1 - 5 - 2 , cycle1 : 3 - 6 - 8 - 9 - 7

    auto [child1, child2] = cycleCrossoverImpl(parent1, parent2);

    REQUIRE(child1.chromosome == Chromosome<int>{ { 0, 1, 2, 6, 4, 5, 8, 3, 9, 7 } });
    REQUIRE(child2.chromosome == Chromosome<int>{ { 4, 5, 0, 3, 1, 2, 6, 7, 8, 9 } });
}

TEMPLATE_TEST_CASE("edge_crossover", "[crossover]", int, unsigned)
{
    Candidate<TestType> parent1{ { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 } };
    Candidate<TestType> parent2{ { 4, 5, 0, 6, 1, 2, 8, 3, 9, 7 } };

    auto child1 = edgeCrossoverImpl(parent1, parent2);
    auto child2 = edgeCrossoverImpl(parent2, parent1);

    REQUIRE(child1.chromosome == Chromosome<TestType>{ { 0, 5, 4, 1, 6, 7, 9, 8, 3, 2 } });
    REQUIRE(child2.chromosome == Chromosome<TestType>{ { 4, 5, 0, 6, 1, 2, 3, 9, 8, 7 } });
}

TEMPLATE_TEST_CASE("pmx_crossover", "[crossover]", int, unsigned)
{
    Candidate<TestType> parent1{ { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 } };
    Candidate<TestType> parent2{ { 4, 5, 0, 6, 1, 2, 8, 3, 9, 7 } };

    auto child1 = pmxCrossoverImpl(parent1, parent2, 4, 8);
    auto child2 = pmxCrossoverImpl(parent2, parent1, 4, 8);

    REQUIRE(child1.chromosome == Chromosome<TestType>{ { 1, 2, 0, 8, 4, 5, 6, 7, 9, 3 } });
    REQUIRE(child2.chromosome == Chromosome<TestType>{ { 0, 4, 5, 7, 1, 2, 8, 3, 6, 9 } });
}

TEMPLATE_TEST_CASE("real_crossover", "[crossover]", real::Arithmetic, real::BLXa, real::SimulatedBinary, real::Wright)
{
    using Crossover = TestType;

    constexpr size_t chrom_len = 10;
    constexpr Bounds<RealGene> bounds = { 0.0, 1.0 };

    RCGA context;
    context.solve(DummyFitnessFunction<RealGene>(chrom_len), bounds, 1);

    constexpr Crossover crossover{ 0.8 };

    Candidate<RealGene> parent1{ { 0.0, 0.12, 0.48, 0.19, 1.0, 1.0, 0.0, 0.72, 0.81, 0.03 } };
    Candidate<RealGene> parent2{ { 1.0, 0.34, 0.97, 0.36, 1.0, 0.0, 0.0, 0.28, 0.49, 0.79 } };
    parent1.fitness = { 0.0 }; parent1.is_evaluated = true;
    parent2.fitness = { 0.0 }; parent2.is_evaluated = true;

    auto [child1, child2] = crossover(context, parent1, parent2);

    REQUIRE(child1.chromosome.size() == chrom_len);
    REQUIRE(child2.chromosome.size() == chrom_len);
    REQUIRE(std::all_of(child1.chromosome.begin(), child1.chromosome.end(), detail::between(bounds.lower(), bounds.upper())));
    REQUIRE(std::all_of(child2.chromosome.begin(), child2.chromosome.end(), detail::between(bounds.lower(), bounds.upper())));
}

TEST_CASE("crossover_fitness_eval", "[crossover]")
{
    BinaryGA context;
    context.solve(DummyFitnessFunction<BinaryGene>(10), 1);

    binary::SinglePoint crossover;

    Candidate<BinaryGene> parent1{ { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } };
    Candidate<BinaryGene> parent2{ { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 } };
    parent1.fitness = { 0.0 }; parent1.is_evaluated = true;
    parent2.fitness = { 0.0 }; parent2.is_evaluated = true;

    SECTION("unchanged chromosomes")
    {
        crossover.crossover_rate(0.0);

        auto [child1, child2] = crossover(context, parent1, parent2);

        REQUIRE(child1.is_evaluated);
        REQUIRE(child2.is_evaluated);
        REQUIRE(child1.fitness == parent1.fitness);
        REQUIRE(child2.fitness == parent2.fitness);
    }

    SECTION("changed chromosomes")
    {
        crossover.crossover_rate(1.0);

        auto [child1, child2] = crossover(context, parent1, parent2);

        REQUIRE((!child1.is_evaluated || child1.chromosome == parent1.chromosome || child1.chromosome == parent2.chromosome));
        REQUIRE((!child2.is_evaluated || child2.chromosome == parent1.chromosome || child2.chromosome == parent2.chromosome));
        REQUIRE(child1.fitness == parent1.fitness);
        REQUIRE(child2.fitness == parent1.fitness);
    }
}