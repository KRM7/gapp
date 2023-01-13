/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include "crossover/crossover_dtl.hpp"
#include "population/candidate.hpp"


using namespace genetic_algorithm;
using namespace genetic_algorithm::crossover::dtl;

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

    std::vector<size_t> cx_points{ 1, 3, 7, 11 };

    auto [child1, child2] = nPointCrossoverImpl(parent1, parent2, cx_points);

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

    auto child1 = positionCrossoverImpl(parent1, parent2, { 0, 3, 4, 7 });
    auto child2 = positionCrossoverImpl(parent2, parent1, { 0, 3, 4, 7 });

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

TEST_CASE("edge_crossover", "[crossover]")
{
    Candidate<unsigned> parent1{ { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 } };
    Candidate<unsigned> parent2{ { 4, 5, 0, 6, 1, 2, 8, 3, 9, 7 } };

    auto child1 = edgeCrossoverImpl(parent1, parent2);
    auto child2 = edgeCrossoverImpl(parent2, parent1);

    REQUIRE(child1.chromosome == Chromosome<unsigned>{ { 0, 5, 6, 1, 2, 8, 9, 7, 4, 3 } });
    REQUIRE(child2.chromosome == Chromosome<unsigned>{ { 4, 5, 6, 7, 0, 1, 2, 8, 9, 3 } });
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