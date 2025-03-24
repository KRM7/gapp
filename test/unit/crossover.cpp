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
    const Candidate<int> parent1{ { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } };
    const Candidate<int> parent2{ { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 } };
                                                  //

    auto [child1, child2] = singlePointCrossoverImpl(parent1, parent2, 5);

    REQUIRE(child1.chromosome == Chromosome<int>{ { 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0 } });
    REQUIRE(child2.chromosome == Chromosome<int>{ { 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1 } });


    auto [child3, child4] = singlePointCrossoverImpl(parent2, parent1, 5);

    REQUIRE(child3.chromosome == child2.chromosome);
    REQUIRE(child4.chromosome == child1.chromosome);


    auto [child5, child6] = singlePointCrossoverImpl(parent1, parent2, 0);

    REQUIRE(child5.chromosome == parent1.chromosome);
    REQUIRE(child6.chromosome == parent2.chromosome);


    auto [child7, child8] = singlePointCrossoverImpl(parent1, parent2, 12);

    REQUIRE(child7.chromosome == parent2.chromosome);
    REQUIRE(child8.chromosome == parent1.chromosome);


    auto [child9, child10] = singlePointCrossoverImpl(parent1, parent1, rng::randomInt(0, 12));

    REQUIRE(child9.chromosome == parent1.chromosome);
    REQUIRE(child10.chromosome == parent1.chromosome);
}

TEST_CASE("two_point_crossover", "[crossover]")
{
    const Candidate<int> parent1{ { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } };
    const Candidate<int> parent2{ { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 } };
                                            //                //

    auto [child1, child2] = twoPointCrossoverImpl(parent1, parent2, { 9, 3 });

    REQUIRE(child1.chromosome == Chromosome<int>{ { 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0 } });
    REQUIRE(child2.chromosome == Chromosome<int>{ { 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1 } });


    auto [child3, child4] = twoPointCrossoverImpl(parent2, parent1, { 9, 3 });

    REQUIRE(child3.chromosome == child2.chromosome);
    REQUIRE(child4.chromosome == child1.chromosome);


    auto [child5, child6] = twoPointCrossoverImpl(parent1, parent2, { 3, 9 });

    REQUIRE(child5.chromosome == child1.chromosome);
    REQUIRE(child6.chromosome == child2.chromosome);


    auto [child7, child8] = twoPointCrossoverImpl(parent1, parent2, { 0, 12 });

    REQUIRE(child7.chromosome == parent2.chromosome);
    REQUIRE(child8.chromosome == parent1.chromosome);


    auto [child9, child10] = twoPointCrossoverImpl(parent1, parent1, { rng::randomInt(0, 12), rng::randomInt(0, 12) });

    REQUIRE(child9.chromosome == parent1.chromosome);
    REQUIRE(child10.chromosome == parent1.chromosome);
}

TEST_CASE("npoint_crossover", "[crossover]")
{
    const Candidate<char> parent1{ { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } };
    const Candidate<char> parent2{ { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 } };
                                       //    //          //          //

    auto [child1, child2] = nPointCrossoverImpl(parent1, parent2, { 1, 3, 7, 11 });

    REQUIRE(child1.chromosome == Chromosome<char>{ { 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1 } });
    REQUIRE(child2.chromosome == Chromosome<char>{ { 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0 } });


    auto [child3, child4] = nPointCrossoverImpl(parent2, parent1, { 1, 3, 7, 11 });

    REQUIRE(child3.chromosome == child2.chromosome);
    REQUIRE(child4.chromosome == child1.chromosome);


    auto [child5, child6] = nPointCrossoverImpl(parent1, parent2, { 3, 11, 1, 7 });

    REQUIRE(child5.chromosome == child1.chromosome);
    REQUIRE(child6.chromosome == child2.chromosome);


    auto [child7, child8] = nPointCrossoverImpl(parent1, parent2, { 0 });

    REQUIRE(child7.chromosome == parent1.chromosome);
    REQUIRE(child8.chromosome == parent2.chromosome);


    auto [child9, child10] = nPointCrossoverImpl(parent1, parent2, { 12 });

    REQUIRE(child9.chromosome == parent2.chromosome);
    REQUIRE(child10.chromosome == parent1.chromosome);


    auto [child11, child12] = nPointCrossoverImpl(parent1, parent2, { 0, 12 });

    REQUIRE(child11.chromosome == parent1.chromosome);
    REQUIRE(child12.chromosome == parent2.chromosome);


    auto [child13, child14] = nPointCrossoverImpl(parent1, parent1, rng::sampleUnique(0_sz, 12_sz, 4));

    REQUIRE(child13.chromosome == parent1.chromosome);
    REQUIRE(child14.chromosome == parent1.chromosome);
}


TEMPLATE_TEST_CASE("order1_crossover", "[crossover]", int, unsigned)
{
    const Candidate<TestType> parent1{ { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 } };
    const Candidate<TestType> parent2{ { 4, 5, 0, 6, 1, 2, 8, 3, 9, 7 } };

    const auto child1 = order1CrossoverImpl(parent1, parent2, 4, 8);
    const auto child2 = order1CrossoverImpl(parent2, parent1, 4, 8);

    REQUIRE(child1.chromosome == Chromosome<TestType>{ { 1, 2, 8, 3, 4, 5, 6, 7, 9, 0 } });
    REQUIRE(child2.chromosome == Chromosome<TestType>{ { 4, 5, 6, 7, 1, 2, 8, 3, 9, 0 } });

    const auto child3 = order1CrossoverImpl(parent1, parent1, 4, 8);

    REQUIRE(child3.chromosome == parent1.chromosome);
}

TEMPLATE_TEST_CASE("order2_crossover", "[crossover]", int, unsigned)
{
    const Candidate<TestType> parent1{ { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 } };
    const Candidate<TestType> parent2{ { 4, 5, 0, 6, 1, 2, 8, 3, 9, 7 } };

    const auto child1 = order2CrossoverImpl(parent1, parent2, 4, 8);
    const auto child2 = order2CrossoverImpl(parent2, parent1, 4, 8);

    REQUIRE(child1.chromosome == Chromosome<TestType>{ { 0, 1, 2, 8, 4, 5, 6, 7, 3, 9 } });
    REQUIRE(child2.chromosome == Chromosome<TestType>{ { 0, 4, 5, 6, 1, 2, 8, 3, 7, 9 } });

    const auto child3 = order2CrossoverImpl(parent1, parent1, 4, 8);

    REQUIRE(child3.chromosome == parent1.chromosome);
}

TEMPLATE_TEST_CASE("position_crossover", "[crossover]", int, unsigned)
{
    const Candidate<TestType> parent1{ { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 } };
    const Candidate<TestType> parent2{ { 4, 5, 0, 6, 1, 2, 8, 3, 9, 7 } };

    const auto child1 = positionCrossoverImpl(parent1, parent2, { { 0, 3, 4, 7 } });
    const auto child2 = positionCrossoverImpl(parent2, parent1, { { 0, 3, 4, 7 } });

    REQUIRE(child1.chromosome == Chromosome<TestType>{ { 0, 5, 6, 3, 4, 1, 2, 7, 8, 9 } });
    REQUIRE(child2.chromosome == Chromosome<TestType>{ { 4, 0, 2, 6, 1, 5, 7, 3, 8, 9 } });

    auto child3 = positionCrossoverImpl(parent1, parent1, { { 0, 3, 4, 7 } });

    REQUIRE(child3.chromosome == parent1.chromosome);
}

TEST_CASE("cycle_crossover", "[crossover]")
{
    const Candidate<int> parent1{ { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 } };
    const Candidate<int> parent2{ { 4, 5, 0, 6, 1, 2, 8, 3, 9, 7 } };

    // cycle0 : 0 - 4 - 1 - 5 - 2 , cycle1 : 3 - 6 - 8 - 9 - 7

    auto [child1, child2] = cycleCrossoverImpl(parent1, parent2);

    REQUIRE(child1.chromosome == Chromosome<int>{ { 0, 1, 2, 6, 4, 5, 8, 3, 9, 7 } });
    REQUIRE(child2.chromosome == Chromosome<int>{ { 4, 5, 0, 3, 1, 2, 6, 7, 8, 9 } });

    auto [child3, child4] = cycleCrossoverImpl(parent1, parent1);

    REQUIRE(child3.chromosome == parent1.chromosome);
    REQUIRE(child4.chromosome == parent1.chromosome);
}

TEMPLATE_TEST_CASE("edge_crossover", "[crossover]", int, unsigned)
{
    const Candidate<TestType> parent1{ { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 } };
    const Candidate<TestType> parent2{ { 4, 5, 0, 6, 1, 2, 8, 3, 9, 7 } };

    auto child1 = edgeCrossoverImpl(parent1, parent2);
    auto child2 = edgeCrossoverImpl(parent2, parent1);

    REQUIRE(child1.chromosome == Chromosome<TestType>{ { 0, 5, 4, 1, 6, 7, 9, 8, 3, 2 } });
    REQUIRE(child2.chromosome == Chromosome<TestType>{ { 4, 5, 0, 6, 1, 2, 3, 9, 8, 7 } });

    auto child3 = edgeCrossoverImpl(parent1, parent1);

    REQUIRE(child3.chromosome == parent1.chromosome);
}

TEMPLATE_TEST_CASE("pmx_crossover", "[crossover]", int, unsigned)
{
    const Candidate<TestType> parent1{ { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 } };
    const Candidate<TestType> parent2{ { 4, 5, 0, 6, 1, 2, 8, 3, 9, 7 } };

    auto child1 = pmxCrossoverImpl(parent1, parent2, 4, 8);
    auto child2 = pmxCrossoverImpl(parent2, parent1, 4, 8);

    REQUIRE(child1.chromosome == Chromosome<TestType>{ { 1, 2, 0, 8, 4, 5, 6, 7, 9, 3 } });
    REQUIRE(child2.chromosome == Chromosome<TestType>{ { 0, 4, 5, 7, 1, 2, 8, 3, 6, 9 } });

    auto child3 = pmxCrossoverImpl(parent1, parent1, 4, 8);

    REQUIRE(child3.chromosome == parent1.chromosome);
}

TEMPLATE_TEST_CASE("real_crossover", "[crossover]", real::Arithmetic, real::BLXa, real::SimulatedBinary, real::Wright)
{
    using Crossover = TestType;

    constexpr size_t chrom_len = 10;
    const BoundsVector<RealGene> bounds(chrom_len, { 0.0, 1.0 });

    RCGA context;
    context.solve(DummyFitnessFunction<RealGene>(chrom_len), bounds, 1);

    Crossover crossover{ 0.8 };

    Candidate<RealGene> parent1{ { 0.0, 0.12, 0.48, 0.19, 1.0, 1.0, 0.0, 0.72, 0.81, 0.03 }, bounds };
    Candidate<RealGene> parent2{ { 1.0, 0.34, 0.97, 0.36, 1.0, 0.0, 0.0, 0.28, 0.49, 0.79 }, bounds };
    parent1.fitness = { 0.0 };
    parent2.fitness = { 0.0 };

    auto [child1, child2] = crossover(context, parent1, parent2);

    REQUIRE(child1.chromosome.size() == chrom_len);
    REQUIRE(child2.chromosome.size() == chrom_len);
    REQUIRE(std::all_of(child1.chromosome.begin(), child1.chromosome.end(), detail::between(0.0, 1.0)));
    REQUIRE(std::all_of(child2.chromosome.begin(), child2.chromosome.end(), detail::between(0.0, 1.0)));


    Candidate<RealGene> parent3{ GaTraits<RealGene>::randomChromosome(chrom_len, bounds), bounds };
    Candidate<RealGene> parent4{ GaTraits<RealGene>::randomChromosome(chrom_len, bounds), bounds };
    parent3.fitness = { 0.0 };
    parent4.fitness = { 0.0 };

    auto [child3, child4] = crossover(context, parent3, parent4);

    REQUIRE(child1.chromosome.size() == chrom_len);
    REQUIRE(child2.chromosome.size() == chrom_len);
    REQUIRE(std::all_of(child1.chromosome.begin(), child1.chromosome.end(), detail::between(0.0, 1.0)));
    REQUIRE(std::all_of(child2.chromosome.begin(), child2.chromosome.end(), detail::between(0.0, 1.0)));


    auto [child5, child6] = crossover(context, parent3, parent3);

    REQUIRE(child5.chromosome == parent3.chromosome);
    REQUIRE(child6.chromosome == parent3.chromosome);
}

TEST_CASE("mixed_crossover", "[crossover]")
{
    Mixed crossover{ binary::SinglePoint{}, real::Arithmetic{} };

    crossover.crossover_rates({ 0.3, 0.7 });
    REQUIRE(crossover.crossover_rates() == std::array<Probability, 2>{ 0.3, 0.7 });

    crossover.crossover_rates(1.0);
    REQUIRE(crossover.crossover_rates() == std::array<Probability, 2>{ 1.0, 1.0 });

    crossover.crossover_rate<RealGene>(0.5);
    REQUIRE(crossover.crossover_rate<RealGene>() == 0.5);

    REQUIRE(!crossover.allow_variable_chrom_length<RealGene>());
    REQUIRE(!crossover.allow_variable_chrom_length<BinaryGene>());

    REQUIRE(crossover.component<BinaryGene>().crossover_rate() == 1.0);
    REQUIRE(crossover.component<RealGene>().crossover_rate() == 0.5);


    const std::array chrom_lens = { 3_sz, 4_sz };
    const BoundsVector<RealGene> bounds(chrom_lens[1], { 0.0,1.0 });

    MixedGA<BinaryGene, RealGene> context;
    context.solve(DummyFitnessFunction<MixedGene<BinaryGene, RealGene>>(chrom_lens), bounds, 1);

    Candidate<MixedGene<BinaryGene, RealGene>> parent1{ { Chromosome<BinaryGene>{ 0, 0, 0 }, Chromosome<RealGene>{ 0.0, 0.0, 0.0, 0.0 } }, bounds };
    Candidate<MixedGene<BinaryGene, RealGene>> parent2{ { Chromosome<BinaryGene>{ 1, 1, 1 }, Chromosome<RealGene>{ 1.0, 1.0, 1.0, 1.0 } }, bounds };
    parent1.fitness = { 0.0 };
    parent2.fitness = { 0.0 };

    auto [child1, child2] = crossover(context, parent1, parent2);

    REQUIRE(child1.chrom_len<BinaryGene>() == chrom_lens[0]);
    REQUIRE(child2.chrom_len<BinaryGene>() == chrom_lens[0]);
    REQUIRE(child1.chrom_len<RealGene>() == chrom_lens[1]);
    REQUIRE(child2.chrom_len<RealGene>() == chrom_lens[1]);

    REQUIRE(std::all_of(child1.chrom<BinaryGene>().begin(), child1.chrom<BinaryGene>().end(), detail::between(0, 1)));
    REQUIRE(std::all_of(child2.chrom<BinaryGene>().begin(), child2.chrom<BinaryGene>().end(), detail::between(0, 1)));
    REQUIRE(std::all_of(child1.chrom<RealGene>().begin(), child1.chrom<RealGene>().end(), detail::between(0.0, 1.0)));
    REQUIRE(std::all_of(child2.chrom<RealGene>().begin(), child2.chrom<RealGene>().end(), detail::between(0.0, 1.0)));

    REQUIRE(detail::equal(child1.bounds<RealGene>(), BoundsView<RealGene>(bounds)));
    REQUIRE(detail::equal(child1.bounds<RealGene>(), child2.bounds<RealGene>()));
}

TEST_CASE("crossover_bounds", "[crossover]")
{
    const BoundsVector<RealGene> bounds(10, { 0.0, 1.0 });

    RCGA context;
    context.solve(DummyFitnessFunction<RealGene>(10), bounds, 1);

    real::Arithmetic crossover;

    Candidate<RealGene> parent1{ { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, bounds };
    Candidate<RealGene> parent2{ { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 }, bounds };
    parent1.fitness = { 0.0 };
    parent2.fitness = { 0.0 };

    auto [child1, child2] = crossover(context, parent1, parent2);

    REQUIRE(detail::equal(child1.gene_bounds, parent1.gene_bounds));
    REQUIRE(detail::equal(child2.gene_bounds, parent2.gene_bounds));
}

TEST_CASE("crossover_lambda", "[crossover]")
{
    crossover::Lambda<BinaryGene> crossover([](const GaInfo&, const auto& parent1, const auto& parent2)
    {
        return CandidatePair<BinaryGene>{ parent1, parent2 };
    });

    crossover.crossover_rate(0.1);
    REQUIRE(crossover.crossover_rate() == 0.1);
}

TEST_CASE("crossover_callable", "[mutation]")
{
    RCGA ga;
    ga.crossover_method([](const GaInfo&, const auto& p1, const auto& p2)
    {
        return CandidatePair<RealGene>{ p1, p2 };
    });

    REQUIRE(!ga.crossover_method().allow_variable_chrom_length());
}
