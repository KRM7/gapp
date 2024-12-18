/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include "core/population.hpp"
#include "utility/rng.hpp"
#include "utility/functional.hpp"
#include <algorithm>
#include <numeric>
#include <cmath>

using namespace gapp;
using namespace gapp::detail;
using namespace Catch;


static FitnessVector randomFitness(size_t nobj, double radius)
{
    FitnessVector fvec(nobj, 0.0);

    std::generate(fvec.begin(), fvec.end(), []{ return std::abs(rng::randomNormal()); });
    const double sum = std::inner_product(fvec.begin(), fvec.end(), fvec.begin(), 0.0);
    std::transform(fvec.begin(), fvec.end(), fvec.begin(), detail::divide_by(std::sqrt(sum) / radius));

    return fvec;
}

static Population<RealGene> randomPopulation(size_t popsize, size_t nobj, double radius = 1.0)
{
    Population<RealGene> pop(popsize, Candidate<RealGene>());
    for (auto& candidate : pop) { candidate.fitness = randomFitness(nobj, radius); }

    return pop;
}


TEST_CASE("merge_pareto_sets_size", "[benchmark]")
{
    constexpr size_t num_obj = 3;
    const size_t popsize = GENERATE(40, 200, 1500);

    WARN("Population/set size: " << popsize);

    BENCHMARK_ADVANCED("merge")(Benchmark::Chronometer meter)
    {
        Population<RealGene> lhs = randomPopulation(popsize, num_obj);
        Population<RealGene> rhs = randomPopulation(popsize, num_obj);

        meter.measure([&] { return mergeParetoSets(std::move(lhs), std::move(rhs)); });
    };

    BENCHMARK_ADVANCED("append/naive")(Benchmark::Chronometer meter)
    {
        Population<RealGene> lhs = randomPopulation(popsize, num_obj);
        Population<RealGene> rhs = randomPopulation(popsize, num_obj);

        meter.measure([&]
        {
            lhs.reserve(lhs.size() + rhs.size());
            lhs.insert(lhs.end(), rhs.begin(), rhs.end());
            return findParetoFront(lhs);
        });
    };
}


TEST_CASE("merge_pareto_sets_relative_sizes", "[benchmark]")
{
    constexpr size_t num_obj = 3;
    constexpr size_t large_popsize = 2500;
    constexpr size_t small_popsize = 100;

    BENCHMARK_ADVANCED("merge_left_greater")(Benchmark::Chronometer meter)
    {
        Population<RealGene> lhs = randomPopulation(large_popsize, num_obj);
        Population<RealGene> rhs = randomPopulation(small_popsize, num_obj);

        meter.measure([&] { return mergeParetoSets(std::move(lhs), std::move(rhs)); });
    };

    BENCHMARK_ADVANCED("merge_right_greater")(Benchmark::Chronometer meter)
    {
        Population<RealGene> lhs = randomPopulation(small_popsize, num_obj);
        Population<RealGene> rhs = randomPopulation(large_popsize, num_obj);

        meter.measure([&] { return mergeParetoSets(std::move(lhs), std::move(rhs)); });
    };
}


TEST_CASE("merge_pareto_sets_dominated", "[benchmark]")
{
    constexpr size_t num_obj = 3;
    constexpr size_t popsize = 500;

    BENCHMARK_ADVANCED("merge_equal")(Benchmark::Chronometer meter)
    {
        Population<RealGene> lhs = randomPopulation(popsize, num_obj);
        Population<RealGene> rhs = randomPopulation(popsize, num_obj);

        meter.measure([&] { return mergeParetoSets(std::move(lhs), std::move(rhs)); });
    };

    BENCHMARK_ADVANCED("merge_left_dominated")(Benchmark::Chronometer meter)
    {
        Population<RealGene> lhs = randomPopulation(popsize, num_obj, 0.5);
        Population<RealGene> rhs = randomPopulation(popsize, num_obj, 1.0);

        meter.measure([&] { return mergeParetoSets(std::move(lhs), std::move(rhs)); });
    };

    BENCHMARK_ADVANCED("merge_right_dominated")(Benchmark::Chronometer meter)
    {
        Population<RealGene> lhs = randomPopulation(popsize, num_obj, 1.0);
        Population<RealGene> rhs = randomPopulation(popsize, num_obj, 0.5);

        meter.measure([&] { return mergeParetoSets(std::move(lhs), std::move(rhs)); });
    };
}


TEST_CASE("merge_pareto_sets_objectives", "[benchmark]")
{
    constexpr size_t small_obj = 3;
    constexpr size_t medium_obj = 15;
    constexpr size_t large_obj = 100;
    constexpr size_t popsize = 100;

    BENCHMARK_ADVANCED("merge_small_obj")(Benchmark::Chronometer meter)
    {
        Population<RealGene> lhs = randomPopulation(popsize, small_obj);
        Population<RealGene> rhs = randomPopulation(popsize, small_obj);

        meter.measure([&] { return mergeParetoSets(std::move(lhs), std::move(rhs)); });
    };

    BENCHMARK_ADVANCED("merge_medium_obj")(Benchmark::Chronometer meter)
    {
        Population<RealGene> lhs = randomPopulation(popsize, medium_obj);
        Population<RealGene> rhs = randomPopulation(popsize, medium_obj);

        meter.measure([&] { return mergeParetoSets(std::move(lhs), std::move(rhs)); });
    };

    BENCHMARK_ADVANCED("merge_large_obj")(Benchmark::Chronometer meter)
    {
        Population<RealGene> lhs = randomPopulation(popsize, large_obj);
        Population<RealGene> rhs = randomPopulation(popsize, large_obj);

        meter.measure([&] { return mergeParetoSets(std::move(lhs), std::move(rhs)); });
    };
}
