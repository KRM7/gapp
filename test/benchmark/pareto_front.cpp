/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include "core/population.hpp"
#include "utility/rng.hpp"
#include "utility/functional.hpp"

using namespace gapp;
using namespace gapp::detail;
using namespace Catch;


static FitnessMatrix randomFitnessMatrix(size_t pop_size, size_t num_obj)
{
    FitnessMatrix fmat(pop_size, num_obj);

    for (const auto& row : fmat)
    {
        for (double& val : row) { val = rng::randomReal(); }
    }

    return fmat;
}


TEST_CASE("find_pareto_front_size", "[benchmark]")
{
    constexpr size_t num_obj = 3;
    const size_t popsize = GENERATE(40, 200, 1500);

    WARN("Population size: " << popsize);

    BENCHMARK_ADVANCED("sort")(Benchmark::Chronometer meter)
    {
        FitnessMatrix fmat = randomFitnessMatrix(popsize, num_obj);

        meter.measure([&] { return findParetoFrontSort(fmat); });
    };

    BENCHMARK_ADVANCED("best")(Benchmark::Chronometer meter)
    {
        FitnessMatrix fmat = randomFitnessMatrix(popsize, num_obj);

        meter.measure([&] { return findParetoFrontBest(fmat); });
    };

    BENCHMARK_ADVANCED("kung")(Benchmark::Chronometer meter)
    {
        FitnessMatrix fmat = randomFitnessMatrix(popsize, num_obj);

        meter.measure([&] { return findParetoFrontKung(fmat); });
    };
}


TEST_CASE("find_pareto_front_nobj", "[benchmark]")
{
    constexpr size_t popsize = 200;
    const size_t num_obj = GENERATE(3, 15, 100);

    WARN("Number of objectives: " << num_obj);

    BENCHMARK_ADVANCED("sort")(Benchmark::Chronometer meter)
    {
        FitnessMatrix fmat = randomFitnessMatrix(popsize, num_obj);

        meter.measure([&] { return findParetoFrontSort(fmat); });
    };

    BENCHMARK_ADVANCED("best")(Benchmark::Chronometer meter)
    {
        FitnessMatrix fmat = randomFitnessMatrix(popsize, num_obj);

        meter.measure([&] { return findParetoFrontBest(fmat); });
    };

    BENCHMARK_ADVANCED("kung")(Benchmark::Chronometer meter)
    {
        FitnessMatrix fmat = randomFitnessMatrix(popsize, num_obj);

        meter.measure([&] { return findParetoFrontKung(fmat); });
    };
}


TEST_CASE("find_pareto_front_1D", "[benchmark]")
{
    const size_t popsize = GENERATE(40, 200, 1500);

    WARN("Population size: " << popsize);

    BENCHMARK_ADVANCED("1D_spec")(Benchmark::Chronometer meter)
    {
        FitnessMatrix fmat = randomFitnessMatrix(popsize, 1);

        meter.measure([&] { return findParetoFront1D(fmat); });
    };

    BENCHMARK_ADVANCED("sort")(Benchmark::Chronometer meter)
    {
        FitnessMatrix fmat = randomFitnessMatrix(popsize, 1);

        meter.measure([&] { return findParetoFrontSort(fmat); });
    };

    BENCHMARK_ADVANCED("best")(Benchmark::Chronometer meter)
    {
        FitnessMatrix fmat = randomFitnessMatrix(popsize, 1);

        meter.measure([&] { return findParetoFrontBest(fmat); });
    };
}
