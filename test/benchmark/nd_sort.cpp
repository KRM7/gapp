/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include "algorithm/nd_sort.hpp"
#include "core/population.hpp"
#include "utility/rng.hpp"

using namespace gapp;
using namespace gapp::algorithm::dtl;
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


TEST_CASE("nd_sort_popsize", "[benchmark]")
{
    constexpr size_t num_obj = 3;
    const size_t popsize = GENERATE(40, 200, 1500);

    WARN("Population size: " << popsize);

    BENCHMARK_ADVANCED("FNDS")(Benchmark::Chronometer meter)
    {
        FitnessMatrix fmat = randomFitnessMatrix(popsize, num_obj);

        meter.measure([&] { return fastNonDominatedSort(fmat.begin(), fmat.end()); });
    };

    BENCHMARK_ADVANCED("DDS")(Benchmark::Chronometer meter)
    {
        FitnessMatrix fmat = randomFitnessMatrix(popsize, num_obj);

        meter.measure([&] { return dominanceDegreeSort(fmat.begin(), fmat.end()); });
    };
}


TEST_CASE("nd_sort_num_objectives", "[benchmark]")
{
    constexpr size_t popsize = 200;
    const size_t num_obj = GENERATE(2, 15, 100);

    WARN("Number of objectives: " << num_obj);

    BENCHMARK_ADVANCED("FNDS")(Benchmark::Chronometer meter)
    {
        FitnessMatrix fmat = randomFitnessMatrix(popsize, num_obj);

        meter.measure([&] { return fastNonDominatedSort(fmat.begin(), fmat.end()); });
    };

    BENCHMARK_ADVANCED("DDS")(Benchmark::Chronometer meter)
    {
        FitnessMatrix fmat = randomFitnessMatrix(popsize, num_obj);

        meter.measure([&] { return dominanceDegreeSort(fmat.begin(), fmat.end()); });
    };
}
