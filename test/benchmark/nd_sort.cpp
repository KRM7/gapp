/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include "algorithm/nd_sort.hpp"
#include "population/population.hpp"
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
    constexpr size_t small_size = 40;
    constexpr size_t medium_size = 200;
    constexpr size_t large_size = 1500;

    BENCHMARK_ADVANCED("fnds_small")(Benchmark::Chronometer meter)
    {
        FitnessMatrix fmat = randomFitnessMatrix(small_size, num_obj);

        meter.measure([&] { return fastNonDominatedSort(fmat.begin(), fmat.end()); });
    };

    BENCHMARK_ADVANCED("dds_small")(Benchmark::Chronometer meter)
    {
        FitnessMatrix fmat = randomFitnessMatrix(small_size, num_obj);

        meter.measure([&] { return dominanceDegreeSort(fmat.begin(), fmat.end()); });
    };

    BENCHMARK_ADVANCED("fnds_medium")(Benchmark::Chronometer meter)
    {
        FitnessMatrix fmat = randomFitnessMatrix(medium_size, num_obj);

        meter.measure([&] { return fastNonDominatedSort(fmat.begin(), fmat.end()); });
    };

    BENCHMARK_ADVANCED("dds_medium")(Benchmark::Chronometer meter)
    {
        FitnessMatrix fmat = randomFitnessMatrix(medium_size, num_obj);

        meter.measure([&] { return dominanceDegreeSort(fmat.begin(), fmat.end()); });
    };

    BENCHMARK_ADVANCED("fnds_large")(Benchmark::Chronometer meter)
    {
        FitnessMatrix fmat = randomFitnessMatrix(large_size, num_obj);

        meter.measure([&] { return fastNonDominatedSort(fmat.begin(), fmat.end()); });
    };

    BENCHMARK_ADVANCED("dds_large")(Benchmark::Chronometer meter)
    {
        FitnessMatrix fmat = randomFitnessMatrix(large_size, num_obj);

        meter.measure([&] { return dominanceDegreeSort(fmat.begin(), fmat.end()); });
    };
}


TEST_CASE("nd_sort_num_objectives", "[benchmark]")
{
    constexpr size_t popsize = 200;
    constexpr size_t small_obj = 2;
    constexpr size_t medium_obj = 15;
    constexpr size_t large_obj = 100;

    BENCHMARK_ADVANCED("fnds_small")(Benchmark::Chronometer meter)
    {
        FitnessMatrix fmat = randomFitnessMatrix(popsize, small_obj);

        meter.measure([&] { return fastNonDominatedSort(fmat.begin(), fmat.end()); });
    };

    BENCHMARK_ADVANCED("dds_small")(Benchmark::Chronometer meter)
    {
        FitnessMatrix fmat = randomFitnessMatrix(popsize, small_obj);

        meter.measure([&] { return dominanceDegreeSort(fmat.begin(), fmat.end()); });
    };

    BENCHMARK_ADVANCED("fnds_medium")(Benchmark::Chronometer meter)
    {
        FitnessMatrix fmat = randomFitnessMatrix(popsize, medium_obj);

        meter.measure([&] { return fastNonDominatedSort(fmat.begin(), fmat.end()); });
    };

    BENCHMARK_ADVANCED("dds_medium")(Benchmark::Chronometer meter)
    {
        FitnessMatrix fmat = randomFitnessMatrix(popsize, medium_obj);

        meter.measure([&] { return dominanceDegreeSort(fmat.begin(), fmat.end()); });
    };

    BENCHMARK_ADVANCED("fnds_large")(Benchmark::Chronometer meter)
    {
        FitnessMatrix fmat = randomFitnessMatrix(popsize, large_obj);

        meter.measure([&] { return fastNonDominatedSort(fmat.begin(), fmat.end()); });
    };

    BENCHMARK_ADVANCED("dds_large")(Benchmark::Chronometer meter)
    {
        FitnessMatrix fmat = randomFitnessMatrix(popsize, large_obj);

        meter.measure([&] { return dominanceDegreeSort(fmat.begin(), fmat.end()); });
    };
}
