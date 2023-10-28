/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include "metrics/pop_stats.hpp"
#include "core/candidate.hpp"
#include "utility/rng.hpp"

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


TEST_CASE("hypervolume_size", "[benchmark]")
{
    constexpr size_t ndim = 3;
    constexpr size_t small_size = 40;
    constexpr size_t medium_size = 200;
    constexpr size_t large_size = 1500;

    BENCHMARK_ADVANCED("small")(Benchmark::Chronometer meter)
    {
        FitnessMatrix fmat = randomFitnessMatrix(small_size, ndim);
        FitnessVector ref_point(ndim, 0.0);

        meter.measure([&] { return hypervolume(fmat, ref_point); });
    };

    BENCHMARK_ADVANCED("medium")(Benchmark::Chronometer meter)
    {
        FitnessMatrix fmat = randomFitnessMatrix(medium_size, ndim);
        FitnessVector ref_point(ndim, 0.0);

        meter.measure([&] { return hypervolume(fmat, ref_point); });
    };

    BENCHMARK_ADVANCED("large")(Benchmark::Chronometer meter)
    {
        FitnessMatrix fmat = randomFitnessMatrix(large_size, ndim);
        FitnessVector ref_point(ndim, 0.0);

        meter.measure([&] { return hypervolume(fmat, ref_point); });
    };
}


TEST_CASE("hypervolume_dimensions", "[benchmark]")
{
    constexpr size_t popsize = 30;
    constexpr size_t small_dim = 2;
    constexpr size_t medium_dim = 5;
    constexpr size_t large_dim = 10;

    BENCHMARK_ADVANCED("small")(Benchmark::Chronometer meter)
    {
        FitnessMatrix fmat = randomFitnessMatrix(popsize, small_dim);
        FitnessVector ref_point(small_dim, 0.0);

        meter.measure([&] { return hypervolume(fmat, ref_point); });
    };

    BENCHMARK_ADVANCED("medium")(Benchmark::Chronometer meter)
    {
        FitnessMatrix fmat = randomFitnessMatrix(popsize, medium_dim);
        FitnessVector ref_point(medium_dim, 0.0);

        meter.measure([&] { return hypervolume(fmat, ref_point); });
    };

    BENCHMARK_ADVANCED("large")(Benchmark::Chronometer meter)
    {
        FitnessMatrix fmat = randomFitnessMatrix(popsize, large_dim);
        FitnessVector ref_point(large_dim, 0.0);

        meter.measure([&] { return hypervolume(fmat, ref_point); });
    };
}
