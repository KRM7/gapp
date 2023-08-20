/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "encoding/binary.hpp"
#include "metrics/metrics.hpp"
#include "utility/functional.hpp"
#include "test_utils.hpp"
#include <algorithm>

using namespace gapp;
using namespace gapp::metrics;

constexpr static size_t num_obj = 3;
constexpr static size_t num_gen = 10;
constexpr static size_t popsize = 100;


TEMPLATE_TEST_CASE("fitness_metrics", "[metrics]", FitnessMin, FitnessMax, FitnessMean, FitnessVariance, FitnessStdDev, NadirPoint)
{
    using Metric = TestType;

    BinaryGA GA{ popsize };
    GA.track(Metric{});
    GA.solve(DummyFitnessFunction<BinaryGene>{ 10, num_obj }, num_gen);

    const auto& metric = GA.get_metric<Metric>();

    REQUIRE(metric.size() == num_gen);
    REQUIRE(metric.data().size() == num_gen);

    REQUIRE(std::all_of(metric.begin(), metric.end(), detail::is_size(num_obj)));
    REQUIRE(metric[4].size() == num_obj);

    const auto& val = metric[7];
    REQUIRE(std::all_of(val.begin(), val.end(), detail::equal_to(0.0)));
}

TEST_CASE("nadir_point_metric", "[metrics]")
{
    BinaryGA GA{ popsize };
    GA.track(NadirPoint{});
    GA.solve(DummyFitnessFunction<BinaryGene>{ 10, num_obj }, num_gen);

    const auto& metric = GA.get_metric<NadirPoint>();

    REQUIRE(metric.size() == num_gen);
    REQUIRE(std::all_of(metric.begin(), metric.end(), detail::is_size(num_obj)));

    const auto& val = metric[5];
    REQUIRE(std::all_of(val.begin(), val.end(), detail::equal_to(0.0)));
}

TEST_CASE("hypervolume_metric", "[metrics]")
{
    BinaryGA GA{ popsize };
    GA.track(Hypervolume{ FitnessVector(num_obj, -10.0) });
    GA.solve(DummyFitnessFunction<BinaryGene>{ 10, num_obj }, num_gen);

    const auto& metric = GA.get_metric<Hypervolume>();

    REQUIRE(metric.size() == num_gen);
    REQUIRE(std::all_of(metric.begin(), metric.end(), detail::equal_to(1000.0)));
}

TEST_CASE("hypervolume_auto", "[metrics]")
{
    BinaryGA GA{ popsize };
    GA.track(AutoHypervolume{ });
    GA.solve(DummyFitnessFunction<BinaryGene>{ 10, num_obj }, num_gen);

    const auto& metric = GA.get_metric<AutoHypervolume>();

    REQUIRE(metric.size() == num_gen);
    REQUIRE(std::all_of(metric.begin(), metric.end(), detail::equal_to(0.0)));
}

TEST_CASE("fitness_evaluations", "[metrics]")
{
    BinaryGA GA{ popsize };
    GA.track(FitnessEvaluations{});
    GA.solve(DummyFitnessFunction<BinaryGene>{ 10, num_obj }, num_gen);

    const auto& metric1 = GA.get_metric<FitnessEvaluations>();

    REQUIRE(metric1.size() == num_gen);
    REQUIRE(std::all_of(metric1.begin(), metric1.end(), detail::between(0_sz, popsize)));

    GA.solve(DummyFitnessFunction<BinaryGene>{ 10, num_obj, false, true }, num_gen);

    const auto& metric2 = GA.get_metric<FitnessEvaluations>();

    REQUIRE(metric2.size() == num_gen);
    REQUIRE(std::all_of(metric2.begin(), metric2.end(), detail::equal_to(popsize)));
}