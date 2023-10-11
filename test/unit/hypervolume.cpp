/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "metrics/pop_stats.hpp"
#include "core/candidate.hpp"

using namespace gapp;
using namespace gapp::detail;
using Catch::Approx;


TEST_CASE("hypervolume", "[metrics]")
{
    SECTION("1D")
    {
        const FitnessMatrix fmat = {
            {  1.0 },
            {  1.0 },
            {  3.2 },
            { -2.0 },
            { 10.0 },
            {  5.6 },
            {-10.0 },
            { 10.0 },
            {  0.2 }
        };
        const FitnessVector ref_point = { -10.0 };

        REQUIRE(hypervolume(fmat, ref_point) == Approx(20.0).margin(1E-8));
    }

    SECTION("2D")
    {
        const FitnessMatrix fmat = {
            { 2.0,  12.0 },
            { 10.0, 3.0 },
            { 6.0,  10.0 },
            { 10.0, 10.0 },
            { 10.0, 10.0 },
            { 13.0, 3.0 },
            { 0.0,  0.0 },
            { 12.0, 6.0 },
            { 5.0,  7.0 },
            { 1.0,  2.0 },
            { 20.0, 0.0 }
        };
        const FitnessVector ref_point = { 0.0, 0.0 };

        REQUIRE(hypervolume(fmat, ref_point) == Approx(119.0).margin(1E-8));
    }

    SECTION("3D")
    {
        const FitnessMatrix fmat = {
            { 10.0, 10.0, 10.0 },
            { 11.0,  8.0,  3.0 },
            {  4.0,  4.0, 18.0 },
            {  0.0,  0.0,  0.0 },
            { 12.0,  2.0,  6.0 },
            { 10.0,  8.0, 10.0 },
            { 11.0,  8.0,  3.0 },
            { 11.0,  8.0,  3.0 },
            {  8.0, 13.0,  8.0 },
            {  1.0,  1.0,  9.0 },
            { 40.0,  0.0,  0.0 }
        };
        const FitnessVector ref_point = { 0.0, 0.0, 0.0 };

        REQUIRE(hypervolume(fmat, ref_point) == Approx(1362.0).margin(1E-8));
    }

    SECTION("zero")
    {
        const FitnessMatrix fmat = { { 1.0, 1.0, 1.0, 1.0, 1.0 } };
        const FitnessVector ref_point{ fmat[0] };

        REQUIRE(hypervolume(fmat, ref_point) == 0.0);
    }

    SECTION("inf")
    {
        const FitnessMatrix fmat = {
            { 0.0 },
            { math::inf<double> },
            { 1E+105 }
        };
        const FitnessVector ref_point = { 0.0 };

        REQUIRE(hypervolume(fmat, ref_point) == math::inf<double>);
    }

    SECTION("empty")
    {
        const FitnessMatrix fmat;
        const FitnessVector ref_point = { 1.0, 2.0 };

        REQUIRE(hypervolume(fmat, ref_point) == 0.0);
    }
}
