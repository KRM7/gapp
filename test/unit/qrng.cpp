/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/catch_approx.hpp>
#include "utility/qrng.hpp"
#include "utility/functional.hpp"
#include <algorithm>

using namespace gapp;
using namespace gapp::rng;


TEST_CASE("qrng", "[rng][qrng]")
{
    const size_t dim = GENERATE(0, 1, 2, 3, 10);

    QuasiRandom qrng{ dim };

    REQUIRE(qrng.dim() == dim);

    auto point = qrng();

    REQUIRE(point.size() == dim);
    REQUIRE(std::all_of(point.begin(), point.end(), detail::between(0.0, 1.0)));
}
