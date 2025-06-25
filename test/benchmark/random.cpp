/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include "utility/rng.hpp"
#include <random>

#if __has_include(<boost/random.hpp>)
#  include <boost/random.hpp>
#endif

using namespace gapp;

TEST_CASE("uniform_bool_distribution", "[benchmark]")
{
    rng::uniform_bool_distribution dist1;
    std::uniform_int_distribution<uint64_t> dist2(0, 1);

    BENCHMARK("rng::uniform_bool_distribution") { return dist1(rng::prng); };
    BENCHMARK("std::uniform_int_distribution") { return dist2(rng::prng); };
    BENCHMARK("uniform_generator_last_bit") { return rng::prng() & 1; };
}

TEST_CASE("uniform_int_distribution", "[benchmark]")
{
    rng::uniform_int_distribution<uint64_t> dist1(0, 100);
    std::uniform_int_distribution<uint64_t> dist2(0, 100);

    BENCHMARK("rng::uniform_int_distribution") { return dist1(rng::prng); };
    BENCHMARK("std::uniform_int_distribution") { return dist2(rng::prng); };

#if __has_include(<boost/random.hpp>)
    boost::uniform_int<uint64_t> dist3(0, 100);
    BENCHMARK("boost::uniform_int_distribution") { return dist3(rng::prng); };
#endif
}

TEST_CASE("uniform_real_distribution", "[benchmark]")
{
    rng::uniform_real_distribution<double> dist1(0.0, 100.0);
    std::uniform_real_distribution<double> dist2(0.0, 100.0);

    BENCHMARK("rng::uniform_real_distribution") { return dist1(rng::prng); };
    BENCHMARK("std::uniform_real_distribution") { return dist2(rng::prng); };

#if __has_include(<boost/random.hpp>)
    boost::uniform_real<double> dist3(0.0, 100.0);
    BENCHMARK("boost::uniform_real_distribution") { return dist3(rng::prng); };
#endif
}

TEST_CASE("exponential_distribution", "[benchmark]")
{
    rng::exponential_distribution<double> dist1(5.0);
    std::exponential_distribution<double> dist2(5.0);

    BENCHMARK("rng::exponential_distribution") { return dist1(rng::prng); };
    BENCHMARK("std::exponential_distribution") { return dist2(rng::prng); };

#if __has_include(<boost/random.hpp>)
    boost::exponential_distribution<double> dist3(5.0);
    BENCHMARK("boost::exponential_distribution") { return dist3(rng::prng); };
#endif
}

TEST_CASE("normal_distribution", "[benchmark]")
{
    rng::normal_distribution<double> dist1;
    std::normal_distribution<double> dist3;

    BENCHMARK("rng::normal_distribution") { return dist1(rng::prng); };
    BENCHMARK("std::normal_distribution") { return dist3(rng::prng); };

#if __has_include(<boost/random.hpp>)
    boost::normal_distribution<double> dist4;
    BENCHMARK("boost::normal_distribution") { return dist4(rng::prng); };
#endif
}

TEST_CASE("poisson_distribution", "[benchmark]")
{
    rng::small_poisson_distribution<uint64_t> dist1(6.0);
    std::poisson_distribution<uint64_t> dist2(6.0);

    BENCHMARK("rng::small_poisson_distribution") { return dist1(rng::prng); };
    BENCHMARK("std::poisson_distribution") { return dist2(rng::prng); };

#if __has_include(<boost/random.hpp>)
    boost::poisson_distribution<uint64_t> dist3(6.0);
    BENCHMARK("boost::poisson_distribution") { return dist3(rng::prng); };
#endif
}

TEST_CASE("symmetric_binomial_distribution, n = 100", "[benchmark]")
{
    rng::symmetric_binomial_distribution<uint64_t> dist1(100);
    rng::binomial_distribution<uint64_t> dist2(100, 0.5);
    std::binomial_distribution<uint64_t> dist3(100, 0.5);

    BENCHMARK("rng::symmetric_binomial_distribution") { return dist1(rng::prng); };
    BENCHMARK("rng::binomial_distribution, p = 0.5") { return dist2(rng::prng); };
    BENCHMARK("std::binomial_distribution, p = 0.5") { return dist3(rng::prng); };

#if __has_include(<boost/random.hpp>)
    boost::binomial_distribution<uint64_t> dist4(100, 0.5);
    BENCHMARK("boost::binomial_distribution") { return dist4(rng::prng); };
#endif
}

TEST_CASE("symmetric_binomial_distribution, n = 500", "[benchmark]")
{
    rng::symmetric_binomial_distribution<uint64_t> dist1(500);
    rng::binomial_distribution<uint64_t> dist2(500, 0.5);
    std::binomial_distribution<uint64_t> dist3(500, 0.5);

    BENCHMARK("rng::symmetric_binomial_distribution") { return dist1(rng::prng); };
    BENCHMARK("rng::binomial_distribution, p = 0.5") { return dist2(rng::prng); };
    BENCHMARK("std::binomial_distribution, p = 0.5") { return dist3(rng::prng); };

#if __has_include(<boost/random.hpp>)
    boost::binomial_distribution<uint64_t> dist4(500, 0.5);
    BENCHMARK("boost::binomial_distribution") { return dist4(rng::prng); };
#endif
}

TEST_CASE("binomial_distribution, mean=1.0", "[benchmark]")
{
    rng::binomial_distribution<uint64_t> dist1(100, 0.01);
    std::binomial_distribution<uint64_t> dist2(100, 0.01);

    BENCHMARK("rng::binomial_distribution") { return dist1(rng::prng); };
    BENCHMARK("std::binomial_distribution") { return dist2(rng::prng); };

#if __has_include(<boost/random.hpp>)
    boost::binomial_distribution<uint64_t> dist3(100, 0.01);
    BENCHMARK("boost::binomial_distribution") { return dist3(rng::prng); };
#endif
}

TEST_CASE("binomial_distribution, mean=5.0", "[benchmark]")
{
    rng::binomial_distribution<uint64_t> dist1(100, 0.05);
    std::binomial_distribution<uint64_t> dist2(100, 0.05);

    BENCHMARK("rng::binomial_distribution") { return dist1(rng::prng); };
    BENCHMARK("std::binomial_distribution") { return dist2(rng::prng); };

#if __has_include(<boost/random.hpp>)
    boost::binomial_distribution<uint64_t> dist3(100, 0.05);
    BENCHMARK("boost::binomial_distribution") { return dist3(rng::prng); };
#endif
}

TEST_CASE("binomial_distribution, mean=10.0", "[benchmark]")
{
    rng::binomial_distribution<uint64_t> dist1(100, 0.1);
    std::binomial_distribution<uint64_t> dist2(100, 0.1);

    BENCHMARK("rng::binomial_distribution") { return dist1(rng::prng); };
    BENCHMARK("std::binomial_distribution") { return dist2(rng::prng); };

#if __has_include(<boost/random.hpp>)
    boost::binomial_distribution<uint64_t> dist3(100, 0.1);
    BENCHMARK("boost::binomial_distribution") { return dist3(rng::prng); };
#endif
}

TEST_CASE("binomial_distribution, mean=20.0", "[benchmark]")
{
    rng::binomial_distribution<uint64_t> dist1(100, 0.2);
    std::binomial_distribution<uint64_t> dist2(100, 0.2);

    BENCHMARK("rng::binomial_distribution") { return dist1(rng::prng); };
    BENCHMARK("std::binomial_distribution") { return dist2(rng::prng); };

#if __has_include(<boost/random.hpp>)
    boost::binomial_distribution<uint64_t> dist3(100, 0.2);
    BENCHMARK("boost::binomial_distribution") { return dist3(rng::prng); };
#endif
}

TEST_CASE("binomial_distribution, mean=50.0", "[benchmark]")
{
    rng::binomial_distribution<uint64_t> dist1(1000, 0.05);
    std::binomial_distribution<uint64_t> dist2(1000, 0.05);

    BENCHMARK("rng::binomial_distribution") { return dist1(rng::prng); };
    BENCHMARK("std::binomial_distribution") { return dist2(rng::prng); };

#if __has_include(<boost/random.hpp>)
    boost::binomial_distribution<uint64_t> dist3(1000, 0.05);
    BENCHMARK("boost::binomial_distribution") { return dist3(rng::prng); };
#endif
}
