/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include "utility/thread_pool.hpp"
#include "utility/iterators.hpp"
#include <algorithm>
#include <numeric>
#include <execution>
#include <atomic>

using namespace gapp::detail;


TEST_CASE("parallel_for", "[benchmark]")
{
    const size_t work_size = GENERATE(10, 100, 1000, 10000);

    std::atomic<double> n = 0.0;

    std::vector v(work_size, 0.0);
    std::iota(v.begin(), v.end(), 0.0);

    auto work = [&](int)
    {
        n += std::inner_product(v.begin(), v.end(), v.begin(), 0.0) / std::reduce(v.begin(), v.end(), 0.0);
    };


    BENCHMARK("single_parallel_for")
    {
        parallel_for(iota_iterator(0), iota_iterator(1000), work);
        return n.load();
    };

    BENCHMARK("single_std_for_each")
    {
        std::for_each(std::execution::par, iota_iterator(0), iota_iterator(1000), work);
        return n.load();
    };


    BENCHMARK("nested_parallel_for")
    {
        parallel_for(iota_iterator(0), iota_iterator(10), [&](int)
        {
            parallel_for(iota_iterator(0), iota_iterator(10), [&](int)
            {
                parallel_for(iota_iterator(0), iota_iterator(100), work);
            });
        });
        return n.load();
    };

    BENCHMARK("nested_std_for_each")
    {
        std::for_each(std::execution::par, iota_iterator(0), iota_iterator(10), [&](int)
        {
            std::for_each(std::execution::par, iota_iterator(0), iota_iterator(10), [&](int)
            {
                std::for_each(std::execution::par, iota_iterator(0), iota_iterator(100), work);
            });
        });
        return n.load();
    };
}
