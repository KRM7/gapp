/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include "utility/thread_pool.hpp"
#include "utility/iterators.hpp"
#include <atomic>
#include <vector>
#include <numeric>
#include <thread>
#include <tuple>

using namespace gapp::detail;
using namespace Catch::Matchers;


TEST_CASE("concurrent_queue", "[thread-pool]")
{
    concurrent_queue<int> queue;

    std::vector<int> output;
    std::vector<int> input(1000);
    std::iota(input.begin(), input.end(), 0);

    auto input_task = [&]
    {
        for (int n : input) { std::ignore = queue.emplace(n); }
        queue.close();
    };

    auto output_task = [&]
    {
        for (auto n = queue.take(); n.has_value(); n = queue.take()) { output.push_back(*n); }
    };

    {
        std::jthread t1{ input_task };
        std::jthread t2{ output_task };
    }

    REQUIRE_THAT(output, Equals(input));
}

TEST_CASE("parallel_for", "[thread-pool]")
{
    int n = 0;
    auto increment_n = [&](int) { std::atomic_ref{ n }.fetch_add(1, std::memory_order_relaxed); };

    parallel_for(iota_iterator(0), iota_iterator(100), increment_n);
    REQUIRE(n == 100);

    parallel_for(iota_iterator(0), iota_iterator(100), increment_n);
    REQUIRE(n == 200);
}

TEST_CASE("nested_parallel_for", "[thread-pool]")
{
    int n = 0;

    parallel_for(iota_iterator(0), iota_iterator(10), [&](int)
    {
        parallel_for(iota_iterator(0), iota_iterator(10), [&](int)
        {
            parallel_for(iota_iterator(0), iota_iterator(100), [&](int)
            {
                std::atomic_ref{ n }.fetch_add(1, std::memory_order_relaxed);
            });
        });
    });

    REQUIRE(n == 10000);
}
