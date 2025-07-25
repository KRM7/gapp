/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "utility/thread_pool.hpp"
#include "utility/iterators.hpp"
#include <atomic>
#include <vector>
#include <numeric>
#include <thread>
#include <tuple>

using namespace gapp;
using namespace gapp::detail;

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

            parallel_for(iota_iterator(0), iota_iterator(100), [&](int)
            {
                std::atomic_ref{ n }.fetch_add(1, std::memory_order_relaxed);
            });
        });
    });

    REQUIRE(n == 20000);
}

TEST_CASE("thread_count", "[thread-pool]")
{
    const size_t thread_count = GENERATE(1, 8, 123);

    execution_threads(thread_count);

    REQUIRE(execution_threads() == thread_count);

    int n = 0;

    parallel_for(iota_iterator(0), iota_iterator(10), [&](int)
    {
        parallel_for(iota_iterator(0), iota_iterator(100), [&](int)
        {
            std::atomic_ref{ n }.fetch_add(1, std::memory_order_relaxed);
        });
    });

    REQUIRE(n == 1000);

    execution_threads(std::thread::hardware_concurrency());
}

TEST_CASE("task_exceptions", "[thread_pool]")
{
    REQUIRE_THROWS(
        parallel_for(iota_iterator(0), iota_iterator(10), [&](int)
        {
            throw std::logic_error{ "" };
        })
    );

    REQUIRE_THROWS(
        parallel_for(iota_iterator(0), iota_iterator(10), [&](int)
        {
            parallel_for(iota_iterator(0), iota_iterator(10), [&](int)
            {
                throw std::logic_error{ "" };
            });
        })
    );
}
