/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include "utility/rcu.hpp"
#include <atomic>
#include <shared_mutex>

using namespace gapp::detail;
using namespace std::chrono_literals;

std::shared_mutex rwlock;

volatile size_t number = 0;
volatile std::atomic_size_t atomic_number = 0;
rcu_obj<volatile size_t> rcu_number = 0;

TEST_CASE("rcu_lock", "[benchmark]")
{
    BENCHMARK("read") { return number; };
    BENCHMARK("atomic_fetch_add") { return atomic_number.fetch_add(1); };
    BENCHMARK("rwlock_read") { std::shared_lock _{ rwlock }; return number; };
    BENCHMARK("rcu_read") { std::scoped_lock _{ rcu_number }; return rcu_number.get(); };
}
