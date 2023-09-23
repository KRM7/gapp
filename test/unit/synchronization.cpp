/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include "utility/spinlock.hpp"
#include "utility/shared_spinlock.hpp"
#include <mutex>
#include <shared_mutex>
#include <thread>

using namespace gapp::detail;


TEST_CASE("spinlock", "[synchronization]")
{
    int n = 0;
    spinlock lock;

    auto increment_n = [&]
    {
        for (size_t i = 0; i < 1000; i++) { std::scoped_lock _{ lock }; n++; }
    };

    {
        std::jthread t1{ increment_n };
        std::jthread t2{ increment_n };
        std::jthread t3{ increment_n };
    }

    REQUIRE(n == 3000);
}

TEST_CASE("shared_spinlock", "[synchronization]")
{
    int n = 0;
    int read = 0;
    shared_spinlock lock;

    auto increment_n = [&]
    {
        for (size_t i = 0; i < 1000; i++) { std::scoped_lock _{ lock }; n++; }
    };

    auto reader_func = [&]
    {
        for (size_t i = 0; i < 1000; i++) { std::shared_lock _{ lock }; read = n; }
    };

    {
        std::jthread t1{ increment_n };
        std::jthread t2{ increment_n };
        std::jthread t3{ increment_n };
        std::jthread t4{ reader_func };
    }

    REQUIRE(n == 3000);
    REQUIRE(0 <= read);
    REQUIRE(read <= 3000);
}
