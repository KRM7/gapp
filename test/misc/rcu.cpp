/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include "utility/rcu.hpp"
#include <thread>
#include <shared_mutex>
#include <chrono>
#include <iostream>

using namespace gapp::detail;
using namespace std::chrono_literals;


rcu_obj<int> number = 0;

static const auto reader_func = []
{
    while (true)
    {
        std::shared_lock _{ number };
        [[maybe_unused]] const int& n = number.get();
        std::this_thread::sleep_for(2ms);
        assert(0 <= n && n <= 100);
    }
};

static const auto writer_func = [i = 0]() mutable
{
    while (true) { number = i++ % 100; }
};

static const auto status_func = []
{
    while (true)
    {
        std::cout << "Running RCU tests...\n";
        std::this_thread::sleep_for(5s);
    }
};


int main()
{
    std::jthread reader1{ reader_func };
    std::jthread reader2{ reader_func };
    std::jthread reader3{ reader_func };
    std::jthread reader4{ reader_func };
    std::jthread reader5{ reader_func };

    std::jthread writer1{ writer_func };
    std::jthread writer2{ writer_func };

    std::jthread status{ status_func };
}
