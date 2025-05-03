/* Copyright (c) 2025 Krisztián Rugási. Subject to the MIT License. */

#include "thread_pool.hpp"
#include "rng.hpp"
#include <atomic>
#include <cstdint>

namespace gapp::detail
{
    thread_local std::atomic<std::uint64_t> thread_pool::this_thread_id_ = 0;

} // namespace gapp::detail

namespace gapp::rng
{
    alignas(128) thread_local ConcurrentXoroshiro128p::RegisteredGenerator ConcurrentXoroshiro128p::generator_;

} // namespace gapp::rng
