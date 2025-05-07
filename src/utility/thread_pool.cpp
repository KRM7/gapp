/* Copyright (c) 2025 Krisztián Rugási. Subject to the MIT License. */

#include "thread_pool.hpp"
#include <cstdint>

namespace gapp::detail
{
    thread_local std::atomic<std::uint64_t> thread_pool::this_thread_id_;

} // namespace gapp::detail
