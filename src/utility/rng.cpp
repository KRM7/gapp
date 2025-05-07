/* Copyright (c) 2025 Krisztián Rugási. Subject to the MIT License. */

#include "rng.hpp"

namespace gapp::rng
{
    alignas(128) thread_local ConcurrentXoroshiro128p::RegisteredGenerator ConcurrentXoroshiro128p::generator_;

} // namespace gapp::rng
