/* Copyright (c) 2025 Krisztián Rugási. Subject to the MIT License. */

#include "rcu.hpp"

namespace gapp::detail
{
    alignas(128) thread_local rcu_domain::registered_reader rcu_domain::reader;

} // namespace gapp::detail
