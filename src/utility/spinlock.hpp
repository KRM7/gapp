/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_UTILITY_SPINLOCK_HPP
#define GAPP_UTILITY_SPINLOCK_HPP

#include "utility.hpp"
#include <atomic>

namespace gapp::detail
{
    class spinlock
    {
    public:
        void lock() noexcept
        {
            while (true)
            {
                if (!locked_.test_and_set(std::memory_order_acquire)) break;
                while (locked_.test(std::memory_order_relaxed)) GAPP_PAUSE();
            }
        }

        bool try_lock() noexcept
        {
            return !locked_.test(std::memory_order_relaxed) &&
                   !locked_.test_and_set(std::memory_order_acquire);
        }

        void unlock() noexcept
        {
            locked_.clear(std::memory_order_release);
        }

    private:
        std::atomic_flag locked_;
    };

} // namespace gapp::detail

#endif // !GAPP_UTILITY_SPINLOCK_HPP
