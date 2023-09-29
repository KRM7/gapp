/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_SPINLOCK_HPP
#define GA_UTILITY_SPINLOCK_HPP

#include "utility.hpp"
#include <atomic>
#include <thread>

namespace gapp::detail
{
    class spinlock
    {
    public:
        void lock() noexcept
        {
            while (true)
            {
                while (locked_.test(std::memory_order_relaxed)) GAPP_PAUSE();
                if (!locked_.test_and_set(std::memory_order_acquire)) break;
                std::this_thread::yield();
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

#endif // !GA_UTILITY_SPINLOCK_HPP
