/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_SHARED_SPINLOCK_HPP
#define GA_UTILITY_SHARED_SPINLOCK_HPP

#include "utility.hpp"
#include <atomic>
#include <limits>
#include <cstdint>

namespace gapp::detail
{
    class shared_spinlock
    {
    public:
        void lock() noexcept
        {
            while (true)
            {
                if (try_lock()) break;
                while (cntr_.load(std::memory_order_relaxed)) GAPP_PAUSE();
            }
        }

        bool try_lock() noexcept
        {
            std::uint32_t expected = 0;
            return cntr_.compare_exchange_strong(expected, WRITER, std::memory_order_acquire, std::memory_order_relaxed);
        }

        void unlock() noexcept
        {
            cntr_.fetch_sub(WRITER, std::memory_order_release);
        }

        void lock_shared() noexcept
        {
            cntr_.fetch_add(1, std::memory_order_relaxed);
            while (cntr_.load(std::memory_order_relaxed) >= WRITER) GAPP_PAUSE();
            std::atomic_thread_fence(std::memory_order_acquire);
            GAPP_ANNOTATE_TSAN_ACQUIRE(&cntr_);
        }

        bool try_lock_shared() noexcept
        {
            if (cntr_.fetch_add(1, std::memory_order_acquire) < WRITER) return true;
            cntr_.fetch_sub(1, std::memory_order_relaxed);
            return false;
        }

        void unlock_shared() noexcept
        {
            cntr_.fetch_sub(1, std::memory_order_release);
        }

    private:
        constexpr inline static std::uint32_t WRITER = std::numeric_limits<std::uint32_t>::max() >> 1;
        std::atomic<std::uint32_t> cntr_;
    };

} // namespace gapp::detail

#endif // !GA_UTILITY_SHARED_SPINLOCK_HPP
