/* Copyright (c) 2024 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_LATCH_HPP
#define GA_UTILITY_LATCH_HPP

#include "utility.hpp"
#include <atomic>
#include <thread>
#include <cstdint>

namespace gapp::detail
{
    class latch
    {
    public:
        explicit latch(std::uint32_t n) noexcept :
            count_(n)
        {}

        void count_down(std::uint32_t n = 1) noexcept
        {
            count_.fetch_sub(n, std::memory_order_release);
        }

        void wait() const noexcept
        {
            while (count_.load(std::memory_order_relaxed)) std::this_thread::yield();
            std::atomic_thread_fence(std::memory_order_acquire);
            GAPP_ANNOTATE_TSAN_ACQUIRE(&count_);
        }

        bool try_wait() const noexcept
        {
            return count_.load(std::memory_order_acquire) == 0;
        }

    private:
        alignas(64) std::atomic<std::uint32_t> count_;
    };

} // namespace gapp::detail

#endif // !GA_UTILITY_LATCH_HPP
