/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_PARALLEL_FOR_HPP
#define GA_UTILITY_PARALLEL_FOR_HPP

#include "utility.hpp"
#include <atomic>
#include <algorithm>
#include <execution>
#include <exception>
#include <functional>
#include <iterator>
#include <concepts>
#include <tuple>

namespace gapp::detail
{
    template<std::forward_iterator Iter, typename F>
    requires std::invocable<F&, std::iter_reference_t<Iter>>
    void parallel_for(Iter first, Iter last, F&& f)
    {
        constinit static std::atomic<bool> barrier;
        constexpr auto tsan_memory_barrier = [] { std::ignore = barrier.exchange(false, std::memory_order_acq_rel); };

        std::atomic<bool> has_exception;
        std::exception_ptr exception;

        [[maybe_unused]] auto thread_guard = [&]() noexcept
        {
            if (!has_exception.exchange(true, std::memory_order_relaxed))
                exception = std::current_exception();
        };

        auto thread_func = [&](auto&& elem) noexcept
        {
            tsan_memory_barrier();
            GAPP_TRY { std::invoke(f, std::forward<decltype(elem)>(elem)); }
            GAPP_CATCH (...) { std::invoke(thread_guard); }
            tsan_memory_barrier();
        };

        tsan_memory_barrier();
        std::for_each(std::execution::par, first, last, std::move(thread_func));
        tsan_memory_barrier();

        if (has_exception.load(std::memory_order_relaxed))
        {
            std::rethrow_exception(exception);
        }
    }

} // namespace gapp::detail

#endif // !GA_UTILITY_PARALLEL_FOR_HPP
