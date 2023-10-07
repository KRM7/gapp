﻿/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_THREAD_POOL_HPP
#define GA_UTILITY_THREAD_POOL_HPP

#include "concurrent_queue.hpp"
#include "algorithm.hpp"
#include "functional.hpp"
#include "iterators.hpp"
#include "utility.hpp"
#include <type_traits>
#include <concepts>
#include <thread>
#include <atomic>
#include <future>
#include <latch>
#include <functional>
#include <iterator>
#include <exception>
#include <stdexcept>
#include <utility>
#include <cstddef>

namespace gapp::detail
{
    class thread_pool
    {
    public:
        template<typename F, typename... Args>
        auto execute_task(F f, Args... args) -> std::future<std::invoke_result_t<F, Args...>>
        {
            using R = std::invoke_result_t<F, Args...>;

            std::packaged_task<R()> task{ [f = std::move(f), ...args = std::move(args)]() mutable { std::invoke(std::move(f), std::move(args)...); }};
            std::future<R> result = task.get_future();

            bool success = scheduled_worker_queue().emplace([task = std::move(task)]() mutable { task(); });
            if (!success) GAPP_THROW(std::runtime_error, "Attempting to submit a task to a stopped thread pool.");

            return result;
        }

        template<typename F, typename Iter>
        void execute_loop(Iter first, Iter last, size_t block_size, F&& unary_op)
        {
            if (first == last) return;

            std::exception_ptr exception;
            std::atomic<bool> has_exception;

            [[maybe_unused]] auto thread_guard = [&]() noexcept
            {
                if (!has_exception.exchange(true, std::memory_order_relaxed))
                    exception = std::current_exception();
            };

            const size_t iterations  = std::distance(first, last);
            const size_t block_count = iterations / block_size + bool(iterations % block_size);
            const size_t task_count  = detail::min(workers_.size() + workers_.empty(), iterations, block_count);
            const size_t step_size   = iterations / task_count;
            const size_t remainder   = iterations % task_count;

            std::latch remaining_tasks(task_count - 1);

            for (size_t i = 0; i < task_count - 1; i++)
            {
                Iter block_last  = std::next(first, step_size + (i < remainder));
                Iter block_first = std::exchange(first, block_last);

                auto task = [&, block_first, block_last]() mutable noexcept
                {
                    GAPP_TRY { while (block_first != block_last) { std::invoke(unary_op, *block_first++); } }
                    GAPP_CATCH (...) { thread_guard(); }
                    remaining_tasks.count_down();
                };

                bool success = scheduled_worker_queue().emplace(std::move(task));
                if (!success) GAPP_THROW(std::runtime_error, "Attempting to submit a task to a stopped thread pool.");
            }

            while (first != last) { std::invoke(unary_op, *first++); }

            if (auto* local_queue = local_worker_queue(); !local_queue)
            {
                remaining_tasks.wait();
            }
            else while (true)
            {
                while (!local_queue->empty()) { std::invoke(*local_queue->take()); }
                if (remaining_tasks.try_wait()) break;
                std::this_thread::yield();
            }

            if (has_exception.load(std::memory_order_relaxed)) std::rethrow_exception(exception);
        }

        ~thread_pool() noexcept { stop(); }

    private:
        using task_t = detail::move_only_function<void()>;

        struct alignas(128) worker_t
        {
            static void worker_main(detail::concurrent_queue<task_t>& task_queue) noexcept
            {
                while (true)
                {
                    auto task = task_queue.take();
                    if (!task.has_value()) return;
                    std::invoke(std::move(*task));
                }
            }

            detail::concurrent_queue<task_t> task_queue;
            std::jthread thread{ worker_main, std::ref(task_queue) };
        };

        auto scheduled_worker_queue() noexcept -> detail::concurrent_queue<task_t>&
        {
            const size_t current_turn   = turn_.fetch_add(1, std::memory_order_relaxed);
            const size_t current_worker = current_turn % workers_.size();

            return workers_[current_worker].task_queue;
        }

        auto local_worker_queue() noexcept -> detail::concurrent_queue<task_t>*
        {
            for (worker_t& worker : workers_)
            {
                if (worker.thread.get_id() == std::this_thread::get_id()) { return &worker.task_queue; }
            }
            return nullptr;
        }

        void stop() noexcept
        {
            for (worker_t& worker : workers_) { worker.task_queue.close(); }
        }

        std::vector<worker_t> workers_{ std::max(std::thread::hardware_concurrency(), 1u) - 1u };
        std::atomic<size_t> turn_;
    };

    struct execution_context
    {
        GAPP_API inline static thread_pool global_thread_pool;
    };


    template<typename F, typename Iter>
    requires std::invocable<F, std::iter_reference_t<Iter>>
    void parallel_for(Iter first, Iter last, F&& f)
    {
        execution_context::global_thread_pool.execute_loop(first, last, 1, std::forward<F>(f));
    }

    template<typename F, typename Iter>
    requires std::invocable<F, std::iter_reference_t<Iter>>
    void parallel_for(Iter first, Iter last, size_t block_size, F&& f)
    {
        execution_context::global_thread_pool.execute_loop(first, last, block_size, std::forward<F>(f));
    }

} // namespace gapp::detail

#endif // !GA_UTILITY_THREAD_POOL_HPP
