/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_THREAD_POOL_HPP
#define GA_UTILITY_THREAD_POOL_HPP

#include "concurrent_queue.hpp"
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

            std::packaged_task<R()> task{ [f = std::move(f), ...args = std::move(args)]
            {
                return std::invoke(std::move(f), std::move(args)...);
            }};

            std::future<R> result = task.get_future();

            bool success = scheduled_worker_queue().emplace([task = std::move(task)]() mutable { task(); });
            if (!success) GAPP_THROW(std::runtime_error, "Attempting to submit a task to a stopped thread pool.");

            return result;
        }

        template<typename F, typename Iter>
        void execute_loop(Iter first, Iter last, F&& unary_op)
        {
            if (first == last) return;

            std::exception_ptr exception;
            std::atomic<bool> has_exception;

            [[maybe_unused]] auto thread_guard = [&]() noexcept
            {
                if (!has_exception.exchange(true, std::memory_order_relaxed))
                    exception = std::current_exception();
            };

            const ptrdiff_t iterations = std::distance(first, last);
            const ptrdiff_t remainder  = iterations % workers_.size();
            const ptrdiff_t step_size  = iterations / workers_.size() + bool(remainder);
            const ptrdiff_t rem_task   = iterations / step_size && iterations % step_size;
            const ptrdiff_t task_count = iterations / step_size + rem_task;

            std::latch remaining_tasks(task_count);

            for (auto partition_last = first; first != last; first = partition_last)
            {
                detail::advance_in_range(partition_last, last, step_size);

                auto task = [&, first, partition_last]() mutable noexcept
                {
                    GAPP_TRY { while (first != partition_last) { std::invoke(unary_op, *first++); } }
                    GAPP_CATCH (...) { thread_guard(); }
                    remaining_tasks.count_down();
                };

                bool success = scheduled_worker_queue().emplace(std::move(task));
                if (!success) GAPP_THROW(std::runtime_error, "Attempting to submit a task to a stopped thread pool.");
            }

            remaining_tasks.wait();

            if (has_exception.load(std::memory_order_relaxed)) std::rethrow_exception(exception);
        }

        ~thread_pool() noexcept { stop(); }

    private:
        using task_t = detail::move_only_function<void()>;

        struct worker_t
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

        void stop() noexcept
        {
            for (worker_t& worker : workers_) { worker.task_queue.close(); }
        }

        std::vector<worker_t> workers_{ std::max(std::thread::hardware_concurrency(), 1u) };
        std::atomic<size_t> turn_;
    };

    struct execution
    {
        GAPP_API inline static thread_pool global_thread_pool;
    };

    template<typename F, typename Iter>
    requires std::invocable<F, std::iter_reference_t<Iter>>
    void parallel_for(Iter first, Iter last, F&& f)
    {
        execution::global_thread_pool.execute_loop(first, last, std::forward<F>(f));
    }

} // namespace gapp::detail

#endif // !GA_UTILITY_THREAD_POOL_HPP
