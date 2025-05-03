/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_UTILITY_THREAD_POOL_HPP
#define GAPP_UTILITY_THREAD_POOL_HPP

#include "concurrent_queue.hpp"
#include "algorithm.hpp"
#include "functional.hpp"
#include "iterators.hpp"
#include "latch.hpp"
#include "utility.hpp"
#include <algorithm>
#include <type_traits>
#include <concepts>
#include <thread>
#include <atomic>
#include <future>
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
        template<typename F, typename Iter>
        void execute_loop(Iter first, Iter last, size_t block_size, F&& unary_op)
        {
            if (first == last) return;

            std::exception_ptr exception;
            std::atomic<bool> has_exception;

            [[maybe_unused]] const auto thread_guard = [&]() noexcept
            {
                if (!has_exception.exchange(true, std::memory_order_relaxed))
                    exception = std::current_exception();
            };

            const size_t iterations  = std::distance(first, last);
            const size_t block_count = iterations / block_size + bool(iterations % block_size);
            const size_t task_count  = detail::min(workers_.size() + 1, iterations, block_count);
            const size_t step_size   = iterations / task_count;
            const size_t remainder   = iterations % task_count;

            detail::latch remaining_tasks(task_count - 1);

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

        void reset_scheduler()
        {
            turn_.store(0, std::memory_order_relaxed);
        }

        void thread_count(size_t count)
        {
            reset_scheduler();
            if (count == thread_count()) return;
            stop();
            workers_ = std::vector<worker_t>(count - 1);
        }

        size_t thread_count() const noexcept
        {
            return workers_.size() + 1;
        }

        ~thread_pool() noexcept { stop(); }

        thread_pool() = default;

        thread_pool(const thread_pool&)            = delete;
        thread_pool(thread_pool&&)                 = delete;

        thread_pool& operator=(const thread_pool&) = delete;
        thread_pool& operator=(thread_pool&&)      = delete;

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

namespace gapp
{
    /**
    * Set the number of threads that will be used by the library to run the genetic algorithms.
    * 
    * The value should be between 1 and the number of hardware threads. The default number of threads
    * used will be the value returned by std::thread::hardware_concurrency.
    * 
    * @note This function is not thread-safe and shouldn't be called while a genetic algorithm
    *   is running.
    * 
    * @param count The desired number of threads. Must be at least 1.
    */
    inline void execution_threads(size_t count)
    {
        detail::execution_context::global_thread_pool.thread_count(std::max(count, 1_sz));
    }

    /** @returns The number of threads used by the library. */
    inline size_t execution_threads() noexcept
    {
        return detail::execution_context::global_thread_pool.thread_count();
    }

} // namespace gapp

#endif // !GAPP_UTILITY_THREAD_POOL_HPP
