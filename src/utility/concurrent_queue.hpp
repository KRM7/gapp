/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_CONCURRENT_QUEUE_HPP
#define GA_UTILITY_CONCURRENT_QUEUE_HPP

#include <condition_variable>
#include <mutex>
#include <deque>
#include <optional>
#include <type_traits>

namespace gapp::detail
{
    template<typename T>
    class concurrent_queue
    {
    public:
        template<typename... Args>
        [[nodiscard]] bool emplace(Args&&... args)
        {
            std::scoped_lock lock{ queue_lock_ };
            if (is_closed_) return false;
            queue_.emplace_back(std::forward<Args>(args)...);
            queue_cv_.notify_one();
            return true;
        }

        [[nodiscard]] std::optional<T> take() noexcept(std::is_nothrow_move_constructible_v<T>)
        {
            std::unique_lock lock{ queue_lock_ };
            queue_cv_.wait(lock, [&] { return !queue_.empty() || is_closed_; });

            if (is_closed_ && queue_.empty()) return {};

            T elem = std::move(queue_.front());
            queue_.pop_front();
            return elem;
        }

        void close() noexcept
        {
            std::scoped_lock lock{ queue_lock_ };
            is_closed_ = true;
            queue_cv_.notify_all();
        }

        [[nodiscard]] bool closed() const noexcept
        {
            std::scoped_lock lock{ queue_lock_ };
            return is_closed_;
        }

        [[nodiscard]] bool empty() const noexcept
        {
            std::scoped_lock lock{ queue_lock_ };
            return queue_.empty();
        }

    private:
        std::deque<T> queue_;
        mutable std::mutex queue_lock_;
        std::condition_variable queue_cv_;
        bool is_closed_ = false;
    };

} // namespace gapp::detail

#endif // !GA_UTILITY_CONCURRENT_QUEUE_HPP
