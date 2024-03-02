/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_RCU_HPP
#define GA_UTILITY_RCU_HPP

#include "shared_spinlock.hpp"
#include "indestructible.hpp"
#include "utility.hpp"
#include <atomic>
#include <memory>
#include <mutex>
#include <shared_mutex>
#include <array>
#include <vector>
#include <functional>
#include <cstdint>
#include <cstddef>

namespace gapp::detail
{
    class disposal_queue
    {
    public:
        template<typename T>
        void emplace(T* ptr) noexcept requires(!std::is_array_v<T>)
        {
            constexpr auto deleter = [](void* ptr) { delete static_cast<T*>(ptr); };
            deleters_[size_++] = { ptr, deleter };
        }

        bool is_full() const noexcept
        {
            return size_ == deleters_.size();
        }

        void collect() noexcept
        {
            while (size_) std::invoke(deleters_[--size_]);
        }

        ~disposal_queue() noexcept { collect(); }

    private:
        struct deleter_t
        {
            void operator()() const noexcept { std::invoke(deleter, ptr); }

            void* ptr             = nullptr;
            void(*deleter)(void*) = nullptr;
        };

        std::array<deleter_t, 16> deleters_;
        std::uint8_t size_ = 0;
    };

    class rcu_domain
    {
    public:
        static void read_lock() noexcept
        {
            const std::uint64_t read_counter = reader.counter.load(std::memory_order_relaxed);

            if (nesting_depth(read_counter) == 0)
            {
                GAPP_ASSERT(nesting_depth(writer_counter.load(std::memory_order_relaxed)) == 1);
                reader.counter.store(writer_counter.load(std::memory_order_relaxed), std::memory_order_release);
                std::ignore = reader.counter.load(std::memory_order_acquire);
            }
            else
            {
                GAPP_ASSERT(nesting_depth(read_counter + 1) > nesting_depth(read_counter));
                reader.counter.store(read_counter + 1, std::memory_order_relaxed);
            }
        }

        static void read_unlock() noexcept
        {
            GAPP_ASSERT(nesting_depth(reader.counter.load(std::memory_order_relaxed)));
            reader.counter.store(reader.counter.load(std::memory_order_relaxed) - 1, std::memory_order_release);
        }

        static void synchronize() noexcept
        {
            GAPP_ASSERT(!nesting_depth(reader.counter.load(std::memory_order_relaxed)));

            std::uint64_t current = writer_counter.load(std::memory_order_acquire);
            std::uint64_t target  = current + 0x100;
            writer_counter.compare_exchange_strong(current, target, std::memory_order_release);

            std::shared_lock _{ tls_readers->lock };

            for (const registered_reader* tls_reader : tls_readers->list)
            {
                while (true)
                {
                    const std::uint64_t counter = tls_reader->counter.load(std::memory_order_acquire);
                    if (!nesting_depth(counter) || counter >= target) break;
                    std::this_thread::yield();
                }
            }
        }

        template<typename T>
        static void retire(T* ptr) noexcept
        {
            GAPP_ASSERT(!garbage_queue.is_full());

            garbage_queue.emplace(ptr);

            if (garbage_queue.is_full())
            {
                rcu_domain::synchronize();
                garbage_queue.collect();
            }
        }

    private:
        struct registered_reader
        {
            registered_reader() noexcept // NOLINT(*exception-escape)
            {
                std::unique_lock _{ tls_readers->lock };
                tls_readers->list.push_back(this);
            }

            ~registered_reader() noexcept
            {
                std::unique_lock _{ tls_readers->lock };
                std::erase(tls_readers->list, this);
            }

            std::atomic<std::uint64_t> counter = 0;
        };

        struct tls_reader_list
        {
            detail::shared_spinlock lock;
            std::vector<registered_reader*> list;
        };

        static std::uint64_t nesting_depth(std::uint64_t counter) noexcept
        {
            return counter & 0xFF;
        }

        GAPP_API inline static detail::Indestructible<tls_reader_list> tls_readers;
        GAPP_API inline static constinit std::atomic<std::uint64_t> writer_counter = 1;
        alignas(128) inline static thread_local registered_reader reader;
        inline static thread_local constinit disposal_queue garbage_queue;
    };

    // NOLINTBEGIN(*owning-memory)

    template<typename T>
    class rcu_obj
    {
    public:
        template<typename... Args>
        constexpr rcu_obj(Args&&... args) :
            data_(new T(std::forward<Args>(args)...))
        {}

        ~rcu_obj() noexcept
        {
            delete data_.load(std::memory_order_consume);
        }

        template<typename U>
        rcu_obj& operator=(U&& value)
        {
            T* new_ptr = new T(std::forward<U>(value));
            T* old_ptr = data_.exchange(new_ptr, std::memory_order_acq_rel);
            rcu_domain::retire(old_ptr);

            return *this;
        }

        rcu_obj(const rcu_obj&)            = delete;
        rcu_obj(rcu_obj&&)                 = delete;

        rcu_obj& operator=(const rcu_obj&) = delete;
        rcu_obj& operator=(rcu_obj&&)      = delete;

        T& get() noexcept { return *data_.load(std::memory_order_consume); }
        const T& get() const noexcept { return *data_.load(std::memory_order_consume); }

        T& operator*() noexcept { return get(); }
        const T& operator*() const noexcept { return get(); }

        T* operator->() noexcept { return std::addressof(get()); }
        const T* operator->() const noexcept { return std::addressof(get()); }

        void lock_shared() const noexcept { rcu_domain::read_lock(); }
        void unlock_shared() const noexcept { rcu_domain::read_unlock(); }
        bool try_lock_shared() const noexcept { rcu_domain::read_lock(); return true; }

    private:
        std::atomic<T*> data_;
    };

    // NOLINTEND(*owning-memory)

} // namespace gapp::detail

#endif // !GA_UTILITY_RCU_HPP
