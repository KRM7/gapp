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
#include <vector>
#include <cstdint>
#include <cstddef>

namespace gapp::detail
{
    class rcu_domain
    {
    public:
        static void read_lock() noexcept
        {
            reader.epoch.store(writer_epoch.load(std::memory_order_relaxed), std::memory_order_release);
            std::ignore = reader.epoch.load(std::memory_order_acquire);
        }

        static void read_unlock() noexcept
        {
            reader.epoch.store(NOT_READING, std::memory_order_release);
        }

        static void synchronize() noexcept
        {
            uint64_t current = writer_epoch.load(std::memory_order_acquire);
            uint64_t target  = current + 1;
            writer_epoch.compare_exchange_strong(current, target, std::memory_order_release);

            std::shared_lock _{ tls_readers->lock };

            for (const registered_reader* tls_reader : tls_readers->list)
            {
                while (tls_reader->epoch.load(std::memory_order_acquire) < target) GAPP_PAUSE();
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

            std::atomic<uint64_t> epoch = NOT_READING;
        };

        struct tls_reader_list
        {
            detail::shared_spinlock lock;
            std::vector<registered_reader*> list;
        };

        inline static constexpr uint64_t NOT_READING = std::numeric_limits<uint64_t>::max();

        GAPP_API inline static detail::Indestructible<tls_reader_list> tls_readers;
        GAPP_API inline static constinit std::atomic<uint64_t> writer_epoch = 0;
        alignas(128) inline static thread_local registered_reader reader;
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
            rcu_domain::synchronize();
            delete old_ptr;

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

        void wait_for_readers() const noexcept { rcu_domain::synchronize(); }

    private:
        std::atomic<T*> data_;
    };

    // NOLINTEND(*owning-memory)

} // namespace gapp::detail

#endif // !GA_UTILITY_RCU_HPP
