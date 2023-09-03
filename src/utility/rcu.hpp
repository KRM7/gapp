/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_RCU_HPP
#define GA_UTILITY_RCU_HPP

#include "utility.hpp"
#include <atomic>
#include <mutex>
#include <shared_mutex>
#include <vector>
#include <tuple>
#include <cstdint>
#include <cstddef>

namespace gapp::detail
{
    struct default_rcu_domain_tag {};

    template<typename = default_rcu_domain_tag>
    class rcu_domain
    {
    public:
        inline static void read_lock() noexcept
        {
            reader.epoch.store(writer_epoch.load(std::memory_order_relaxed), std::memory_order_release);
            std::ignore = reader.epoch.load(std::memory_order_acquire);
        }

        inline static void read_unlock() noexcept
        {
            reader.epoch.store(NOT_READING, std::memory_order_release);
        }

        inline static void synchronize() noexcept
        {
            uint64_t current = writer_epoch.load(std::memory_order_acquire);
            uint64_t target  = current + 1;
            writer_epoch.compare_exchange_strong(current, target, std::memory_order_acq_rel);

            std::shared_lock _{ reader_list_mtx };

            for (const registered_reader* reader_ : reader_list)
            {
                while (reader_->epoch.load(std::memory_order_acquire) < target) { GAPP_PAUSE(); }
            }
        }

    private:
        struct registered_reader
        {
            registered_reader() noexcept
            {
                std::unique_lock _{ reader_list_mtx };
                reader_list.push_back(this);
            }

            ~registered_reader() noexcept
            {
                std::unique_lock _{ reader_list_mtx };
                std::erase(reader_list, this);
            }

            std::atomic<uint64_t> epoch = NOT_READING;
        };

        inline static constexpr uint64_t NOT_READING = std::numeric_limits<uint64_t>::max();

        GAPP_API inline static std::vector<registered_reader*> reader_list;
        GAPP_API inline static std::shared_mutex reader_list_mtx;

        alignas(128) GAPP_API inline static constinit std::atomic<uint64_t> writer_epoch = 0;
        alignas(128) inline static thread_local registered_reader reader;
    };


    template<typename T, typename D = default_rcu_domain_tag>
    class rcu_obj
    {
    public:
        template<typename... Args>
        constexpr rcu_obj(Args&&... args) :
            data_(new T(std::forward<Args>(args)...))
        {}

        ~rcu_obj() noexcept
        {
            T* ptr = data_.load(std::memory_order_consume);
            rcu_domain<D>::synchronize();
            delete ptr;
        }

        template<typename U>
        rcu_obj& operator=(U&& value)
        {
            T* new_ptr = new T(std::forward<U>(value));
            T* old_ptr = data_.exchange(new_ptr, std::memory_order_acq_rel);
            rcu_domain<D>::synchronize();
            delete old_ptr;

            return *this;
        }

        T& get() const noexcept
        {
            return *data_.load(std::memory_order_consume);
        }

        T& operator*() const noexcept { return get(); }
        T* operator->() const noexcept { return std::addressof(get()); }

        void lock() const noexcept { rcu_domain<D>::read_lock(); }
        void unlock() const noexcept { rcu_domain<D>::read_unlock(); }

    private:
        std::atomic<T*> data_;
    };

} // namespace gapp::detail

#endif // !GA_UTILITY_RCU_HPP
