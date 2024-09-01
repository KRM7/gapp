/* Copyright (c) 2024 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_UTILITY_CACHE_HPP
#define GAPP_UTILITY_CACHE_HPP

#include "circular_buffer.hpp"
#include "concepts.hpp"
#include "utility.hpp"
#include <unordered_map>
#include <type_traits>
#include <memory>
#include <utility>
#include <cstddef>

namespace gapp::detail
{
    template<typename Key, typename Value>
    class fifo_cache
    {
    public:
        using key_type        = Key;
        using value_type      = Value;
        using reference       = Value&;
        using const_reference = const Value&;
        using pointer         = Value*;
        using const_pointer   = const Value*;
        using size_type       = std::size_t;
        using difference_type = std::ptrdiff_t;

        static_assert(detail::hashable<Key> && std::equality_comparable<Key>);

        //-----------------------------------//
        //            CONSTRUCTORS           //
        //-----------------------------------//

        constexpr fifo_cache() noexcept = default;

        constexpr fifo_cache(size_type capacity) :
            cache_(capacity + 1),
            order_(capacity)
        {}

        constexpr fifo_cache(const fifo_cache& other) :
            cache_(other.capacity() + 1),
            order_(other.capacity())
        {
            for (auto it : other.order_)
            {
                order_.push_back(cache_.insert(*it).first);
            }
        }

        constexpr fifo_cache(fifo_cache&& other) noexcept
        {
            swap(other);
        }

        constexpr fifo_cache& operator=(fifo_cache other) noexcept
        {
            swap(other);
            return *this;
        }

        //-----------------------------------//
        //              CAPACITY             //
        //-----------------------------------//

        constexpr size_type size() const noexcept { return cache_.size(); }
        constexpr size_type capacity() const noexcept { return order_.capacity(); }

        constexpr bool empty() const noexcept { return cache_.empty(); }
        constexpr bool full() const noexcept { return size() == capacity(); }

        //-----------------------------------//
        //               LOOKUP              //
        //-----------------------------------//

        constexpr bool contains(const key_type& key) const
        {
            return cache_.contains(key);
        }

        constexpr pointer get(const key_type& key)
        {
            if (empty()) return nullptr;
            const auto it = cache_.find(key);
            if (it == cache_.end()) return nullptr;
            return std::addressof(it->second);
        }

        constexpr const_pointer get(const key_type& key) const
        {
            if (empty()) return nullptr;
            const auto it = cache_.find(key);
            if (it == cache_.end()) return nullptr;
            return std::addressof(it->second);
        }

        //-----------------------------------//
        //             MODIFIERS             //
        //-----------------------------------//

        template<typename K, typename V>
        constexpr void insert(K&& key, V&& value)
        {
            GAPP_ASSERT(free_capacity() >= 1);

            if (capacity() == 0) [[unlikely]] return;

            const auto [it, inserted] = cache_.insert_or_assign(std::forward<K>(key), std::forward<V>(value));
            if (!inserted) return;

            if (order_.full()) cache_.erase(order_.front());
            order_.push_back(it);
        }

        template<typename K, typename... Args>
        constexpr void try_insert(K&& key, Args&&... vargs)
        {
            GAPP_ASSERT(free_capacity() >= 1);

            if (capacity() == 0) [[unlikely]] return;

            const auto [it, inserted] = cache_.try_emplace(std::forward<K>(key), std::forward<Args>(vargs)...);
            if (!inserted) return;

            if (order_.full()) cache_.erase(order_.front());
            order_.push_back(it);
        }

        template<std::forward_iterator Iter, std::invocable<std::iter_reference_t<Iter>> F>
        constexpr void insert(Iter first, Iter last, F&& f)
        {
            GAPP_ASSERT(free_capacity() >= 1);

            if (capacity() == 0) [[unlikely]] return;

            const size_type range_len = std::distance(first, last);
            const size_type insert_count = std::min(capacity(), range_len);
            first = std::prev(last, insert_count);

            for (; first != last; ++first)
            {
                const auto [it, inserted] = cache_.insert_or_assign(*first, std::invoke(f, *first));
                if (!inserted) continue;

                if (order_.full()) cache_.erase(order_.front());
                order_.push_back(it);
            }
        }

        constexpr void clear() noexcept
        {
            cache_.clear();
            order_.clear();
        }

        constexpr void reset(size_type new_capacity)
        {
            clear();
            cache_.reserve(new_capacity + 1);
            order_.reset(new_capacity);
        }

        constexpr void swap(fifo_cache& other) noexcept
        {
            cache_.swap(other.cache_);
            order_.swap(other.order_);
        }

        //-----------------------------------//
        //               OTHER               //
        //-----------------------------------//

        constexpr friend bool operator==(const fifo_cache& lhs, const fifo_cache& rhs)
        {
            if (lhs.size() != rhs.size()) return false;

            for (size_type i = 0; i < lhs.size(); i++)
            {
                if (*lhs.order_[i] != *rhs.order_[i]) return false;
            }
            return true;
        }

    private:
        using Iter = typename std::unordered_map<Key, Value>::iterator;

        std::unordered_map<Key, Value> cache_;
        detail::circular_buffer<Iter> order_;

        constexpr size_type free_capacity() const noexcept
        {
            GAPP_ASSERT(cache_.max_load_factor() == 1.0);
            GAPP_ASSERT(cache_.size() <= cache_.bucket_count());

            return cache_.bucket_count() - cache_.size();
        }
    };

    template<typename K, typename V>
    constexpr void swap(fifo_cache<K, V>& lhs, fifo_cache<K, V>& rhs) noexcept
    {
        lhs.swap(rhs);
    }

} // namespace gapp::detail

#endif // !GAPP_UTILITY_CACHE_HPP
