/* Copyright (c) 2024 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_UTILITY_CIRCULAR_BUFFER_HPP
#define GAPP_UTILITY_CIRCULAR_BUFFER_HPP

#include "small_vector.hpp"
#include "iterators.hpp"
#include "scope_exit.hpp"
#include "utility.hpp"
#include <algorithm>
#include <iterator>
#include <initializer_list>
#include <type_traits>
#include <compare>
#include <memory>
#include <stdexcept>
#include <utility>
#include <cstddef>

namespace gapp::detail
{
    template<typename T, typename A = std::allocator<T>>
    class circular_buffer
    {
    public:
        using value_type      = T;
        using allocator_type  = A;
        using reference       = T&;
        using const_reference = const T&;
        using pointer         = T*;
        using const_pointer   = const T*;
        using size_type       = std::size_t;
        using difference_type = std::ptrdiff_t;

        using iterator               = detail::stable_iterator<circular_buffer>;
        using const_iterator         = detail::const_stable_iterator<circular_buffer>;
        using reverse_iterator       = std::reverse_iterator<iterator>;
        using const_reverse_iterator = std::reverse_iterator<const_iterator>;

        static_assert(std::allocator_traits<A>::is_always_equal::value);

        //-----------------------------------//
        //            CONSTRUCTORS           //
        //-----------------------------------//

        constexpr explicit circular_buffer(A alloc = A()) noexcept :
            alloc_(std::move(alloc))
        {}

        constexpr explicit circular_buffer(size_type capacity, A alloc = A()) :
            buffer_(detail::allocate<T>(alloc, capacity).data),
            capacity_(capacity),
            alloc_(std::move(alloc))
        {}

        template<std::forward_iterator Iter>
        constexpr circular_buffer(size_type capacity, Iter first, Iter last, A alloc = A()) :
            capacity_(capacity),
            alloc_(std::move(alloc))
        {
            buffer_ = detail::allocate<T>(alloc_, capacity_).data;
            detail::scope_exit guard{ [&] { detail::deallocate(alloc_, buffer_, capacity_); } };
            const size_type size = std::min<size_type>(capacity, std::distance(first, last));
            detail::construct_range(alloc_, buffer_, buffer_ + size, std::prev(last, size));
            size_ = size;
            guard.release();
        }

        constexpr circular_buffer(size_type capacity, std::initializer_list<T> ilist, A alloc = A()) :
            circular_buffer(capacity, ilist.begin(), ilist.end(), std::move(alloc))
        {}

        constexpr circular_buffer(const circular_buffer& other) :
            circular_buffer(other.capacity(), other.begin(), other.end(), std::allocator_traits<A>::select_on_container_copy_construction(other.alloc_))
        {}

        constexpr circular_buffer(circular_buffer&& other) noexcept
        {
            swap(other);
        }

        //-----------------------------------//
        //             DESTRUCTOR            //
        //-----------------------------------//

        constexpr ~circular_buffer() noexcept
        { 
            destroy();
        }

        //-----------------------------------//
        //             ASSIGNMENT            //
        //-----------------------------------//

        constexpr circular_buffer& operator=(circular_buffer other) noexcept
        {
            swap(other);
            return *this;
        }

        //-----------------------------------//
        //             ITERATORS             //
        //-----------------------------------//

        constexpr iterator begin() noexcept { return iterator(this, 0); }
        constexpr iterator end() noexcept { return iterator(this, size()); }

        constexpr const_iterator begin() const noexcept { return const_iterator(this, 0); }
        constexpr const_iterator end() const noexcept { return const_iterator(this, size()); }

        constexpr const_iterator cbegin() const noexcept { return begin(); }
        constexpr const_iterator cend() const noexcept { return end(); }

        constexpr reverse_iterator rbegin() noexcept { return std::reverse_iterator(end()); }
        constexpr reverse_iterator rend() noexcept { return std::reverse_iterator(begin()); }

        constexpr const_reverse_iterator rbegin() const noexcept { return std::reverse_iterator(end()); }
        constexpr const_reverse_iterator rend() const noexcept { return std::reverse_iterator(begin()); }

        constexpr const_reverse_iterator crbegin() const noexcept { return rbegin(); }
        constexpr const_reverse_iterator crend() const noexcept { return rend(); }

        //-----------------------------------//
        //              CAPACITY             //
        //-----------------------------------//

        constexpr size_type size() const noexcept { return size_; }
        constexpr size_type capacity() const noexcept { return capacity_; }
        constexpr size_type max_size() const noexcept { return std::allocator_traits<A>::max_size(alloc_); }

        constexpr bool empty() const noexcept { return size() == 0; }
        constexpr bool full() const noexcept { return size() == capacity(); }

        //-----------------------------------//
        //           ELEMENT ACCESS          //
        //-----------------------------------//

        constexpr reference operator[](size_type pos)
        {
            GAPP_ASSERT(pos < size());
            return *pointer_to(pos);
        }

        constexpr const_reference operator[](size_type pos) const
        {
            GAPP_ASSERT(pos < size());
            return *pointer_to(pos);
        }

        constexpr reference at(size_type pos)
        {
            if (pos >= size()) GAPP_THROW(std::out_of_range, "Bad buffer index.");
            return (*this)[pos];
        }

        constexpr const_reference at(size_type pos) const
        {
            if (pos >= size()) GAPP_THROW(std::out_of_range, "Bad buffer index.");
            return (*this)[pos];
        }

        constexpr reference front()
        {
            GAPP_ASSERT(!empty());
            return buffer_[first_];
        }

        constexpr const_reference front() const
        {
            GAPP_ASSERT(!empty());
            return buffer_[first_];
        }

        constexpr reference back()
        {
            GAPP_ASSERT(!empty());
            return (*this)[size() - 1];
        }

        constexpr const_reference back() const
        {
            GAPP_ASSERT(!empty());
            return (*this)[size() - 1];
        }

        //-----------------------------------//
        //             MODIFIERS             //
        //-----------------------------------//

        constexpr void push_back(const value_type& value)
        noexcept(std::is_nothrow_copy_constructible_v<T>)
        {
            GAPP_ASSERT(capacity() > 0);
            emplace_back(value);
        }

        constexpr void push_back(value_type&& value)
        noexcept(std::is_nothrow_move_constructible_v<T>)
        {
            GAPP_ASSERT(capacity() > 0);
            emplace_back(std::move(value));
        }

        constexpr void push_front(const value_type& value)
        noexcept(std::is_nothrow_move_constructible_v<T>)
        {
            GAPP_ASSERT(capacity() > 0);
            emplace_front(value);
        }

        constexpr void push_front(value_type&& value)
        noexcept(std::is_nothrow_move_constructible_v<T>)
        {
            GAPP_ASSERT(capacity() > 0);
            emplace_front(std::move(value));
        }

        template<typename... Args>
        constexpr reference emplace_back(Args&&... args)
        noexcept(std::is_nothrow_constructible_v<T, Args&&...>)
        {
            GAPP_ASSERT(capacity() > 0);

            if (full())
            {
                pointer last = buffer_ + first_;
                *last = T(std::forward<Args>(args)...);
                first_ = next_idx(first_);
                return *last;
            }

            pointer last = pointer_to(size_);
            detail::construct(alloc_, last, std::forward<Args>(args)...);
            size_++;
            return *last;
        }

        template<typename... Args>
        constexpr reference emplace_front(Args&&... args)
        noexcept(std::is_nothrow_constructible_v<T, Args&&...>)
        {
            GAPP_ASSERT(capacity() > 0);

            const size_type new_first = prev_idx(first_);

            if (full())
            {
                buffer_[new_first] = T(std::forward<Args>(args)...);
                first_ = new_first;
            }
            else
            {
                detail::construct(alloc_, buffer_ + new_first, std::forward<Args>(args)...);
                first_ = new_first;
                size_++;
            }

            return buffer_[first_];
        }

        constexpr void pop_front() noexcept
        {
            GAPP_ASSERT(!empty());
            detail::destroy(alloc_, std::addressof(front()));
            first_ = next_idx(first_);
            size_--;
        }

        constexpr void pop_back() noexcept
        {
            GAPP_ASSERT(!empty());
            detail::destroy(alloc_, std::addressof(back()));
            size_--;
        }

        constexpr void set_capacity(size_type new_capacity)
        {
            if (new_capacity == capacity()) return;

            pointer new_buffer = detail::allocate<T>(alloc_, new_capacity).data;
            detail::scope_exit guard{ [&] { detail::deallocate(alloc_, buffer_, new_capacity); } };
            const size_type new_size = std::min<size_type>(new_capacity, size());
            std::uninitialized_copy_n(detail::make_move_iterator_if_noexcept(begin()), new_size, new_buffer);
            destroy();
            buffer_ = new_buffer;
            size_ = new_size;
            capacity_ = new_capacity;
            guard.release();
        }

        constexpr void clear() noexcept
        {
            for (T& elem : *this) detail::destroy(alloc_, std::addressof(elem));
            first_ = 0;
            size_  = 0;
        }

        constexpr void reset(size_type new_capacity)
        {
            clear();
            if (new_capacity == capacity_) return;
            if (buffer_) detail::deallocate(alloc_, buffer_, capacity_);
            buffer_ = detail::allocate<T>(alloc_, new_capacity).data;
            capacity_ = new_capacity;
        }

        constexpr void swap(circular_buffer& other) noexcept
        {
            using std::swap;
            swap(buffer_, other.buffer_);
            swap(first_, other.first_);
            swap(size_, other.size_);
            swap(capacity_, other.capacity_);
            swap(alloc_, other.alloc_);
        }

        //-----------------------------------//
        //               OTHER               //
        //-----------------------------------//

        constexpr allocator_type get_allocator() const noexcept(std::is_nothrow_copy_constructible_v<A>)
        {
            return alloc_;
        }

        constexpr friend bool operator==(const circular_buffer& lhs, const circular_buffer& rhs) noexcept
        {
            return std::equal(lhs.begin(), lhs.end(), rhs.begin(), rhs.end());
        }

        constexpr friend auto operator<=>(const circular_buffer& lhs, const circular_buffer& rhs) noexcept
        {
            return std::lexicographical_compare_three_way(lhs.begin(), lhs.end(), rhs.begin(), rhs.end());
        }

    private:
        pointer buffer_     = nullptr;
        size_type first_    = 0;
        size_type size_     = 0;
        size_type capacity_ = 0;
        GAPP_NO_UNIQUE_ADDRESS A alloc_;

        constexpr size_type next_idx(size_type idx) const noexcept
        {
            return detail::next_mod(idx, capacity_);
        }

        constexpr size_type prev_idx(size_type idx) const noexcept
        {
            return detail::prev_mod(idx, capacity_);
        }

        constexpr pointer pointer_to(size_type pos) const noexcept
        {
            const size_type idx = first_ + pos;
            return idx >= capacity_ ? (buffer_ + idx - capacity_) : (buffer_ + idx);
        }

        constexpr void destroy() noexcept
        {
            if (!buffer_) return;
            clear();
            detail::deallocate(alloc_, buffer_, capacity_);
        }
    };

    template<typename T, typename A>
    constexpr void swap(circular_buffer<T, A>& lhs, circular_buffer<T, A>& rhs) noexcept
    {
        lhs.swap(rhs);
    }

    template<std::forward_iterator Iter, typename Alloc = std::allocator<std::iter_value_t<Iter>>>
    circular_buffer(Iter, Iter, Alloc = Alloc()) -> circular_buffer<std::iter_value_t<Iter>, Alloc>;

} // namespace gapp::detail

#endif // !GAPP_UTILITY_CIRCULAR_BUFFER_HPP
