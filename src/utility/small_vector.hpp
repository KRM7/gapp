/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_UTILITY_SMALL_VECTOR_HPP
#define GAPP_UTILITY_SMALL_VECTOR_HPP

#include "scope_exit.hpp"
#include "utility.hpp"
#include <algorithm>
#include <iterator>
#include <type_traits>
#include <new>
#include <memory>
#include <initializer_list>
#include <utility>
#include <stdexcept>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <cassert>

// NOLINTBEGIN(*pointer-arithmetic, *no-malloc, *reinterpret-cast)

namespace gapp::detail
{
    //----------------------------------------- HELPER TYPE TRAITS ---------------------------------------------------

    template<typename Allocator>
    using alloc_pointer_t = typename std::allocator_traits<Allocator>::pointer;

    template<typename Allocator>
    using alloc_size_t = typename std::allocator_traits<Allocator>::size_type;


    // 'Allocator' has a trivial construct method for the type 'T' with the arguments 'TArgs'
    // if calling allocator_traits<Allocator>::construct is equivalent to directly calling the
    // constructor of T
    template<typename Allocator, typename T, typename... Args>
    concept has_construct_method = requires { std::declval<Allocator>().construct(std::declval<T*>(), std::declval<Args>()...); };

    template<typename Allocator, typename T, typename... Args>
    inline constexpr bool has_trivial_construct_v = !has_construct_method<Allocator, T, Args...>;

    // 'Allocator' has a trivial destroy method for the type 'T' if calling
    // allocator_traits<Allocator>::destroy is equivalent to directly calling the
    // destructor of T
    template<typename Allocator, typename T>
    concept has_destroy_method = requires { std::declval<Allocator>().destroy(std::declval<T*>()); };

    template<typename Allocator, typename T>
    inline constexpr bool has_trivial_destroy_v = !has_destroy_method<Allocator, T>;


    // The type 'T' is considered trivially relocatable if its move construction and move assignment
    // can be replaced with memcpy()/memmove().
    template<typename T>
    struct is_trivially_relocatable :
        std::conjunction<std::is_trivially_move_constructible<T>, std::is_trivially_move_assignable<T>, std::is_trivially_destructible<T>>
    {};

    template<typename T>
    inline constexpr bool is_trivially_relocatable_v = is_trivially_relocatable<T>::value;


    //----------------------------------------- ALLOCATE / DEALLOCATE ---------------------------------------------------

    template<typename T, typename A>
    constexpr alloc_pointer_t<A> allocate(A& allocator, alloc_size_t<A> count)
    {
        return std::allocator_traits<A>::allocate(allocator, count);
    }

    template<typename A>
    constexpr void deallocate(A& allocator, alloc_pointer_t<A> data, alloc_size_t<A> count) noexcept
    {
        std::allocator_traits<A>::deallocate(allocator, data, count);
    }

    //--------------------------- CONSTRUCT / DESTROY ONE IN UNINITIALIZED MEMORY ---------------------------------------

    template<typename T, typename A>
    constexpr void destroy(A& allocator, T* at) noexcept
    {
        if constexpr (!std::is_trivially_destructible_v<T> || !has_trivial_destroy_v<A, T>)
        {
            std::allocator_traits<A>::destroy(allocator, at);
        }
    }

    // default construct
    template<typename T, typename A>
    constexpr void construct(A& allocator, T* at)
    noexcept(std::is_nothrow_default_constructible_v<T> && has_trivial_construct_v<A, T>)
    {
        if constexpr (std::is_trivially_default_constructible_v<T> && has_trivial_construct_v<A, T>)
        {
            std::memset(at, 0, sizeof(T));
        }
        else { std::allocator_traits<A>::construct(allocator, at); }
    }

    // copy consruct
    template<typename T, typename A>
    constexpr void construct(A& allocator, T* at, const std::type_identity_t<T>& from)
    noexcept(std::is_nothrow_copy_constructible_v<T> && has_trivial_construct_v<A, T, const T&>)
    {
        std::allocator_traits<A>::construct(allocator, at, from);
    }

    // move construct
    template<typename T, typename A>
    constexpr void construct(A& allocator, T* at, std::type_identity_t<T>&& from)
    noexcept(std::is_nothrow_move_constructible_v<T> && has_trivial_construct_v<A, T, T&&>)
    {
        std::allocator_traits<A>::construct(allocator, at, std::move(from));
    }

    // construct from args
    template<typename T, typename A, typename... Args>
    constexpr void construct(A& allocator, T* at, Args&&... args)
    noexcept(std::is_nothrow_constructible_v<T, Args...> && has_trivial_construct_v<A, Args...>)
    {
        std::allocator_traits<A>::construct(allocator, at, std::forward<Args>(args)...);
    }

    //--------------------------- CONSTRUCT / DESTROY RANGE IN UNINITIALIZED MEMORY ---------------------------------------

    template<typename T, typename A>
    constexpr void destroy_range(A& allocator, T* first, T* last) noexcept
    {
        for (; first != last; ++first) detail::destroy(allocator, first);
    }

    // default construct
    template<typename T, typename A>
    constexpr void construct_range(A& allocator, T* first, T* last)
    noexcept(std::is_nothrow_default_constructible_v<T> && has_trivial_construct_v<A, T>)
    {
        if constexpr (std::is_trivially_default_constructible_v<T> && has_trivial_construct_v<A, T>)
        {
            std::memset(first, 0, sizeof(T) * (last - first));
        }
        else
        {
            T* next = first;
            scope_exit guard{ [&] { detail::destroy_range(allocator, first, next); } };
            for (; next != last; ++next) detail::construct(allocator, next);
            guard.release();
        }
    }

    // copy construct
    template<typename T, typename A>
    constexpr void construct_range(A& allocator, T* first, T* last, const T& val)
    noexcept(noexcept(detail::construct(allocator, first, val)))
    {
        T* next = first;
        scope_exit guard{ [&] { detail::destroy_range(allocator, first, next); } };
        for (; next != last; ++next) detail::construct(allocator, next, val);
        guard.release();
    }

    // construct from range
    template<typename T, typename A, std::forward_iterator Iter>
    constexpr void construct_range(A& allocator, T* first, T* last, Iter src_first)
    noexcept(noexcept(detail::construct(allocator, first, *src_first)))
    {
        using R = std::iter_reference_t<Iter>;

        if constexpr (std::contiguous_iterator<Iter> && std::is_trivially_constructible_v<T, R> && has_trivial_construct_v<A, T, R>)
        {
            std::memcpy(first, std::addressof(*src_first), sizeof(T) * (last - first));
        }
        else
        {
            T* next = first;
            scope_exit guard{ [&] { detail::destroy_range(allocator, first, next); } };
            for (; next != last; ++next, ++src_first) detail::construct(allocator, next, *src_first);
            guard.release();
        }
    }

    // move construct from another range if noexcept
    template<typename T, typename A>
    constexpr void relocate_range_strong(A& allocator, T* first, T* last, T* dest)
    noexcept(std::is_nothrow_move_constructible_v<T> && has_trivial_construct_v<A, T, T&&>)
    {
        if constexpr (is_trivially_relocatable_v<T> && has_trivial_construct_v<A, T, T&&>)
        {
            std::memcpy(dest, first, sizeof(T) * (last - first));
        }
        else
        {
            T* next = dest;
            scope_exit guard{ [&] { detail::destroy_range(allocator, dest, next); } };
            for (; first != last; ++first, ++next) detail::construct(allocator, next, std::move_if_noexcept(*first));
            guard.release();
        }
    }

    // move construct from another range
    template<typename T, typename A>
    constexpr void relocate_range_weak(A& allocator, T* first, T* last, T* dest)
    noexcept(std::is_nothrow_move_constructible_v<T> && has_trivial_construct_v<A, T, T&&>)
    {
        if constexpr (is_trivially_relocatable_v<T> && has_trivial_construct_v<A, T, T&&>)
        {
            std::memcpy(dest, first, sizeof(T) * (last - first));
        }
        else
        {
            T* next = dest;
            scope_exit guard{ [&] { detail::destroy_range(allocator, dest, next); } };
            for (; first != last; ++first, ++next) detail::construct(allocator, next, std::move(*first));
            guard.release();
        }
    }

    //---------------------------------------- ASSIGNMENT METHODS -------------------------------------------------------

    template<typename T, std::forward_iterator Iter>
    constexpr void assign_range(T* first, T* last, Iter src_first)
    noexcept(std::is_nothrow_assignable_v<T&, decltype(*src_first)>)
    {
        std::copy(src_first, src_first + std::distance(first, last), first);
    }


    //--------------------------------------- SMALL VECTOR BUFFER -------------------------------------------------------

    inline constexpr size_t cache_line_size = 64;

    template<typename T, size_t Size>
    struct small_vector_buffer
    {
    public:
        auto begin() noexcept { return std::launder(reinterpret_cast<T*>(std::addressof(data_[0]))); }
        auto begin() const noexcept { return std::launder(reinterpret_cast<const T*>(std::addressof(data_[0]))); }

        auto end() noexcept { return reinterpret_cast<T*>(std::addressof(data_[0])) + Size; }
        auto end() const noexcept { return reinterpret_cast<const T*>(std::addressof(data_[0])) + Size; }

        constexpr size_t size() const noexcept { return Size; }

    private:
        inline constexpr static size_t buffer_size = sizeof(T) * Size;

        using storage_type = unsigned char[buffer_size]; // NOLINT(*avoid-c-arrays)

        alignas(T) storage_type data_;
    };

    template<typename T>
    struct default_small_size
    {
    private:
        inline constexpr static size_t overall_size = cache_line_size;
        inline constexpr static size_t buffer_size = overall_size - 3 * sizeof(T*);
        inline constexpr static size_t buffer_min_count = 4;
    public:
        inline constexpr static size_t value = std::max(buffer_min_count, buffer_size / sizeof(T));
    };

    template<typename T>
    inline constexpr size_t default_small_size_v = default_small_size<T>::value;

} // namespace gapp::detail

namespace gapp
{
    template<typename T, size_t Size = detail::default_small_size_v<T>, typename A = std::allocator<T>>
    class small_vector
    {
    public:
        using value_type      = T;
        using allocator_type  = A;
        using reference       = T&;
        using const_reference = const T&;
        using pointer         = typename std::allocator_traits<A>::pointer;
        using const_pointer   = typename std::allocator_traits<A>::const_pointer;
        using size_type       = std::size_t;
        using difference_type = std::ptrdiff_t;

        using iterator               = T*;
        using const_iterator         = const T*;
        using reverse_iterator       = std::reverse_iterator<iterator>;
        using const_reverse_iterator = std::reverse_iterator<const_iterator>;

        //-----------------------------------//
        //            CONSTRUCTORS           //
        //-----------------------------------//

        small_vector() noexcept(std::is_nothrow_default_constructible_v<A>) :
            first_(buffer_.begin()),
            last_(buffer_.begin()),
            last_alloc_(buffer_.end())
        {}

        explicit small_vector(const A& allocator) noexcept(std::is_nothrow_copy_constructible_v<A>) :
            first_(buffer_.begin()),
            last_(buffer_.begin()),
            last_alloc_(buffer_.end()),
            alloc_(allocator)
        {}

        explicit small_vector(size_type count, const A& allocator = {}) :
            alloc_(allocator)
        {
            allocate_n(count);
            detail::scope_exit guard{ [&] { deallocate(); } };
            detail::construct_range(alloc_, first_, first_ + count);
            last_ = first_ + count;
            guard.release();
        }

        small_vector(size_type count, const T& value, const A& allocator = {}) :
            alloc_(allocator)
        {
            allocate_n(count);
            detail::scope_exit guard{ [&] { deallocate(); } };
            detail::construct_range(alloc_, first_, first_ + count, value);
            last_ = first_ + count;
            guard.release();
        }

        template<std::forward_iterator Iter>
        small_vector(Iter src_first, Iter src_last, const A& allocator = {}) :
            alloc_(allocator)
        {
            const auto src_len = std::distance(src_first, src_last);
            allocate_n(src_len);
            detail::scope_exit guard{ [&] { deallocate(); } };
            detail::construct_range(alloc_, first_, first_ + src_len, src_first);
            last_ = first_ + src_len;
            guard.release();
        }

        small_vector(std::initializer_list<T> init, const A& allocator = {}) :
            small_vector(init.begin(), init.end(), allocator)
        {}

        small_vector(const small_vector& other) :
            small_vector(other.begin(), other.end(), std::allocator_traits<A>::select_on_container_copy_construction(other.alloc_))
        {}

        small_vector(const small_vector& other, const A& allocator) :
            small_vector(other.begin(), other.end(), allocator)
        {}

        small_vector(small_vector&& other) noexcept(std::is_nothrow_move_constructible_v<T> && detail::has_trivial_construct_v<A, T, T&&>) :
            alloc_(std::move(other.alloc_))
        {
            if (other.is_small())
            {
                set_buffer_storage(0);
                detail::relocate_range_weak(alloc_, other.first_, other.last_, first_);
                detail::destroy_range(alloc_, other.first_, other.last_);
                set_buffer_storage(other.size());
                other.last_ = other.first_;
            }
            else
            {
                std::swap(first_, other.first_);
                std::swap(last_, other.last_);
                std::swap(last_alloc_, other.last_alloc_);
            }
        }

        //-----------------------------------//
        //             DESTRUCTOR            //
        //-----------------------------------//

        ~small_vector() noexcept
        {
            detail::destroy_range(alloc_, first_, last_);
            deallocate();
        }

        //-----------------------------------//
        //             ASSIGNMENT            //
        //-----------------------------------//

        template<std::forward_iterator Iter>
        void assign(Iter src_first, Iter src_last)
        {
            const auto src_size = std::distance(src_first, src_last);
            const auto old_size = std::distance(first_, last_);

            if (capacity() >= size_type(src_size))
            {
                detail::assign_range(first_, first_ + std::min(old_size, src_size), src_first);
                if (old_size < src_size) detail::construct_range(alloc_, first_ + old_size, first_ + src_size, src_first + old_size);
                else if (old_size > src_size) detail::destroy_range(alloc_, first_ + src_size, last_);
                last_ = first_ + src_size;
            }
            else
            {
                pointer new_first = detail::allocate<T>(alloc_, src_size);
                detail::scope_exit guard{ [&] { detail::deallocate(alloc_, new_first, src_size); } };
                detail::construct_range(alloc_, new_first, new_first + src_size, src_first);
                detail::destroy_range(alloc_, first_, last_);
                guard.release();
                deallocate();
                set_storage(new_first, src_size, src_size);
            }
        }

        small_vector& operator=(const small_vector& other)
        {
            if (std::addressof(other) == this) [[unlikely]] return *this;

            if constexpr (std::allocator_traits<A>::propagate_on_container_copy_assignment::value)
            {
                if (alloc_ != other.alloc_)
                {
                    reset();
                    alloc_ = other.alloc_;
                    allocate_n(other.size());
                    detail::construct_range(alloc_, first_, first_ + other.size(), other.first_);
                    last_ = first_ + other.size();
                    return *this;
                }
            }
            assign(other.begin(), other.end());
            return *this;
        }

        small_vector& operator=(small_vector&& other)
        {
            if (std::addressof(other) == this) [[unlikely]] return *this;

            if (!this->is_small() && !other.is_small())
            {
                swap(other);
                return *this;
            }
            if (this->is_small() && !other.is_small())
            {
                detail::destroy_range(alloc_, first_, last_);
                if constexpr (std::allocator_traits<A>::propagate_on_container_move_assignment::value) { alloc_ = std::move(other.alloc_); }
                this->set_storage(other.first_, other.last_, other.last_alloc_);
                other.set_storage(nullptr, nullptr, nullptr);
                return *this;
            }
            if constexpr (std::allocator_traits<A>::propagate_on_container_move_assignment::value)
            {
                if (alloc_ != other.alloc_)
                {
                    reset();
                    alloc_ = std::move(other.alloc_);
                    allocate_n(other.size());
                    detail::relocate_range_weak(alloc_, other.first_, other.last_, first_);
                    last_ = first_ + other.size();
                    detail::destroy_range(alloc_, other.first_, other.last_);
                    other.last_ = other.first_;
                    return *this;
                }
            }
            assign(std::make_move_iterator(other.first_), std::make_move_iterator(other.last_));
            return *this;
        }

        small_vector& operator=(std::initializer_list<T> list)
        {
            assign(list.begin(), list.end());
            return *this;
        }

        //-----------------------------------//
        //             ITERATORS             //
        //-----------------------------------//

        iterator begin() noexcept { return first_; }
        const_iterator begin() const noexcept { return first_; }
        const_iterator cbegin() const noexcept { return first_; }

        iterator end() noexcept { return last_; }
        const_iterator end() const noexcept { return last_; }
        const_iterator cend() const noexcept { return last_; }

        reverse_iterator rbegin() noexcept { return std::make_reverse_iterator(last_); }
        const_reverse_iterator rbegin() const noexcept { return std::make_reverse_iterator(last_); }
        const_reverse_iterator crbegin() const noexcept { return std::make_reverse_iterator(last_); }

        reverse_iterator rend() noexcept { return std::make_reverse_iterator(first_); }
        const_reverse_iterator rend() const noexcept { return std::make_reverse_iterator(first_); }
        const_reverse_iterator crend() const noexcept { return std::make_reverse_iterator(first_); }

        //-----------------------------------//
        //           ELEMENT ACCESS          //
        //-----------------------------------//

        reference operator[](size_type pos) noexcept
        {
            GAPP_ASSERT(pos < size());
            return first_[pos];
        }

        const_reference operator[](size_type pos) const noexcept
        {
            GAPP_ASSERT(pos < size());
            return first_[pos];
        }

        reference at(size_type pos)
        {
            if (pos >= size()) GAPP_THROW(std::out_of_range, "Bad vector index.");
            return first_[pos];
        }

        const_reference at(size_type pos) const
        {
            if (pos >= size()) GAPP_THROW(std::out_of_range, "Bad vector index.");
            return first_[pos];
        }

        reference front() noexcept { GAPP_ASSERT(!empty()); return *first_; }
        const_reference front() const noexcept { GAPP_ASSERT(!empty()); return *first_; }

        reference back() noexcept { GAPP_ASSERT(!empty()); return *(last_ - 1); }
        const_reference back() const noexcept { GAPP_ASSERT(!empty()); return *(last_ - 1); }

        pointer data() noexcept { return first_; }
        const_pointer data() const noexcept { return first_; }

        //-----------------------------------//
        //              CAPACITY             //
        //-----------------------------------//

        bool empty() const noexcept { return first_ == last_; }

        size_type size() const noexcept { return size_type(last_ - first_); }
        size_type capacity() const noexcept { return size_type(last_alloc_ - first_); }
        size_type max_size() const noexcept { return std::allocator_traits<A>::max_size(alloc_); }

        bool is_small() const noexcept { return first_ == buffer_.begin(); }
        size_type small_capacity() const noexcept { return Size; }

        void reserve(size_type new_capacity) { if (new_capacity > capacity()) reallocate_n(new_capacity); }
        void shrink_to_fit() {}

        //-----------------------------------//
        //             MODIFIERS             //
        //-----------------------------------//

        void clear() noexcept
        {
            detail::destroy_range(alloc_, first_, last_);
            last_ = first_;
        }

        void reset() noexcept
        {
            detail::destroy_range(alloc_, first_, last_);
            deallocate();
            set_buffer_storage(0);
        }

        void swap(small_vector& other)
        noexcept(std::is_nothrow_swappable_v<T> && std::is_nothrow_move_constructible_v<T> && detail::has_trivial_construct_v<A, T, T&&>)
        {
            if (std::addressof(other) == this) [[unlikely]] return;

            if (!this->is_small() && !other.is_small())
            {
                std::swap(first_, other.first_);
                std::swap(last_, other.last_);
                std::swap(last_alloc_, other.last_alloc_);
            }
            else if (this->is_small() && other.is_small())
            {
                small_vector& big   = (this->size() < other.size()) ? other : *this;
                small_vector& small = (this->size() < other.size()) ? *this : other;
                const auto small_size = small.size();

                std::swap_ranges(small.first_, small.last_, big.first_);
                detail::relocate_range_strong(small.alloc_, big.first_ + small.size(), big.last_, small.last_);
                detail::destroy_range(big.alloc_, big.first_ + small.size(), big.last_);

                small.set_buffer_storage(big.size());
                big.set_buffer_storage(small_size);
            }
            else
            {
                small_vector& big   = this->is_small() ? other : *this;
                small_vector& small = this->is_small() ? *this : other;
                const auto small_size = small.size();

                detail::relocate_range_strong(big.alloc_, small.first_, small.last_, big.buffer_.begin());
                detail::destroy_range(small.alloc_, small.first_, small.last_);

                small.set_storage(big.first_, big.last_, big.last_alloc_);
                big.set_buffer_storage(small_size);
            }

            if constexpr (std::allocator_traits<A>::propagate_on_container_swap::value)
            {
                using std::swap;
                swap(alloc_, other.alloc_);
            }
        }

        void push_back(const T& value) { emplace_back(value); }
        void push_back(T&& value) { emplace_back(std::move(value)); }

        template<typename... Args>
        reference emplace_back(Args&&... args)
        {
            if (size() == capacity())
                reallocate_append(next_capacity(), std::forward<Args>(args)...);
            else
                detail::construct(alloc_, last_++, std::forward<Args>(args)...);
            return back();
        }

        void pop_back() noexcept { GAPP_ASSERT(!empty()); detail::destroy(alloc_, last_--); }

        void resize(size_type count) { resize_impl(count); }
        void resize(size_type count, const T& value) { resize_impl(count, value); }

        iterator insert(const_iterator pos, const T& value) { return emplace(pos, value); }
        iterator insert(const_iterator pos, T&& value) { return emplace(pos, std::move(value)); }
        iterator insert(const_iterator pos, std::initializer_list<T> list) { return insert(pos, list.begin(), list.end()); }
        iterator erase(const_iterator pos) noexcept(std::is_nothrow_move_assignable_v<T>) { GAPP_ASSERT(pos != end()); return erase(pos, pos + 1); }

        template<typename... Args>
        iterator emplace(const_iterator pos, Args&&... args)
        {
            if (pos == cend()) return std::addressof(emplace_back(std::forward<Args>(args)...));

            T new_elem{ std::forward<Args>(args)... };
            const auto offset = std::distance(cbegin(), pos);

            if (size() == capacity()) reallocate_n(next_capacity());
            detail::construct(alloc_, last_, std::move(back()));
            std::shift_right(first_ + offset, last_++, 1);
            *(first_ + offset) = std::move(new_elem);
            return first_ + offset;
        }

        template<std::forward_iterator Iter>
        iterator insert(const_iterator pos, Iter src_first, Iter src_last)
        {
            const auto offset = std::distance(cbegin(), pos);
            const auto src_size = std::distance(src_first, src_last);

            reserve(size() + src_size);

            const auto middle = std::max(last_ - src_size, first_ + offset);
            const auto moved_size = last_ - middle;
            const auto new_last = last_ + src_size;
            const auto new_middle = last_ + src_size - moved_size;

            detail::relocate_range_weak(alloc_, middle, last_, new_middle);
            detail::scope_exit guard{ [&] { detail::destroy_range(alloc_, new_middle, new_last); } };
            detail::assign_range(middle, last_, src_first);
            detail::construct_range(alloc_, last_, new_middle, src_first + moved_size);
            last_ = new_last;
            guard.release();

            return first_ + offset;
        }

        iterator erase(const_iterator first, const_iterator last) noexcept(std::is_nothrow_move_assignable_v<T>)
        {
            const auto erase_first = first_ + std::distance(cbegin(), first);
            const auto erase_count = std::distance(first, last);
            const auto new_last = std::shift_left(erase_first, last_, erase_count);
            detail::destroy_range(alloc_, new_last, last_);
            last_ = new_last;
            return erase_first;
        }

        //-----------------------------------//
        //               OTHER               //
        //-----------------------------------//

        allocator_type get_allocator() const noexcept(std::is_nothrow_copy_constructible_v<A>) { return alloc_; }

        friend bool operator==(const small_vector& lhs, const small_vector& rhs) noexcept
        {
            return std::equal(lhs.begin(), lhs.end(), rhs.begin(), rhs.end());
        }

        friend auto operator<=>(const small_vector& lhs, const small_vector& rhs) noexcept
        {
            return std::lexicographical_compare_three_way(lhs.begin(), lhs.end(), rhs.begin(), rhs.end());
        }

    private:

        detail::small_vector_buffer<T, Size> buffer_;
        pointer first_ = nullptr;
        pointer last_ = nullptr;
        pointer last_alloc_ = nullptr;
        GAPP_NO_UNIQUE_ADDRESS allocator_type alloc_;

        static inline constexpr double growth_factor_ = 1.618;


        void allocate_n(size_type count)
        {
            if (count <= buffer_.size())
            {
                set_buffer_storage(0);
            }
            else
            {
                first_ = detail::allocate<T>(alloc_, count);
                last_alloc_ = first_ + count;
            }
        }

        void reallocate_n(size_type new_capacity)
        {
            GAPP_ASSERT(new_capacity > capacity());

            const size_type old_size = size();

            pointer new_first = detail::allocate<T>(alloc_, new_capacity);
            detail::scope_exit guard{ [&] { detail::deallocate(alloc_, new_first, new_capacity); } };
            detail::relocate_range_strong(alloc_, first_, last_, new_first);
            detail::destroy_range(alloc_, first_, last_);
            deallocate();
            guard.release();
            set_storage(new_first, old_size, new_capacity);
        }

        template<typename... Args>
        void reallocate_append(size_type new_capacity, Args&&... args)
        {
            GAPP_ASSERT(new_capacity > capacity());

            const size_type old_size = size();

            pointer new_first = detail::allocate<T>(alloc_, new_capacity);
            detail::scope_exit guard1{ [&] { detail::deallocate(alloc_, new_first, new_capacity); } };
            detail::construct(alloc_, new_first + old_size, std::forward<Args>(args)...);
            detail::scope_exit guard2{ [&] { detail::destroy(alloc_, new_first + old_size); } };
            detail::relocate_range_strong(alloc_, first_, last_, new_first);
            { guard1.release(); guard2.release(); }
            detail::destroy_range(alloc_, first_, last_);
            deallocate();
            set_storage(new_first, old_size + 1, new_capacity);
        }

        void deallocate() noexcept
        {
            if (!is_small() && data()) detail::deallocate(alloc_, first_, capacity());
        }

        template<typename... Args>
        void resize_impl(size_type count, Args&&... args)
        {
            if (count < size())
            {
                detail::destroy_range(alloc_, first_ + count, last_);
                last_ = first_ + count;
            }
            else if (count > size())
            {
                if (count > capacity()) reallocate_n(count);
                detail::construct_range(alloc_, last_, first_ + count, std::forward<Args>(args)...);
                last_ = first_ + count;
            }
        }

        void set_storage(pointer first, pointer last, pointer last_alloc) noexcept
        {
            first_ = first;
            last_ = last;
            last_alloc_ = last_alloc;
        }

        void set_storage(pointer first, size_type size, size_type capacity) noexcept
        {
            set_storage(first, first + size, first + capacity);
        }

        void set_buffer_storage(size_type size) noexcept
        {
            set_storage(buffer_.begin(), buffer_.begin() + size, buffer_.end());
        }

        size_type next_capacity() const noexcept
        {
            return size_type(growth_factor_ * capacity()) + 1;
        }

    }; // class small_vector

    template<std::forward_iterator Iter, std::size_t Size = detail::default_small_size_v<std::iter_value_t<Iter>>, typename Alloc = std::allocator<std::iter_value_t<Iter>>>
    small_vector(Iter, Iter, Alloc = {}) -> small_vector<std::iter_value_t<Iter>, Size, Alloc>;

    template<typename T, std::size_t Size, typename A>
    void swap(small_vector<T, Size, A>& lhs, small_vector<T, Size, A>& rhs)
    noexcept(noexcept(lhs.swap(rhs)))
    {
        lhs.swap(rhs);
    }

} // namespace gapp

// NOLINTEND(*pointer-arithmetic, *no-malloc, *reinterpret-cast)

#endif // !GAPP_UTILITY_SMALL_VECTOR_HPP
