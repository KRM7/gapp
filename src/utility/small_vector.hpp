/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_UTILITY_SMALL_VECTOR_HPP
#define GAPP_UTILITY_SMALL_VECTOR_HPP

#include "scope_exit.hpp"
#include "utility.hpp"
#include <array>
#include <algorithm>
#include <iterator>
#include <initializer_list>
#include <memory>
#include <new>
#include <type_traits>
#include <utility>
#include <stdexcept>
#include <cstring>
#include <cstdlib>
#include <cstddef>

// NOLINTBEGIN(*pointer-arithmetic, *union-access)

namespace gapp::detail
{
    #if __cpp_lib_allocate_at_least
        inline constexpr bool has_allocate_at_least = true;
    #else
        inline constexpr bool has_allocate_at_least = false;
    #endif

    //----------------------------------------- HELPER TYPE TRAITS ------------------------------------------------------

    template<typename Allocator>
    using alloc_pointer_t = typename std::allocator_traits<Allocator>::pointer;

    template<typename Allocator>
    using alloc_size_t = typename std::allocator_traits<Allocator>::size_type;


    // 'Allocator' has a trivial construct method for the type 'T' with the arguments 'Args'
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


    template<typename Allocator>
    inline constexpr bool copy_allocators_v = std::allocator_traits<Allocator>::propagate_on_container_copy_assignment::value &&
        !std::allocator_traits<Allocator>::is_always_equal::value;

    template<typename Allocator>
    inline constexpr bool move_allocators_v = std::allocator_traits<Allocator>::propagate_on_container_move_assignment::value &&
        !std::allocator_traits<Allocator>::is_always_equal::value;

    template<typename Allocator>
    inline constexpr bool swap_allocators_v = std::allocator_traits<Allocator>::propagate_on_container_swap::value &&
        !std::allocator_traits<Allocator>::is_always_equal::value;

    template<typename Allocator>
    inline constexpr bool steal_pointers_v = std::allocator_traits<Allocator>::propagate_on_container_move_assignment::value ||
        std::allocator_traits<Allocator>::is_always_equal::value;


    //----------------------------------------- ALLOCATE / DEALLOCATE ---------------------------------------------------

    template<typename A>
    struct alloc_result_t { alloc_pointer_t<A> data; std::size_t size; };

    template<typename T, typename A>
    constexpr alloc_result_t<A> allocate(A& allocator, alloc_size_t<A> count) requires(!has_allocate_at_least)
    {
        return { std::allocator_traits<A>::allocate(allocator, count), count };
    }

    template<typename T, typename A>
    constexpr alloc_result_t<A> allocate(A& allocator, alloc_size_t<A> count) requires(has_allocate_at_least)
    {
        auto [data, size] = std::allocator_traits<A>::allocate_at_least(allocator, count);
        return { data, size };
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
        if constexpr (!std::is_trivially_destructible_v<T> || !has_trivial_destroy_v<A&, T>)
        {
            std::allocator_traits<A>::destroy(allocator, at);
        }
    }

    // default construct
    template<typename T, typename A>
    constexpr void construct(A& allocator, T* at)
    noexcept(std::is_nothrow_default_constructible_v<T> && has_trivial_construct_v<A&, T>)
    {
        std::allocator_traits<A>::construct(allocator, at);
    }

    // copy construct
    template<typename T, typename A>
    constexpr void construct(A& allocator, T* at, const std::type_identity_t<T>& from)
    noexcept(std::is_nothrow_copy_constructible_v<T> && has_trivial_construct_v<A&, T, const T&>)
    {
        std::allocator_traits<A>::construct(allocator, at, from);
    }

    // move construct
    template<typename T, typename A>
    constexpr void construct(A& allocator, T* at, std::type_identity_t<T>&& from)
    noexcept(std::is_nothrow_move_constructible_v<T> && has_trivial_construct_v<A&, T, T&&>)
    {
        std::allocator_traits<A>::construct(allocator, at, std::move(from));
    }

    // construct from args
    template<typename T, typename A, typename... Args>
    constexpr void construct(A& allocator, T* at, Args&&... args)
    noexcept(std::is_nothrow_constructible_v<T, Args...> && has_trivial_construct_v<A&, Args...>)
    {
        std::allocator_traits<A>::construct(allocator, at, std::forward<Args>(args)...);
    }

    //--------------------------- CONSTRUCT / DESTROY RANGE IN UNINITIALIZED MEMORY -------------------------------------

    template<typename T, typename A>
    constexpr void destroy_range(A& allocator, T* first, T* last) noexcept
    {
        for (; first != last; ++first) detail::destroy(allocator, first);
    }

    // default construct
    template<typename T, typename A>
    constexpr void construct_range(A& allocator, T* first, T* last)
    noexcept(std::is_nothrow_default_constructible_v<T> && has_trivial_construct_v<A&, T>)
    {
        if constexpr (std::is_trivially_default_constructible_v<T> && has_trivial_construct_v<A&, T>)
        {
            if (!std::is_constant_evaluated())
            {
                std::memset(first, 0, sizeof(T) * (last - first));
                return;
            }
        }

        T* next = first;
        scope_exit guard{ [&] { detail::destroy_range(allocator, first, next); } };
        for (; next != last; ++next) detail::construct(allocator, next);
        guard.release();
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
    template<typename T, typename A, std::input_iterator Iter>
    constexpr void construct_range(A& allocator, T* first, T* last, Iter src_first)
    noexcept(noexcept(detail::construct(allocator, first, *src_first)))
    {
        using R = std::iter_reference_t<Iter>;

        if constexpr (std::contiguous_iterator<Iter> && std::is_trivially_constructible_v<T, R> && has_trivial_construct_v<A&, T, R>)
        {
            if (!std::is_constant_evaluated())
            {
                std::memcpy(first, std::to_address(src_first), sizeof(T) * (last - first));
                return;
            }
        }

        T* next = first;
        scope_exit guard{ [&] { detail::destroy_range(allocator, first, next); } };
        for (; next != last; ++next, ++src_first) detail::construct(allocator, next, *src_first);
        guard.release();
    }

    // move construct from another range if noexcept
    template<typename T, typename A>
    constexpr void relocate_range_strong(A& allocator, T* first, T* last, T* dest)
    noexcept(std::is_nothrow_move_constructible_v<T> && has_trivial_construct_v<A&, T, T&&>)
    {
        if constexpr (is_trivially_relocatable_v<T> && has_trivial_construct_v<A&, T, T&&>)
        {
            if (!std::is_constant_evaluated())
            {
                std::memcpy(dest, first, sizeof(T) * (last - first));
                return;
            }
        }

        T* next = dest;
        scope_exit guard{ [&] { detail::destroy_range(allocator, dest, next); } };
        for (; first != last; ++first, ++next) detail::construct(allocator, next, std::move_if_noexcept(*first));
        guard.release();
    }

    // move construct from another range
    template<typename T, typename A>
    constexpr void relocate_range_weak(A& allocator, T* first, T* last, T* dest)
    noexcept(std::is_nothrow_move_constructible_v<T> && has_trivial_construct_v<A&, T, T&&>)
    {
        if constexpr (is_trivially_relocatable_v<T> && has_trivial_construct_v<A&, T, T&&>)
        {
            if (!std::is_constant_evaluated())
            {
                std::memcpy(dest, first, sizeof(T) * (last - first));
                return;
            }
        }

        T* next = dest;
        scope_exit guard{ [&] { detail::destroy_range(allocator, dest, next); } };
        for (; first != last; ++first, ++next) detail::construct(allocator, next, std::move(*first));
        guard.release();
    }

    //---------------------------------------- ASSIGNMENT METHODS -------------------------------------------------------

    template<typename T, std::forward_iterator Iter>
    constexpr void assign_range(T* first, T* last, Iter src_first)
    noexcept(std::is_nothrow_assignable_v<T&, std::iter_reference_t<Iter>>)
    {
        std::copy(src_first, src_first + std::distance(first, last), first);
    }

    template<typename T>
    constexpr void assign_range(T* first, T* last, const T& value)
    noexcept(std::is_nothrow_copy_assignable_v<T>)
    {
        std::fill(first, last, value);
    }

    //------------------------------------ ALLOCATOR MANGAGED OBJECT ----------------------------------------------------

    template<typename T, typename Allocator>
    class allocator_managed
    {
    public:
        template<typename... Args>
        constexpr allocator_managed(Allocator& alloc, Args&&... args) :
            alloc_(alloc)
        {
            detail::construct(alloc_, std::addressof(**this), std::forward<Args>(args)...);
        }

        constexpr ~allocator_managed() noexcept
        {
            detail::destroy(alloc_, std::addressof(**this));
        }

        constexpr T& operator*() noexcept { return data_; }
        constexpr const T& operator*() const noexcept { return data_; }

    private:
        union { T data_; };
        Allocator& alloc_;
    };

    //--------------------------------------- SMALL VECTOR BUFFER -------------------------------------------------------

    inline constexpr std::size_t cache_line_size = 64;

    template<typename T, std::size_t Size>
    struct small_vector_buffer
    {
    public:
        constexpr small_vector_buffer() noexcept {};    // NOLINT(*default)
        constexpr ~small_vector_buffer() noexcept {};   // NOLINT(*default)

        constexpr auto begin() noexcept { return data_.data(); }
        constexpr auto begin() const noexcept { return data_.data(); }

        constexpr auto end() noexcept { return data_.data() + Size; }
        constexpr auto end() const noexcept { return data_.data() + Size; }

        constexpr std::size_t size() const noexcept { return Size; }
    private:
        union { std::array<T, Size> data_; };
    };


    template<typename T>
    struct default_small_size
    {
    private:
        inline constexpr static std::size_t overall_size = cache_line_size;
        inline constexpr static std::size_t buffer_size = overall_size - 3 * sizeof(T*);
        inline constexpr static std::size_t buffer_min_count = 4;
    public:
        inline constexpr static std::size_t value = std::max(buffer_min_count, buffer_size / sizeof(T));
    };

    template<typename T>
    inline constexpr std::size_t default_small_size_v = default_small_size<T>::value;

} // namespace gapp::detail

namespace gapp
{
    template<typename T, std::size_t Size = detail::default_small_size_v<T>, typename A = std::allocator<T>>
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

        static_assert(Size, "The size of the inline buffer must be at least 1.");

        //-----------------------------------//
        //            CONSTRUCTORS           //
        //-----------------------------------//

        constexpr small_vector() noexcept(std::is_nothrow_default_constructible_v<A>) :
            first_(buffer_.begin()),
            last_(buffer_.begin()),
            last_alloc_(buffer_.end())
        {}

        constexpr explicit small_vector(const A& allocator) noexcept(std::is_nothrow_copy_constructible_v<A>) :
            first_(buffer_.begin()),
            last_(buffer_.begin()),
            last_alloc_(buffer_.end()),
            alloc_(allocator)
        {}

        constexpr explicit small_vector(size_type count, const A& allocator = A()) :
            alloc_(allocator)
        {
            allocate_n(count);
            detail::scope_exit guard{ [&] { deallocate(); } };
            detail::construct_range(alloc_, first_, first_ + count);
            last_ = first_ + count;
            guard.release();
        }

        constexpr small_vector(size_type count, const T& value, const A& allocator = A()) :
            alloc_(allocator)
        {
            allocate_n(count);
            detail::scope_exit guard{ [&] { deallocate(); } };
            detail::construct_range(alloc_, first_, first_ + count, value);
            last_ = first_ + count;
            guard.release();
        }

        template<std::forward_iterator Iter>
        constexpr small_vector(Iter src_first, Iter src_last, const A& allocator = A()) :
            alloc_(allocator)
        {
            if (src_first == src_last) return;

            const auto src_len = std::distance(src_first, src_last);
            allocate_n(src_len);
            detail::scope_exit guard{ [&] { deallocate(); } };
            detail::construct_range(alloc_, first_, first_ + src_len, src_first);
            last_ = first_ + src_len;
            guard.release();
        }

        template<std::input_iterator Iter>
        constexpr small_vector(Iter src_first, Iter src_last, const A& allocator = A()) :
            first_(buffer_.begin()),
            last_(buffer_.begin()),
            last_alloc_(buffer_.end()),
            alloc_(allocator)
        {
            while (src_first != src_last) emplace_back(*src_first++);
        }

        constexpr small_vector(std::initializer_list<T> init, const A& allocator = A()) :
            small_vector(init.begin(), init.end(), allocator)
        {}

        constexpr small_vector(const small_vector& other) :
            small_vector(other.begin(), other.end(), std::allocator_traits<A>::select_on_container_copy_construction(other.alloc_))
        {}

        constexpr small_vector(const small_vector& other, const A& allocator) :
            small_vector(other.begin(), other.end(), allocator)
        {}

        constexpr small_vector(small_vector&& other) noexcept(std::is_nothrow_move_constructible_v<T> && detail::has_trivial_construct_v<A&, T, T&&>) :
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

        constexpr ~small_vector() noexcept
        {
            detail::destroy_range(alloc_, first_, last_);
            deallocate();
        }

        //-----------------------------------//
        //             ASSIGNMENT            //
        //-----------------------------------//

        constexpr void assign(size_type count, const T& value)
        {
            const auto src_size = difference_type(count);
            const auto old_size = std::distance(first_, last_);
            const auto com_size = std::min(old_size, src_size);

            if (last_alloc_ - first_ >= difference_type(count))
            {
                detail::assign_range(first_, first_ + com_size, value);
                detail::construct_range(alloc_, first_ + com_size, first_ + src_size, value);
                detail::destroy_range(alloc_, first_ + com_size, last_);
                last_ = first_ + src_size;
            }
            else
            {
                size_type new_cap = next_capacity(src_size - old_size);
                detail::alloc_result_t<A> alloc_result = detail::allocate<T>(alloc_, new_cap);
                detail::scope_exit guard{ [&] { detail::deallocate(alloc_, alloc_result.data, alloc_result.size); } };
                detail::construct_range(alloc_, alloc_result.data, alloc_result.data + src_size, value);
                detail::destroy_range(alloc_, first_, last_);
                guard.release();
                deallocate();
                set_storage(alloc_result.data, count, alloc_result.size);
            }
        }

        template<std::forward_iterator Iter>
        constexpr void assign(Iter src_first, Iter src_last)
        {
            const auto src_size = std::distance(src_first, src_last);
            const auto old_size = std::distance(first_, last_);
            const auto com_size = std::min(old_size, src_size);

            if (last_alloc_ - first_ >= src_size)
            {
                detail::assign_range(first_, first_ + com_size, src_first);
                detail::construct_range(alloc_, first_ + com_size, first_ + src_size, src_first + com_size);
                detail::destroy_range(alloc_, first_ + com_size, last_);
                last_ = first_ + src_size;
            }
            else
            {
                size_type new_cap = next_capacity(size_type(src_size - old_size));
                detail::alloc_result_t<A> alloc_result = detail::allocate<T>(alloc_, new_cap);
                detail::scope_exit guard{ [&] { detail::deallocate(alloc_, alloc_result.data, alloc_result.size); } };
                detail::construct_range(alloc_, alloc_result.data, alloc_result.data + src_size, src_first);
                detail::destroy_range(alloc_, first_, last_);
                guard.release();
                deallocate();
                set_storage(alloc_result.data, src_size, alloc_result.size);
            }
        }

        template<std::input_iterator Iter>
        constexpr void assign(Iter src_first, Iter src_last)
        {
            pointer next = first_;
            while (next != last_ && src_first != src_last) { *next++ = *src_first++; }

            detail::destroy_range(alloc_, next, last_);
            last_ = next;

            while (src_first != src_last) { emplace_back(*src_first++); }
        }

        constexpr void assign(std::initializer_list<T> list)
        {
            assign(list.begin(), list.end());
        }

        constexpr small_vector& operator=(const small_vector& other)
        {
            if (std::addressof(other) == this) [[unlikely]] return *this;

            if constexpr (detail::copy_allocators_v<A>)
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

        constexpr small_vector& operator=(small_vector&& other)
        noexcept(std::is_nothrow_move_constructible_v<T> && detail::has_trivial_construct_v<A&, T, T&&> &&
                 std::is_nothrow_move_assignable_v<T> &&
                 detail::steal_pointers_v<A>)
        {
            if (std::addressof(other) == this) [[unlikely]] return *this;

            if (!this->is_small() && !other.is_small())
            {
                using std::swap;

                if constexpr (detail::move_allocators_v<A>) swap(alloc_, other.alloc_);
                swap(first_, other.first_);
                swap(last_, other.last_);
                swap(last_alloc_, other.last_alloc_);

                return *this;
            }

            if (!other.is_small())
            {
                if constexpr (detail::steal_pointers_v<A>)
                {
                    detail::destroy_range(alloc_, first_, last_);
                    this->set_storage(other.first_, other.last_, other.last_alloc_);
                    other.set_buffer_storage(0);
                    alloc_ = std::move(other.alloc_);
                    return *this;
                }
                else if (alloc_ == other.alloc_)
                {
                    detail::destroy_range(alloc_, first_, last_);
                    this->set_storage(other.first_, other.last_, other.last_alloc_);
                    other.set_buffer_storage(0);
                    return *this;
                }
            }

            if constexpr (detail::move_allocators_v<A>)
            {
                if (alloc_ != other.alloc_)
                {
                    reset();
                    detail::relocate_range_weak(other.alloc_, other.first_, other.last_, first_);
                    last_ = first_ + other.size();
                    detail::destroy_range(other.alloc_, other.first_, other.last_);
                    other.set_buffer_storage(0);
                    alloc_ = std::move(other.alloc_);
                    return *this;
                }
            }

            assign(std::make_move_iterator(other.first_), std::make_move_iterator(other.last_));
            return *this;
        }

        constexpr small_vector& operator=(std::initializer_list<T> list)
        {
            assign(list.begin(), list.end());
            return *this;
        }

        //-----------------------------------//
        //             ITERATORS             //
        //-----------------------------------//

        constexpr iterator begin() noexcept { return first_; }
        constexpr const_iterator begin() const noexcept { return first_; }
        constexpr const_iterator cbegin() const noexcept { return first_; }

        constexpr iterator end() noexcept { return last_; }
        constexpr const_iterator end() const noexcept { return last_; }
        constexpr const_iterator cend() const noexcept { return last_; }

        constexpr reverse_iterator rbegin() noexcept { return std::make_reverse_iterator(last_); }
        constexpr const_reverse_iterator rbegin() const noexcept { return std::make_reverse_iterator(last_); }
        constexpr const_reverse_iterator crbegin() const noexcept { return std::make_reverse_iterator(last_); }

        constexpr reverse_iterator rend() noexcept { return std::make_reverse_iterator(first_); }
        constexpr const_reverse_iterator rend() const noexcept { return std::make_reverse_iterator(first_); }
        constexpr const_reverse_iterator crend() const noexcept { return std::make_reverse_iterator(first_); }

        //-----------------------------------//
        //           ELEMENT ACCESS          //
        //-----------------------------------//

        constexpr reference operator[](size_type pos) noexcept
        {
            GAPP_ASSERT(pos < size());
            return first_[pos];
        }

        constexpr const_reference operator[](size_type pos) const noexcept
        {
            GAPP_ASSERT(pos < size());
            return first_[pos];
        }

        constexpr reference at(size_type pos)
        {
            if (pos >= size()) GAPP_THROW(std::out_of_range, "Bad vector index.");
            return first_[pos];
        }

        constexpr const_reference at(size_type pos) const
        {
            if (pos >= size()) GAPP_THROW(std::out_of_range, "Bad vector index.");
            return first_[pos];
        }

        constexpr reference front() noexcept { GAPP_ASSERT(!empty()); return *first_; }
        constexpr const_reference front() const noexcept { GAPP_ASSERT(!empty()); return *first_; }

        constexpr reference back() noexcept { GAPP_ASSERT(!empty()); return *(last_ - 1); }
        constexpr const_reference back() const noexcept { GAPP_ASSERT(!empty()); return *(last_ - 1); }

        constexpr pointer data() noexcept { return first_; }
        constexpr const_pointer data() const noexcept { return first_; }

        //-----------------------------------//
        //              CAPACITY             //
        //-----------------------------------//

        constexpr bool empty() const noexcept { return first_ == last_; }

        constexpr size_type size() const noexcept { return size_type(last_ - first_); }
        constexpr difference_type ssize() const noexcept { return last_ - first_; }
        constexpr size_type capacity() const noexcept { return size_type(last_alloc_ - first_); }
        constexpr size_type max_size() const noexcept { return std::allocator_traits<A>::max_size(alloc_); }

        constexpr bool is_small() const noexcept { return first_ == buffer_.begin(); }
        static constexpr size_type inline_capacity() noexcept { return Size; }

        constexpr void reserve(size_type new_capacity) { if (new_capacity > capacity()) reallocate_n(next_capacity(new_capacity - capacity())); }
        constexpr void shrink_to_fit() { if (size() > inline_capacity()) reallocate_n(size()); }

        //-----------------------------------//
        //             MODIFIERS             //
        //-----------------------------------//

        constexpr void clear() noexcept
        {
            detail::destroy_range(alloc_, first_, last_);
            last_ = first_;
        }

        constexpr void reset() noexcept
        {
            detail::destroy_range(alloc_, first_, last_);
            deallocate();
            set_buffer_storage(0);
        }

        constexpr void swap(small_vector& other)
        noexcept(std::is_nothrow_swappable_v<T> && std::is_nothrow_move_constructible_v<T> && detail::has_trivial_construct_v<A&, T, T&&>)
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
                detail::relocate_range_strong(big.alloc_, big.first_ + small.size(), big.last_, small.last_);
                detail::destroy_range(big.alloc_, big.first_ + small.size(), big.last_);

                small.set_buffer_storage(big.size());
                big.set_buffer_storage(small_size);
            }
            else
            {
                small_vector& big   = this->is_small() ? other : *this;
                small_vector& small = this->is_small() ? *this : other;
                const auto small_size = small.size();

                detail::relocate_range_strong(small.alloc_, small.first_, small.last_, big.buffer_.begin());
                detail::destroy_range(small.alloc_, small.first_, small.last_);

                small.set_storage(big.first_, big.last_, big.last_alloc_);
                big.set_buffer_storage(small_size);
            }

            if constexpr (detail::swap_allocators_v<A>)
            {
                using std::swap;
                swap(alloc_, other.alloc_);
            }
        }

        constexpr void push_back(const T& value) { emplace_back(value); }
        constexpr void push_back(T&& value) { emplace_back(std::move(value)); }

        constexpr void push_back_unchecked(const T& value) noexcept(std::is_nothrow_copy_constructible_v<T>) { emplace_back_unchecked(value); }
        constexpr void push_back_unchecked(T&& value) noexcept(std::is_nothrow_move_constructible_v<T>) { emplace_back_unchecked(std::move(value)); }

        template<typename... Args>
        constexpr reference emplace_back(Args&&... args)
        {
            if (last_ != last_alloc_) return emplace_back_unchecked(std::forward<Args>(args)...);
            return *reallocate_append(next_capacity(), std::forward<Args>(args)...);
        }

        template<typename... Args>
        constexpr reference emplace_back_unchecked(Args&&... args) noexcept(std::is_nothrow_constructible_v<T, Args...>)
        {
            detail::construct(alloc_, last_, std::forward<Args>(args)...);
            return *last_++;
        }

        constexpr void pop_back() noexcept { GAPP_ASSERT(!empty()); detail::destroy(alloc_, last_--); }

        constexpr void resize(size_type count) { resize_impl(count); }
        constexpr void resize(size_type count, const T& value) { resize_impl(count, value); }

        constexpr iterator insert(const_iterator pos, const T& value) { return emplace(pos, value); }
        constexpr iterator insert(const_iterator pos, T&& value) { return emplace(pos, std::move(value)); }
        constexpr iterator insert(const_iterator pos, std::initializer_list<T> list) { return insert(pos, list.begin(), list.end()); }
        constexpr iterator erase(const_iterator pos) noexcept(std::is_nothrow_move_assignable_v<T>) { return erase(pos, pos + 1); }

        template<typename... Args>
        constexpr iterator emplace(const_iterator pos, Args&&... args)
        {
            if (last_ != last_alloc_)
            {
                if (pos == cend()) return std::addressof(emplace_back_unchecked(std::forward<Args>(args)...));

                detail::allocator_managed<T, A> new_elem(alloc_, std::forward<Args>(args)...);

                const difference_type offset = std::distance(cbegin(), pos);

                detail::construct(alloc_, last_, std::move(back()));
                std::shift_right(first_ + offset, last_++, 1);
                *(first_ + offset) = std::move(*new_elem);
                return first_ + offset;
            }

            return reallocate_emplace(next_capacity(), pos, std::forward<Args>(args)...);
        }

        constexpr iterator insert(const_iterator pos, size_type count, const T& value)
        {
            if (last_alloc_ - last_ >= count)
            {
                const difference_type offset = std::distance(cbegin(), pos);
                const difference_type src_size = difference_type(count);

                const auto middle = first_ + std::max(ssize() - src_size, offset);
                const auto moved_size = last_ - middle;
                const auto old_last   = last_;
                const auto new_last   = last_ + src_size;
                const auto new_middle = middle + src_size;

                detail::construct_range(alloc_, last_, new_middle, value);
                last_ = new_middle;
                detail::relocate_range_weak(alloc_, middle, old_last, new_middle);
                last_ = new_last;
                detail::assign_range(middle, old_last, value);

                return first_ + offset;
            }

            return reallocate_insert(next_capacity(count), pos, count, value);
        }

        template<std::forward_iterator Iter>
        constexpr iterator insert(const_iterator pos, Iter src_first, Iter src_last)
        {
            const difference_type src_size = std::distance(src_first, src_last);

            if (last_alloc_ - last_ >= src_size)
            {
                const difference_type offset = std::distance(cbegin(), pos);

                const auto middle = first_ + std::max(ssize() - src_size, offset);
                const auto moved_size = last_ - middle;
                const auto old_last   = last_;
                const auto new_last   = last_ + src_size;
                const auto new_middle = middle + src_size;

                detail::construct_range(alloc_, last_, new_middle, std::next(src_first, moved_size));
                last_ = new_middle;
                detail::relocate_range_weak(alloc_, middle, old_last, new_middle);
                last_ = new_last;
                detail::assign_range(middle, old_last, src_first);

                return first_ + offset;
            }

            return reallocate_insert(next_capacity(size_type(src_size)), pos, src_first, src_last);
        }

        template<std::input_iterator Iter>
        constexpr iterator insert(const_iterator pos, Iter src_first, Iter src_last)
        {
            const auto offset = std::distance(cbegin(), pos);
            const auto old_size = std::distance(first_, last_);

            while (src_first != src_last) emplace_back(*src_first++);
            std::rotate(first_ + offset, first_ + old_size, last_);

            return first_ + offset;
        }

        constexpr iterator erase(const_iterator first, const_iterator last) noexcept(std::is_nothrow_move_assignable_v<T>)
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

        constexpr allocator_type get_allocator() const noexcept(std::is_nothrow_copy_constructible_v<A>) { return alloc_; }

        constexpr friend bool operator==(const small_vector& lhs, const small_vector& rhs) noexcept
        {
            return std::equal(lhs.begin(), lhs.end(), rhs.begin(), rhs.end());
        }

        constexpr friend auto operator<=>(const small_vector& lhs, const small_vector& rhs) noexcept
        {
            return std::lexicographical_compare_three_way(lhs.begin(), lhs.end(), rhs.begin(), rhs.end());
        }

    private:
        static constexpr std::size_t alignment = std::max(alignof(T), detail::cache_line_size);

        alignas(alignment)
        GAPP_NO_UNIQUE_ADDRESS detail::small_vector_buffer<T, Size> buffer_;
        pointer first_      = nullptr;
        pointer last_       = nullptr;
        pointer last_alloc_ = nullptr;
        GAPP_NO_UNIQUE_ADDRESS allocator_type alloc_;


        constexpr void allocate_n(size_type count)
        {
            if (count <= inline_capacity()) return set_buffer_storage(0);

            detail::alloc_result_t<A> alloc_result = detail::allocate<T>(alloc_, count);
            first_      = alloc_result.data;
            last_       = alloc_result.data;
            last_alloc_ = alloc_result.data + alloc_result.size;
        }

        constexpr void reallocate_n(size_type new_capacity)
        {
            const size_type old_size = size();

            detail::alloc_result_t<A> alloc_result = detail::allocate<T>(alloc_, new_capacity);
            detail::scope_exit guard{ [&] { detail::deallocate(alloc_, alloc_result.data, alloc_result.size); } };
            detail::relocate_range_strong(alloc_, first_, last_, alloc_result.data);
            detail::destroy_range(alloc_, first_, last_);
            guard.release();
            deallocate();
            set_storage(alloc_result.data, old_size, alloc_result.size);
        }

        template<typename... Args>
        constexpr iterator reallocate_append(size_type new_capacity, Args&&... args)
        {
            const size_type old_size = size();

            detail::alloc_result_t<A> alloc_result = detail::allocate<T>(alloc_, new_capacity);
            detail::scope_exit guard1{ [&] { detail::deallocate(alloc_, alloc_result.data, alloc_result.size); } };
            detail::construct(alloc_, alloc_result.data + old_size, std::forward<Args>(args)...);
            detail::scope_exit guard2{ [&] { detail::destroy(alloc_, alloc_result.data + old_size); } };
            detail::relocate_range_strong(alloc_, first_, last_, alloc_result.data);
            { guard1.release(); guard2.release(); }
            detail::destroy_range(alloc_, first_, last_);
            deallocate();
            set_storage(alloc_result.data, old_size + 1, alloc_result.size);

            return alloc_result.data + old_size;
        }

        template<typename... Args>
        constexpr iterator reallocate_emplace(size_t new_capacity, const_iterator pos, Args&&... args)
        {
            const size_type old_size = size();
            const difference_type offset = std::distance(cbegin(), pos);

            detail::alloc_result_t<A> alloc_result = detail::allocate<T>(alloc_, new_capacity);
            detail::scope_exit guard1{ [&] { detail::deallocate(alloc_, alloc_result.data, alloc_result.size); } };
            detail::construct(alloc_, alloc_result.data + offset, std::forward<Args>(args)...);
            detail::scope_exit guard2{ [&] { detail::destroy(alloc_, alloc_result.data + offset); } };
            detail::relocate_range_strong(alloc_, first_, first_ + offset, alloc_result.data);
            detail::scope_exit guard3{ [&] { detail::destroy_range(alloc_, alloc_result.data, alloc_result.data + offset); } };
            detail::relocate_range_strong(alloc_, first_ + offset, last_, alloc_result.data + offset + 1);
            { guard1.release(); guard2.release(); guard3.release(); }
            detail::destroy_range(alloc_, first_, last_);
            deallocate();
            set_storage(alloc_result.data, old_size + 1, alloc_result.size);

            return alloc_result.data + offset;
        }

        constexpr iterator reallocate_insert(size_t new_capacity, const_iterator pos, size_type count, const T& value)
        {
            const size_type old_size = size();
            const difference_type src_size = difference_type(count);
            const difference_type offset = std::distance(cbegin(), pos);

            detail::alloc_result_t<A> alloc_result = detail::allocate<T>(alloc_, new_capacity);
            detail::scope_exit guard1{ [&] { detail::deallocate(alloc_, alloc_result.data, alloc_result.size); } };
            detail::relocate_range_weak(alloc_, first_, first_ + offset, alloc_result.data);
            detail::scope_exit guard2{ [&] { detail::destroy_range(alloc_, alloc_result.data, alloc_result.data + offset); } };
            detail::construct_range(alloc_, alloc_result.data + offset, alloc_result.data + offset + src_size, value);
            detail::scope_exit guard3{ [&] { detail::destroy_range(alloc_, alloc_result.data + offset, alloc_result.data + offset + src_size); } };
            detail::relocate_range_weak(alloc_, first_ + offset, last_, alloc_result.data + offset + src_size);
            { guard1.release(); guard2.release(); guard3.release(); }
            detail::destroy_range(alloc_, first_, last_);
            deallocate();
            set_storage(alloc_result.data, old_size + count, alloc_result.size);

            return alloc_result.data + offset;
        }

        template<std::forward_iterator Iter>
        constexpr iterator reallocate_insert(size_t new_capacity, const_iterator pos, Iter src_first, Iter src_last)
        {
            const size_type old_size = size();
            const difference_type src_size = std::distance(src_first, src_last);
            const difference_type offset = std::distance(cbegin(), pos);

            detail::alloc_result_t<A> alloc_result = detail::allocate<T>(alloc_, new_capacity);
            detail::scope_exit guard1{ [&] { detail::deallocate(alloc_, alloc_result.data, alloc_result.size); } };
            detail::relocate_range_weak(alloc_, first_, first_ + offset, alloc_result.data);
            detail::scope_exit guard2{ [&] { detail::destroy_range(alloc_, alloc_result.data, alloc_result.data + offset); } };
            detail::construct_range(alloc_, alloc_result.data + offset, alloc_result.data + offset + src_size, src_first);
            detail::scope_exit guard3{ [&] { detail::destroy_range(alloc_, alloc_result.data + offset, alloc_result.data + offset + src_size); } };
            detail::relocate_range_weak(alloc_, first_ + offset, last_, alloc_result.data + offset + src_size);
            { guard1.release(); guard2.release(); guard3.release(); }
            detail::destroy_range(alloc_, first_, last_);
            deallocate();
            set_storage(alloc_result.data, old_size + size_type(src_size), alloc_result.size);

            return alloc_result.data + offset;
        }

        constexpr void deallocate() noexcept
        {
            if (!is_small() && data()) detail::deallocate(alloc_, first_, capacity());
        }

        template<typename... Args>
        constexpr void resize_impl(size_type count, Args&&... args)
        {
            if (count <= size())
            {
                detail::destroy_range(alloc_, first_ + count, last_);
                last_ = first_ + count;
            }
            else
            {
                reserve(count);
                detail::construct_range(alloc_, last_, first_ + count, std::forward<Args>(args)...);
                last_ = first_ + count;
            }
        }

        constexpr void set_storage(pointer first, pointer last, pointer last_alloc) noexcept
        {
            first_ = first;
            last_ = last;
            last_alloc_ = last_alloc;
        }

        constexpr void set_storage(pointer first, size_type size, size_type capacity) noexcept
        {
            set_storage(first, first + size, first + capacity);
        }

        constexpr void set_buffer_storage(size_type size) noexcept
        {
            set_storage(buffer_.begin(), buffer_.begin() + size, buffer_.end());
        }

        constexpr size_type next_capacity(size_type min_growth = 1) const
        {
            const size_type current_capacity = capacity();
            const size_type growth_capacity  = current_capacity >> 1;
            const size_type max_capacity     = max_size();

            if (min_growth > max_capacity - current_capacity)
            {
                GAPP_THROW(std::length_error, "Too big vector.");
            }

            if (current_capacity > max_capacity - growth_capacity)
            {
                return max_capacity;
            }

            return std::max(current_capacity + min_growth, current_capacity + growth_capacity);
        }

    }; // class small_vector

    template<std::input_iterator Iter, std::size_t Size = detail::default_small_size_v<std::iter_value_t<Iter>>, typename Alloc = std::allocator<std::iter_value_t<Iter>>>
    small_vector(Iter, Iter, Alloc = Alloc()) -> small_vector<std::iter_value_t<Iter>, Size, Alloc>;

    template<typename T, std::size_t Size, typename A>
    constexpr void swap(small_vector<T, Size, A>& lhs, small_vector<T, Size, A>& rhs)
    noexcept(noexcept(lhs.swap(rhs)))
    {
        lhs.swap(rhs);
    }

} // namespace gapp

// NOLINTEND(*pointer-arithmetic, *union-access)

#endif // !GAPP_UTILITY_SMALL_VECTOR_HPP
