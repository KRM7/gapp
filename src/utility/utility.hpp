/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_UTILITY_HPP
#define GA_UTILITY_UTILITY_HPP

#include <execution>
#include <vector>
#include <concepts>
#include <type_traits>
#include <cassert>
#include <cstddef>


#if __has_cpp_attribute(assume)
#   define GAPP_ASSUME(expr) [[assume(expr)]]
#elif defined(_MSC_VER) && !defined(__clang__)
#   define GAPP_ASSUME(expr) //__assume(expr)
#elif defined(__clang__)
#   define GAPP_ASSUME(expr) //__builtin_assume(expr)
#else
#   define GAPP_ASSUME(expr)
#endif


#ifdef __GNUC__
#   define GAPP_UNREACHABLE() __builtin_unreachable()
#elif defined(_MSC_VER)
#   define GAPP_UNREACHABLE() __assume(false)
#else
#   define GAPP_UNREACHABLE() assert(false);
#endif


#if defined(_MSC_VER) && !defined(__clang__)
#   define GAPP_NO_UNIQUE_ADDRESS [[msvc::no_unique_address]]
#elif __has_cpp_attribute(no_unique_address)
#   define GAPP_NO_UNIQUE_ADDRESS [[no_unique_address]]
#else
#   define GAPP_NO_UNIQUE_ADDRESS
#endif


#if defined(__GNUC__) || defined(__clang__)
#   define GAPP_CONST __attribute__((const))
#elif defined(_MSC_VER)
#   define GAPP_CONST __declspec(noalias)
#else
#   define GAPP_CONST
#endif


#if defined(__GNUC__) || defined(__clang__)
#   define GAPP_PURE __attribute__((pure))
#else
#   define GAPP_PURE
#endif


#if defined(__GNUC__) || defined(__clang__)
#   define GAPP_NON_NULL __attribute__((nonnull))
#else
#   define GAPP_NON_NULL
#endif


#if defined(__GNUC__) || defined(__clang__)
#   define GAPP_NOINLINE __attribute((noinline))
#elif defined(_MSC_VER)
#   define GAPP_NOINLINE __declspec(noinline)
#else
#   define GAPP_NOINLINE
#endif


#if defined(__GNUC__) || defined(__clang__)
#   define GAPP_PAUSE() __builtin_ia32_pause()
#elif defined(_MSC_VER)
#   define GAPP_PAUSE() _mm_pause()
#else
#   define GAPP_PAUSE()
#endif


#if defined(_WIN32) && !defined(GAPP_BUILD_STATIC)
#   if defined(gapp_EXPORTS)
#       define GAPP_API __declspec(dllexport)
#   else
#       define GAPP_API __declspec(dllimport)
#   endif
#else
#   define GAPP_API
#endif


#ifndef NDEBUG
#   define GAPP_ASSERT(condition, ...) assert( (condition) __VA_OPT__(&& (__VA_ARGS__)) )
#else
#   define GAPP_ASSERT(condition, ...) GAPP_ASSUME(condition)
#endif


#ifndef GAPP_NO_EXCEPTIONS
#   define GAPP_THROW(exception_type, msg) throw exception_type(msg)
#else
#   define GAPP_THROW(exception_type, msg) GAPP_UNREACHABLE()
#endif


#ifndef GAPP_EXCUTION_UNSEQ
#   define GAPP_EXEC_UNSEQ std::execution::par_unseq
#endif

#ifndef GAPP_EXEC_SEQ
#   define GAPP_EXEC_SEQ std::execution::par
#endif


namespace gapp
{
    constexpr std::size_t operator ""_sz(unsigned long long arg) noexcept
    {
        return static_cast<std::size_t>(arg);
    }

    constexpr std::ptrdiff_t operator ""_pd(unsigned long long arg) noexcept
    {
        return static_cast<std::ptrdiff_t>(arg);
    }

} // namespace gapp


namespace gapp::detail
{
    template<typename T>
    inline void clear_reserve(std::vector<T>& vec, size_t new_capacity = 0)
    {
        std::vector<T> temp;
        temp.reserve(new_capacity);
        temp.swap(vec);
    }

    // returns true if the signs of the parameters are the same
    template<std::signed_integral T>
    constexpr bool same_sign(T left, T right) noexcept
    {
        return (left ^ right) >= 0;
    }

    // returns the length of the range [low, high)
    template<std::integral T>
    constexpr size_t range_length(T low, T high) noexcept
    {
        GAPP_ASSERT(low <= high);

        if constexpr (std::is_unsigned_v<T>)
        {
            return high - low;
        }
        else if (same_sign(low, high))
        {
            return high - low;
        }
        else
        {
            return static_cast<size_t>(-(low + 1)) + high + 1;
        }
    }

} // namespace gapp::detail

#endif // !GA_UTILITY_UTILITY_HPP