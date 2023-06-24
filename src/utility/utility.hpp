/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_UTILITY_HPP
#define GA_UTILITY_UTILITY_HPP

#include <execution>
#include <vector>
#include <concepts>
#include <type_traits>
#include <cassert>
#include <cstddef>


#define GA_ASSERT(condition, ...) assert( (condition) __VA_OPT__(&& (__VA_ARGS__)) )
#define GA_THROW(exception_type, msg) throw exception_type(msg)


#ifdef __GNUC__
#   define GA_UNREACHABLE() __builtin_unreachable()
#elif defined(_MSC_VER)
#   define GA_UNREACHABLE() __assume(false)
#else
#   define GA_UNREACHABLE() assert(false);
#endif


#if defined(__GNUC__) || defined(__clang__)
#   define GA_CONST __attribute__((const))
#elif defined(_MSC_VER)
#   define GA_CONST __declspec(noalias)
#else
#   define GA_CONST
#endif


#if defined(__GNUC__) || defined(__clang__)
#   define GA_PURE __attribute__((pure))
#else
#   define GA_PURE
#endif


#if defined(__GNUC__) || defined(__clang__)
#   define GA_NON_NULL __attribute__((nonnull))
#else
#   define GA_NON_NULL
#endif


#if defined(__GNUC__) || defined(__clang__)
#   define GA_NOINLINE __attribute((noinline))
#elif defined(_MSC_VER)
#   define GA_NOINLINE __declspec(noinline)
#else
#   define GA_NOINLINE
#endif


#if defined(__GNUC__) || defined(__clang__)
#   define GA_PAUSE() __builtin_ia32_pause()
#elif defined(_MSC_VER)
#   define GA_PAUSE() _mm_pause()
#else
#   define GA_PAUSE()
#endif


#if defined(_WIN32) && !defined(GAPP_BUILD_STATIC)
#   if defined(gapp_EXPORTS)
#       define GA_API __declspec(dllexport)
#   else
#       define GA_API __declspec(dllimport)
#   endif
#else
#   define GA_API
#endif


#ifndef GA_EXCUTION_UNSEQ
#   define GA_EXECUTION_UNSEQ std::execution::par_unseq
#endif

#ifndef GA_EXECUTION_SEQ
#   define GA_EXECUTION_SEQ std::execution::par
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
        GA_ASSERT(low <= high);

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