/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_UTILITY_HPP
#define GA_UTILITY_UTILITY_HPP

#include <execution>
#include <vector>
#include <concepts>
#include <type_traits>
#include <cassert>
#include <cstddef>


#define GAPP_GCC_COMPILER ( __GNUC__ && !__clang__ )
#define GAPP_CLANG_COMPILER ( __GNUC__ && __clang__ )
#define GAPP_MSVC_COMPILER ( _MSC_VER && !__clang__ )
#define GAPP_CLANG_CL_COMPILER ( _MSC_VER && __clang__ )


#if defined(_MSC_VER)
#   define GAPP_ASSUME(expr) __assume(expr)
#elif defined(__clang__)
#   define GAPP_ASSUME(expr) __builtin_assume(expr)
#else
#   define GAPP_ASSUME(expr)
#endif


#if defined(__GNUC__)
#   define GAPP_UNREACHABLE() __builtin_unreachable()
#elif defined(_MSC_VER)
#   define GAPP_UNREACHABLE() __assume(false)
#else
#   define GAPP_UNREACHABLE() assert(false);
#endif


#if __has_cpp_attribute(no_unique_address)
#   if GAPP_MSVC_COMPILER
#       define GAPP_NO_UNIQUE_ADDRESS [[msvc::no_unique_address]]
#   else
#       define GAPP_NO_UNIQUE_ADDRESS [[no_unique_address]]
#   endif
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
#   define GAPP_ASSERT_1(condition) assert(condition)
#   define GAPP_ASSERT_2(condition, msg) assert( (condition) && (msg) )
#   define GAPP_ASSERT_(_1, _2, NAME, ...) NAME
#   define GAPP_ASSERT(...) GAPP_ASSERT_(__VA_ARGS__, GAPP_ASSERT_2, GAPP_ASSERT_1, _0)(__VA_ARGS__)
#else
#   define GAPP_ASSERT(...)
#endif


#if ((__GNUC__ && __cpp_exceptions) || (_MSC_VER && _HAS_EXCEPTIONS))
#   define GAPP_THROW(exception_type, msg) throw exception_type(msg)
#   define GAPP_TRY try
#   define GAPP_CATCH(exception) catch (exception)
#else
#   define GAPP_THROW(exception_type, msg) GAPP_UNREACHABLE()
#   define GAPP_TRY if constexpr (true)
#   define GAPP_CATCH(exception) if constexpr (false)
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

    // returns the length of the range [low, high) without overflow
    template<std::integral T>
    constexpr size_t range_length(T low, T high) noexcept
    {
        GAPP_ASSERT(low <= high);

        if constexpr (std::is_unsigned_v<T>)
        {
            return high - low;
        }
        else
        {
            return same_sign(low, high) ? high - low : -(low + 1) + high + 1;
        }
    }

} // namespace gapp::detail

#endif // !GA_UTILITY_UTILITY_HPP