/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_UTILITY_UTILITY_HPP
#define GAPP_UTILITY_UTILITY_HPP

#include <vector>
#include <concepts>
#include <type_traits>
#include <cassert>
#include <cstdint>
#include <cstddef>


#if defined(__GNUC__) && !defined(__clang__)
#   define GAPP_GCC_COMPILER
#elif defined(__GNUC__) && defined(__clang__)
#   define GAPP_CLANG_COMPILER
#elif defined(_MSC_VER) && !defined(__clang__)
#   define GAPP_MSVC_COMPILER
#elif defined(_MSC_VER) && defined(__clang__)
#   define GAPP_CLANG_CL_COMPILER
#endif


#if defined(_M_IX86) || defined(_M_X64) || defined(__x86_64__) || defined(__i386__)
#   define GAPP_X86_ARCH
#endif

#if defined(_M_ARM) || defined(_M_ARM64) || defined(__arm__) || defined(__aarch64__)
#   define GAPP_ARM_ARCH
#endif


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
#   if defined(GAPP_MSVC_COMPILER)
#       define GAPP_NO_UNIQUE_ADDRESS [[msvc::no_unique_address]]
#   else
#       define GAPP_NO_UNIQUE_ADDRESS [[no_unique_address]]
#   endif
#else
#   define GAPP_NO_UNIQUE_ADDRESS
#endif


#if defined(__GNUC__) || defined(__clang__)
#   define GAPP_NOINLINE __attribute((noinline))
#elif defined(_MSC_VER)
#   define GAPP_NOINLINE __declspec(noinline)
#else
#   define GAPP_NOINLINE
#endif


#if (defined(__GNUC__) || defined(__clang__)) && defined(GAPP_X86_ARCH)
#   define GAPP_PAUSE() __builtin_ia32_pause()
#elif defined(GAPP_MSVC_COMPILER) && defined(GAPP_X86_ARCH)
#   define GAPP_PAUSE() _mm_pause()
#elif defined(GAPP_MSVC_COMPILER) && defined(GAPP_ARM_ARCH)
#   define GAPP_PAUSE() __yield()
#else
#   define GAPP_PAUSE()
#endif


#if defined(__has_feature)
#   if __has_feature(thread_sanitizer) && __has_include(<sanitizer/tsan_interface.h>)
#       include <sanitizer/tsan_interface.h>
#       define GAPP_ANNOTATE_TSAN_ACQUIRE(p) __tsan_acquire((void*)p)
#       define GAPP_ANNOTATE_TSAN_RELEASE(p) __tsan_release((void*)p)
#   else
#       define GAPP_ANNOTATE_TSAN_ACQUIRE(p)
#       define GAPP_ANNOTATE_TSAN_RELEASE(p)
#   endif
#elif defined(__SANITIZE_THREAD__) && __has_include(<sanitizer/tsan_interface.h>)
#   include <sanitizer/tsan_interface.h>
#   define GAPP_ANNOTATE_TSAN_ACQUIRE(p) __tsan_acquire((void*)p)
#   define GAPP_ANNOTATE_TSAN_RELEASE(p) __tsan_release((void*)p)
#else
#   define GAPP_ANNOTATE_TSAN_ACQUIRE(p)
#   define GAPP_ANNOTATE_TSAN_RELEASE(p)
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


#if !defined(NDEBUG) && !defined(GAPP_DISABLE_ASSERTS)
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

    // returns the length of the range [low, high] without overflow
    template<std::integral T>
    constexpr std::uint64_t range_length(T low, T high) noexcept
    {
        GAPP_ASSERT(low <= high);

        return std::uint64_t(high) - std::uint64_t(low);
    }

    template<std::integral T>
    constexpr T next_mod(T value, T mod) noexcept
    {
        GAPP_ASSERT(mod > 0);
        GAPP_ASSERT(0 <= value && value < mod);

        return (value + 1 == mod) ? T(0) : (value + 1);
    }

    template<std::integral T>
    constexpr T prev_mod(T value, T mod) noexcept
    {
        GAPP_ASSERT(mod > 0);
        GAPP_ASSERT(0 <= value && value < mod);

        return (value == 0) ? (mod - 1) : (value - 1);
    }

    template<std::integral T>
    constexpr void increment_mod(T& value, T mod) noexcept
    {
        value = next_mod(value, mod);
    }

    template<std::integral T>
    constexpr void decrement_mod(T& value, T mod) noexcept
    {
        value = prev_mod(value, mod);
    }

} // namespace gapp::detail

#endif // !GAPP_UTILITY_UTILITY_HPP
