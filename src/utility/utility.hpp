/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_UTILITY_HPP
#define GA_UTILITY_UTILITY_HPP

#include <execution>
#include <cassert>
#include <cstddef>


#define GA_ASSERT(condition, msg) assert((condition) && (msg))
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


#if defined(_WIN32) && !defined(GA_BUILD_STATIC)
#   if defined(GeneticAlgorithm_EXPORTS)
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


namespace genetic_algorithm
{
    constexpr std::size_t operator ""_sz(unsigned long long arg) noexcept
    {
        return static_cast<std::size_t>(arg);
    }

    constexpr std::ptrdiff_t operator ""_pd(unsigned long long arg) noexcept
    {
        return static_cast<std::ptrdiff_t>(arg);
    }

} // namespace genetic_algorithm

#endif // !GA_UTILITY_UTILITY_HPP