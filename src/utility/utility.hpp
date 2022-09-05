/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_UTILITY_HPP
#define GA_UTILITY_UTILITY_HPP

#include <limits>
#include <execution>
#include <cassert>
#include <cstddef>


#ifndef GA_SEED
#define GA_SEED 0x3da99432ab975d26LL
#endif


#define GA_ASSERT(condition, msg) assert((condition) && msg)
#define GA_THROW(exception_type, msg) throw exception_type(msg)

#define GA_UNUSED(...) (void)(sizeof(__VA_ARGS__))


#ifndef GA_EXCUTION_UNSEQ
#define GA_EXECUTION_UNSEQ std::execution::par_unseq
#endif
#ifndef GA_EXECUTION_SEQ
#define GA_EXECUTION_SEQ std::execution::par
#endif


#ifndef GA_EPSILON
#define GA_EPSILON 4
#endif

namespace genetic_algorithm::math
{
    template<typename T>
    inline constexpr T epsilon = GA_EPSILON * std::numeric_limits<T>::epsilon();

    template<typename T>
    inline constexpr T inf = std::numeric_limits<T>::infinity();

} // namespace genetic_algorithm::math

namespace genetic_algorithm
{
    constexpr std::size_t operator ""_sz(unsigned long long arg)
    {
        return static_cast<std::size_t>(arg);
    }

    constexpr std::ptrdiff_t operator ""_pd(unsigned long long arg)
    {
        return static_cast<std::ptrdiff_t>(arg);
    }

} // namespace genetic_algorithm

#endif // !GA_UTILITY_UTILITY_HPP