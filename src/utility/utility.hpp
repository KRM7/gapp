/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILS_HPP
#define GA_UTILS_HPP

#include <limits>
#include <execution>
#include <cstddef>

#ifndef GA_EPSILON
#define GA_EPSILON (4 * std::numeric_limits<double>::epsilon())
#endif

#ifndef GA_EXCUTION_UNSEQ
#define GA_EXECUTION_UNSEQ std::execution::par_unseq
#endif
#ifndef GA_EXECUTION_SEQ
#define GA_EXECUTION_SEQ std::execution::par
#endif

#ifndef GA_SEED
#define GA_SEED 0x3da99432ab975d26LL
#endif

namespace genetic_algorithm
{
    constexpr std::size_t operator ""_sz(unsigned long long n)
    {
        return size_t{ n };
    }

} // namespace genetic_algorithm

#endif // !GA_UTILS_HPP