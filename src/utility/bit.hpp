/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_BIT_HPP
#define GA_UTILITY_BIT_HPP

#include "utility.hpp"
#include <concepts>
#include <limits>
#include <climits>

namespace gapp::detail
{
    template<typename T>
    inline constexpr size_t bitsizeof = CHAR_BIT * sizeof(T);

    template<typename T>
    inline constexpr T lsb_mask = T{ 1 };

    template<typename T>
    inline constexpr T msb_mask = T{ 1 } << (bitsizeof<T> - 1);


    template<std::floating_point T>
    inline constexpr size_t mantissa_bits = std::numeric_limits<T>::digits - 1;

    template<std::floating_point T>
    inline constexpr size_t exponent_bits = bitsizeof<T> - mantissa_bits<T> - 1;


    template<typename T>
    constexpr bool is_nth_bit_set(T value, size_t n) noexcept
    {
        GAPP_ASSERT(n < bitsizeof<T>);

        return value & (T{ 1 } << n);
    }

    template<typename T>
    constexpr bool first_bit(T value) noexcept
    {
        return value & detail::msb_mask<T>;
    }

    template<typename T>
    constexpr bool last_bit(T value) noexcept
    {
        return value & detail::lsb_mask<T>;
    }

} // namespace gapp::detail

#endif // !GA_UTILITY_BIT_HPP
