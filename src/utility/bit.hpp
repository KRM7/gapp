/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_BIT_HPP
#define GA_UTILITY_BIT_HPP

#include "utility.hpp"
#include <concepts>
#include <limits>
#include <climits>
#include <cstddef>

namespace gapp::detail
{
    template<typename T>
    inline constexpr size_t bitsizeof = CHAR_BIT * sizeof(T);

    template<std::integral T>
    inline constexpr T lsb_mask = T{ 1 };

    template<std::integral T>
    inline constexpr T msb_mask = T{ 1 } << (bitsizeof<T> - 1);


    template<std::integral T>
    inline constexpr T ones = ~T{ 0 };

    template<std::integral T>
    inline constexpr T zeros = T{ 0 };


    template<std::floating_point T>
    inline constexpr size_t mantissa_bits = std::numeric_limits<T>::digits - 1;

    template<std::floating_point T>
    inline constexpr size_t exponent_bits = bitsizeof<T> - mantissa_bits<T> - 1;


    template<std::integral T>
    constexpr bool is_nth_bit_set(T value, size_t n) noexcept
    {
        GAPP_ASSERT(n < bitsizeof<T>);

        return value & (T{ 1 } << n);
    }

    template<std::integral T>
    constexpr bool msb(T value) noexcept
    {
        return value & detail::msb_mask<T>;
    }

    template<std::integral T>
    constexpr bool lsb(T value) noexcept
    {
        return value & detail::lsb_mask<T>;
    }

    template<std::integral T>
    constexpr T mask_right_n(std::size_t n) noexcept
    {
        GAPP_ASSERT(n <= bitsizeof<T>);

        return n ? ones<T> >> (bitsizeof<T> - n) : zeros<T>;
    }

    template<std::integral T>
    constexpr T mask_left_n(std::size_t n) noexcept
    {
        GAPP_ASSERT(n <= bitsizeof<T>);

        return n ? ones<T> << (bitsizeof<T> - n) : zeros<T>;
    }

    template<std::integral T>
    constexpr T block_of(bool value) noexcept
    {
        return static_cast<T>(-1 * static_cast<int>(value)); // NOLINT(*widening-cast)
    }

} // namespace gapp::detail

#endif // !GA_UTILITY_BIT_HPP
