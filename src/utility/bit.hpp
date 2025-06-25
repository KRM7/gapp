/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_BIT_HPP
#define GA_UTILITY_BIT_HPP

#include "utility.hpp"
#include <bit>
#include <type_traits>
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

    template<std::floating_point T>
    inline constexpr size_t implicit_mantissa_bits = std::numeric_limits<T>::digits;


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

        return n ? ones<std::make_unsigned_t<T>> >> (bitsizeof<T> - n) : zeros<std::make_unsigned_t<T>>;
    }

    template<std::integral T>
    constexpr T mask_left_n(std::size_t n) noexcept
    {
        GAPP_ASSERT(n <= bitsizeof<T>);

        return n ? ones<std::make_unsigned_t<T>> << (bitsizeof<T> - n) : zeros<std::make_unsigned_t<T>>;
    }

    template<std::size_t first, std::size_t last, std::integral T>
    constexpr T extract_bits(T value) noexcept
    {
        static_assert(first < last);
        static_assert(first < detail::bitsizeof<T>);
        static_assert(last < detail::bitsizeof<T> + 1);

        return (std::make_unsigned_t<T>(value) >> first) & detail::mask_right_n<std::make_unsigned_t<T>>(last - first);
    }

    template<std::integral T>
    constexpr T block_of(bool value) noexcept
    {
        return static_cast<T>(-1 * static_cast<int>(value)); // NOLINT(*widening-cast)
    }

    template<std::same_as<float> T>
    constexpr T set_sign_bit(T value, bool sign_bit) noexcept
    {
        return std::bit_cast<T>((std::bit_cast<uint32_t>(value) & ~(uint32_t(1) << 31u)) | (uint32_t(sign_bit) << 31u));
    }

    template<std::same_as<double> T>
    constexpr T set_sign_bit(T value, bool sign_bit) noexcept
    {
        return std::bit_cast<T>((std::bit_cast<uint64_t>(value) & ~(uint64_t(1) << 63u)) | (uint64_t(sign_bit) << 63u));
    }

} // namespace gapp::detail

#endif // !GA_UTILITY_BIT_HPP
