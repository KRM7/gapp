/* Copyright (c) 2025 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_UTILITY_HASH_HPP
#define GAPP_UTILITY_HASH_HPP

#include "utility.hpp"
#include <concepts>
#include <functional>
#include <iterator>
#include <ranges>
#include <cstddef>
#include <cstdint>

namespace gapp::detail
{
    template<typename T>
    concept hashable = requires(T arg)
    {
        { std::hash<T>{}(arg) } -> std::convertible_to<std::size_t>;
    };

    template<detail::hashable T>
    constexpr std::uint64_t hash(const T& value) noexcept
    {
        return std::hash<T>{}(value);
    }

    constexpr std::uint64_t hash_combine(std::uint64_t value) noexcept
    {
        return value;
    }

    constexpr std::uint64_t hash_combine(std::uint64_t first, std::uint64_t second) noexcept
    {
        std::uint64_t hash = first + second + 0x9e3779b97f4a7c15;
        hash = 0xbf58476d1ce4e5b9 * (hash ^ (hash >> 30));
        hash = 0x94d049bb133111eb * (hash ^ (hash >> 27));
        return hash ^ (hash >> 31);
    }

    template<std::convertible_to<std::uint64_t>... Ts>
    constexpr std::uint64_t hash_combine(std::uint64_t first, std::uint64_t second, Ts... rest) noexcept
    {
        return detail::hash_combine(detail::hash_combine(first, second), rest...);
    }

    template<std::forward_iterator Iter>
    constexpr std::uint64_t hash_range(std::uint64_t seed, Iter first, Iter last) noexcept
    {
        GAPP_ASSERT(first <= last);

        for (; first != last; ++first) seed = hash_combine(seed, detail::hash(*first));
        return seed;
    }

    template<std::forward_iterator Iter>
    constexpr std::uint64_t hash_range(Iter first, Iter last) noexcept
    {
        return detail::hash_range(std::distance(first, last), first, last);
    }

    template<std::ranges::forward_range R>
    constexpr std::uint64_t hash_range(R&& range) noexcept
    {
        return detail::hash_range(range.begin(), range.end());
    }

} // namespace gapp::detail

#endif // !GAPP_UTILITY_HASH_HPP
