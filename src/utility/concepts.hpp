/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_CONCEPTS_HPP
#define GA_UTILITY_CONCEPTS_HPP

#include <concepts>
#include <functional>

namespace gapp::detail
{
    template<typename T>
    concept hashable = requires(T arg)
    {
        { std::hash<T>{}(arg) } -> std::convertible_to<size_t>;
    };

} // namespace gapp::detail

#endif // !GA_UTILITY_CONCEPTS_HPP