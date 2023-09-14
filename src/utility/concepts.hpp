﻿/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_CONCEPTS_HPP
#define GA_UTILITY_CONCEPTS_HPP

#include "type_traits.hpp"
#include <concepts>
#include <functional>
#include <type_traits>
#include <limits>

namespace gapp::detail
{
    template<typename T>
    concept Hashable = requires(T arg)
    {
        { std::hash<T>{}(arg) } -> std::convertible_to<size_t>;
    };

    template<typename C>
    concept Container = requires(C c, const C cc, C&& rc)
    {
        typename C::value_type;
        typename C::reference;
        typename C::const_reference;
        typename C::iterator;
        typename C::const_iterator;
        typename C::difference_type;
        typename C::size_type;

        requires std::forward_iterator<typename C::iterator>;
        requires std::forward_iterator<typename C::const_iterator>;

        requires std::unsigned_integral<typename C::size_type>;
        requires std::signed_integral<typename C::difference_type>;
        requires std::same_as<typename C::difference_type, typename std::iterator_traits<typename C::iterator>::difference_type>;
        requires std::same_as<typename C::difference_type, typename std::iterator_traits<typename C::const_iterator>::difference_type>;
        requires std::numeric_limits<typename C::size_type>::max() >= std::numeric_limits<typename C::difference_type>::max();

        { c.begin() }   -> std::same_as<typename C::iterator>;
        { cc.begin() }  -> std::same_as<typename C::const_iterator>;
        { c.end() }     -> std::same_as<typename C::iterator>;
        { cc.end() }    -> std::same_as<typename C::const_iterator>;
        { c.cbegin() }  -> std::same_as<typename C::const_iterator>;
        { c.cbegin() }  -> std::same_as<typename C::const_iterator>;
        { cc.cbegin() } -> std::same_as<typename C::const_iterator>;
        { cc.cbegin() } -> std::same_as<typename C::const_iterator>;

        { c.size() }  -> std::same_as<typename C::size_type>;
        { cc.size() } -> std::same_as<typename C::size_type>;

        { c.empty() }  -> std::convertible_to<bool>;
        { cc.empty() } -> std::convertible_to<bool>;
    };

    template<typename C>
    concept IndexableContainer = requires(C c, const C cc, size_t idx)
    {
        requires Container<C>;

        { c.operator[](idx)  } -> std::same_as<typename C::reference>;
        { cc.operator[](idx) } -> std::same_as<typename C::const_reference>;
    };

} // namespace gapp::detail

#endif // !GA_UTILITY_CONCEPTS_HPP