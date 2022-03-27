/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CONCEPTS_HPP
#define GA_CONCEPTS_HPP

#include <concepts>
#include <functional>
#include <type_traits>
#include <limits>
#include <iterator>
#include <compare>

namespace genetic_algorithm
{
    template<typename geneType>
    class GA;

    namespace detail
    {
        /* Types that std::hash is specialized for. */
        template<typename T>
        concept Hashable = requires(T arg)
        {
            { std::hash<T>{}(arg) } -> std::unsigned_integral;
        };

        template<typename T, template<typename...> class Templ>
        struct SpecializationOfImpl : std::false_type {};

        template<template<typename...> class Templ, typename... TArgs>
        struct SpecializationOfImpl<Templ<TArgs...>, Templ> : std::true_type {};

        /* Types that are specializations of a class template. */
        template<typename S, template<typename...> class Templ>
        concept SpecializationOf = SpecializationOfImpl<S, Templ>::value;

        template<template<typename...> class BaseTempl, typename... TArgs>
        void derivedFromSpecializationOfImpl(const BaseTempl<TArgs...>&);

        /* Types that are derived from a specialization of a class template. */
        template<typename Derived, template<typename...> class BaseTempl>
        concept DerivedFromSpecializationOf = requires(const Derived& arg)
        {
            derivedFromSpecializationOfImpl<BaseTempl>(arg);
        };

        /* based on named requirements */
        template<typename C>
        concept Container = requires(C c, const C cc, C && rc)
        {
            requires C::value_type;
            requires C::reference;
            requires C::const_reference;
            requires C::iterator;
            requires C::const_iterator;
            requires C::difference_type;
            requires C::size_type;

            requires std::same_as<typename C::reference, typename C::value_type&>;
            requires std::same_as<typename C::const_reference, const typename C::value_type&>;

            requires std::forward_iterator<typename C::iterator>;
            requires std::same_as<typename std::iterator_traits<typename C::iterator>::value_type, typename C::value_type>;
            requires std::forward_iterator<typename C::const_iterator>;
            requires std::same_as<typename std::iterator_traits<typename C::const_iterator>::value_type, const typename C::value_type>;

            requires std::signed_integral<typename C::difference_type>;
            requires std::same_as<typename C::difference_type, typename std::iterator_traits<typename C::iterator>::difference_type>;
            requires std::same_as<typename C::difference_type, typename std::iterator_traits<typename C::const_iterator>::difference_type>;
            requires std::unsigned_integral<typename C::size_type>;
            requires std::numeric_limits<typename C::size_type>::max() >= std::numeric_limits<typename C::difference_type>::max();

            requires std::destructible<typename C::value_type>;
            requires std::regular<C>;
            requires std::destructible<C>;

            { c.operator=(c) } -> std::same_as<C&>;
            { c.operator=(cc) } -> std::same_as<C&>;
            { c.operator=(rc) } -> std::same_as<C&>;

            { c.begin() } -> std::same_as<typename C::iterator>;
            { cc.begin() } -> std::same_as<typename C::const_iterator>;
            { c.end() } -> std::same_as<typename C::iterator>;
            { cc.end() } -> std::same_as<typename C::const_iterator>;
            { c.cbegin() } -> std::same_as<typename C::const_iterator>;
            { c.cbegin() } -> std::same_as<typename C::const_iterator>;
            { cc.cbegin() } -> std::same_as<typename C::const_iterator>;
            { cc.cbegin() } -> std::same_as<typename C::const_iterator>;

            { c.size() } -> std::same_as<typename C::size_type>;
            { c.max_size() } -> std::same_as<typename C::size_type>;
            { cc.size() } -> std::same_as<typename C::size_type>;
            { cc.max_size() } -> std::same_as<typename C::size_type>;

            { c.empty() } -> std::convertible_to<bool>;
            { cc.empty() } -> std::convertible_to<bool>;
        };

        template<typename C>
        concept IndexableContainer = requires(C c, const C cc, size_t idx)
        {
            requires Container<C>;
            { c.operator[](idx) } -> std::same_as<typename C::reference>;
            { cc.operator[](idx) } -> std::same_as<typename C::const_reference>;
        };

    } // namespace detail

    /** Valid gene types in the genetic algorithms. */
    template<typename T>
    concept Gene = requires
    {
        requires detail::Hashable<T>;
        requires std::regular<T>;
        requires std::destructible<T>;
        requires std::three_way_comparable<T>;
    };

    /** Genetic algorithm types. */
    template<typename T>
    concept GeneticAlgorithm = detail::DerivedFromSpecializationOf<T, GA>;

} // namespace genetic_algorithm

#endif // !GA_CONCEPTS_HPP