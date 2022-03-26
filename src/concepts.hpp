/*
*  MIT License
*
*  Copyright (c) 2022 Krisztián Rugási
*
*  Permission is hereby granted, free of charge, to any person obtaining a copy
*  of this software and associated documentation files (the "Software"), to deal
*  in the Software without restriction, including without limitation the rights
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*  copies of the Software, and to permit persons to whom the Software is
*  furnished to do so, subject to the following conditions:
*
*  The above copyright notice and this permission notice shall be included in all
*  copies or substantial portions of the Software.
*
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
*  SOFTWARE.
*/

#ifndef GA_CONCEPTS_HPP
#define GA_CONCEPTS_HPP

#include <concepts>
#include <functional>
#include <type_traits>
#include <limits>
#include <iterator>

namespace genetic_algorithm
{
    template<typename geneType>
    class GA;

    namespace detail
    {
        /* Types that std::hash is specialized for. */
        template<typename T>
        concept hashable = requires(T arg)
        {
            { std::hash<T>{}(arg) } -> std::convertible_to<size_t>;
        };

        /* Types that are regular and std::hash is specialized for. */
        template<typename T>
        concept regular_hashable = std::regular<T> && detail::hashable<T>;

        template<typename T, template<typename...> class Templ>
        struct SpecializationOfImpl : std::false_type {};

        template<template<typename...> class Templ, typename... TArgs>
        struct SpecializationOfImpl<Templ<TArgs...>, Templ> : std::true_type {};

        /* Types that are specializations of a class template. */
        template<typename S, template<typename...> class Templ>
        concept specialization_of = SpecializationOfImpl<S, Templ>::value;

        template<template<typename...> class BaseTempl, typename... TArgs>
        void derivedFromSpecializationOfImpl(const BaseTempl<TArgs...>&);

        /* Types that are derived from a specialization of a class template. */
        template<typename Derived, template<typename...> class BaseTempl>
        concept derived_from_specialization_of = requires(const Derived& arg)
        {
            derivedFromSpecializationOfImpl<BaseTempl>(arg);
        };

    } // namespace detail

    /** Valid gene types in the genetic algorithms. */
    template<typename T>
    concept gene = detail::regular_hashable<T>;

    /** Specializations of the base GA class template. */
    //template<typename T>
    //concept genetic_algorithm = detail::specialization_of<T, GA>;

    /** Genetic algorithm types. (Types that are derived from the GA class.) */
    template<typename T>
    concept genetic_algorithm = detail::derived_from_specialization_of<T, GA>;

    /* based on named requirements */
    template<typename C>
    concept Container = requires(C c, const C cc, C&& rc)
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

} // namespace genetic_algorithm

#endif // !GA_CONCEPTS_HPP