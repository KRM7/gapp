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
        concept derived_from_specialization_of = requires(const Derived & arg)
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

} // namespace genetic_algorithm

#endif // !GA_CONCEPTS_HPP