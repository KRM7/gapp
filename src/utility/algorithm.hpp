/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_HPP
#define GA_ALGORITHM_HPP

#include "concepts.hpp"
#include <vector>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <functional>
#include <utility>
#include <type_traits>
#include <concepts>
#include <cstddef>
#include <cassert>

namespace genetic_algorithm::detail
{
    template<typename T>
    constexpr auto lforward(std::remove_reference_t<T>& t) noexcept
    {
        return std::ref<std::remove_reference_t<T>>(t);
    }

    template<typename T>
    requires(!std::is_lvalue_reference_v<T>)
    constexpr T&& lforward(std::remove_reference_t<T>&& t) noexcept
    {
        return static_cast<T&&>(t);
    }

    template<typename F>
    constexpr auto compose(F&& f) noexcept
    {
        return [f = lforward<F>(f)] <typename... Args>
        (Args&&... args) requires std::invocable<F, Args...>
        {
            return std::invoke(f, std::forward<Args>(args)...);
        };
    }

    template<typename F, typename... Fs>
    constexpr auto compose(F&& f, Fs&&... fs) noexcept
    {
        return [f = lforward<F>(f), ...fs = lforward<Fs>(fs)] <typename... Args>
        (Args&&... args) requires std::invocable<F, Args...>
        {
            return compose(fs...)(std::invoke(f, std::forward<Args>(args)...));
        };
    }

    template<template<typename...> class ContainerType, typename ValueType, typename... Rest, typename F>
    requires Container<ContainerType<ValueType, Rest...>> && std::invocable<F, ValueType>
    auto map(const ContainerType<ValueType, Rest...>& cont, F&& f)
    {
        using ResultType = ContainerType<std::invoke_result_t<F, ValueType>>;
        
        ResultType result;
        if constexpr (std::is_same_v<ResultType, std::vector<std::invoke_result_t<F, ValueType>>>)
        {
            result.reserve(cont.size());
        }

        std::transform(std::begin(cont), std::end(cont), std::back_inserter(result), // only works if the container has push_back
        [f = lforward<F>(f)](const ValueType& elem)
        {
            return std::invoke(f, elem);
        });

        return result;
    }

    template<std::random_access_iterator Iter,
             typename Comp = std::less<typename std::iterator_traits<Iter>::value_type>>
    requires std::strict_weak_order<Comp, typename std::iterator_traits<Iter>::value_type,
                                          typename std::iterator_traits<Iter>::value_type>
    auto argsort(Iter first, Iter last, Comp&& comp = std::less<typename std::iterator_traits<Iter>::value_type>{})
    {
        assert(std::distance(first, last) >= 0);

        std::vector<size_t> indices(std::distance(first, last));
        std::iota(indices.begin(), indices.end(), size_t{ 0 });

        std::sort(indices.begin(), indices.end(),
        [first, comp = lforward<Comp>(comp)](size_t lidx, size_t ridx)
        {
            return std::invoke(comp, *std::next(first, lidx), *std::next(first, ridx));
        });

        return indices;
    }

    template<std::random_access_iterator Iter,
             typename Comp = std::less<typename std::iterator_traits<Iter>::value_type>>
    requires std::strict_weak_order<Comp, typename std::iterator_traits<Iter>::value_type,
                                          typename std::iterator_traits<Iter>::value_type>
    auto partial_argsort(Iter first, Iter middle, Iter last, Comp&& comp = std::less<typename std::iterator_traits<Iter>::value_type>{})
    {
        assert(std::distance(first, middle) >= 0);
        assert(std::distance(middle, last) >= 0);

        std::vector<size_t> indices(std::distance(first, last));
        std::iota(indices.begin(), indices.end(), size_t{ 0 });

        std::partial_sort(indices.begin(), indices.begin() + (middle - first), indices.end(),
        [first, comp = lforward<Comp>(comp)](size_t lidx, size_t ridx)
        {
            return std::invoke(comp, *std::next(first, lidx), *std::next(first, ridx));
        });

        return indices;
    }

    template<std::input_iterator Iter, typename Pred>
    requires std::predicate<Pred, typename std::iterator_traits<Iter>::value_type>
    auto find_all(Iter first, Iter last, Pred&& pred)
    {
        assert(std::distance(first, last) >= 0);

        std::vector<Iter> result;
        while (first != last)
        {
            if (std::invoke(pred, *first))
            {
                result.push_back(first);
            }
            first++;
        }

        return result;
    }

    template<std::input_iterator Iter, typename Pred>
    requires std::predicate<Pred, typename std::iterator_traits<Iter>::value_type>
    auto find_all_v(Iter first, Iter last, Pred&& pred)
    {
        assert(std::distance(first, last) >= 0);

        std::vector<typename std::iterator_traits<Iter>::value_type> result;
        while (first != last)
        {
            if (std::invoke(pred, *first))
            {
                result.push_back(*first);
            }
            first++;
        }

        return result;
    }

    template<std::input_iterator Iter>
    bool contains(Iter first, Iter last, const typename std::iterator_traits<Iter>::value_type& val)
    {
        return std::find(first, last, val) != last;
    }

    template<std::input_iterator Iter, typename Pred>
    requires std::predicate<Pred, typename std::iterator_traits<Iter>::value_type>
    bool contains(Iter first, Iter last, Pred&& pred)
    {
        return std::find_if(first, last, std::forward<Pred>(pred)) != last;
    }

    template<typename T>
    size_t index_of(const std::vector<T>& cont, const T& val)
    {
        auto pos = std::find(cont.begin(), cont.end(), val);
        assert(pos != cont.end());

        return size_t(pos - cont.begin());
    }

    template<typename T, std::predicate<T> Pred>
    size_t index_of(const std::vector<T>& cont, Pred&& pred)
    {
        auto pos = std::find_if(cont.begin(), cont.end(), std::forward<Pred>(pred));
        assert(pos != cont.end());

        return size_t(pos - cont.begin());
    }

    template<std::random_access_iterator Iter, typename Comp = std::less<typename std::iterator_traits<Iter>::value_type>>
    requires std::strict_weak_order<Comp, typename std::iterator_traits<Iter>::value_type,
                                          typename std::iterator_traits<Iter>::value_type>
    size_t argmax(Iter first, Iter last, Comp&& comp = std::less<typename std::iterator_traits<Iter>::value_type>{})
    {
        return size_t(std::max_element(first, last, std::forward<Comp>(comp)) - first);
    }

    template<std::random_access_iterator Iter, typename Comp = std::less<typename std::iterator_traits<Iter>::value_type>>
    requires std::strict_weak_order<Comp, typename std::iterator_traits<Iter>::value_type,
                                          typename std::iterator_traits<Iter>::value_type>
    size_t argmin(Iter first, Iter last, Comp&& comp = std::less<typename std::iterator_traits<Iter>::value_type>{})
    {
        return size_t(std::min_element(first, last, std::forward<Comp>(comp)) - first);
    }

    template<std::random_access_iterator Iter, typename Comp = std::less<typename std::iterator_traits<Iter>::value_type>>
    requires std::strict_weak_order<Comp, typename std::iterator_traits<Iter>::value_type,
                                          typename std::iterator_traits<Iter>::value_type>
    auto argmax_v(Iter first, Iter last, Comp&& comp = std::less<typename std::iterator_traits<Iter>::value_type>{})
    {
        auto it = std::max_element(first, last, std::forward<Comp>(comp));

        return std::make_pair(it - first, *it);
    }

    template<std::random_access_iterator Iter, typename Comp = std::less<typename std::iterator_traits<Iter>::value_type>>
    requires std::strict_weak_order<Comp, typename std::iterator_traits<Iter>::value_type,
                                          typename std::iterator_traits<Iter>::value_type>
    auto argmin_v(Iter first, Iter last, Comp&& comp = std::less<typename std::iterator_traits<Iter>::value_type>{})
    {
        auto it = std::min_element(first, last, std::forward<Comp>(comp));

        return std::make_pair(it - first, *it);
    }

}

#endif // !GA_ALGORITHM_HPP