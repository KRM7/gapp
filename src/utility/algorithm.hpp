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
    constexpr auto lforward(T&& t) noexcept
    {
        return std::conditional_t<std::is_lvalue_reference_v<T>,
                                  std::reference_wrapper<std::remove_reference_t<T>>,
                                  T>
               { std::forward<T>(t) };
    }


    template<typename F>
    constexpr auto compose(F&& f) noexcept
    {
        return [f = lforward(f)] <typename... Args>
        (Args&&... args) requires std::invocable<F, Args...>
        {
            return std::invoke(f, std::forward<Args>(args)...);
        };
    }

    template<typename F, typename... Fs>
    constexpr auto compose(F&& f, Fs&&... fs) noexcept
    {
        return [f = lforward(f), ...fs = lforward(fs)] <typename... Args>
        (Args&&... args) requires std::invocable<F, Args...>
        {
            return compose(fs...)(std::invoke(f, std::forward<Args>(args)...));
        };
    }

    template<Container C, typename F>
    requires std::invocable<F, const typename C::value_type>
    auto map(const C& cont, F&& f)
    {
        using result_t = std::vector<std::invoke_result_t<F, const typename C::value_type>>;
        
        result_t result;
        result.reserve(cont.size());

        std::transform(std::begin(cont), std::end(cont), std::back_inserter(result),
        [f = lforward(f)](const C::value_type& elem)
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
        [first, last, comp = lforward(comp)](size_t lidx, size_t ridx)
        {
            return std::invoke(comp, *std::next(first, lidx), *std::next(first, ridx));
        });

        return indices;
    }

    template<std::input_iterator Iter, typename Pred>
    requires std::predicate<Pred, typename std::iterator_traits<Iter>::value_type>
    auto find_all(Iter first, Iter last, Pred&& pred) -> std::vector<Iter>
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
    auto find_all_v(Iter first, Iter last, Pred&& pred) -> std::vector<typename std::iterator_traits<Iter>::value_type>
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

    template<std::random_access_iterator Iter, typename Comp = std::less<typename std::iterator_traits<Iter>::value_type>>
    requires std::strict_weak_order<Comp, typename std::iterator_traits<Iter>::value_type,
                                          typename std::iterator_traits<Iter>::value_type>
    auto argmax(Iter first, Iter last, Comp&& comp = std::less<typename std::iterator_traits<Iter>::value_type>{}) -> size_t
    {
        assert(first < last);

        return size_t(std::max_element(first, last, std::forward<Comp>(comp)) - first);
    }

    template<std::random_access_iterator Iter, typename Comp = std::less<typename std::iterator_traits<Iter>::value_type>>
    requires std::strict_weak_order<Comp, typename std::iterator_traits<Iter>::value_type,
                                          typename std::iterator_traits<Iter>::value_type>
    auto argmin(Iter first, Iter last, Comp&& comp = std::less<typename std::iterator_traits<Iter>::value_type>{}) -> size_t
    {
        assert(first < last);

        return size_t(std::min_element(first, last, std::forward<Comp>(comp)) - first);
    }

    template<std::random_access_iterator Iter, typename Comp = std::less<typename std::iterator_traits<Iter>::value_type>>
    requires std::strict_weak_order<Comp, typename std::iterator_traits<Iter>::value_type,
                                          typename std::iterator_traits<Iter>::value_type>
    auto argmax_v(Iter first, Iter last, Comp&& comp = std::less<typename std::iterator_traits<Iter>::value_type>{})
    {
        assert(first < last);

        auto it = std::max_element(first, last, std::forward<Comp>(comp));

        return std::make_pair(it - first, *it);
    }

    template<std::random_access_iterator Iter, typename Comp = std::less<typename std::iterator_traits<Iter>::value_type>>
    requires std::strict_weak_order<Comp, typename std::iterator_traits<Iter>::value_type,
                                          typename std::iterator_traits<Iter>::value_type>
    auto argmin_v(Iter first, Iter last, Comp&& comp = std::less<typename std::iterator_traits<Iter>::value_type>{})
    {
        assert(first < last);

        auto it = std::min_element(first, last, std::forward<Comp>(comp));

        return std::make_pair(it - first, *it);
    }

}

#endif // !GA_ALGORITHM_HPP