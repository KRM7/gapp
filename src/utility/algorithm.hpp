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
    /* Forward function for lambda captures, same as std::forward, but uses std::ref_wrapper for lvalues. */
    template<typename T>
    constexpr auto lforward(std::remove_reference_t<T>& t) noexcept
    {
        return std::ref<std::remove_reference_t<T>>(t);
    }

    /* Forward function for lambda captures, same as std::forward, but uses std::ref_wrapper for lvalues. */
    template<typename T>
    requires(!std::is_lvalue_reference_v<T>)
    constexpr T&& lforward(std::remove_reference_t<T>&& t) noexcept
    {
        return static_cast<T&&>(t);
    }

    /* Returns a function that is the result of the composition of the argument functions. */
    template<typename F>
    constexpr auto compose(F&& f) noexcept
    {
        return [f = lforward<F>(f)] <typename... Args>
        (Args&&... args) noexcept(std::is_nothrow_invocable_v<F, Args...>)
        requires std::invocable<F, Args...>
        {
            return std::invoke(f, std::forward<Args>(args)...);
        };
    }

    /* Returns a function that is the result of the composition of the argument functions. */
    template<typename F, typename... Fs>
    constexpr auto compose(F&& f, Fs&&... fs) noexcept
    {
        return [f = lforward<F>(f), ...fs = lforward<Fs>(fs)] <typename... Args>
        (Args&&... args) noexcept(std::is_nothrow_invocable_v<F, Args...>)
        requires std::invocable<F, Args...>
        {
            return compose(fs...)(std::invoke(f, std::forward<Args>(args)...));
        };
    }

    /* Returns a vector of the sorted indices (first = 0). */
    template<std::random_access_iterator Iter,
             typename Comp = std::less<typename std::iterator_traits<Iter>::value_type>>
    requires std::strict_weak_order<Comp, typename std::iterator_traits<Iter>::value_type,
                                          typename std::iterator_traits<Iter>::value_type>
    auto argsort(Iter first, Iter last, Comp&& comp = std::less<typename std::iterator_traits<Iter>::value_type>{})
    {
        assert(std::distance(first, last) >= 0);

        std::vector<size_t> indices(std::distance(first, last));
        std::iota(indices.begin(), indices.end(), 0);

        std::sort(indices.begin(), indices.end(),
        [first, comp = lforward<Comp>(comp)](size_t lidx, size_t ridx)
        {
            return std::invoke(comp, *(first + lidx), *(first + ridx));
        });

        return indices;
    }

    /* Returns a vector of the partially sorted (up to middle) indices (first = 0). */
    template<std::random_access_iterator Iter,
             typename Comp = std::less<typename std::iterator_traits<Iter>::value_type>>
    requires std::strict_weak_order<Comp, typename std::iterator_traits<Iter>::value_type,
                                          typename std::iterator_traits<Iter>::value_type>
    auto partial_argsort(Iter first, Iter middle, Iter last, Comp&& comp = std::less<typename std::iterator_traits<Iter>::value_type>{})
    {
        assert(std::distance(first, middle) >= 0);
        assert(std::distance(middle, last) >= 0);

        std::vector<size_t> indices(std::distance(first, last));
        std::iota(indices.begin(), indices.end(), 0);

        std::partial_sort(indices.begin(), indices.begin() + std::distance(first, middle), indices.end(),
        [first, comp = lforward<Comp>(comp)](size_t lidx, size_t ridx)
        {
            return std::invoke(comp, *(first + lidx), *(first + ridx));
        });

        return indices;
    }

    /* Returns iterators to all elements in the range [first, last) that satisfy a predicate. */
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
            ++first;
        }

        return result;
    }

    /* Returns all values in the range [first, last) that satisfy a predicate. */
    template<std::input_iterator Iter, typename Pred>
    requires std::predicate<Pred, typename std::iterator_traits<Iter>::value_type>
    auto find_all_v(Iter first, Iter last, Pred&& pred)
    {
        assert(std::distance(first, last) >= 0);

        using ValueType = typename std::iterator_traits<Iter>::value_type;

        std::vector<ValueType> result;
        while (first != last)
        {
            if (std::invoke(pred, *first))
            {
                result.push_back(*first);
            }
            ++first;
        }

        return result;
    }

    /* Returns true if the range [first, last) contains an element that is equal to val. */
    template<std::input_iterator Iter>
    constexpr bool contains(Iter first, Iter last, const typename std::iterator_traits<Iter>::value_type& val)
    {
        assert(std::distance(first, last) >= 0);

        while (first != last)
        {
            if (*first++ == val) return true;
        }
        return false;
    }

    /* Returns the index of the first element that is equal to val, assuming the container contains this element. */
    template<typename T>
    constexpr size_t index_of(const std::vector<T>& container, const T& val)
    {
        assert(!container.empty());

        auto pos = std::find(container.begin(), container.end(), val);
        size_t idx = size_t(pos - container.begin());

        assert(idx < container.size());
        return idx;  
    }

    /* Returns the index of the first element that satisfies the predicate, assuming the container contains this element. */
    template<typename T, std::predicate<T> Pred>
    constexpr size_t index_of(const std::vector<T>& container, Pred&& pred)
    {
        assert(!container.empty());

        auto pos = std::find_if(container.begin(), container.end(), std::forward<Pred>(pred));
        size_t idx = size_t(pos - container.begin());

        assert(idx < container.size());
        return idx;
    }

    /* Erases the first element in container the is equal to val, without changing the order of the other elements. */
    template<typename T>
    constexpr bool erase_first_v(std::vector<T>& container, const T& val)
    {
        auto pos = std::find(container.cbegin(), container.cend(), val);
        if (pos != container.cend())
        {
            container.erase(pos);
            return true;
        }
        return false;
    }

    /* Erases all duplicate elements from the container. */
    template<typename T,
             std::predicate<T, T> Pred = std::equal_to<T>,
             std::strict_weak_order<T, T> Comp = std::less<T>>
    constexpr void erase_duplicates(std::vector<T>& container, Pred&& pred = std::equal_to<T>{}, Comp&& comp = std::less<T>{})
    {
        std::sort(container.begin(), container.end(), std::forward<Comp>(comp));
        auto last = std::unique(container.begin(), container.end(), std::forward<Pred>(pred));
        container.erase(last, container.end());
    }

    /* Returns the index of the max element in the range [first, last). */
    template<std::random_access_iterator Iter, typename Comp = std::less<typename std::iterator_traits<Iter>::value_type>>
    requires std::strict_weak_order<Comp, typename std::iterator_traits<Iter>::value_type,
                                          typename std::iterator_traits<Iter>::value_type>
    constexpr size_t argmax(Iter first, Iter last, Comp&& comp = std::less<typename std::iterator_traits<Iter>::value_type>{})
    {
        assert(first != last);
        return size_t(std::max_element(first, last, std::forward<Comp>(comp)) - first);
    }

    /* Returns the index of the min element in the range [first, last). */
    template<std::random_access_iterator Iter, typename Comp = std::less<typename std::iterator_traits<Iter>::value_type>>
    requires std::strict_weak_order<Comp, typename std::iterator_traits<Iter>::value_type,
                                          typename std::iterator_traits<Iter>::value_type>
    constexpr size_t argmin(Iter first, Iter last, Comp&& comp = std::less<typename std::iterator_traits<Iter>::value_type>{})
    {
        assert(first != last);
        return size_t(std::min_element(first, last, std::forward<Comp>(comp)) - first);
    }

}

#endif // !GA_ALGORITHM_HPP