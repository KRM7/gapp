/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_HPP
#define GA_ALGORITHM_HPP

#include "concepts.hpp"
#include <vector>
#include <tuple>
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
        (Args&&... args) noexcept(std::is_nothrow_invocable_v<F, Args...>)
        requires std::invocable<F, Args...>
        {
            return std::invoke(f, std::forward<Args>(args)...);
        };
    }

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

    template<typename T>
    constexpr size_t index_of(const std::vector<T>& container, const T& val)
    {
        assert(!container.empty());

        for (size_t i = 0; i < container.size(); i++)
        {
            if (container[i] == val) return i;
        }
        assert(false);
        return 0;
    }

    template<typename T, std::predicate<T> Pred>
    constexpr size_t find_index(const std::vector<T>& container, Pred&& pred)
    {
        assert(!container.empty());

        for (size_t i = 0; i < container.size(); i++)
        {
            if (std::invoke(pred, container[i])) return i;
        }
        assert(false);
        return 0;
    }

    template<typename T, std::predicate<T> Pred>
    std::vector<size_t> find_indices(const std::vector<T>& container, Pred&& pred)
    {
        std::vector<size_t> indices;

        for (size_t i = 0; i < container.size(); i++)
        {
            if (std::invoke(pred, container[i])) indices.push_back(i);
        }
        return indices;
    }

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

    template<typename T,
             std::predicate<T, T> Pred = std::equal_to<T>,
             std::strict_weak_order<T, T> Comp = std::less<T>>
    constexpr void erase_duplicates(std::vector<T>& container, Pred&& pred = std::equal_to<T>{}, Comp&& comp = std::less<T>{})
    {
        std::sort(container.begin(), container.end(), std::forward<Comp>(comp));
        auto last = std::unique(container.begin(), container.end(), std::forward<Pred>(pred));
        container.erase(last, container.end());
    }

    namespace _
    {
        template<typename... Ts, typename Acc, typename TransformOp, typename ReduceOp>
        Acc transform_reduce_impl(TransformOp&&, ReduceOp&&, Acc&& acc, Ts&&...)
        {
            return acc;
        }

        template<typename T, typename... Ts, typename Acc, typename TransformOp, typename ReduceOp>
        Acc transform_reduce_impl(TransformOp&& tr, ReduceOp&& rd, Acc&& acc, T&& arg, Ts&&... args)
        {
            auto transform_result = std::invoke(tr, std::forward<T>(arg));

            acc = static_cast<std::remove_reference_t<Acc>>(std::invoke(rd, std::forward<Acc>(acc), std::move(transform_result)));

            return transform_reduce_impl(std::forward<TransformOp>(tr), std::forward<ReduceOp>(rd), std::forward<Acc>(acc), std::forward<Ts>(args)...);
        }
    }

    template<typename Tuple, typename Acc, typename TransformOp, typename ReduceOp>
    requires is_specialization_of_v<std::remove_reference_t<Tuple>, std::tuple>
    Acc transform_reduce(Tuple&& tup, Acc&& init, TransformOp&& transform, ReduceOp&& reduce)
    {
        auto transform_reduce_ =
        [acc = lforward<Acc>(init),
         tr = lforward<TransformOp>(transform),
         rd = lforward<ReduceOp>(reduce)]
        (auto&&... args) mutable -> Acc
        {
            return _::transform_reduce_impl(
                std::forward<TransformOp>(tr),
                std::forward<ReduceOp>(rd),
                std::forward<Acc>(acc),
                std::forward<decltype(args)>(args)...);
        };

        return std::apply(std::move(transform_reduce_), std::forward<Tuple>(tup));
    }


    template<std::random_access_iterator Iter, typename Comp = std::less<typename std::iterator_traits<Iter>::value_type>>
    requires std::strict_weak_order<Comp, typename std::iterator_traits<Iter>::value_type,
                                          typename std::iterator_traits<Iter>::value_type>
    constexpr size_t argmax(Iter first, Iter last, Comp&& comp = std::less<typename std::iterator_traits<Iter>::value_type>{})
    {
        assert(first != last);
        return size_t(std::max_element(first, last, std::forward<Comp>(comp)) - first);
    }

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