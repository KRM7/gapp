/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_ALGORITHM_HPP
#define GA_UTILITY_ALGORITHM_HPP

#include "concepts.hpp"
#include "type_traits.hpp"
#include "utility.hpp"
#include <vector>
#include <tuple>
#include <optional>
#include <algorithm>
#include <numeric>
#include <functional>
#include <iterator>
#include <random>
#include <type_traits>
#include <concepts>
#include <utility>
#include <cstddef>
#include <cassert>

namespace genetic_algorithm::detail
{
    template<std::integral T>
    constexpr void increment_mod(T& value, T mod)
    {
        GA_ASSERT(mod > 0);
        GA_ASSERT(0 <= value && value < mod);

        value = (value + 1 == mod) ? T(0) : value + 1;
    }

    inline std::vector<size_t> index_vector(size_t n, size_t first = 0)
    {
        std::vector<size_t> indices(n);
        std::iota(indices.begin(), indices.end(), first);

        return indices;
    }

    template<std::random_access_iterator Iter, typename Comp = std::less<typename std::iterator_traits<Iter>::value_type>>
    requires std::strict_weak_order<Comp, typename std::iterator_traits<Iter>::value_type,
                                          typename std::iterator_traits<Iter>::value_type>
    std::vector<size_t> argsort(Iter first, Iter last, Comp&& comp = std::less<typename std::iterator_traits<Iter>::value_type>{})
    {
        GA_ASSERT(std::distance(first, last) >= 0);

        auto indices = detail::index_vector(last - first);

        if constexpr (detail::is_reverse_iterator_v<Iter>)
        {
            const size_t last_idx = std::distance(first, last) - 1; // wraparound is ok
            std::sort(indices.begin(), indices.end(), [&](size_t lidx, size_t ridx)
            {
                return std::invoke(comp, *(first + last_idx - ridx), *(first + last_idx - lidx));
            });
        }
        else
        {
            std::sort(indices.begin(), indices.end(), [&](size_t lidx, size_t ridx)
            {
                return std::invoke(comp, *(first + lidx), *(first + ridx));
            });
        }

        return indices;
    }

    template<std::random_access_iterator Iter, typename Comp = std::less<typename std::iterator_traits<Iter>::value_type>>
    requires std::strict_weak_order<Comp, typename std::iterator_traits<Iter>::value_type,
                                          typename std::iterator_traits<Iter>::value_type>
    std::vector<size_t> partial_argsort(Iter first, Iter middle, Iter last, Comp&& comp = std::less<typename std::iterator_traits<Iter>::value_type>{})
    {
        GA_ASSERT(std::distance(first, middle) >= 0);
        GA_ASSERT(std::distance(middle, last) >= 0);

        auto indices = detail::index_vector(last - first);

        if constexpr (detail::is_reverse_iterator_v<Iter>)
        {
            const size_t last_idx = std::distance(first, last) - 1; // wraparound is ok
            std::partial_sort(indices.begin(), indices.begin() + std::distance(first, middle), indices.end(),
            [&](size_t lidx, size_t ridx)
            {
                return std::invoke(comp, *(first + last_idx - ridx), *(first + last_idx - lidx));
            });
        }
        else
        {
            std::partial_sort(indices.begin(), indices.begin() + std::distance(first, middle), indices.end(),
            [&](size_t lidx, size_t ridx)
            {
                return std::invoke(comp, *(first + lidx), *(first + ridx));
            });
        }

        return indices;
    }

    template<std::random_access_iterator Iter, typename Comp = std::less<typename std::iterator_traits<Iter>::value_type>>
    requires std::strict_weak_order<Comp, typename std::iterator_traits<Iter>::value_type,
                                          typename std::iterator_traits<Iter>::value_type>
    constexpr size_t argmax(Iter first, Iter last, Comp&& comp = std::less<typename std::iterator_traits<Iter>::value_type>{})
    {
        GA_ASSERT(std::distance(first, last) > 0);

        const auto it = std::max_element(first, last, std::forward<Comp>(comp));
        const size_t idx = std::distance(first, it);

        if constexpr (detail::is_reverse_iterator_v<Iter>)
        {
            const size_t last_idx = std::distance(first, last) - 1;
            return last_idx - idx;
        }
        else
        {
            return idx;
        }
    }

    template<std::random_access_iterator Iter, typename Comp = std::less<typename std::iterator_traits<Iter>::value_type>>
    requires std::strict_weak_order<Comp, typename std::iterator_traits<Iter>::value_type,
                                          typename std::iterator_traits<Iter>::value_type>
    constexpr size_t argmin(Iter first, Iter last, Comp&& comp = std::less<typename std::iterator_traits<Iter>::value_type>{})
    {
        GA_ASSERT(std::distance(first, last) > 0);

        const auto it = std::min_element(first, last, std::forward<Comp>(comp));
        const size_t idx = std::distance(first, it);

        if constexpr (detail::is_reverse_iterator_v<Iter>)
        {
            const size_t last_idx = std::distance(first, last) - 1;
            return last_idx - idx;
        }
        else
        {
            return idx;
        }
    }

    template<std::random_access_iterator Iter, typename URBG>
    void partial_shuffle(Iter first, Iter middle, Iter last, URBG&& gen)
    {
        GA_ASSERT(std::distance(first, middle) >= 0);
        GA_ASSERT(std::distance(middle, last) >= 0);

        for (; first != middle; ++first)
        {
            const auto max_offset = std::distance(first, last) - 1;
            const auto offset = std::uniform_int_distribution{ 0_pd, max_offset }(gen);
            const auto new_pos = std::next(first, offset);

            std::iter_swap(first, new_pos);
        }
    }

    template<std::input_iterator Iter>
    constexpr bool contains(Iter first, Iter last, const typename std::iterator_traits<Iter>::value_type& val)
    {
        GA_ASSERT(std::distance(first, last) >= 0);

        return std::any_of(first, last, [&](const auto& elem) { return elem == val; });
    }

    template<std::input_iterator Iter, typename Pred>
    requires std::predicate<Pred, typename std::iterator_traits<Iter>::value_type>
    std::vector<Iter> find_all(Iter first, Iter last, Pred&& pred)
    {
        GA_ASSERT(std::distance(first, last) >= 0);

        std::vector<Iter> result;
        result.reserve(last - first);

        for (; first != last; ++first)
        {
            if (std::invoke(pred, *first)) result.push_back(first);
        }

        return result;
    }

    template<std::input_iterator Iter, typename Pred>
    requires std::predicate<Pred, typename std::iterator_traits<Iter>::value_type>
    auto find_all_v(Iter first, Iter last, Pred&& pred)
    {
        GA_ASSERT(std::distance(first, last) >= 0);

        using ValueType = typename std::iterator_traits<Iter>::value_type;

        std::vector<ValueType> result;
        result.reserve(last - first);

        std::copy_if(first, last, std::back_inserter(result), std::forward<Pred>(pred));

        return result;
    }


    template<typename T, std::predicate<T> Pred>
    std::vector<size_t> find_indices(const std::vector<T>& container, Pred&& pred)
    {
        std::vector<size_t> indices;
        indices.reserve(indices.size());

        for (size_t i = 0; i < container.size(); i++)
        {
            if (std::invoke(pred, container[i])) indices.push_back(i);
        }

        return indices;
    }

    template<typename T>
    std::optional<size_t> index_of(const std::vector<T>& container, const T& val)
    {
        const auto found = std::find(container.begin(), container.end(), val);
        const size_t idx = std::distance(container.begin(), found);
        
        return (idx == container.size()) ? std::optional<size_t>{} : idx;
    }

    template<typename T, std::predicate<T> Pred>
    std::optional<size_t> find_index(const std::vector<T>& container, Pred&& pred)
    {
        const auto found = std::find_if(container.begin(), container.end(), std::forward<Pred>(pred));
        const size_t idx = std::distance(container.begin(), found);

        return (idx == container.size()) ? std::optional<size_t>{} : idx;
    }

    template<typename T>
    std::vector<T> elementwise_min(std::vector<T> left, const std::vector<T>& right)
    {
        GA_ASSERT(left.size() == right.size());

        std::transform(left.begin(), left.end(), right.begin(), left.begin(), [](const T& lhs, const T& rhs) { return lhs < rhs ? lhs : rhs; });

        return left;
    }

    template<typename T>
    std::vector<T> elementwise_max(std::vector<T> left, const std::vector<T>& right)
    {
        GA_ASSERT(left.size() == right.size());

        std::transform(left.begin(), left.end(), right.begin(), left.begin(), [](const T& lhs, const T& rhs) { return lhs < rhs ? rhs : lhs; });

        return left;
    }

    template<typename T>
    bool erase_first_stable(std::vector<T>& container, const T& value)
    {
        const auto found = std::find(container.cbegin(), container.cend(), value);
        if (found != container.cend())
        {
            container.erase(found);
            return true;
        }
        return false;
    }

    template<typename ValueType>
    std::vector<ValueType> select(const std::vector<ValueType>& cont, const std::vector<size_t>& indices)
    {
        GA_ASSERT(std::all_of(indices.begin(), indices.end(), [&](size_t idx) { return idx < cont.size(); }));

        std::vector<ValueType> selected;
        selected.reserve(indices.size());

        std::transform(indices.begin(), indices.end(), std::back_inserter(selected), [&](size_t idx) { return cont[idx]; });

        return selected;
    }

    template<typename ValueType>
    std::vector<ValueType> select(std::vector<ValueType>&& cont, const std::vector<size_t>& indices)
    {
        GA_ASSERT(std::all_of(indices.begin(), indices.end(), [&](size_t idx) { return idx < cont.size(); }));

        std::vector<ValueType> selected;
        selected.reserve(indices.size());

        std::transform(indices.begin(), indices.end(), std::back_inserter(selected), [&](size_t idx) { return std::move(cont[idx]); });

        return selected;
    }

    template<typename T,
             std::predicate<T, T> Pred = std::equal_to<T>,
             std::strict_weak_order<T, T> Comp = std::less<T>>
    constexpr void erase_duplicates(std::vector<T>& container, Pred&& pred = std::equal_to<T>{}, Comp&& comp = std::less<T>{})
    {
        std::sort(container.begin(), container.end(), std::forward<Comp>(comp));
        const auto last = std::unique(container.begin(), container.end(), std::forward<Pred>(pred));
        container.erase(last, container.end());
    }

    namespace tr
    {
        template<typename... Ts, typename R, typename Tr, typename Rd>
        constexpr R transform_reduce_impl(Tr&&, Rd&&, R&& acc, Ts&&...)
        {
            return acc;
        }

        template<typename T, typename... Ts, typename R, typename Tr, typename Rd>
        constexpr R transform_reduce_impl(Tr&& tr, Rd&& rd, R&& acc, T&& arg, Ts&&... args)
        {
            auto transform_result = std::invoke(tr, std::forward<T>(arg));

            acc = static_cast<std::remove_reference_t<R>>(std::invoke(rd, std::forward<R>(acc), std::move(transform_result)));

            return transform_reduce_impl(std::forward<Tr>(tr), std::forward<Rd>(rd), std::forward<R>(acc), std::forward<Ts>(args)...);
        }

    } // namespace tr

    template<typename Tuple, typename Acc, typename TransformOp, typename ReduceOp>
    requires is_specialization_of_v<std::remove_cvref_t<Tuple>, std::tuple>
    constexpr Acc transform_reduce(Tuple&& tup, Acc&& init, TransformOp&& transform, ReduceOp&& reduce)
    {
        auto transform_reduce_ =
        [&] (auto&&... args) mutable -> Acc
        {
            return tr::transform_reduce_impl(std::forward<TransformOp>(transform),
                                             std::forward<ReduceOp>(reduce), 
                                             std::forward<Acc>(init),
                                             std::forward<decltype(args)>(args)...);
        };

        return std::apply(std::move(transform_reduce_), std::forward<Tuple>(tup));
    }

} // namespace genetic_algorithm::detail

#endif // !GA_UTILITY_ALGORITHM_HPP