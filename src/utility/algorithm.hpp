/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_ALGORITHM_HPP
#define GA_UTILITY_ALGORITHM_HPP

#include "concepts.hpp"
#include "type_traits.hpp"
#include "utility.hpp"
#include <vector>
#include <span>
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

namespace gapp::detail
{
    template<typename T, typename U>
    constexpr auto&& forward_like(U&& u) noexcept // from P2445R1
    {
        using NorefU = std::remove_reference_t<U>;

        using CastType = std::conditional_t<std::is_lvalue_reference_v<T>,
            detail::copy_const_t<T, NorefU&>,
            detail::copy_const_t<T, NorefU&&>>;

        return static_cast<CastType>(u);
    }

    inline std::vector<size_t> index_vector(size_t n, size_t first = 0)
    {
        std::vector<size_t> indices(n);
        std::iota(indices.begin(), indices.end(), first);

        return indices;
    }

    template<std::random_access_iterator Iter, typename Comp = std::less<iterator_value_t<Iter>>>
    requires std::strict_weak_order<Comp, iterator_value_t<Iter>, iterator_value_t<Iter>>
    std::vector<size_t> argsort(Iter first, Iter last, Comp&& comp = {})
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

    template<std::random_access_iterator Iter, typename Comp = std::less<iterator_value_t<Iter>>>
    requires std::strict_weak_order<Comp, iterator_value_t<Iter>, iterator_value_t<Iter>>
    std::vector<size_t> partial_argsort(Iter first, Iter middle, Iter last, Comp&& comp = {})
    {
        GA_ASSERT(std::distance(first, middle) >= 0);
        GA_ASSERT(std::distance(middle, last) >= 0);

        if (std::distance(first, middle) >= 0.2 * std::distance(first, last))
        {
            return detail::argsort(first, last, std::forward<Comp>(comp));
        }

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

    template<std::random_access_iterator Iter, typename Comp = std::less<iterator_value_t<Iter>>>
    requires std::strict_weak_order<Comp, iterator_value_t<Iter>, iterator_value_t<Iter>>
    constexpr size_t argmax(Iter first, Iter last, Comp&& comp = {})
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

    template<std::random_access_iterator Iter, typename Comp = std::less<iterator_value_t<Iter>>>
    requires std::strict_weak_order<Comp, iterator_value_t<Iter>, iterator_value_t<Iter>>
    constexpr size_t argmin(Iter first, Iter last, Comp&& comp = {})
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
    constexpr void partial_shuffle(Iter first, Iter middle, Iter last, URBG&& gen)
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
    constexpr bool contains(Iter first, Iter last, const iterator_value_t<Iter>& val)
    {
        GA_ASSERT(std::distance(first, last) >= 0);
        
        return std::find(first, last, val) != last;
    }

    template<std::input_iterator Iter, std::predicate<iterator_value_t<Iter>> Pred>
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

    template<std::input_iterator Iter, std::predicate<iterator_value_t<Iter>> Pred>
    auto find_all_v(Iter first, Iter last, Pred&& pred)
    {
        GA_ASSERT(std::distance(first, last) >= 0);

        std::vector<iterator_value_t<Iter>> result;
        result.reserve(last - first);

        for (; first != last; ++first)
        {
            if (std::invoke(pred, *first)) result.push_back(*first);
        }

        return result;
    }


    template<detail::IndexableContainer Container, std::predicate<typename Container::value_type> Pred>
    std::vector<size_t> find_indices(const Container& container, Pred&& pred)
    {
        std::vector<size_t> indices;
        indices.reserve(container.size() / 2);

        for (size_t i = 0; i < container.size(); i++)
        {
            if (std::invoke(pred, container[i])) indices.push_back(i);
        }

        return indices;
    }

    template<detail::IndexableContainer Container>
    constexpr std::optional<size_t> index_of(const Container& container, const typename Container::value_type& val)
    {
        const auto found = std::find(container.begin(), container.end(), val);
        const size_t idx = std::distance(container.begin(), found);
        
        return (found == container.end()) ? std::optional<size_t>{} : idx;
    }

    template<detail::IndexableContainer Container, std::predicate<typename Container::value_type> Pred>
    constexpr std::optional<size_t> find_index(const Container& container, Pred&& pred)
    {
        const auto found = std::find_if(container.begin(), container.end(), std::forward<Pred>(pred));
        const size_t idx = std::distance(container.begin(), found);

        return (found == container.end()) ? std::optional<size_t>{} : idx;
    }

    template<typename T>
    std::vector<T> elementwise_min(std::vector<T> left, std::span<const std::type_identity_t<T>> right)
    {
        GA_ASSERT(left.size() == right.size());

        std::transform(left.begin(), left.end(), right.begin(), left.begin(),
        [](const T& lhs, const T& rhs)
        {
            return lhs < rhs ? lhs : rhs;
        });

        return left;
    }

    template<typename T>
    std::vector<T> elementwise_max(std::vector<T> left, std::span<const std::type_identity_t<T>> right)
    {
        GA_ASSERT(left.size() == right.size());

        std::transform(left.begin(), left.end(), right.begin(), left.begin(),
        [](const T& lhs, const T& rhs)
        {
            return lhs < rhs ? rhs : lhs;
        });

        return left;
    }

    template<detail::Container T>
    constexpr bool erase_first_stable(T& container, const typename T::value_type& value)
    {
        const auto found = std::find(container.cbegin(), container.cend(), value);
        if (found != container.cend())
        {
            container.erase(found);
            return true;
        }
        return false;
    }

    template<typename Container>
    requires detail::IndexableContainer<std::remove_cvref_t<Container>>
    auto select(Container&& container, const std::vector<size_t>& indices)
    {
        GA_ASSERT(std::all_of(indices.begin(), indices.end(), [&](size_t idx) { return idx < container.size(); }));

        using ValueType = typename std::remove_cvref_t<Container>::value_type;

        std::vector<ValueType> selected;
        selected.reserve(indices.size());

        for (size_t idx : indices)
        {
            selected.push_back(detail::forward_like<Container>(container[idx]));
        }

        return selected;
    }

    template<std::totally_ordered T>
    constexpr void erase_duplicates(std::vector<T>& container)
    {
        std::sort(container.begin(), container.end());
        const auto last = std::unique(container.begin(), container.end());
        container.erase(last, container.end());
    }

    template<std::integral T>
    constexpr void increment_mod(T& value, T mod)
    {
        GA_ASSERT(mod > 0);
        GA_ASSERT(0 <= value && value < mod);

        value = (value + 1 == mod) ? T(0) : value + 1;
    }

} // namespace gapp::detail

#endif // !GA_UTILITY_ALGORITHM_HPP