#ifndef GA_FUNCTIONAL_HPP
#define GA_FUNCTIONAL_HPP

#include "concepts.hpp"
#include "type_traits.hpp"
#include <array>
#include <vector>
#include <list>
#include <forward_list>
#include <deque>
#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <utility>
#include <limits>
#include <type_traits>
#include <cassert>

namespace genetic_algorithm::detail
{
    namespace _
    {
        template<template<typename...> class T>
        inline constexpr bool is_set_container =
            is_one_of_templ_v<T,
                std::set,
                std::multiset,
                std::unordered_set,
                std::unordered_multiset>;

        template<template<typename...> class T>
        inline constexpr bool is_map_container =
            is_one_of_templ_v<T,
                std::map,
                std::multimap,
                std::unordered_map,
                std::unordered_multimap>;

        template<template<typename...> class T>
        inline constexpr bool is_unordered_container =
            is_one_of_templ_v<T,
                std::unordered_set,
                std::unordered_multiset,
                std::unordered_map,
                std::unordered_multimap>;

        template<template<typename...> class T>
        concept MapContainer = requires { is_map_container<T>; };
    }

    template<typename ValueType, typename... Rest, typename F>
    requires std::invocable<F, ValueType>
    auto map(const std::vector<ValueType, Rest...>& cont, F&& f)
    {
        using MappedType = std::invoke_result_t<F, ValueType>;
        using ResultType = std::vector<MappedType>;

        if constexpr (std::is_scalar_v<ResultType>)
        {
            ResultType result(cont.size());
            for (size_t i = 0; i < cont.size(); i++)
            {
                result[i] = std::invoke(f, cont[i]);
            }
            return result;
        }
        else
        {
            ResultType result;
            result.reserve(cont.size());
            for (const auto& elem : cont)
            {
                result.push_back(std::invoke(f, elem));
            }
            return result;
        }
    }

    template<typename ValueType, size_t N, typename F>
    requires std::invocable<F, ValueType>
    constexpr auto map(const std::array<ValueType, N>& cont, F&& f)
        noexcept(std::is_nothrow_invocable_v<F, ValueType> &&
                 std::is_nothrow_move_assignable_v<std::invoke_result_t<F, ValueType>> &&
                 std::is_nothrow_default_constructible_v<std::invoke_result_t<F, ValueType>>)
    {
        using MappedType = std::invoke_result_t<F, ValueType>;
        using ResultType = std::array<MappedType, N>;

        ResultType result;
        for (size_t i = 0; i < N; i++)
        {
            result[i] = std::invoke(f, cont[i]);
        }
        return result;
    }

    template<template<typename...> class ContainerType, typename ValueType, typename... Rest, typename F>
    requires Container<ContainerType<ValueType, Rest...>> && std::invocable<F, ValueType>
    auto map(const ContainerType<ValueType, Rest...>& cont, F&& f)
    {
        using MappedType = std::invoke_result_t<F, ValueType>;
        using ResultType = ContainerType<MappedType>;

        ResultType result;

        if constexpr (is_one_of_templ_v<ContainerType, std::deque, std::list>)
        {
            for (const auto& elem : cont)
            {
                result.push_back(std::invoke(f, elem));
            }
            return result;
        }
        else if constexpr (_::is_set_container<ContainerType>)
        {
            if constexpr (_::is_unordered_container<ContainerType>)
            {
                result.reserve(cont.size());
            }
            for (const auto& elem : cont)
            {
                result.insert(std::invoke(f, elem));
            }
            return result;
        }
        else if constexpr (is_same_template_v<ContainerType, std::forward_list>)
        {
            for (const auto& elem : cont)
            {
                result.push_front(std::invoke(f, elem));
            }
            result.reverse();
            return result;
        }
    }

    template<template<typename...> class ContainerType, typename KeyType, typename ValueType, typename... Rest, typename F>
    requires Container<ContainerType<KeyType, ValueType, Rest...>> &&
             _::MapContainer<ContainerType> &&
             std::invocable<F, KeyType, ValueType>
    auto map(const ContainerType<KeyType, ValueType, Rest...>& cont, F&& f)
    {
        using MappedType = std::invoke_result_t<F, KeyType, ValueType>;
        static_assert(is_specialization_of_v<MappedType, std::pair>);
        using ResultType = ContainerType<typename MappedType::first_type, typename MappedType::second_type>;

        ResultType result;
        if constexpr (_::is_unordered_container<ContainerType>)
        {
            result.reserve(cont.size());
        }

        for (const auto& [key, value] : cont)
        {
            result.insert(std::invoke(f, key, value));
        }

        return result;
    }

    template<typename T>
    std::vector<T> flatten(std::vector<std::vector<T>>&& in)
    {
        size_t out_size = 0;
        for (const auto& vec : in)
        {
            assert( out_size <= (std::numeric_limits<size_t>::max() - vec.size()) );
            out_size += vec.size();
        }

        std::vector<T> out;
        out.reserve(out_size);

        for (size_t i = 0; i < in.size(); i++)
        {
            for (size_t j = 0; j < in[i].size(); j++)
            {
                out.push_back(std::move(in[i][j]));
            }
        }

        return out;
    }

    template<typename T>
    std::vector<T> flatten(std::vector<std::pair<T, T>>&& in)
    {
        assert( in.size() <= (std::numeric_limits<size_t>::max() / 2) );

        std::vector<T> out;
        out.reserve(2 * in.size());

        for (size_t i = 0; i < in.size(); i++)
        {
            out.push_back(std::move(in[i].first));
            out.push_back(std::move(in[i].second));
        }

        return out;
    }

} // namespace genetic_algorithm::detail

#endif // !GA_FUNCTIONAL_HPP