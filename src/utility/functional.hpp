#ifndef GA_FUNCTIONAL_HPP
#define GA_FUNCTIONAL_HPP

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
#include <type_traits>
#include "concepts.hpp"
#include "type_traits.hpp"

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

    template<template<typename...> class ContainerType, typename ValueType, typename... Rest, typename F>
    requires Container<ContainerType<ValueType, Rest...>> && std::invocable<F, ValueType>
    auto map(const ContainerType<ValueType, Rest...>& cont, F&& f)
    {
        using MappedType = std::invoke_result_t<F, ValueType>;
        using ResultType = ContainerType<MappedType>;

        ResultType result;

        if constexpr (is_same_template_v<ContainerType, std::vector> || _::is_unordered_container<ContainerType>)
        {
            result.reserve(cont.size());
        }

        if constexpr (is_one_of_templ_v<ContainerType, std::vector, std::deque, std::list>)
        {
            for (const auto& elem : cont)
            {
                result.push_back(std::invoke(f, elem));
            }
            return result;
        }
        else if constexpr (_::is_set_container<ContainerType>)
        {
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
        else static_assert(false);
    }

    template<typename ValueType, size_t N, typename F>
    requires std::invocable<F, ValueType>
    auto map(const std::array<ValueType, N>& cont, F&& f)
    {
        using MappedType = std::invoke_result_t<F, ValueType>;
        using ResultType = std::array<MappedType, N>;

        ResultType result;

        for (size_t i = 0; i < cont.size(); i++)
        {
            result[i] = std::invoke(f, cont[i]);
        }

        return result;
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

} // namespace genetic_algorithm::detail

#endif // !GA_FUNCTIONAL_HPP