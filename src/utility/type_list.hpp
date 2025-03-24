/* Copyright (c) 2025 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_UTILITY_TYPE_LIST_HPP
#define GAPP_UTILITY_TYPE_LIST_HPP

#include "type_traits.hpp"
#include "type_id.hpp"
#include <type_traits>
#include <tuple>
#include <optional>
#include <cstddef>

namespace gapp::detail
{
    template<typename... Ts>
    struct type_list;


    template<typename T>
    struct args_to_list : std::type_identity<type_list<>> {};

    template<template<typename...> class Templ, typename... Ts>
    struct args_to_list<Templ<Ts...>> : std::type_identity<type_list<Ts...>> {};

    template<typename T>
    using args_to_list_t = typename args_to_list<T>::type;


    template<typename T>
    struct tuple_to_list {};

    template<typename... Ts>
    struct tuple_to_list<std::tuple<Ts...>> : std::type_identity<type_list<Ts...>> {};

    template<typename T>
    using tuple_to_list_t = typename tuple_to_list<T>::type;


    template<typename... Ts>
    struct type_list
    {
        static constexpr std::size_t size = sizeof...(Ts);

        template<typename T>
        static constexpr bool contains = (std::is_same_v<Ts, T> || ...);

        template<typename T>
        static constexpr std::size_t index_of = detail::index_of_type_v<T, Ts...>;

        using to_tuple = std::tuple<Ts...>;

        template<template<typename> class Pred>
        using filter_types_t = detail::tuple_to_list_t<detail::filter_types_t<Pred, Ts...>>;

        template<typename F>
        constexpr static auto apply(F&& f)
        {
            return f.template operator()<Ts...>();
        }

        template<typename F>
        constexpr static void for_each(F&& f)
        {
            size_t i = 0;
            (f.template operator()<Ts>(i++), ...);
        }

        template<typename F>
        constexpr static std::optional<size_t> find_index(F&& f)
        {
            size_t i = 0;
            std::optional<size_t> idx;

            ((f.template operator()<Ts>() ? (idx = i, true) : (i++, false)) || ...);

            return idx;
        }

        constexpr static std::optional<size_t> index_of_id(size_t type_id)
        {
            return find_index([=]<typename U>() { return detail::type_id<U>() == type_id; });
        }
    };

} // namespace gapp::detail

#endif // !GAPP_UTILITY_TYPE_LIST_HPP
