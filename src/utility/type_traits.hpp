/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_UTILITY_TYPE_TRAITS_HPP
#define GAPP_UTILITY_TYPE_TRAITS_HPP

#include <type_traits>
#include <tuple>
#include <optional>
#include <iterator>
#include <cstddef>

namespace gapp::detail
{
    struct empty_t {};

    struct inplace_t {};


    template<typename... Ts>
    inline constexpr bool always_false = false;

    template<typename... Ts>
    inline constexpr bool always_true = false;


    template<typename T>
    using iterator_t = decltype(std::declval<T&>().begin());

    template<typename T>
    using const_iterator_t = decltype(std::declval<T&>().cbegin());

    template<typename T>
    using value_t = std::iter_value_t<detail::iterator_t<T>>;

    template<typename T>
    using reference_t = std::iter_reference_t<detail::iterator_t<T>>;

    template<typename T>
    using const_reference_t = decltype(*std::declval<detail::const_iterator_t<T>&>());

    template<typename T>
    using size_type = decltype(std::declval<T&>().size());

    template<typename T>
    using difference_t = std::iter_difference_t<detail::iterator_t<T>>;



    template<typename... Ts>
    struct concat_tup { static_assert(always_false<Ts...>, "Can't concatenate non tuple types."); };

    template<typename... Ts, typename... Us>
    struct concat_tup<std::tuple<Ts...>, std::tuple<Us...>> : std::type_identity<std::tuple<Ts..., Us...>> {};

    template<typename T, typename... Us>
    struct concat_tup<T, std::tuple<Us...>> : std::type_identity<std::tuple<T, Us...>> {};

    template<typename... Ts, typename U>
    struct concat_tup<std::tuple<Ts...>, U> : std::type_identity<std::tuple<Ts..., U>> {};

    template<typename... Ts>
    using concat_tup_t = typename concat_tup<Ts...>::type;



    template<template<typename> class Pred, typename... Ts>
    struct filter_types {};

    template<template<typename> class Pred>
    struct filter_types<Pred> : std::type_identity<std::tuple<>> {};

    template<template<typename> class Pred, typename T, typename... Ts>
    struct filter_types<Pred, T, Ts...>
    {
        using type = std::conditional_t<Pred<T>::value, concat_tup_t<T, typename filter_types<Pred, Ts...>::type>, typename filter_types<Pred, Ts...>::type>;
    };

    template<template<typename> class Pred, typename... Ts>
    using filter_types_t = typename filter_types<Pred, Ts...>::type;



    template<template<typename...> class F, typename... Ts>
    struct map_types {};

    template<template<typename...> class F, typename... Ts>
    struct map_types<F, std::tuple<Ts...>> : std::type_identity<std::tuple<F<Ts>...>> {};

    template<template<typename...> class F, typename... Ts>
    using map_types_t = typename map_types<F, Ts...>::type;



    template<template<typename...> class T1, template<typename...> class T2>
    struct is_same_template : std::false_type {};

    template<template<typename...> class T>
    struct is_same_template<T, T> : std::true_type {};

    template<template<typename...> class T1, template<typename...> class T2>
    inline constexpr bool is_same_template_v = is_same_template<T1, T2>::value;



    template<template<typename...> class...>
    struct is_one_of_templates : std::false_type {};

    template<template<typename...> class Match,
             template<typename...> class Target,
             template<typename...> class... Targets>
    struct is_one_of_templates<Match, Target, Targets...>
        : std::disjunction<is_same_template<Match, Target>, is_one_of_templates<Match, Targets...>>
    {};

    template<template<typename...> class T, template<typename...> class... Ts>
    inline constexpr bool is_one_of_templates_v = is_one_of_templates<T, Ts...>::value;



    template<typename... Ts>
    struct unique_types : std::true_type {};

    template<typename T, typename... Ts>
    struct unique_types<T, Ts...> : std::conjunction<std::negation<std::is_same<T, Ts>>..., unique_types<Ts...>> {};

    template<typename... Ts>
    inline constexpr bool unique_types_v = unique_types<Ts...>::value;



    template<size_t N, typename... Args>
    struct nth_type { static_assert(always_false<Args...>, "No such index."); };

    template<typename Arg, typename... Args>
    struct nth_type<0, Arg, Args...> : std::type_identity<Arg> {};

    template<size_t N, typename Arg, typename... Args>
    struct nth_type<N, Arg, Args...>
    {
        using type = typename nth_type<N - 1, Args...>::type;
    };

    template<size_t N, typename... Args>
    using nth_type_t = typename nth_type<N, Args...>::type;



    template<typename T, typename... Ts>
    struct index_of_type { static_assert(always_false<T>, "No such type."); };

    template<typename T, typename... Ts>
    struct index_of_type<T, T, Ts...> : std::integral_constant<std::size_t, 0> {};

    template<typename T, typename U, typename... Ts>
    struct index_of_type<T, U, Ts...> : std::integral_constant<std::size_t, index_of_type<T, Ts...>::value + 1> {};

    template<typename T, typename... Ts>
    inline constexpr std::size_t index_of_type_v = index_of_type<T, Ts...>::value;



    template<typename Base, typename Derived>
    struct is_proper_base_of :
        std::conjunction<std::is_base_of<std::remove_cv_t<Base>, std::remove_cv_t<Derived>>,
                         std::negation<std::is_same<std::remove_cv_t<Base>, std::remove_cv_t<Derived>>>>
    {};
    
    template<typename Base, typename Derived>
    inline constexpr bool is_proper_base_of_v = is_proper_base_of<Base, Derived>::value;



    template<typename Derived, template<typename...> class BaseTempl>
    struct is_derived_from_spec_of
    {
    private:
        template<typename... TArgs>
        static std::true_type f(BaseTempl<TArgs...>*);
        static std::false_type f(...);
    public:
        constexpr static bool value = decltype( f(static_cast<Derived*>(nullptr)) )::value;
    };

    template<typename Derived, template<typename...> class BaseTempl>
    inline constexpr bool is_derived_from_spec_of_v = is_derived_from_spec_of<Derived, BaseTempl>::value;



    template<typename S, template<typename...> class Templ>
    struct is_specialization_of : std::false_type {};

    template<template<typename...> class Templ, typename... TArgs>
    struct is_specialization_of<Templ<TArgs...>, Templ> : std::true_type {};

    template<typename S, template<typename...> class Templ>
    inline constexpr bool is_specialization_of_v = is_specialization_of<S, Templ>::value;



    template<typename T>
    struct is_reverse_iterator : std::false_type {};

    template<typename T>
    struct is_reverse_iterator<std::reverse_iterator<T>> : std::bool_constant<!is_reverse_iterator<T>::value> {};

    template<typename T>
    inline constexpr bool is_reverse_iterator_v = is_reverse_iterator<T>::value;



    template<typename T>
    struct dereference : std::type_identity<decltype(*std::declval<T>())> {};

    template<typename T>
    using dereference_t = typename dereference<T>::type;



    template<typename T>
    struct remove_rvalue_ref : std::type_identity<T> {};

    template<typename T>
    struct remove_rvalue_ref<T&&> : std::type_identity<T> {};

    template<typename T>
    using remove_rvalue_ref_t = typename remove_rvalue_ref<T>::type;


    
    template<typename From, typename To>
    using copy_const_t = std::conditional_t<std::is_const_v<std::remove_reference_t<From>>, const To, To>;

    template<typename From, typename To>
    using copy_volatile_t = std::conditional_t<std::is_volatile_v<std::remove_reference_t<From>>, volatile To, To>;

    template<typename From, typename To>
    using copy_cv_t = copy_const_t<From, copy_volatile_t<From, To>>;



    template<typename T>
    struct promoted : std::conditional<std::is_signed_v<T>, std::ptrdiff_t, std::size_t> {};

    template<typename T>
    using promoted_t = typename promoted<T>::type;



    template<typename T>
    struct is_nothrow_dereferenceable : std::bool_constant<noexcept(*std::declval<T>())> {};

    template<typename T>
    inline constexpr bool is_nothrow_dereferenceable_v = is_nothrow_dereferenceable<T>::value;

} // namespace gapp::detail

#endif // !GAPP_UTILITY_TYPE_TRAITS_HPP
