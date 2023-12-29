/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_SCOPE_EXIT_HPP
#define GA_UTILITY_SCOPE_EXIT_HPP

#include <functional>
#include <type_traits>
#include <tuple>
#include <utility>

namespace gapp::detail
{
    template<typename F>
    class [[nodiscard]] scope_exit
    {
    public:
        constexpr explicit scope_exit(F on_exit)
        noexcept(std::is_nothrow_move_constructible_v<F>) :
            on_exit_(std::move(on_exit))
        {}

        scope_exit(const scope_exit&)            = delete;
        scope_exit(scope_exit&&)                 = delete;
        scope_exit& operator=(const scope_exit&) = delete;
        scope_exit& operator=(scope_exit&&)      = delete;

        constexpr ~scope_exit() noexcept
        {
            if (active_) std::invoke(std::move(on_exit_));
        }

        constexpr void release() noexcept { active_ = false; }

    private:
        F on_exit_;
        bool active_ = true;
    };

    template<typename... Ts>
    class [[nodiscard]] restore_on_exit
    {
    public:
        template<typename... Us>
        constexpr explicit restore_on_exit(Us&&... vars)
        noexcept(( std::is_nothrow_constructible_v<std::remove_reference_t<Us>, Us&&> && ... )) :
            vars_(vars...), old_values_(std::forward<Us>(vars)...)
        {}

        restore_on_exit(const restore_on_exit&)            = delete;
        restore_on_exit(restore_on_exit&&)                 = delete;
        restore_on_exit& operator=(const restore_on_exit&) = delete;
        restore_on_exit& operator=(restore_on_exit&&)      = delete;

        constexpr ~restore_on_exit() noexcept { vars_ = std::move(old_values_); }

    private:
        std::tuple<Ts&...> vars_;
        std::tuple<Ts...> old_values_;
    };

    template<typename... Ts>
    restore_on_exit(Ts&&...) -> restore_on_exit<std::remove_reference_t<Ts>...>;

} // namespace gapp::detail

#endif // !GA_UTILITY_SCOPE_EXIT_HPP
