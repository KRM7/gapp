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
    class [[nodiscard]] ScopeExit
    {
    public:
        constexpr explicit ScopeExit(F on_exit)
        noexcept(std::is_nothrow_move_constructible_v<F>) :
            on_exit_(std::move(on_exit))
        {}

        ScopeExit(const ScopeExit&)            = delete;
        ScopeExit(ScopeExit&&)                 = delete;
        ScopeExit& operator=(const ScopeExit&) = delete;
        ScopeExit& operator=(ScopeExit&&)      = delete;

        constexpr ~ScopeExit() noexcept
        {
            if (active_) std::invoke(std::move(on_exit_));
        }

        constexpr void release() noexcept { active_ = false; }

    private:
        F on_exit_;
        bool active_ = true;
    };

    template<typename... Ts>
    class [[nodiscard]] RestoreOnExit
    {
    public:
        template<typename... Us>
        constexpr explicit RestoreOnExit(Us&&... vars)
        noexcept(( std::is_nothrow_constructible_v<std::remove_reference_t<Us>, Us&&> && ... )) :
            vars_(vars...), old_values_(std::forward<Us>(vars)...)
        {}

        constexpr ~RestoreOnExit() noexcept { vars_ = std::move(old_values_); }

    private:
        std::tuple<Ts&...> vars_;
        std::tuple<Ts...> old_values_;
    };

    template<typename... Ts>
    RestoreOnExit(Ts&&...) -> RestoreOnExit<std::remove_reference_t<Ts>...>;

} // namespace gapp::detail

#endif // !GA_UTILITY_SCOPE_EXIT_HPP