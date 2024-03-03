/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_INDESTRUCTIBLE_HPP
#define GA_UTILITY_INDESTRUCTIBLE_HPP

#include <memory>
#include <type_traits>
#include <utility>

namespace gapp::detail
{
    // NOLINTBEGIN(*union-access)

    template<typename T>
    class Indestructible
    {
    public:
        template<typename... Args>
        constexpr Indestructible(Args&&... args) noexcept(std::is_nothrow_constructible_v<T, Args...>)
        {
            std::construct_at(std::addressof(data_), std::forward<Args>(args)...);
        }

        constexpr ~Indestructible() noexcept {} // NOLINT(*default)

        Indestructible(const Indestructible&)            = delete;
        Indestructible(Indestructible&&)                 = delete;
        Indestructible& operator=(const Indestructible&) = delete;
        Indestructible& operator=(Indestructible&&)      = delete;

        constexpr T& get() noexcept { return data_; }
        constexpr const T& get() const noexcept { return data_; }

        constexpr T& operator*() noexcept { return get(); }
        constexpr const T& operator*() const noexcept { return get(); }

        constexpr T* operator->() noexcept { return std::addressof(get()); }
        constexpr const T* operator->() const noexcept { return std::addressof(get()); }

        constexpr explicit(false) operator T&() & noexcept { return get(); }
        constexpr explicit(false) operator const T&() const& noexcept { return get(); }

    private:
        union { T data_; };
    };
    
    // NOLINTEND(*union-access)

} // namespace gapp::detail

#endif // !GA_UTILITY_INDESTRUCTIBLE_HPP
