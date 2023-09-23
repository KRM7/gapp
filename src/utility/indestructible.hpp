/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_INDESTRUCTIBLE_HPP
#define GA_UTILITY_INDESTRUCTIBLE_HPP

#include <bit>
#include <new>
#include <memory>
#include <type_traits>
#include <utility>

namespace gapp::detail
{
    template<typename T>
    class Indestructible
    {
    public:
        template<typename... Args>
        constexpr Indestructible(Args&&... args)
        noexcept(std::is_nothrow_constructible_v<T, Args...>)
        {
            ::new(std::addressof(data_)) T(std::forward<Args>(args)...);
        }

        Indestructible(const Indestructible&)            = delete;
        Indestructible(Indestructible&&)                 = delete;
        Indestructible& operator=(const Indestructible&) = delete;
        Indestructible& operator=(Indestructible&&)      = delete;

        ~Indestructible() = default;

        // These can't be constexpr because bit_cast between ptr types isn't constexpr
        T& get() noexcept { return *std::bit_cast<T*>(std::addressof(data_)); }
        const T& get() const noexcept { return *std::bit_cast<const T*>(std::addressof(data_)); }

        T& operator*() noexcept { return get(); }
        const T& operator*() const noexcept { return get(); }

        T* operator->() noexcept { return std::addressof(get()); }
        const T* operator->() const noexcept { return std::addressof(get()); }

        /* implicit */ operator T&() & noexcept { return get(); }
        /* implicit */ operator const T&() const& noexcept { return get(); }

    private:
        using storage_type = unsigned char[sizeof(T)]; // NOLINT(*avoid-c-arrays)
        alignas(T) storage_type data_{};
    };

} // namespace gapp::detail

#endif // !GA_UTILITY_INDESTRUCTIBLE_HPP
