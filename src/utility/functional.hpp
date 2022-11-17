#ifndef GA_UTILITY_FUNCTIONAL_HPP
#define GA_UTILITY_FUNCTIONAL_HPP

#include <functional>
#include <array>
#include <vector>
#include <cmath>
#include <utility>
#include <limits>
#include <type_traits>
#include <concepts>
#include <cassert>

namespace genetic_algorithm
{
    template<auto F>
    constexpr auto Fn = [](auto&&... args) noexcept(std::is_nothrow_invocable_v<decltype(F), decltype(args)...>) -> decltype(auto)
    {
        return std::invoke(F, std::forward<decltype(args)>(args)...);
    };

} // namespace genetic_algorithm

namespace genetic_algorithm::detail
{
    template<typename T>
    constexpr auto lforward(std::remove_reference_t<T>& t) noexcept
    {
        return std::ref<std::remove_reference_t<T>>(t);
    }

    template<typename T>
    requires(!std::is_lvalue_reference_v<T>)
    constexpr T&& lforward(std::remove_reference_t<T>&& t) noexcept
    {
        return static_cast<T&&>(t);
    }

    template<typename F>
    constexpr auto compose(F&& f) noexcept
    {
        return [f = lforward<F>(f)] <typename... Args> (Args&&... args)
        noexcept(std::is_nothrow_invocable_v<F, Args...>) -> decltype(auto)
        requires std::invocable<F, Args...>
        {
            return std::invoke(f, std::forward<Args>(args)...);
        };
    }

    template<typename F, typename... Fs>
    constexpr auto compose(F&& f, Fs&&... fs) noexcept
    {
        return [f = lforward<F>(f), ...fs = lforward<Fs>(fs)] <typename... Args> (Args&&... args)
        noexcept(std::is_nothrow_invocable_v<F, Args...> &&
                 std::is_nothrow_invocable_v<decltype(compose(fs...)), std::invoke_result_t<F, Args...>>) -> decltype(auto)
        requires (std::invocable<F, Args...> &&
                  std::invocable<decltype(compose(fs...)), std::invoke_result_t<F, Args...>>)
        {
            return compose(fs...)(std::invoke(f, std::forward<Args>(args)...));
        };
    }

    template<typename ValueType, std::invocable<ValueType> F>
    auto map(const std::vector<ValueType>& cont, F&& f)
    {
        using MappedType = std::invoke_result_t<std::remove_reference_t<F>, ValueType>;
        using ResultType = std::vector<std::decay_t<MappedType>>;

        ResultType result;
        result.reserve(cont.size());
        for (const auto& elem : cont)
        {
            result.push_back(std::invoke(f, elem));
        }

        return result;
    }

    template<typename ValueType, std::invocable<ValueType> F>
    auto map(const std::vector<ValueType>& cont, F&& f) requires std::is_scalar_v<ValueType>
    {
        using MappedType = std::invoke_result_t<std::remove_reference_t<F>, ValueType>;
        using ResultType = std::vector<std::decay_t<MappedType>>;

        ResultType result(cont.size());
        for (size_t i = 0; i < cont.size(); i++)
        {
            result[i] = std::invoke(f, cont[i]);
        }

        return result;
    }

    template<typename ValueType, size_t N, typename F>
    requires std::invocable<std::remove_reference_t<F>, ValueType>
    constexpr auto map(const std::array<ValueType, N>& cont, F&& f)
    {
        using MappedType = std::invoke_result_t<std::remove_reference_t<F>, ValueType>;
        using ResultType = std::array<std::decay_t<MappedType>, N>;

        ResultType result;
        for (size_t i = 0; i < N; i++)
        {
            result[i] = std::invoke(f, cont[i]);
        }
        return result;
    }

    template<typename T>
    std::vector<T> flatten(const std::vector<std::pair<T, T>>& pairs)
    {
        std::vector<T> flat;
        flat.reserve(2 * pairs.size());

        for (size_t i = 0; i < pairs.size(); i++)
        {
            flat.push_back(pairs[i].first);
            flat.push_back(pairs[i].second);
        }

        return flat;
    }

    template<typename T>
    std::vector<T> flatten(std::vector<std::pair<T, T>>&& pairs)
    {
        std::vector<T> flat;
        flat.reserve(2 * pairs.size());

        for (size_t i = 0; i < pairs.size(); i++)
        {
            flat.push_back(std::move(pairs[i].first));
            flat.push_back(std::move(pairs[i].second));
        }

        return flat;
    }

    template<typename T>
    constexpr auto multiply_by(const T& multiplier)
    noexcept(std::is_nothrow_copy_constructible_v<T>)
    {
        return [=](const auto& value) noexcept(noexcept(value * multiplier))
        {
            return value * multiplier;
        };
    }

    template<typename T>
    constexpr auto divide_by(const T& divisor)
    noexcept(std::is_nothrow_copy_constructible_v<T>)
    {
        return [=](const auto& value) noexcept(noexcept(value / divisor))
        {
            return value / divisor;
        };
    }

    template<typename T>
    constexpr auto add(const T& increment)
    noexcept(std::is_nothrow_copy_constructible_v<T>)
    {
        return [=](const auto& value) noexcept(noexcept(value + increment))
        {
            return value + increment;
        };
    }

    template<typename T>
    constexpr auto subtract(const T& decrement)
    noexcept(std::is_nothrow_copy_constructible_v<T>)
    {
        return [=](const auto& value) noexcept(noexcept(value - decrement))
        {
            return value - decrement;
        };
    }

    template<typename T, typename U>
    constexpr auto multiply_add(const T& multiplier, const U& increment)
    noexcept(std::is_nothrow_copy_constructible_v<T> && std::is_nothrow_copy_constructible_v<U>)
    {
        return [=](const auto& value) noexcept(noexcept(multiplier * value + increment))
        {
            return multiplier * value + increment;
        };
    }

    template<typename T>
    constexpr auto equal_to(const T& rhs)
    noexcept(std::is_nothrow_copy_constructible_v<T>)
    {
        return [=](const auto& lhs) noexcept(noexcept(lhs == rhs))
        {
            return lhs == rhs;
        };
    }

    template<typename T>
    constexpr auto not_equal_to(const T& rhs)
    noexcept(std::is_nothrow_copy_constructible_v<T>)
    {
        return [=](const auto& lhs) noexcept(noexcept(lhs != rhs))
        {
            return lhs != rhs;
        };
    }

    template<typename T>
    constexpr auto greater_than(const T& rhs)
    noexcept(std::is_nothrow_copy_constructible_v<T>)
    {
        return [=](const auto& lhs) noexcept(noexcept(lhs > rhs))
        {
            return lhs > rhs;
        };
    }

    template<typename T>
    constexpr auto greater_eq_than(const T& rhs)
    noexcept(std::is_nothrow_copy_constructible_v<T>)
    {
        return [=](const auto& lhs) noexcept(noexcept(lhs >= rhs))
        {
            return lhs >= rhs;
        };
    }

    template<typename T>
    constexpr auto less_than(const T& rhs)
    noexcept(std::is_nothrow_copy_constructible_v<T>)
    {
        return [=](const auto& lhs) noexcept(noexcept(lhs < rhs))
        {
            return lhs < rhs;
        };
    }

    template<typename T>
    constexpr auto less_eq_than(const T& rhs)
    noexcept(std::is_nothrow_copy_constructible_v<T>)
    {
        return [=](const auto& lhs) noexcept(noexcept(lhs <= rhs))
        {
            return lhs <= rhs;
        };
    }

    constexpr auto is_size(size_t size) noexcept
    {
        return [=](const auto& container) noexcept
        {
            return container.size() == size;
        };
    }
    
} // namespace genetic_algorithm::detail

#endif // !GA_UTILITY_FUNCTIONAL_HPP