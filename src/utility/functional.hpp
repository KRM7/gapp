/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_UTILITY_FUNCTIONAL_HPP
#define GAPP_UTILITY_FUNCTIONAL_HPP

#include "type_traits.hpp"
#include "utility.hpp"
#include <algorithm>
#include <functional>
#include <array>
#include <vector>
#include <memory>
#include <limits>
#include <type_traits>
#include <concepts>
#include <utility>
#include <cmath>
#include <cstddef>
#include <cassert>

namespace gapp
{
    template<auto F>
    inline constexpr auto Fn = [](auto&&... args) -> decltype(auto)
    {
        return std::invoke(F, std::forward<decltype(args)>(args)...);
    };

} // namespace gapp

namespace gapp::detail
{
    template<typename ValueType, std::invocable<ValueType&> F>
    auto map(const std::vector<ValueType>& cont, F&& f)
    {
        using MappedType = std::invoke_result_t<F&, ValueType>;
        using ResultType = std::vector<std::decay_t<MappedType>>;

        ResultType result;
        result.reserve(cont.size());
        for (const auto& elem : cont)
        {
            result.push_back(std::invoke(f, elem));
        }

        return result;
    }

    template<typename ValueType, std::invocable<ValueType&> F>
    auto map(const std::vector<ValueType>& cont, F&& f) requires std::is_scalar_v<ValueType>
    {
        using MappedType = std::invoke_result_t<F&, ValueType>;
        using ResultType = std::vector<std::decay_t<MappedType>>;

        ResultType result(cont.size());
        for (size_t i = 0; i < cont.size(); i++)
        {
            result[i] = std::invoke(f, cont[i]);
        }

        return result;
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

    template<typename T>
    constexpr auto between(const T& low, const T& high)
    noexcept(std::is_nothrow_copy_constructible_v<T>)
    {
        return [=](const auto& val) noexcept(noexcept(val < high))
        {
            return (low <= val && val <= high);
        };
    }

    constexpr auto is_size(size_t size) noexcept
    {
        return [=](const auto& container) noexcept
        {
            return container.size() == size;
        };
    }

    constexpr auto element_at(size_t idx) noexcept
    {
        return [=](const auto& container) noexcept(noexcept(container[idx]))
        {
            return container[idx];
        };
    }

    constexpr auto reference_to(const auto& target) noexcept
    {
        return [&](const auto& value) noexcept
        {
            return std::addressof(value) == std::addressof(target);
        };
    }

    constexpr auto element_of(const auto& container) noexcept
    {
        return [&](const auto& elem) noexcept(noexcept(elem == elem))
        {
            return std::any_of(container.begin(), container.end(), detail::equal_to(elem));
        };
    }

    constexpr auto points_into(const auto& container) noexcept
    {
        return [&](const auto* pointer) noexcept
        {
            return pointer && std::any_of(container.begin(), container.end(), detail::reference_to(*pointer));
        };
    }


    template<typename...>
    class function_ref;

    template<typename Ret, typename... Args>
    class function_ref<Ret(Args...)>
    {
    public:
        function_ref() = default;
        function_ref(std::nullptr_t) noexcept {}

        template<typename Callable>
        requires(!std::is_same_v<std::remove_const_t<Callable>, function_ref> && std::is_invocable_r_v<Ret, Callable&, Args...>)
        function_ref(Callable& f) noexcept :
            callable_(reinterpret_cast<void*>(std::addressof(f))),
            invoke_(invoke_fn<Callable>)
        {}

        template<typename Callable>
        requires(!std::is_same_v<std::remove_const_t<Callable>, function_ref> && std::is_invocable_r_v<Ret, Callable&, Args...>)
        function_ref& operator=(Callable& f) noexcept
        {
            callable_ = reinterpret_cast<void*>(std::addressof(f));
            invoke_ = invoke_fn<Callable>;
            return *this;
        }

        Ret operator()(Args... args)
        {
            GAPP_ASSERT(callable_, "Attempting to invoke an empty function_ref.");
            return invoke_(callable_, std::forward<Args>(args)...);
        }

        explicit operator bool() const noexcept { return static_cast<bool>(callable_); }

    private:
        using InvokeFn = Ret(void*, Args...);

        template<typename Callable>
        static Ret invoke_fn(void* f, Args... args) noexcept(std::is_nothrow_invocable_r_v<Ret, Callable&, Args...>)
        {
            return std::invoke(*reinterpret_cast<Callable*>(f), std::forward<Args>(args)...);
        }

        void* callable_ = nullptr;
        InvokeFn* invoke_ = nullptr;
    };


    template<typename...>
    class move_only_function;

    template<typename Ret, typename... Args>
    class move_only_function<Ret(Args...)>
    {
    public:
        move_only_function() noexcept = default;
        move_only_function(std::nullptr_t) noexcept {}

        template<typename F>
        requires(!std::is_same_v<std::remove_reference_t<F>, move_only_function> && std::is_invocable_r_v<Ret, F&, Args...>)
        move_only_function(F&& f) :
            fptr_(std::make_unique<Impl<std::decay_t<F>>>(std::forward<F>(f)))
        {}

        template<typename F>
        requires(!std::is_same_v<std::remove_reference_t<F>, move_only_function> && std::is_invocable_r_v<Ret, F&, Args...>)
        move_only_function& operator=(F&& f)
        {
            fptr_ = std::make_unique<Impl<std::decay_t<F>>>(std::forward<F>(f));
            return *this;
        }

        move_only_function(move_only_function&&)            = default;
        move_only_function& operator=(move_only_function&&) = default;

        Ret operator()(Args... args)
        {
            GAPP_ASSERT(fptr_, "Attempting to invoke an empty move_only_function.");
            return fptr_->invoke(std::forward<Args>(args)...);
        }

        void swap(move_only_function& other) noexcept
        {
            fptr_.swap(other.fptr_);
        }

        explicit operator bool() const noexcept { return bool(fptr_); }

    private:
        struct ImplBase
        {
            virtual Ret invoke(Args...) = 0;
            virtual ~ImplBase() = default;
        };

        template<typename Callable>
        struct Impl : public ImplBase
        {
            Impl(Callable func) :
                func_(std::move(func))
            {}

            Ret invoke(Args... args) override
            {
                return std::invoke(func_, std::forward<Args>(args)...);
            }

            Callable func_;
        };

        std::unique_ptr<ImplBase> fptr_ = nullptr;
    };
    
} // namespace gapp::detail

#endif // !GAPP_UTILITY_FUNCTIONAL_HPP
