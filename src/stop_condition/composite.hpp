/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_STOP_CONDITION_COMPOSITE
#define GA_STOP_CONDITION_COMPOSITE

#include "stop_condition_base.hpp"
#include "../core/ga_info.hpp"
#include "../utility/concepts.hpp"
#include <functional>
#include <concepts>
#include <utility>

namespace gapp::stopping::dtl
{
    template<std::derived_from<StopCondition> Left, std::derived_from<StopCondition> Right>
    class OR final : public StopCondition
    {
    public:
        constexpr OR(Left left, Right right) noexcept :
            left_(std::move(left)), right_(std::move(right))
        {}

    private:
        void initialize(const GaInfo& ga) override
        {
            left_.initialize(ga);
            right_.initialize(ga);
        }

        bool stop_condition(const GaInfo& ga) override
        {
            // short-circuiting is fine here, since the GA will be stopped anyway if left_ returns true.
            return std::invoke(left_, ga) || std::invoke(right_, ga);   
        }

        Left left_;
        Right right_;
    };

    template<std::derived_from<StopCondition> Left, std::derived_from<StopCondition> Right>
    class AND final : public StopCondition
    {
    public:
        constexpr AND(Left left, Right right) noexcept :
            left_(std::move(left)), right_(std::move(right))
        {}

    private:
        void initialize(const GaInfo& ga) override
        {
            left_.initialize(ga);
            right_.initialize(ga);
        }

        bool stop_condition(const GaInfo& ga) override
        {
            // Both of the stop conditions should be evaluated, so avoid short circuits with logical ops
            // (the stop conditions might rely on side effects of operator() to maintain their state).
            return std::invoke(left_, ga) * std::invoke(right_, ga);
        }

        Left left_;
        Right right_;
    };

} // namespace gapp::stopping::dtl


namespace gapp::stopping
{
    constexpr auto operator&&(detail::Leaf<StopCondition> auto lhs, detail::Leaf<StopCondition> auto rhs) noexcept
    {
        return dtl::AND{ std::move(lhs), std::move(rhs) };
    }

    constexpr auto operator||(detail::Leaf<StopCondition> auto lhs, detail::Leaf<StopCondition> auto rhs) noexcept
    {
        return dtl::OR{ std::move(lhs), std::move(rhs) };
    }

} // namespace gapp::stopping

#endif // !GA_STOP_CONDITION_COMPOSITE