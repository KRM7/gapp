/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_STOP_CONDITION_COMPOSITE
#define GA_STOP_CONDITION_COMPOSITE

#include "stop_condition_base.fwd.hpp"
#include "../core/ga_info.hpp"
#include <tuple>

namespace genetic_algorithm::stopping
{
    /**
    * A stop condition that is the combination of other stop conditions, and
    * will evaluate to true if ALL of its members evaluate to true.
    */
    template<StopConditionType T1, StopConditionType T2, StopConditionType... Rest>
    class OR final : public StopCondition
    {
    public:
        /**
        * Create a composite stop condition from other stop conditions. \n
        * This composite stop condition will return true if ANY of its members return true.
        */
        OR(T1 first, T2 second, Rest... rest) noexcept;

    private:
        bool stop_condition(const GaInfo& ga) override;

        std::tuple<T1, T2, Rest...> parts_;
    };

    /**
    * A stop condition that is the combination of other stop conditions,
    * and will evaluate to true if ANY of its members evaluate to true.
    */
    template<StopConditionType T1, StopConditionType T2, StopConditionType... Rest>
    class AND final : public StopCondition
    {
    public:
        /**
        * Create a composite stop condition from other stop conditions. \n
        * This composite stop condition will return true if ALL of its members return true.
        */
        AND(T1 first, T2 second, Rest... rest) noexcept;

    private:
        bool stop_condition(const GaInfo& ga) override;

        std::tuple<T1, T2, Rest...> parts_;
    };

} // namespace genetic_algorithm::stopping


/* IMPLEMENTATION */

#include "../utility/algorithm.hpp"
#include <utility>
#include <functional>

namespace genetic_algorithm::stopping
{
    template<StopConditionType T1, StopConditionType T2, StopConditionType... Rest>
    OR<T1, T2, Rest...>::OR(T1 first, T2 second, Rest... rest) noexcept
        : StopCondition(), parts_(std::move(first), std::move(second), std::move(rest)...)
    {}

    template<StopConditionType T1, StopConditionType T2, StopConditionType... Rest>
    bool OR<T1, T2, Rest...>::stop_condition(const GaInfo& ga)
    {
        /* All of the member stop conditions should be evaluated, so avoid short circuits with logical ops. */
        /* (The stop conditions might rely on side effects of operator() to maintain their state.) */
        return detail::transform_reduce(parts_, false, [&](auto&& part) { return part(ga); }, std::plus{});
    }

    template<StopConditionType T1, StopConditionType T2, StopConditionType... Rest>
    AND<T1, T2, Rest...>::AND(T1 first, T2 second, Rest... rest) noexcept
        : StopCondition(), parts_(std::move(first), std::move(second), std::move(rest)...)
    {}

    template<StopConditionType T1, StopConditionType T2, StopConditionType... Rest>
    bool AND<T1, T2, Rest...>::stop_condition(const GaInfo& ga)
    {
        /* All of the member stop conditions should be evaluated, so avoid short circuits with logical ops. */
        /* (The stop conditions might rely on side effects of operator() to maintain their state.) */
        return detail::transform_reduce(parts_, false, [&](auto&& part) { return part(ga); }, std::multiplies{});
    }

} // namespace genetic_algorithm::stopping

#endif // !GA_STOP_CONDITION_COMPOSITE