/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "stop_condition_base.hpp"
#include <utility>

namespace genetic_algorithm::stopping
{
    bool StopCondition::operator()(const GaInfo& ga)
    {
        return stop_condition(ga);
    }

} // namespace genetic_algorithm::stopping

namespace genetic_algorithm::stopping::dtl
{
    Lambda::Lambda(StopConditionFunction f)
        : f_(std::move(f))
    {
    }

    bool Lambda::stop_condition(const GaInfo& ga)
    {
        return f_(ga);
    };

} // namespace genetic_algorithm::stopping::dtl