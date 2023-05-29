/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "stop_condition_base.hpp"
#include "../core/ga_info.hpp"

namespace gapp::stopping
{
    bool StopCondition::operator()(const GaInfo& ga)
    {
        return stop_condition(ga) || (ga.generation_cntr() + 1 >= ga.max_gen());
    }

} // namespace gapp::stopping