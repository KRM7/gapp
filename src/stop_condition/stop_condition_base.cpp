/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "stop_condition_base.hpp"
#include "../core/ga_info.hpp"

namespace genetic_algorithm::stopping
{
    bool StopCondition::operator()(const GaInfo& ga)
    {
        return stop_condition(ga) || (ga.generation_cntr() > ga.max_gen());
    }

} // namespace genetic_algorithm::stopping