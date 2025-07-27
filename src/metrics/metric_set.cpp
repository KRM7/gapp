/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include "metric_set.hpp"

namespace gapp::detail
{
    void MetricSet::initialize(const GaInfo& ga)
    {
        for (auto& metric : metrics_)
        {
            metric->initialize(ga);
        }
    }

    void MetricSet::update(const GaInfo& ga)
    {
        for (auto& metric : metrics_)
        {
            metric->update(ga);
        }
    }

} // namespace gapp::detail
