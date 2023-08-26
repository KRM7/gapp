﻿/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_METRICS_METRIC_SET_HPP
#define GA_METRICS_METRIC_SET_HPP

#include "monitor_base.hpp"
#include <vector>
#include <memory>
#include <type_traits>
#include <concepts>

namespace gapp
{
    class GaInfo;

} // namespace gapp

namespace gapp::detail
{
    /* A collection of metrics to be tracked in the GA. */
    class MetricSet final
    {
    public:
        template<typename... Metrics>
        requires (std::derived_from<Metrics, metrics::MonitorBase> && ...)
        MetricSet(Metrics... metrics);

        template<typename Metric>
        requires std::derived_from<Metric, metrics::MonitorBase>
        const Metric* get() const noexcept;

        void initialize(const GaInfo& ga);
        void update(const GaInfo& ga);

    private:
        std::vector<std::unique_ptr<metrics::MonitorBase>> metrics_;
    };

} // namespace gapp::detail


/* IMPLEMENTATION */

#include "../utility/type_id.hpp"
#include <algorithm>

namespace gapp::detail
{
    template<typename... Metrics>
    requires (std::derived_from<Metrics, metrics::MonitorBase> && ...)
    MetricSet::MetricSet(Metrics... metrics)
    {
        metrics_.reserve(sizeof...(metrics));
        ( metrics_.push_back(std::make_unique<Metrics>(std::move(metrics))), ... );
    }

    template<typename Metric>
    requires std::derived_from<Metric, metrics::MonitorBase>
    const Metric* MetricSet::get() const noexcept
    {
        auto found = std::find_if(metrics_.begin(), metrics_.end(), [](const auto& metric) { return metric->type_id() == detail::type_id<Metric>; });

        return found != metrics_.end() ? static_cast<Metric*>(found->get()) : nullptr;
    }

} // namespace gapp::detail

#endif // !GA_METRICS_METRIC_SET_HPP
