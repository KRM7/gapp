/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

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
        requires ((std::derived_from<Metrics, metrics::MonitorBase> && std::is_final_v<Metrics>) && ...)
        MetricSet(Metrics... metrics);

        template<typename Metric>
        requires (std::derived_from<Metric, metrics::MonitorBase> && std::is_final_v<Metric>)
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
    requires ((std::derived_from<Metrics, metrics::MonitorBase> && std::is_final_v<Metrics>) && ...)
    MetricSet::MetricSet(Metrics... metrics)
    {
        metrics_.reserve(sizeof...(metrics));
        (metrics_.push_back(std::make_unique<Metrics>(std::move(metrics))), ...);
    }

    template<typename Metric>
    requires (std::derived_from<Metric, metrics::MonitorBase> && std::is_final_v<Metric>)
    const Metric* MetricSet::get() const noexcept
    {
        auto found = std::find_if(metrics_.begin(), metrics_.end(), [](const auto& metric) { return metric->type_id() == detail::type_id<Metric>; });

        if (found == metrics_.end()) return nullptr;
        return static_cast<const Metric*>(found->get());
    }

} // namespace gapp::detail

#endif // !GA_METRICS_METRIC_SET_HPP
