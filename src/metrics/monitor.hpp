/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_METRICS_MONITOR_HPP
#define GA_METRICS_MONITOR_HPP

#include "monitor_base.hpp"
#include "../utility/type_id.hpp"
#include "../utility/type_traits.hpp"
#include "../utility/utility.hpp"
#include <ranges>
#include <vector>
#include <type_traits>
#include <cstddef>

namespace gapp::metrics
{
    /**
    * The base class used for all of the metrics.
    * Metrics can be used to track certain attributes of the %GAs
    * in every generation throughout a run.
    * 
    * New metrics should be derived from this class, and they should implement the following methods:
    *   
    *   - initialize (optional) : Used to initialize the monitor at the start of a run.
    *   - update                : Used to update the monitored attribute once every generation.
    * 
    * @note The monitor doesn't have access to any encoding specific information about a %GA
    *   by default, so no such information can be tracked.
    * 
    * @tparam Derived The type of the derived class.
    * @tparam MetricData The type of the container used to store the monitored metrics (eg. std::vector).
    */
    template<typename Derived, std::ranges::random_access_range MetricData>
    class Monitor : public MonitorBase
    {
    public:
        /** @returns The value of the tracked metric in the specified @p generation. */
        [[nodiscard]]
        constexpr auto operator[](size_t generation) const noexcept { GAPP_ASSERT(generation < data_.size()); return data_[generation]; }

        /** @returns The data collected by the monitor throughout the run. */
        [[nodiscard]]
        constexpr const MetricData& data() const noexcept { return data_; }

        /** @returns The number of data points recorded. Equal to the number of generations the %GA ran for. */
        [[nodiscard]]
        constexpr size_t size() const noexcept { return data_.size(); }

        /** @returns An iterator to the first element of the metric data. */
        constexpr auto begin() const noexcept { return data_.begin(); }

        /** @returns An iterator to one past the last element of the metric data. */
        constexpr auto end() const noexcept { return data_.end(); }

        void initialize(const GaInfo&) override { data_.clear(); }

    protected:
        MetricData data_;

        size_t type_id() const noexcept final { return detail::type_id<Derived>(); }
    };

} // namespace gapp::metrics

#endif // !GA_METRICS_MONITOR_HPP