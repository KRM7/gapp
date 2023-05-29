/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_METRICS_MONITOR_HPP
#define GA_METRICS_MONITOR_HPP

#include "monitor_base.hpp"
#include "../utility/type_id.hpp"
#include "../utility/type_traits.hpp"
#include <vector>
#include <type_traits>
#include <cstddef>

namespace genetic_algorithm::metrics
{
    /**
    * The base class used for all of the metrics.
    * Metrics can be used to track certain attributes of the %GAs
    * in every generation throughout a run.
    * 
    * New metrics should be derived from this class, and they should implement the following 3 methods:
    *   
    *   - initialize : Used to initialize the monitor at the start of a run.
    *   - update     : Used to update the monitored attribute once every generation.
    *   - value_at   : Returns the value of the monitored metric in a specific generation
                        (must be implemented as a public, non-virtual function).
    * 
    * @note The monitor doesn't have access to any encoding specific information about a %GA
    *   by default, so no such information can be tracked.
    * 
    * @tparam Derived The type of the derived class.
    * @tparam MetricData The type of the container used to store the monitored metrics (eg. std::vector).
    */
    template<typename Derived, typename MetricData>
    class Monitor : public MonitorBase
    {
    public:
        /** @returns The value of the tracked metric in the specified @p generation. */
        [[nodiscard]]
        constexpr decltype(auto) operator[](size_t generation) const { return derived().value_at(generation); }

        /** @returns The data collected by the monitor throughout the run. */
        [[nodiscard]]
        constexpr const MetricData& data() const noexcept { return data_; }

        /** @returns The number of data points recorded. Equal to the number of generations the %GA ran for. */
        [[nodiscard]]
        constexpr size_t size() const noexcept { return data_.size(); }

        /** @returns An iterator to the first element of the metric data. */
        constexpr auto begin() const { return data_.begin(); }

        /** @returns An iterator to one past the last element of the metric data. */
        constexpr auto end() const   { return data_.end(); }

    protected:
        MetricData data_;

        constexpr Monitor() noexcept
        {
            static_assert(!std::is_abstract_v<Derived>,
                          "The Derived class should implement all the virtual functions of Monitor.");
            static_assert(detail::is_derived_from_spec_of_v<Derived, Monitor>,
                          "The first type parameter of the Monitor class must be the derived class.");
            static_assert(std::is_invocable_v<decltype(&Derived::value_at), const Derived&, size_t>,
                          "Classes derived from Monitor must implement a public 'value_at(size_t) const' function.");
        }

        size_t type_id() const noexcept final { return detail::type_id<Derived>; }

    private:
        constexpr const Derived& derived() const noexcept { return static_cast<const Derived&>(*this); }
    };

} // namespace genetic_algorithm::metrics

#endif // !GA_METRICS_MONITOR_HPP