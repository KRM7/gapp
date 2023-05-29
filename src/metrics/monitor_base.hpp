/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_METRICS_MONITOR_BASE_HPP
#define GA_METRICS_MONITOR_BASE_HPP

#include <cstddef>

namespace gapp
{
    class GaInfo;

} // namespace gapp

namespace gapp::detail
{
    class MetricSet;

} // namespace gapp::detail

namespace gapp::metrics
{
    /**
    * The base class used for the metrics.
    * 
    * @note The Monitor class should be used as the base class of new metrics
    * instead of directly inheriting from this class.
    */
    class MonitorBase
    {
    public:
        /**
        * Initialize the monitor.
        * This method will be called exactly once at the start of a run.
        * 
        * @param ga The %GA that is being monitored.
        */
        virtual void initialize(const GaInfo& ga) = 0;
        
        /**
        * Update the metric data based on the current generation of the %GA.
        * This method will be called exactly once at the end of every generation.
        * 
        * @param ga The %GA that is being monitored.
        */
        virtual void update(const GaInfo& ga) = 0;

        /** Destructor. */
        virtual ~MonitorBase()                     = default;

    protected:

        MonitorBase()                              = default;
        MonitorBase(const MonitorBase&)            = default;
        MonitorBase(MonitorBase&&)                 = default;
        MonitorBase& operator=(const MonitorBase&) = default;
        MonitorBase& operator=(MonitorBase&&)      = default;

        virtual size_t type_id() const noexcept = 0;

        friend class detail::MetricSet;
    };

} // namespace gapp::metrics

#endif // !GA_METRICS_MONITOR_BASE_HPP
