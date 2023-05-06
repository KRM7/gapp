/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_METRICS_MONITOR_BASE_HPP
#define GA_METRICS_MONITOR_BASE_HPP

#include <cstddef>

namespace genetic_algorithm
{
    class GaInfo;

} // namespace genetic_algorithm

namespace genetic_algorithm::detail
{
    class MetricSet;

} // namespace genetic_algorithm::detail

/** Metrics that can be used to track certain attributes of a GA throughout a run. */
namespace genetic_algorithm::metrics
{
    /**
    * The base class used for the metrics. The Monitor class should be
    * used as the base class of new metrics instead of directly inheriting from this class.
    */
    class MonitorBase
    {
    public:
        /**
        * Initialize the monitor. \n
        * This method will be called exactly once at the start of a run.
        * 
        * @param ga The algorithm being monitored.
        */
        virtual void initialize(const GaInfo& ga) = 0;
        
        /**
        * Update the metric data based on the current generation of the algorithm. \n
        * This method will be called once every generation.
        * 
        * @param ga The algorithm being monitored.
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

} // namespace genetic_algorithm::metrics

#endif // !GA_METRICS_MONITOR_BASE_HPP
