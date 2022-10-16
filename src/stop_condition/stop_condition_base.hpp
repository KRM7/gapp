/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_STOP_CONDITION_BASE_HPP
#define GA_STOP_CONDITION_BASE_HPP

#include "stop_condition_base.fwd.hpp"
#include "../population/candidate.hpp"
#include <functional>
#include <utility>

namespace genetic_algorithm
{
    class GaInfo;
}

/** Stop conditions that can be used in the algorithms. */
namespace genetic_algorithm::stopping
{
    /**
    * Base class used for all of the stop conditions. \n
    * The stop condition will be evaluated only once per generation at the end of the generation,
    * and should return true if the algorithm should be stopped.
    */
    class StopCondition
    {
    public:
        /** Evaluate the stop condition and return true if the genetic algorithm should stop. */
        bool operator()(const GaInfo& ga);

        StopCondition()                                 = default;
        StopCondition(const StopCondition&)             = default;
        StopCondition(StopCondition&&)                  = default;
        StopCondition& operator=(const StopCondition&)  = default;
        StopCondition& operator=(StopCondition&&)       = default;
        virtual ~StopCondition()                        = default;

    private:
        /* Implementation of the stop condition. */
        virtual bool stop_condition(const GaInfo& ga) = 0;
    };

} // namespace genetic_algorithm::stopping


namespace genetic_algorithm::stopping
{
    /*
    * Wraps a callable with the right signature so that it can be used as a
    * stop condition in the GAs.
    */
    class Lambda final : public StopCondition
    {
    public:
        using StopConditionFunction = std::function<bool(const GaInfo&)>;

        explicit Lambda(StopConditionFunction f)
            : stop_condition_(std::move(f))
        {}

    private:
        StopConditionFunction stop_condition_;

        bool stop_condition(const GaInfo& ga) override
        {
            return stop_condition_(ga);
        }
    };

} // namespace genetic_algorithm::stopping

#endif