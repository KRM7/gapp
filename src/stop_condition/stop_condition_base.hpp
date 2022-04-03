/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_STOP_CONDITION_BASE_HPP
#define GA_STOP_CONDITION_BASE_HPP

#include "../population/candidate.hpp"

namespace genetic_algorithm
{
    class GaBase;
}

/** Early stop conditions used in the algorithms. */
namespace genetic_algorithm::stopping
{
    /**
    * Base class used for all of the stop conditions.
    * The stop condition will be evaluated only once per generation and
    * should return true if the algorithm should be stopped.
    */
    class StopCondition
    {
    public:
        /** Evaluate the stop condition and return true if the genetic algorithm should stop. */
        virtual bool operator()(const GaBase& ga) = 0;

        StopCondition()                                 = default;
        StopCondition(const StopCondition&)             = default;
        StopCondition(StopCondition&&)                  = default;
        StopCondition& operator=(const StopCondition&)  = default;
        StopCondition& operator=(StopCondition&&)       = default;
        virtual ~StopCondition()                        = default;
    };

} // namespace genetic_algorithm::stopping

#endif