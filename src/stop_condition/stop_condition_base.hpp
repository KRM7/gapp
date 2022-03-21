/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_STOP_CONDITION_BASE_HPP
#define GA_STOP_CONDITION_BASE_HPP

#include "../candidate.hpp"
#include "../concepts.hpp"

namespace genetic_algorithm
{
    template<typename T>
    class GA;
}

/** Early stop conditions used in the algorithms. */
namespace genetic_algorithm::stopping
{
    /**
    * Base class used for all of the stop conditions.
    * The stop condition will be evaluated only once per generation and
    * should return true if the algorithm should be stopped.
    */
    template<gene T>
    class StopCondition
    {
    public:

        StopCondition() = default;
        virtual ~StopCondition() = default;

        /** Evaluate the stop condition and return true if the genetic algorithm should stop. */
        virtual bool operator()(const GA<T>& ga) = 0;
    };

} // namespace genetic_algorithm::stopping

#endif