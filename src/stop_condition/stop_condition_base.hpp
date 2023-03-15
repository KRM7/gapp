/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_STOP_CONDITION_BASE_HPP
#define GA_STOP_CONDITION_BASE_HPP

#include "../utility/utility.hpp"
#include <functional>
#include <utility>
#include <stdexcept>

namespace genetic_algorithm
{
    class GaInfo;

} // namespace genetic_algorithm

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

        /** Destructor. */
        virtual ~StopCondition()                        = default;

    protected:

        StopCondition()                                 = default;
        StopCondition(const StopCondition&)             = default;
        StopCondition(StopCondition&&)                  = default;
        StopCondition& operator=(const StopCondition&)  = default;
        StopCondition& operator=(StopCondition&&)       = default;

    private:

        /* Implementation of the stop condition. */
        virtual bool stop_condition(const GaInfo& ga) = 0;
    };


    /*
    * Wraps a callable with the right signature so that it can be used as a
    * stop condition in the GAs.
    */
    class Lambda final : public StopCondition
    {
    public:
        using StopConditionCallable = std::function<bool(const GaInfo&)>;

        explicit Lambda(StopConditionCallable f) noexcept
        {
            GA_ASSERT(f, "The stop condition can't be a nullptr.");

            stop_condition_ = std::move(f);
        }

    private:
        StopConditionCallable stop_condition_;

        bool stop_condition(const GaInfo& ga) override
        {
            GA_ASSERT(stop_condition_);

            return stop_condition_(ga);
        }
    };

} // namespace genetic_algorithm::stopping

#endif