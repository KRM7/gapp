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

/** Early-stop conditions that can be used to control when the %GA stops. */
namespace genetic_algorithm::stopping
{
    /**
    * The base class used for all of the early-stop conditions.
    * The stop condition will be evaluated once at the end of every generation,
    * and it should return true when the %GA should stop running.
    * 
    * New stop conditions should be derived from this class, and there are 
    * 2 virtual methods that should be implemented by them:
    * 
    *   - initialize (optional) : Initialize the stop condition at the start of a run. Does nothing by default.
    *   - stop_condition        : Evaluate the stop condition and return true when the %GA should stop running.
    */
    class StopCondition
    {
    public:
        /**
        * Initialize the stop condition.
        * This method will be called exactly once at the start of each run.
        * The default implementation does nothing.
        * 
        * @param ga The genetic algorithm the stop condition is used in.
        */
        virtual void initialize(const GaInfo&) {};

        /**
        * Evaluate the stop condition and return true if the %GA should be stopped.
        * This method will be called exactly once at the end of each generation of
        * a run.
        * 
        * Implemented by stop_condition().
        * 
        * @param ga The genetic algorithm the stop condition is used in.
        * @returns True when the genetic algorithm should stop running.
        */
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

        /**
        * The implementation of the early-stop condition.
        * This method will be called exactly once at the end of each generation
        * of a run, and it should return true when the %GA should be stopped.
        * 
        * @param ga The genetic algorithm the stop condition is used in.
        * @returns True when genetic algorithm should stop running.
        */
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