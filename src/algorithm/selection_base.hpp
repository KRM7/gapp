/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_ALGORITHM_SOGA_SELECTION_BASE_HPP
#define GAPP_ALGORITHM_SOGA_SELECTION_BASE_HPP

#include "../core/population.hpp"
#include <type_traits>
#include <cstddef>

namespace gapp
{
    class GaInfo;

} // namespace gapp

/** %Selection operators for the single-objective algorithm. */
namespace gapp::selection
{
    /**
    * This is the base class used for all of the single-objective selection operators.
    * The selection operator is used to select candidates from the population for
    * the crossovers.
    * 
    * New selection methods for the single-objective algorithm should be derived from
    * this class and there are 3 methods that should be implemented by them:
    * 
    *  - initializeImpl        (optional) : Initializes the selection method at the start of a run.
    *  - prepareSelectionsImpl (optional) : Prepare the operator for the selections if necessary.
    *  - selectImpl                       : Select a candidate from the population for crossover.
    */
    class Selection
    {
    public:
        /**
        * Initialize the selection operator if necessary.
        *
        * This method will be called exactly once at start of the run,
        * after the initial population of the %GA has already been created.
        *
        * The default implementation does nothing.
        *
        * @param ga The %GA that uses the algorithm.
        */
        virtual void initializeImpl(const GaInfo&) {}

        /**
        * Prepare the operator for the selections if necessary.
        *
        * This method will be called exactly once every generation right before the
        * selections are performed.
        *
        * The default implementation does nothing.
        *
        * @param ga The %GA that uses the algorithm.
        * @param pop A view of the population without the encoding dependent parts of the candidates.
        */
        virtual void prepareSelectionsImpl(const GaInfo&, const PopulationView&) {}

        /**
        * Select a single candidate for crossover from the population.
        * 
        * This method will be called exactly (population_size) or (population_size + 1)
        * times in every generation (depending on which number is even).
        * 
        * The method should return the index of the selected candidate from the current
        * population @p pop.
        *
        * The implementation should be thread-safe if parallel execution is enabled for the GAs
        * (which is true by default).
        *
        * @param ga The %GA that uses the algorithm.
        * @param pop A view of the population without the encoding dependent parts of the candidates.
        * @returns The index of the candidate selected from the population @p pop.
        */
        virtual size_t selectImpl(const GaInfo& ga, const PopulationView& pop) const = 0;


        /** Destructor. */
        virtual ~Selection()                    = default;

    protected:

        Selection()                             = default;
        Selection(const Selection&)             = default;
        Selection(Selection&&)                  = default;
        Selection& operator=(const Selection&)  = default;
        Selection& operator=(Selection&&)       = default;

        template<typename, typename>
        friend class AnyObjective;
    };

} // namespace gapp::selection

#endif // !GAPP_ALGORITHM_SOGA_SELECTION_BASE_HPP
