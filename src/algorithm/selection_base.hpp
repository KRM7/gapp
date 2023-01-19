/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_SOGA_SELECTION_BASE_HPP
#define GA_ALGORITHM_SOGA_SELECTION_BASE_HPP

#include "../population/population.hpp"
#include <type_traits>
#include <cstddef>

namespace genetic_algorithm
{
    class GaInfo;

} // namespace genetic_algorithm

namespace genetic_algorithm::selection
{
    using detail::FitnessVector;
    using detail::FitnessMatrix;

    /**
    * Base class used for all of the single-objective selection operators. \n
    * 
    * The selection operator is used to select candidates from the population for
    * the crossovers. The class has 3 methods that can be implemented in the derived
    * classes: \n
    * 
    *  - initializeImpl        (optional) : Initializes the selection method (at the start of a run). \n
    *  - prepareSelectionsImpl (optional) : Prepares the operator for the selections if needed. \n
    *  - selectImpl                       : Selects a candidate from the population for crossover (this should be thread-safe). \n
    */
    class Selection
    {
    public:
        /**
        * Initialize the selection operator if needed. \n
        *
        * This method will be called exactly once at start of the run,
        * after the initial population has already been created. \n
        *
        * The default implementation does nothing.
        *
        * @param ga The GA that uses the algorithm.
        */
        virtual void initializeImpl(const GaInfo&) {}

        /**
        * Prepare the operator for the selections if needed. \n
        *
        * This method will be called exactly once every generation right before the selections are performed. \n
        *
        * The default implementation does nothing.
        *
        * @param ga The GA that uses the algorithm.
        * @param fmat The fitness matrix of the population.
        */
        virtual void prepareSelectionsImpl(const GaInfo&, const FitnessMatrix&) {}

        /**
        * Select a single candidate for crossover from the population. \n
        * The method should return the index of the selected candidate based on the fitness matrix
        * of the current population. (An index from the fitness matrix.) \n
        *
        * This method will be called exactly (population_size) or (population_size + 1)
        * times in every generation (depending on which number is even). \n
        *
        * The implementation should be thread-safe if parallel execution is enabled for the GAs (enabled by default).
        *
        * @param ga The GA that uses the algorithm.
        * @param fmat The fitness matrix of the current population.
        * @returns The index of the candidate selected from the population.
        */
        virtual size_t selectImpl(const GaInfo& ga, const FitnessMatrix& fmat) const = 0;


        virtual ~Selection()                    = default;

    protected:

        Selection()                             = default;
        Selection(const Selection&)             = default;
        Selection(Selection&&)                  = default;
        Selection& operator=(const Selection&)  = default;
        Selection& operator=(Selection&&)       = default;

    };

    /** Single-objective selection method types. */
    template<typename T>
    concept SelectionType = requires
    {
        requires std::derived_from<T, Selection>;
    };

} // namespace genetic_algorithm::selection

#endif // !GA_ALGORITHM_SOGA_SELECTION_BASE_HPP