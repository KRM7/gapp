/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_SELECTION_BASE_HPP
#define GA_SELECTION_BASE_HPP

#include "../population/population.hpp"
#include <vector>
#include <cstddef>

namespace genetic_algorithm
{
    class GaInfo;
}

/** Selection methods for the genetic algorithms. */
namespace genetic_algorithm::selection
{
    /**
    * Base class used for all of the selection methods. \n
    * The selection methods define most parts of a genetic algorithm (eg. single-objective or multi-objective,
    * how to create the next population etc.), and not just the method for selecting a candidate from the population.
    */
    class Selection
    {
    public:
        using FitnessVector = detail::FitnessVector;
        using FitnessMatrix = detail::FitnessMatrix;

        /** 
        * Initialize the selection method if needed. \n
        * Will be called exactly once at start of the genetic algorithm after the initial population
        * has already been created. \n
        * The default implementation of this function does nothing.
        * 
        * @param ga The genetic algorithm that uses the selection method.
        */
        virtual void init(const GaInfo& ga);

        /**
        * Prepare the selection method for the selections beforehand if neccesary. \n
        * Called exactly once every generation before the selections take place. \n
        * The default implementation of this function does nothing.
        * 
        * @param ga The genetic algorithm that uses the selection method.
        * @param population_fmat The fitness matrix of the current population of the algorithm.
        */
        virtual void prepare(const GaInfo& ga, const FitnessMatrix& population_fmat);

        /**
        * Select a single Candidate from the population. \n
        * Called (population_size) number of times in every generation.
        * 
        * @param ga The genetic algorithm that uses the selection method.
        * @param population_fmat The fitness matrix of the current population.
        * @returns The index of the selected Candidate in the fitness matrix.
        */
        virtual size_t select(const GaInfo& ga, const FitnessMatrix& population_fmat) const = 0;

        /**
        * Select the Candidates of the next generation (next population) from the Candidates of the
        * current population and the child population generated from the current population. \n
        * Called once at the end of each generation. \n
        * 
        * The default implementation simply chooses the best (population_size) number of candidates
        * from the combined current and child populations (assuming fitness maximization).
        * 
        * @param ga The genetic algorithm that uses the selection method.
        * @param population_fmat The fitness matrix of the current population and the children of the algorithm.
        * @returns The indices selected from the fitness matrix for the next population of the algorithm.
        */
        virtual std::vector<size_t> nextPopulation(const GaInfo& ga, const FitnessMatrix& population_fmat);

        Selection()                             = default;
        Selection(const Selection&)             = default;
        Selection(Selection&&)                  = default;
        Selection& operator=(const Selection&)  = default;
        Selection& operator=(Selection&&)       = default;
        virtual ~Selection()                    = default;
    };

} // namespace genetic_algorithm::selection

#endif // !GA_SELECTION_BASE_HPP