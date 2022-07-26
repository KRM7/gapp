/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_ALGORITHM_BASE_HPP
#define GA_ALGORITHM_ALGORITHM_BASE_HPP

#include "../population/population.hpp"
#include <vector>
#include <cstddef>
#include <concepts>

namespace genetic_algorithm
{
    class GaInfo;
}

/** Algorithm types that can be used in the genetic algorithms (contains both single- and multi-objective algorithms). */
namespace genetic_algorithm::algorithm
{
    /**
    * Base class used for all of the algorithms. \n
    * 
    * The algorithms define the way the population is evolved over the generations (the selection and
    * population update methods used). They may either be single- or multi-objective (or both), and have 4
    * methods which must be overriden in the derived classes: \n
    * 
    *  - initialize:        Initializes the algorithm. \n
    *  - prepareSelections: Prepares the algorithm for the selections. \n
    *  - select:            Selects a candidate from the population (this should be thread-safe). \n
    *  - nextPopulation:    Select the candidates of the next population from the parents and the children.
    */
    class Algorithm
    {
    public:
        using FitnessVector = detail::FitnessVector;
        using FitnessMatrix = detail::FitnessMatrix;

        /**
        * Initialize the algorithm method if needed. \n
        * This method will be called exactly once at start of the genetic algorithm,
        * after the initial population has already been created. \n
        *
        * @param ga The GA that uses the algorithm.
        */
        virtual void initialize(const GaInfo& ga) = 0;

        /**
        * Prepare the algorithm for the selections beforehand if neccesary. \n
        * This method will be called exactly once every generation right before the selections are performed. \n
        *
        * @param ga The GA that uses the algorithm.
        * @param fmat The fitness matrix of the population of the GA.
        */
        virtual void prepareSelections(const GaInfo& ga, const FitnessMatrix& fmat) = 0;

        /**
        * Select a single candidate from the population. \n
        * This method will be called exactly (population_size) or (population_size + 1)
        * times in every generation (if the population_size is even or odd respectively). \n
        * The implementation should be thread-safe if parallel execution is enabled for the GAs
        * (enabled by default).
        *
        * @param ga The GA that uses the algorithm.
        * @param fmat The fitness matrix of the current population.
        * @returns The selected candidate's index in the fitness matrix.
        */
        virtual size_t select(const GaInfo& ga, const FitnessMatrix& fmat) const = 0;

        /**
        * Select the Candidates of the next generation (next population) from the Candidates of the
        * current population and the child population that was generated from the current population. \n
        * This method will be called exactly once at the end of each generation. \n
        * 
        * The fitness matrix is given as the range [first, last), where
        * the subrange [first, children_first) belongs to the parents, and
        * the subrange [children_first, last) belongs to the children.
        *
        * @param ga The GA that uses the algorithm.
        * @param first The first element of the fitness matrix (first parent).
        * @param children_first The first element of the fitness matrix that belongs to a child.
        * @param last The end of the fitness matrix.
        * @returns The selected candidates' indices in the fitness matrix, assuming that the index of @p first is 0.
        */
        virtual std::vector<size_t> nextPopulation(const GaInfo& ga,
                                                   FitnessMatrix::const_iterator first,
                                                   FitnessMatrix::const_iterator children_first,
                                                   FitnessMatrix::const_iterator last) = 0;


        Algorithm()                             = default;
        Algorithm(const Algorithm&)             = default;
        Algorithm(Algorithm&&)                  = default;
        Algorithm& operator=(const Algorithm&)  = default;
        Algorithm& operator=(Algorithm&&)       = default;
        virtual ~Algorithm()                    = default;
    };

    /** Algorithm types. */
    template<typename T>
    concept AlgorithmType = requires
    {
        requires std::derived_from<T, Algorithm>;
    };

} // namespace genetic_algorithm::algorithm

#endif // !GA_ALGORITHM_ALGORITHM_BASE_HPP