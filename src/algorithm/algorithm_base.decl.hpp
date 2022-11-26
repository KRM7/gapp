/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_ALGORITHM_BASE_DECL_HPP
#define GA_ALGORITHM_ALGORITHM_BASE_DECL_HPP

#include "algorithm_base.fwd.hpp"
#include "../population/population.hpp"
#include <vector>
#include <cstddef>
#include <optional>

namespace genetic_algorithm
{
    class GaInfo;

} // namespace genetic_algorithm

/** Algorithm types that can be used in the genetic algorithms (contains both single- and multi-objective algorithms). */
namespace genetic_algorithm::algorithm
{
    /**
    * Base class used for all of the algorithms. \n
    * 
    * The algorithms define the way the population is evolved over the generations (the selection and
    * population update methods used). They may be single- or multi-objective (or both), and have 5
    * methods that can be implemented in the derived classes: \n
    * 
    *  - initialize        (optional) : Initializes the algorithm (at the start of a run). \n
    *  - prepareSelections (optional) : Prepares the algorithm for the selections if needed. \n
    *  - select                       : Selects a candidate from the population for crossover (this should be thread-safe). \n
    *  - nextPopulation               : Selects the candidates of the next population from the parent and the child populations. \n 
    *  - optimalSolutions  (optional) : Selects the optimal solutions of the population. \n
    */
    class Algorithm
    {
    public:
        using FitnessVector = detail::FitnessVector;
        using FitnessMatrix = detail::FitnessMatrix;

        /**
        * Initialize the algorithm if needed. \n
        * 
        * This method will be called exactly once at start of the genetic algorithm,
        * after the initial population has already been created. \n
        * 
        * Implemented by initializeImpl. \n
        * The default implementation does nothing.
        *
        * @param ga The GA that uses the algorithm.
        */
        void initialize(const GaInfo& ga) { initializeImpl(ga); };

        /**
        * Prepare the algorithm for the selections beforehand if neccesary. \n
        * 
        * This method will be called exactly once every generation right before the selections are performed. \n
        * 
        * Implemented by prepareSelectionsImpl. \n
        * The default implementation does nothing.
        *
        * @param ga The GA that uses the algorithm.
        * @param fmat The fitness matrix of the population.
        */
        void prepareSelections(const GaInfo& ga, const FitnessMatrix& fmat) { prepareSelectionsImpl(ga, fmat); }

        /**
        * Select a single candidate for crossover from the population. \n
        * 
        * This method will be called exactly (population_size) or (population_size + 1)
        * times in every generation (depending on which number is even). \n
        * The population will be the population that was returned by nextPopulation in the previous generation (unchanged). \n
        * 
        * Implemented by selectImpl. \n
        * The implementation should be thread-safe if parallel execution is enabled for the GAs
        * (enabled by default).
        *
        * @param ga The GA that uses the algorithm.
        * @param pop The current population.
        * @param fmat The fitness matrix of the current population.
        * @returns The candidate selected from the population.
        */
        template<Gene T>
        const Candidate<T>& select(const GaInfo& ga, const Population<T>& pop, const FitnessMatrix& fmat) const;

        /**
        * Select the candidates of the next generation from the candidates of the
        * current and the child populations. \n
        * 
        * This method will be called exactly once at the end of each generation
        * (before the call to optimalSolutions). \n
        * 
        * Implemented by nextPopulationImpl. \n
        *
        * @param ga The GA that uses the algorithm.
        * @param parents The parent population (current population of the GA).
        * @param children The child population (created from the parent population).
        * 
        * @param first The first element of the fitness matrix (first parent).
        * @param children_first The first element of the fitness matrix that belongs to a child.
        * @param last The end of the fitness matrix.
        * 
        * @returns The selected candidates' indices in the fitness matrix, assuming that the index of @p first is 0.
        */
        template<Gene T>
        Population<T> nextPopulation(const GaInfo& ga, Population<T>&& parents, Population<T>&& children);

        /**
        * Find the optimal solutions in the unchanged population that was created by nextPopulation. \n
        * 
        * Implemented by optimalSolutionsImpl. \n
        * The default implementation will use detail::findParetoFront to find the optimal
        * solutions in the population.
        * 
        * @param ga The GA that uses the algorithm.
        * @param pop The current population of the GA.
        * @returns The pareto optimal solutions of the population pop.
        */
        template<Gene T>
        Candidates<T> optimalSolutions(const GaInfo& ga, const Population<T>& pop) const;

    protected:

        Algorithm()                             = default;
        Algorithm(const Algorithm&)             = default;
        Algorithm(Algorithm&&)                  = default;
        Algorithm& operator=(const Algorithm&)  = default;
        Algorithm& operator=(Algorithm&&)       = default;

    public:

        virtual ~Algorithm()                    = default;

    private:

        /** Implementation of the initialize function. */
        virtual void initializeImpl(const GaInfo&) {}

        /** Implementation of the prepareSelections function. */
        virtual void prepareSelectionsImpl(const GaInfo&, const FitnessMatrix&) {}

        /** Implementation of the select function. Selects a solution index based on the fitness matrix. */
        virtual size_t selectImpl(const GaInfo& ga, const FitnessMatrix& fmat) const = 0;

        /**
        * Implementation of the nextPopulation function. \n
        * Returns the indices of the fitness vectors that were picked for the next generation
        * (from the fitness matrix). \n
        * 
        * The fitness matrix is given as the range [first, last), where
        * the subrange [first, children_first) belongs to the parents, and
        * the subrange [children_first, last) belongs to the children.
        */
        virtual std::vector<size_t> nextPopulationImpl(const GaInfo& ga,
                                                       FitnessMatrix::const_iterator first,
                                                       FitnessMatrix::const_iterator children_first,
                                                       FitnessMatrix::const_iterator last) = 0;

        /**
        * Implementation of the optimalSolutions function. \n
        * Returns the indices of the optimal solutions in the current population.
        */
        virtual std::optional<std::vector<size_t>> optimalSolutionsImpl(const GaInfo&) const { return {}; }
    };

} // namespace genetic_algorithm::algorithm

#endif // !GA_ALGORITHM_ALGORITHM_BASE_DECL_HPP