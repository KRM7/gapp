/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_ALGORITHM_BASE_DECL_HPP
#define GA_ALGORITHM_ALGORITHM_BASE_DECL_HPP

#include "algorithm_base.fwd.hpp"
#include "selection_base.hpp"
#include "updater_base.hpp"
#include "../population/population.hpp"
#include <vector>
#include <cstddef>
#include <optional>

namespace genetic_algorithm
{
    class GaInfo;

} // namespace genetic_algorithm

/** Algorithm types that can be used in the genetic algorithms (both single- and multi-objective algorithms). */
namespace genetic_algorithm::algorithm
{
    /**
    * Base class used for all of the algorithms. \n
    * 
    * The algorithms define the way the population is evolved over the generations (i.e. the selection and
    * population update methods used). They may be single-, multi-objective, or both, and have 5
    * methods that can be implemented in the derived classes: \n
    * 
    *  - initializeImpl        (optional) : Initializes the algorithm (at the start of a run). \n
    *  - prepareSelectionsImpl (optional) : Prepares the algorithm for the selections if needed. \n
    *  - selectImpl                       : Selects a candidate from the population for crossover (this should be thread-safe). \n
    *  - nextPopulationImpl               : Selects the candidates of the next population from the parent and the child populations. \n 
    *  - optimalSolutionsImpl  (optional) : Selects the optimal solutions of the population. \n
    */
    class Algorithm : protected selection::Selection, protected update::Updater
    {
    public:
        using FitnessVector = detail::FitnessVector;
        using FitnessMatrix = detail::FitnessMatrix;

        /**
        * Initialize the algorithm if needed. \n
        * 
        * This method will be called exactly once at start of the run,
        * after the initial population has already been created. \n
        * 
        * Implemented by initializeImpl. \n
        * The default implementation does nothing.
        *
        * @param ga The GA that uses the algorithm.
        */
        void initialize(const GaInfo& ga) { initializeImpl(ga); };

        /**
        * Prepare the algorithm for the selections if neccesary. \n
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
        * The population will be the unchanged population that was returned by nextPopulation in the previous generation. \n
        * 
        * Implemented by selectImpl. \n
        * The implementation should be thread-safe if parallel execution is enabled for the GAs (enabled by default).
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
        * This method will be called exactly once at the end of each generation before the call to optimalSolutions. \n
        * 
        * Implemented by nextPopulationImpl. \n
        *
        * @param ga The GA that uses the algorithm.
        * @param parents The parent population (current population of the GA).
        * @param children The child population, created from the parent population.
        * 
        * @returns The candidates of the next generation of the algorithm.
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

    private:

        /**
        * Implementation of the optimalSolutions function. \n
        * Returns the indices of the optimal solutions in the current population. \n
        * 
        * If the optimal solutions can't be found trivially, it should just return an empty
        * vector (this is the default behaviour).
        */
        virtual std::vector<size_t> optimalSolutionsImpl(const GaInfo&) const { return {}; }
    };

} // namespace genetic_algorithm::algorithm

#endif // !GA_ALGORITHM_ALGORITHM_BASE_DECL_HPP