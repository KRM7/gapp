/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_ALGORITHM_BASE_DECL_HPP
#define GA_ALGORITHM_ALGORITHM_BASE_DECL_HPP

#include "selection_base.hpp"
#include "replacement_base.hpp"
#include "../population/population.hpp"
#include <vector>
#include <cstddef>
#include <optional>

namespace genetic_algorithm
{
    class GaInfo;

    template<typename T>
    class GA;

} // namespace genetic_algorithm

namespace genetic_algorithm::algorithm
{
    /**
    * This is the base class used for all of the algorithms.
    * 
    * The algorithms define the way the population is evolved over the generations (i.e. the selection and
    * population replacement methods used). They may be single-, multi-objective, or both.
    *
    * New algorithms should be derived from this class, and there are 5 virtual methods that should be
    * implemented by them:
    * 
    *  - initializeImpl        (optional) : Initializes the algorithm at the start of a run.
    *  - prepareSelectionsImpl (optional) : Prepares the algorithm for the selections if needed.
    *  - selectImpl                       : Selects a candidate from the population for crossover.
    *  - nextPopulationImpl               : Selects the candidates of the next population from the parent and the child populations.
    *  - optimalSolutionsImpl  (optional) : Selects the optimal solutions of the population.
    */
    class Algorithm : protected selection::Selection, protected replacement::Replacement
    {
    public:
        /**
        * Initialize the algorithm if needed.
        * 
        * This method will be called exactly once at start of the run,
        * after the initial population has already been created.
        * 
        * Implemented by initializeImpl().
        *
        * @param ga The %GA that uses the algorithm.
        */
        void initialize(const GaInfo& ga) { initializeImpl(ga); };

        /**
        * Prepare the algorithm for the selections if neccesary.
        * 
        * This method will be called exactly once every generation before the selections are performed.
        * The population will be unchanged since the last call to nextPopulation().
        * 
        * Implemented by prepareSelectionsImpl().
        *
        * @param ga The %GA that uses the algorithm.
        * @param fmat The fitness matrix of the population.
        */
        void prepareSelections(const GaInfo& ga, const FitnessMatrix& fmat) { prepareSelectionsImpl(ga, fmat); }

        /**
        * Select a single candidate from the population for crossover.
        * 
        * This method will be called exactly (population_size) or (population_size + 1)
        * times in every generation, depending on which one is even.
        * The population will be the unchanged population that was returned by
        * the last call to nextPopulation() in the previous generation.
        * 
        * Implemented by selectImpl().
        *
        * @param ga The %GA that uses the algorithm.
        * @param pop The current population.
        * @param fmat The fitness matrix of the current population.
        * @returns The candidate selected from the population.
        */
        template<typename T>
        const Candidate<T>& select(const GA<T>& ga, const Population<T>& pop, const FitnessMatrix& fmat) const;

        /**
        * Select the candidates of the next generation from the candidates of the
        * current and the child populations.
        * 
        * This method will be called exactly once at the end of each generation before
        * the call to optimalSolutions().
        * 
        * Implemented by nextPopulationImpl().
        *
        * @param ga The %GA that uses the algorithm.
        * @param parents The parent population (current population of the GA).
        * @param children The child population, created from the parent population.
        * @returns The candidates of the next generation of the algorithm.
        */
        template<typename T>
        Population<T> nextPopulation(const GA<T>& ga, Population<T> parents, Population<T> children);

        /**
        * Find the optimal solutions in the population that was created by nextPopulation().
        * 
        * Implemented by optimalSolutionsImpl().
        * The default implementation will use a generic method to find the pareto-optimal
        * solutions in the population. This default will be used if optimalSolutionsImpl()
        * returns an empty vector.
        * 
        * @param ga The %GA that uses the algorithm.
        * @param pop The current population of the GA.
        * @returns The pareto optimal solutions of the population pop.
        */
        template<typename T>
        Candidates<T> optimalSolutions(const GA<T>& ga, const Population<T>& pop) const;

    protected:

        Algorithm()                             = default;
        Algorithm(const Algorithm&)             = default;
        Algorithm(Algorithm&&)                  = default;
        Algorithm& operator=(const Algorithm&)  = default;
        Algorithm& operator=(Algorithm&&)       = default;

    private:

        /**
        * The implementation of the optimalSolutions() function.
        * 
        * Returns the indices of the optimal solutions in the current population, which is
        * the population that was returned by the last call to nextPopulation().
        * If the optimal solutions can't be found trivially, it should just return an empty
        * vector (this is the default behaviour).
        * 
        * @param ga The %GA that uses the algorithm.
        * @returns The indices of the pareto optimal solutions in the current population.
        */
        virtual std::vector<size_t> optimalSolutionsImpl(const GaInfo&) const { return {}; }
    };

} // namespace genetic_algorithm::algorithm

#endif // !GA_ALGORITHM_ALGORITHM_BASE_DECL_HPP