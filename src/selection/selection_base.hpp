/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_SELECTION_BASE_HPP
#define GA_SELECTION_BASE_HPP

#include "../concepts.hpp"
#include "../population.hpp"

namespace genetic_algorithm
{
    template<typename T>
    class GA;
}

/** Selection methods used in the algorithms. */
namespace genetic_algorithm::selection
{
    /**
    * Base class used for all of the selection methods.
    * The selection methods define most parts of a genetic algorithm (eg. single-objective or multi-objective,
    * how to create the next population etc.), and not just the method for selecting a candidate from the population.
    */
    template<Gene T>
    class Selection
    {
    public:
        Selection() = default;
        Selection(const GA<T>&) {};
        virtual ~Selection() = default;

        /** 
        * Initialize the selection method if needed.
        * Called exactly once at start of the genetic algorithm after the initial population
        * has already been created.
        * 
        * @param ga The genetic algorithm that uses the selection method.
        */
        virtual void init(const GA<T>&) {}

        /**
        * Prepare the selection method for the selections beforehand if neccesary. 
        * Called exactly once every generation before the selections take place.
        * 
        * @param ga The genetic algorithm that uses the selection method.
        * @param pop The current population of the algorithm.
        */
        virtual void prepare(const GA<T>& ga, const Population<T>& pop) = 0;

        /**
        * Select a single Candidate from the population.
        * Called (population_size) number of times in every generation.
        * 
        * @param ga The genetic algorithm that uses the selection method.
        * @param pop The current population of the algorithm.
        * @returns The selected Candidate.
        */
        virtual Candidate<T> select(const GA<T>& ga, const Population<T>& pop) = 0;

        /**
        * Select the Candidates of the next generation (next population) from the Candidates of the
        * current population and the child population generated from the current population.
        * Called once at the end of each generation.
        * 
        * The default implementation simply chooses the best (population_size) number of candidates
        * from the combined current and child populations (assuming fitness maximization).
        * 
        * @param ga The genetic algorithm that uses the selection method.
        * @param pop The current population of the algorithm.
        * @param children The child population generated from the current population.
        * @returns The next population of the algorithm.
        */
        virtual Population<T> nextPopulation(const GA<T>& ga, Population<T>& pop, Population<T>& children);
    };

} // namespace genetic_algorithm::selection


/* IMPLEMENTATION */

#include <algorithm>
#include <iterator>
#include <cstddef>
#include <cassert>

namespace genetic_algorithm::selection
{
    template<Gene T>
    inline Population<T> Selection<T>::nextPopulation(const GA<T>& ga, Population<T>& old_pop, Population<T>& children)
    {
        old_pop.insert(old_pop.end(), std::make_move_iterator(children.begin()),
                                      std::make_move_iterator(children.end()));

        assert(std::all_of(old_pop.begin(), old_pop.end(), [](const Candidate<T>& sol) { return sol.is_evaluated && sol.fitness.size() > 0; }));

        std::partial_sort(old_pop.begin(), old_pop.begin() + ga.population_size(), old_pop.end(),
        [](const Candidate<T>& lhs, const Candidate<T>& rhs)
        {
            return detail::paretoCompareLess(rhs.fitness, lhs.fitness);
        });
        old_pop.resize(ga.population_size());

        return old_pop;
    }

} // namespace genetic_algorithm::selection

#endif // !GA_SELECTION_BASE_HPP