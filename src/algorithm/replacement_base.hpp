/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_SOGA_REPLACEMENT_BASE_HPP
#define GA_ALGORITHM_SOGA_REPLACEMENT_BASE_HPP

#include "../core/population.hpp"
#include <vector>
#include <cstddef>

namespace gapp
{
    class GaInfo;

} // namespace gapp

/** Population replacement policies for the single-objective algorithms. */
namespace gapp::replacement
{
    /**
    * This is the base class used for all of the single-objective population replacement policies.
    * The replacement operator is used to select the candidates of the next population from the
    * combined parent and child populations.
    * 
    * New replacement policies for the single-objective algorithm should be derived from this
    * class, and there is a single method that should be implemented in the derived classes:
    * 
    *   - nextPopulationImpl :  Selects the candidates for the next population.
    */
    class Replacement
    {
    public:
        /**
        * Select the candidates of the next generation from the candidates of the
        * combined current and child populations.
        *
        * The fitness matrix is given as the contiguous range [first, last), where
        * the subrange [first, children_first) belongs to the current population,
        * and the [children_first, last) subrange belongs to the child population.
        * 
        * This method should return (population_size) number of unique indices from this fitness matrix,
        * assuming the index of the first iterator is 0.
        * 
        * @param ga The %GA that uses the update method.
        * @param first The first element of the fitness matrix (first parent).
        * @param children_first The first element of the fitness matrix that belongs to a child.
        * @param last The end iterator of the fitness matrix.
        * 
        * @returns The indices of the candidates selected from the fitness matrix.
        */
        virtual std::vector<size_t> nextPopulationImpl(const GaInfo& ga,
                                                       FitnessMatrix::const_iterator first,
                                                       FitnessMatrix::const_iterator children_first,
                                                       FitnessMatrix::const_iterator last) = 0;
        

        /** Destructor. */
        virtual ~Replacement()                      = default;

    protected:

        Replacement()                               = default;
        Replacement(const Replacement&)             = default;
        Replacement(Replacement&&)                  = default;
        Replacement& operator=(const Replacement&)  = default;
        Replacement& operator=(Replacement&&)       = default;
    };

} // namespace gapp::replacement

#endif // !GA_ALGORITHM_SOGA_REPLACEMENT_BASE_HPP