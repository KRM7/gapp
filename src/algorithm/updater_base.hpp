/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_SOGA_UPDATER_BASE_HPP
#define GA_ALGORITHM_SOGA_UPDATER_BASE_HPP

#include "../population/population.hpp"
#include <vector>
#include <cstddef>

namespace genetic_algorithm
{
    class GaInfo;

} // namespace genetic_algorithm

/** Population update methods for the single-objective algorithms. */
namespace genetic_algorithm::update
{
    using detail::FitnessVector;
    using detail::FitnessMatrix;

    /**
    * The base class used for all of the single-objective population update operators. \n
    *
    * This operator is used to select the candidates of the next population from the candidates
    * of the current population and the children created from this population. The class only has
    * a single method that has to be implemented: \n
    * 
    *   - nextPopulationImpl :  Selects the candidates of the next population.
    */
    class Updater
    {
    public:
        /**
        * Select the candidates of the next generation from the candidates of the
        * combined current and child populations. \n
        *
        * The fitness matrix is given as the contiguous range [first, last), where
        * the subrange [first, children_first) belongs to the current population,
        * and the [children_first, last) subrange belongs to the child population.
        * 
        * The method should return (population_size) number of unique indices from this fitness matrix,
        * assuming the index of the first iterator is 0.
        * 
        * @param ga The GA that uses the update method.
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
        virtual ~Updater()                  = default;

    protected:

        Updater()                           = default;
        Updater(const Updater&)             = default;
        Updater(Updater&&)                  = default;
        Updater& operator=(const Updater&)  = default;
        Updater& operator=(Updater&&)       = default;

    };

    /** Single-objective population update method types. */
    template<typename T>
    concept UpdaterType = requires
    {
        requires std::derived_from<T, Updater>;
    };

} // namespace genetic_algorithm::update

#endif // !GA_ALGORITHM_SOGA_UPDATER_BASE_HPP