/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_ALGORITHM_FWD_HPP
#define GA_ALGORITHM_ALGORITHM_FWD_HPP

#include <concepts>
#include <type_traits>
#include <cstddef>
#include "../population/population.hpp"

namespace genetic_algorithm
{
    class GaInfo;
}

namespace genetic_algorithm::algorithm
{
    class Algorithm;

    template<typename T>
    concept AlgorithmType = requires
    {
        requires std::derived_from<T, Algorithm>;
    };

} // namespace genetic_algorithm::algorithm

namespace genetic_algorithm::selection_
{
    using detail::FitnessVector;
    using detail::FitnessMatrix;

    /**
    * Concept specifying the interface required for the selection methods. \n
    *
    * ...
    */
    template<typename T>
    concept Selection = requires(T selection, const GaInfo & ga, FitnessMatrix fmat)
    {
        requires std::copyable<std::decay_t<T>>;
        requires std::destructible<std::decay_t<T>>;

        { selection.initialize(ga) }              -> std::same_as<void>;
        { selection.prepareSelections(ga, fmat) } -> std::same_as<void>;
        { selection.select(ga, fmat) }            -> std::same_as<size_t>;
    };

} // namespace genetic_algorithm::selection_

namespace genetic_algorithm::pop_update
{
    using detail::FitnessVector;
    using detail::FitnessMatrix;

    /**
     * Concept specifying the interface required for the population update methods. \n
     *
     * The method should be callable with a const& to a GaInfo object, and the FitnessMatrix
     * of the combined parent and child populations specified by 3 iterators:
     * [parents_first, ... , children_first, ... , children_last),
     * and return the indices of the candidates selected for the next generation's population,
     * assuming that parents_first points to the idx 0.
     */
    template<typename T>
    concept Updater = requires(T updater, const GaInfo & ga, FitnessMatrix::const_iterator it)
    {
        requires std::copyable<std::decay_t<T>>;
        requires std::destructible<std::decay_t<T>>;

        { updater(ga, it, it, it) } -> std::same_as<std::vector<size_t>>;
    };

} // namespace genetic_algorithm::pop_update

#endif // !GA_ALGORITHM_ALGORITHM_FWD_HPP