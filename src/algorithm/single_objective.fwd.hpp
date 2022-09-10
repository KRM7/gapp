/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_SINGLE_OBJECTIVE_FWD_HPP
#define GA_ALGORITHM_SINGLE_OBJECTIVE_FWD_HPP

#include "../population/population.hpp"
#include <vector>
#include <concepts>
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

    /** Concept specifying the interface required for the single-objective selection methods. */
    template<typename T>
    concept SelectionType = requires(T selection, const GaInfo& ga, FitnessMatrix fmat)
    {
        requires !std::is_reference_v<T>;
        requires std::copyable<T>;
        requires std::destructible<T>;

        { selection.initialize(ga) }              -> std::same_as<void>;
        { selection.prepareSelections(ga, fmat) } -> std::same_as<void>;
        { selection.select(ga, fmat) }            -> std::same_as<size_t>;
    };

} // namespace genetic_algorithm::selection

namespace genetic_algorithm::update
{
    using detail::FitnessVector;
    using detail::FitnessMatrix;

    /**
     * Concept specifying the interface required for the population update methods. \n
     *
     * The method should be callable with a const& to a GaInfo object, and the FitnessMatrix
     * of the combined parent and child populations specified by 3 iterators: \n
     * [parents_first, ... , children_first, ... , children_last), \n
     * and return the indices of the candidates selected for the next generation's population,
     * assuming that parents_first points to the idx 0.
     */
    template<typename T>
    concept UpdaterType = requires(T updater, const GaInfo& ga, FitnessMatrix::const_iterator it)
    {
        requires !std::is_reference_v<T>;
        requires std::copyable<T>;
        requires std::destructible<T>;

        { updater(ga, it, it, it) } -> std::same_as<std::vector<size_t>>;
    };

} // namespace genetic_algorithm::update

namespace genetic_algorithm::algorithm
{
    template<selection::SelectionType Selection, update::UpdaterType Updater>
    class SingleObjective;

} // namespace genetic_algorithm::algorithm

#endif // !GA_ALGORITHM_SINGLE_OBJECTIVE_FWD_HPP