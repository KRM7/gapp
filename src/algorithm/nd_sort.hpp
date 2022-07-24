/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_NDSORT_HPP
#define GA_ALGORITHM_NDSORT_HPP

#include "../population/population.hpp"
#include <vector>
#include <utility>
#include <cstddef>

namespace genetic_algorithm::algorithm::dtl
{
    using detail::FitnessVector;
    using detail::FitnessMatrix;

    /* [solution idx, solution rank] pairs. */
    using ParetoFronts = std::vector<std::pair<size_t, size_t>>;

    /* Non-dominated sorting for the multi-objective algorithms.
       Returns the pareto fronts (in the form of [idx, rank] pairs) of the population, in non-descending order based on the ranks. */
    ParetoFronts nonDominatedSort(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);

    /* Returns the pareto rank of each candidate based on the pareto fronts (best rank is 0). */
    std::vector<size_t> paretoRanks(const ParetoFronts& pareto_fronts);

    /* Finds the first element of the front following the front which @p current is a part of.
       (Assumes that the range is sorted in non-descending order based on the ranks of the elements.) */
    ParetoFronts::iterator nextFrontBegin(ParetoFronts::iterator current, ParetoFronts::iterator last) noexcept;

    /* Finds the first and last elements of each front in the pareto fronts vector.
       (Assumes that the range is sorted in non-descending order based on the ranks of the elements.) */
    auto paretoFrontBounds(ParetoFronts& pareto_fronts) -> std::vector<std::pair<ParetoFronts::iterator, ParetoFronts::iterator>>;

} // namespace genetic_algorithm::algorithm::dtl

#endif // !GA_ALGORITHM_NDSORT_HPP