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

    /* Sorted (idx, rank) pairs. */
    using ParetoFronts = std::vector<std::pair<size_t, size_t>>;

    /* Non-dominated sorting for the multi-objective algorithms. Returns the pareto fronts (idx, rank pairs) of the population. */
    ParetoFronts nonDominatedSort(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);

    /* Returns the rank of each candidate based on the pareto fronts. */
    std::vector<size_t> paretoRanks(const ParetoFronts& pareto_fronts);

    /* Finds the first element of the front following the front which current is a part of. */
    ParetoFronts::iterator nextFrontBegin(ParetoFronts::iterator current, ParetoFronts::iterator last) noexcept;

    /* Finds the first and last elements of each front in the pareto fronts vector. */
    auto paretoFrontBounds(ParetoFronts& pareto_fronts) -> std::vector<std::pair<ParetoFronts::iterator, ParetoFronts::iterator>>;

} // namespace genetic_algorithm::algorithm::dtl

#endif // !GA_ALGORITHM_NDSORT_HPP