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

    struct FrontInfo
    {
        size_t idx;  // The solution's idx in the fitness matrix
        size_t rank; // The rank of the pareto front the solution belongs to

        friend bool operator==(const FrontInfo&, const FrontInfo&) = default;
    };

    using ParetoFronts      = std::vector<FrontInfo>;
    using ParetoFrontsRange = std::pair<ParetoFronts::iterator, ParetoFronts::iterator>;

    ParetoFronts fastNonDominatedSort(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);

    ParetoFronts dominanceDegreeSort(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);

    /* Non-dominated sorting for the multi-objective algorithms.
       Returns the pareto fronts (in the form of [idx, rank] pairs) of the population, in non-descending order based on the ranks. */
    ParetoFronts nonDominatedSort(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);

    /* Returns the pareto rank of each candidate in the pareto fronts (best rank is 0). */
    std::vector<size_t> paretoRanks(const ParetoFronts& pareto_fronts);

    /* The functions below assume that their ParetoFronts parameters are in the format returned by
       nonDominatedSort (sorted in non-descending order based on the ranks of the elements). */

    /* Find the first element of the front following the front which @p current is a part of. */
    ParetoFronts::iterator nextFrontBegin(ParetoFronts::iterator current, ParetoFronts::iterator last) noexcept;

    /* Find the first and last elements of each front in the pareto fronts. */
    std::vector<ParetoFrontsRange> paretoFrontBounds(ParetoFronts& pareto_fronts);
    
    /* Find the pareto front with the lowest rank that can't be added to the next population in its entirety. */
    ParetoFrontsRange findPartialFront(ParetoFronts::iterator first, ParetoFronts::iterator last, size_t popsize);

} // namespace genetic_algorithm::algorithm::dtl

#endif // !GA_ALGORITHM_NDSORT_HPP