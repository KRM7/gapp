/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "nd_sort.hpp"
#include "../utility/math.hpp"
#include <algorithm>
#include <iterator>
#include <vector>
#include <utility>
#include <cstddef>

namespace genetic_algorithm::algorithm::dtl
{
    ParetoFronts nonDominatedSort(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        assert(std::distance(first, last) >= 0);
        assert(std::all_of(first, last, [first](const FitnessVector& f) { return f.size() == first->size(); }));

        size_t popsize = size_t(last - first);

        /* Find the number of candidates which dominate each candidate, and the indices of the candidates it dominates. */
        std::vector<size_t> better_count(popsize, 0);

        static std::vector<std::vector<size_t>> worse_indices;
        for (auto& vec : worse_indices) vec.clear();
        if (worse_indices.size() != popsize)
        {
            worse_indices.resize(popsize);
            for (auto& vec : worse_indices) vec.shrink_to_fit(), vec.reserve(popsize);
        }

        for (size_t lhs = 0; lhs < popsize; lhs++)
        {
            for (size_t rhs = 1; rhs < lhs; rhs++)
            {
                auto comp = detail::paretoCompare(first[lhs], first[rhs]);
                if (comp < 0)
                {
                    better_count[lhs]++;
                    worse_indices[rhs].push_back(lhs);
                }
                else if (comp > 0)
                {
                    better_count[rhs]++;
                    worse_indices[lhs].push_back(rhs);
                }
            }
        }

        ParetoFronts sorted_indices;
        sorted_indices.reserve(popsize);

        /* Find the indices of all non-dominated candidates (the first/best pareto front). */
        for (size_t i = 0; i < popsize; i++)
        {
            if (better_count[i] == 0) sorted_indices.emplace_back(i, 0);
        }

        /* Find all the other pareto fronts. */
        auto front_first = sorted_indices.cbegin();
        auto front_last  = sorted_indices.cend();

        while (sorted_indices.size() != popsize)
        {
            size_t next_rank = front_first->rank + 1;

            /* Remove the current front from the population and find the next one. */
            for (; front_first != front_last; ++front_first)
            {
                for (size_t worse_idx : worse_indices[front_first->idx])
                {
                    if (--better_count[worse_idx] == 0)
                    {
                        front_last--;
                        sorted_indices.emplace_back(worse_idx, next_rank);
                        front_last++;
                    }
                }
            }
            front_last = sorted_indices.cend();
        }

        return sorted_indices;
    }

    std::vector<size_t> paretoRanks(const ParetoFronts& pfronts)
    {
        std::vector<size_t> ranks(pfronts.size());
        
        for (const auto& [idx, rank] : pfronts) ranks[idx] = rank;

        return ranks;
    }

    ParetoFronts::iterator nextFrontBegin(ParetoFronts::iterator current, ParetoFronts::iterator last) noexcept
    {
        return std::find_if(current, last,
        [current](const FrontInfo& elem) noexcept
        {
            return elem.rank != current->rank;
        });
    }

    std::vector<ParetoFrontsRange> paretoFrontBounds(ParetoFronts& pareto_fronts)
    {
        if (pareto_fronts.empty()) return {};   /* No bounds exist. */

        using Iter = ParetoFronts::iterator;
        using IterPair = std::pair<Iter, Iter>;

        std::vector<IterPair> front_bounds;
        front_bounds.reserve(pareto_fronts.back().rank + 1);

        for (auto first = pareto_fronts.begin(); first != pareto_fronts.end(); )
        {
            auto last = nextFrontBegin(first, pareto_fronts.end());
            front_bounds.emplace_back(first, last);
            first = last;
        }

        return front_bounds;
    }

    ParetoFrontsRange findPartialFront(ParetoFronts::iterator first, ParetoFronts::iterator last, size_t popsize)
    {
        assert(size_t(last - first) >= popsize);

        auto partial_front_first = first + popsize;
        auto partial_front_last  = first + popsize;

        bool has_partial_front = (size_t(last - first) > popsize) && 
                                 (first[popsize - 1].rank == first[popsize].rank);

        if (has_partial_front)
        {
            size_t partial_front_rank = first[popsize - 1].rank;

            partial_front_first = std::find_if(first, last,
            [&](const FrontInfo& sol) noexcept
            {
                return sol.rank == partial_front_rank;
            });

            partial_front_last = nextFrontBegin(partial_front_first, last);
        }

        return { partial_front_first, partial_front_last };
    }

} // namespace genetic_algorithm::algorithm::dtl