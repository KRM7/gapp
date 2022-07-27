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

        size_t pop_size = size_t(last - first);

        /* Find the number of candidates which dominate each candidate, and the indices of the candidates it dominates. */
        std::vector<size_t> better_count(pop_size, 0);

        static std::vector<std::vector<size_t>> worse_indices;
        for (auto& vec : worse_indices) vec.clear();
        if (worse_indices.size() != pop_size)
        {
            worse_indices.resize(pop_size);
            for (auto& vec : worse_indices) vec.shrink_to_fit(), vec.reserve(pop_size);
        }

        for (size_t lhs = 0; lhs < pop_size; lhs++)
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

        /* [idx, rank] */
        ParetoFronts sorted_indices;
        sorted_indices.reserve(pop_size);

        /* Find the indices of all non-dominated candidates (the first/best pareto front). */
        for (size_t i = 0; i < pop_size; i++)
        {
            if (better_count[i] == 0)
            {
                sorted_indices.emplace_back(i, 0);
            }
        }

        /* Find all the other pareto fronts. */
        auto front_first = sorted_indices.cbegin();
        auto front_last  = sorted_indices.cend();

        while (sorted_indices.size() != pop_size)
        {
            size_t next_rank = front_first->second + 1;

            /* Remove the current front from the population and find the next one. */
            for (; front_first != front_last; ++front_first)
            {
                for (size_t worse_idx : worse_indices[front_first->first])
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

    std::vector<size_t> paretoRanks(const ParetoFronts& pareto_fronts)
    {
        std::vector<size_t> ranks(pareto_fronts.size());

        for (const auto& [idx, rank] : pareto_fronts)
        {
            ranks[idx] = rank;
        }

        return ranks;
    }

    ParetoFronts::iterator nextFrontBegin(ParetoFronts::iterator current, ParetoFronts::iterator last) noexcept
    {
        if (current == last) return last;    /* There isn't a next front. */

        return std::find_if(current, last,
        [current_rank = current->second](const std::pair<size_t, size_t>& elem) noexcept
        {
            size_t rank = elem.second;
            return rank != current_rank;
        });
    }

    auto paretoFrontBounds(ParetoFronts& pareto_fronts) -> std::vector<std::pair<ParetoFronts::iterator, ParetoFronts::iterator>>
    {
        if (pareto_fronts.empty()) return {};   /* No bounds exist. */

        using Iter = ParetoFronts::iterator;
        using IterPair = std::pair<Iter, Iter>;

        std::vector<IterPair> front_bounds;
        front_bounds.reserve(pareto_fronts.back().second);

        for (auto first = pareto_fronts.begin(); first != pareto_fronts.end(); )
        {
            auto last = nextFrontBegin(first, pareto_fronts.end());
            front_bounds.emplace_back(first, last);
            first = last;
        }

        return front_bounds;
    }

} // namespace genetic_algorithm::algorithm::dtl