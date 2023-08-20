/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "nd_sort.hpp"
#include "../population/population.hpp"
#include "../utility/algorithm.hpp"
#include "../utility/functional.hpp"
#include "../utility/iterators.hpp"
#include "../utility/math.hpp"
#include "../utility/utility.hpp"
#include "../utility/matrix.hpp"
#include <algorithm>
#include <functional>
#include <execution>
#include <iterator>
#include <atomic>
#include <vector>
#include <utility>
#include <cstddef>

namespace gapp::algorithm::dtl
{
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
        GAPP_ASSERT(std::distance(current, last) >= 0);

        return std::ranges::find_if(current, last, detail::not_equal_to(current->rank), &FrontInfo::rank);
    }

    std::vector<ParetoFrontsRange> paretoFrontBounds(ParetoFronts& pareto_fronts)
    {
        using Iter = ParetoFronts::iterator;
        using IterPair = std::pair<Iter, Iter>;

        if (pareto_fronts.empty()) return {};   /* No bounds exist. */

        std::vector<IterPair> front_bounds;
        front_bounds.reserve(pareto_fronts.back().rank + 1);

        for (auto front_first = pareto_fronts.begin(); front_first != pareto_fronts.end(); )
        {
            auto front_last = nextFrontBegin(front_first, pareto_fronts.end());
            front_bounds.emplace_back(front_first, front_last);
            front_first = front_last;
        }

        return front_bounds;
    }

    ParetoFrontsRange findPartialFront(ParetoFronts::iterator first, ParetoFronts::iterator last, size_t popsize)
    {
        GAPP_ASSERT(0 < popsize && popsize <= size_t(last - first));

        auto last_in = first + popsize;

        if (last_in == last) return { last, last };
        auto partial_front_first = std::ranges::find_if(first, last_in, detail::equal_to(last_in->rank), &FrontInfo::rank);
        auto partial_front_last = nextFrontBegin(std::prev(last_in), last);

        return { partial_front_first, partial_front_last };
    }


    /* Fast non-dominated sorting algorithm (FNDS)
    * 
    * See:
    *  Deb, K., et al. "A fast and elitist multiobjective genetic algorithm: NSGA-II."
    *  IEEE transactions on evolutionary computation 6, no. 2 (2002): 182-197.
    */

    struct DominanceList
    {
        std::vector<size_t> worse_indices;  /* Indices of the solutions dominated by this one. */
        size_t better_count = 0;            /* Number of solutions dominating this one. */
    };
    using DominanceLists = std::vector<DominanceList>;

    static DominanceLists& constructDominanceLists(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        const size_t popsize = std::distance(first, last);

        thread_local DominanceLists dom_lists;

        if (dom_lists.size() != popsize)
        {
            dom_lists.resize(popsize);
            for (auto& list : dom_lists) detail::clear_reserve(list.worse_indices, popsize);
        }

        for (auto& list : dom_lists)
        {
            list.worse_indices.clear();
            list.better_count = 0;
        }

        /* Find the number of candidates which dominate each candidate, and the indices of the candidates it dominates. */
        for (size_t lhs = 0; lhs < popsize; lhs++)
        {
            for (size_t rhs = 0; rhs < lhs; rhs++)
            {
                auto comp = math::paretoCompare(first[lhs], first[rhs]);
                if (comp < 0)
                {
                    dom_lists[lhs].better_count++;
                    dom_lists[rhs].worse_indices.push_back(lhs);
                }
                else if (comp > 0)
                {
                    dom_lists[rhs].better_count++;
                    dom_lists[lhs].worse_indices.push_back(rhs);
                }
            }
        }

        return dom_lists;
    }

    ParetoFronts fastNonDominatedSort(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        const size_t popsize = std::distance(first, last);
        auto& dom_lists = constructDominanceLists(first, last);

        ParetoFronts pfronts;
        pfronts.reserve(popsize);

        /* Find the indices of all non-dominated candidates (the first/best pareto front). */
        for (size_t i = 0; i < popsize; i++)
        {
            if (dom_lists[i].better_count == 0) pfronts.push_back(FrontInfo{ i, 0 });
        }

        /* Find all the other pareto fronts. */
        auto front_first = detail::stable_begin(pfronts);
        auto front_last  = detail::stable_end(pfronts);
        
        while (pfronts.size() != popsize)
        {
            /* Remove the current front from the population and find the next one. */
            const size_t next_front_rank = front_first->rank + 1;

            for (; front_first != front_last; ++front_first)
            {
                const size_t sol_idx = front_first->idx;

                for (size_t worse_idx : dom_lists[sol_idx].worse_indices)
                {
                    if (--dom_lists[worse_idx].better_count == 0)
                    {
                        pfronts.push_back(FrontInfo{ worse_idx, next_front_rank });
                    }
                }
            }
            front_last = detail::stable_end(pfronts);
        }

        return pfronts;
    }


    /* Dominance-degree sorting algorithm (DDS).
     * 
     * See:
     *  Zhou, Yuren, Zefeng Chen, and Jun Zhang. "Ranking vectors by means of the dominance degree matrix."
     *  IEEE Transactions on Evolutionary Computation 21, no. 1 (2016): 34-51.
     * 
     *  Mishra, S., et al. "Time complexity analysis of the dominance degree approach for non-dominated sorting."
     *  Proceedings of the 2020 Genetic and Evolutionary Computation Conference Companion, pp. 169-170. 2020.
     */

    using DominanceMatrix = detail::Matrix<std::atomic_bool>;
    enum : bool { MAXIMAL = false, NONMAXIMAL = true };

    struct Col { size_t idx; size_t sum; }; /* Index and colwise sum of a column in the dominance matrix. */

    static DominanceMatrix constructDominanceMatrix(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        const size_t popsize = std::distance(first, last);
        DominanceMatrix dmat(popsize, popsize /*, MAXIMAL */);

        std::for_each(GAPP_EXEC_UNSEQ, detail::iota_iterator(0_sz), detail::iota_iterator(first->size()), [&](size_t obj)
        {
            FitnessVector fvec(popsize);
            std::transform(first, last, fvec.begin(), detail::element_at(obj));

            const auto indices = detail::argsort(fvec.rbegin(), fvec.rend()); // descending

            std::for_each(indices.rbegin(), indices.rend(), [&, current = indices.rbegin()](size_t row) mutable
            {
                const auto last_idx = std::find_if(++current, indices.rend(), [&](size_t prev_row) noexcept
                {
                    return !math::floatIsEqual(fvec[prev_row], fvec[row]);
                });

                std::for_each(last_idx, indices.rend(), [&](size_t col) noexcept
                {
                    if (dmat(row, col).load(std::memory_order_relaxed) == MAXIMAL) dmat(row, col).store(NONMAXIMAL, std::memory_order_relaxed);
                });
            });
        });

        std::for_each(GAPP_EXEC_UNSEQ, detail::iota_iterator(0_sz), detail::iota_iterator(popsize), [&](size_t row) noexcept
        {
            dmat(row, row).store(NONMAXIMAL, std::memory_order_relaxed); // diagonal is all nonmax

            for (size_t col = row + 1; col < popsize; col++)
            {
                if (dmat(row, col).load(std::memory_order_relaxed) == MAXIMAL && dmat(col, row).load(std::memory_order_relaxed) == MAXIMAL)
                {
                    dmat(row, col).store(NONMAXIMAL, std::memory_order_relaxed);
                    dmat(col, row).store(NONMAXIMAL, std::memory_order_relaxed);
                }
            }
        });

        return dmat;
    }

    static std::vector<Col> colwiseSums(const DominanceMatrix& dmat)
    {
        std::vector<Col> cols(dmat.ncols());
        for (size_t i = 0; i < cols.size(); i++)
        {
            cols[i].idx = i;
            cols[i].sum = dmat.nrows();
        }

        for (size_t row = 0; row < dmat.nrows(); row++)
        {
            for (size_t col = 0; col < dmat.ncols(); col++)
            {
                cols[col].sum -= dmat(row, col).load(std::memory_order_relaxed);
            }
        }

        return cols;
    }

    ParetoFronts dominanceDegreeSort(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        const size_t popsize = std::distance(first, last);
        DominanceMatrix dmat = constructDominanceMatrix(first, last);
        std::vector<Col> cols = colwiseSums(dmat);

        /* Assign each solution to a pareto-front. */
        ParetoFronts pfronts;
        pfronts.reserve(popsize);

        size_t current_rank = 0;
        std::vector<size_t> removed_rows;

        while (pfronts.size() != popsize)
        {
            removed_rows.clear();

            /* Remove cols */
            for (auto& col : cols)
            {
                if (col.sum == 0)
                {
                    pfronts.push_back(FrontInfo{ col.idx, current_rank });
                    removed_rows.push_back(col.idx);
                }
            }
            std::erase_if(cols, [](const Col& col) { return col.sum == 0; });

            /* Remove rows */
            for (size_t row : removed_rows)
            {
                for (auto& [col, sum] : cols)
                {
                    if (dmat(row, col).load(std::memory_order_relaxed) == MAXIMAL)
                    {
                        dmat(row, col).store(NONMAXIMAL, std::memory_order_relaxed);
                        sum--;
                    }
                }
            }

            current_rank++;
        }

        return pfronts;
    }

    ParetoFronts nonDominatedSort(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        GAPP_ASSERT(std::distance(first, last) >= 0);

        return dominanceDegreeSort(first, last);
    }

} // namespace gapp::algorithm::dtl