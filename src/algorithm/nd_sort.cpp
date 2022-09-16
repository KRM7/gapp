/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "nd_sort.hpp"
#include "../utility/algorithm.hpp"
#include "../utility/iterators.hpp"
#include "../utility/math.hpp"
#include "../utility/utility.hpp"
#include "../utility/matrix.hpp"
#include <algorithm>
#include <iterator>
#include <numeric>
#include <atomic>
#include <vector>
#include <utility>
#include <cstddef>

namespace genetic_algorithm::algorithm::dtl
{
    std::vector<size_t> paretoRanks(const ParetoFronts& pfronts)
    {
        std::vector<size_t> ranks(pfronts.size());     
        for (const auto& [idx, rank] : pfronts)
        {
            ranks[idx] = rank;
        }

        return ranks;
    }

    ParetoFronts::iterator nextFrontBegin(ParetoFronts::iterator current, ParetoFronts::iterator last) noexcept
    {
        return std::find_if(current, last, [&](const FrontInfo& sol) noexcept
        {
            return sol.rank != current->rank;
        });
    }

    std::vector<ParetoFrontsRange> paretoFrontBounds(ParetoFronts& pareto_fronts)
    {
        using Iter = ParetoFronts::iterator;
        using IterPair = std::pair<Iter, Iter>;

        if (pareto_fronts.empty()) return {};   /* No bounds exist. */

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

        if (size_t(last - first) == popsize) return { last, last };

        auto partial_front_first = std::find_if(first, first + popsize, [&](const FrontInfo& sol)
        {
            return sol.rank == first[popsize].rank;
        });
        auto partial_front_last = (partial_front_first == first + popsize) ? partial_front_first : nextFrontBegin(partial_front_first, last);

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
        size_t better_count = 0;            /* Number of solutions dominating this one. */
        std::vector<size_t> worse_indices;  /* Indices of the solutions dominated by this one. */
    };
    using DominanceLists = std::vector<DominanceList>;

    static DominanceLists& getDominanceLists(size_t popsize)
    {
        static DominanceLists dlists;

        /* Resize if popsize changes */
        if (dlists.size() != popsize)
        {
            dlists.resize(popsize);
            for (auto& dlist : dlists)
            {
                dlist.worse_indices.shrink_to_fit();
                dlist.worse_indices.reserve(popsize);
            }
        }

        /* Clear */
        for (auto& dlist : dlists)
        {
            dlist.worse_indices.clear();
            dlist.better_count = 0;
        }

        return dlists;
    }

    static DominanceLists& constructDominanceLists(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        const size_t popsize = std::distance(first, last);
        auto& dom_lists = getDominanceLists(popsize);

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

    static ParetoFronts fastNonDominatedSort(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        const size_t popsize = std::distance(first, last);
        auto& dom_lists = constructDominanceLists(first, last);

        ParetoFronts pfronts;
        pfronts.reserve(popsize);

        /* Find the indices of all non-dominated candidates (the first/best pareto front). */
        for (size_t i = 0; i < popsize; i++)
        {
            if (dom_lists[i].better_count == 0) pfronts.emplace_back(i, 0);
        }

        /* Find all the other pareto fronts. */
        auto current_front_first = detail::stable_begin(pfronts);
        auto current_front_last = detail::stable_end(pfronts);

        while (pfronts.size() != popsize)
        {
            /* Remove the current front from the population and find the next one. */
            const size_t next_front_rank = current_front_first->rank + 1;

            for (; current_front_first != current_front_last; ++current_front_first)
            {
                const size_t sol_idx = current_front_first->idx;

                for (size_t worse_idx : dom_lists[sol_idx].worse_indices)
                {
                    if (--dom_lists[worse_idx].better_count == 0)
                    {
                        pfronts.emplace_back(worse_idx, next_front_rank);
                    }
                }
            }
            current_front_last = detail::stable_end(pfronts);
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

    using DominanceMatrix = detail::Matrix<unsigned char>;
    enum : unsigned char { MAXIMAL = true, NONMAXIMAL = false };

    struct Col { size_t idx; size_t sum; }; /* Index and colwise sum of a column in the dominance matrix. */

    static DominanceMatrix constructDominanceMatrix(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        const size_t popsize = std::distance(first, last);
        DominanceMatrix dmat(popsize, popsize, MAXIMAL);

        if (popsize == 0) return dmat;

        const auto objectives = detail::index_vector(first->size());

        std::for_each(GA_EXECUTION_UNSEQ, objectives.cbegin(), objectives.cend(), [&](size_t obj)
        {
            FitnessVector fvec(popsize);
            std::transform(first, last, fvec.begin(), [obj](const FitnessVector& row) { return row[obj]; });

            const auto indices = detail::argsort(fvec.rbegin(), fvec.rend()); // descending

            std::for_each(indices.rbegin(), indices.rend(), [&, current = indices.rbegin()](size_t row) mutable
            {
                auto last_col = std::find_if(++current, indices.rend(), [&](size_t prev_row) noexcept
                {
                    return !math::floatIsEqual(fvec[prev_row], fvec[row]);
                });

                std::for_each(last_col, indices.rend(), [&](size_t col) noexcept
                {
                    std::atomic_ref dmat_entry{ dmat(row, col) };
                    dmat_entry.store(NONMAXIMAL, std::memory_order_relaxed);
                });
            });
        });

        const auto indices = detail::index_vector(popsize);

        std::for_each(GA_EXECUTION_UNSEQ, indices.begin(), indices.end(), [&](size_t row) noexcept
        {
            dmat(row, row) = NONMAXIMAL;
            std::for_each(indices.begin() + row + 1, indices.end(), [&](size_t col) noexcept
            {
                if (dmat(row, col) == MAXIMAL && dmat(col, row) == MAXIMAL)
                {
                    dmat(row, col) = NONMAXIMAL;
                    dmat(col, row) = NONMAXIMAL;
                }
            });
        });

        return dmat;
    }

    static std::vector<Col> colwiseSums(const DominanceMatrix& dmat)
    {
        std::vector<Col> cols(dmat.ncols());
        for (size_t i = 0; i < cols.size(); i++) cols[i].idx = i;

        for (size_t row = 0; row < dmat.nrows(); row++)
            for (size_t col = 0; col < dmat.ncols(); col++)
                cols[col].sum += dmat(row, col);

        return cols;
    }

    static ParetoFronts dominanceDegreeSort(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        const size_t popsize = std::distance(first, last);
        DominanceMatrix dmat = constructDominanceMatrix(first, last);
        std::vector<Col> cols = colwiseSums(dmat);

        /* Assign each solution to a pareto-front. */
        ParetoFronts pfronts;
        pfronts.reserve(popsize);

        size_t current_rank = 0;
        while (pfronts.size() != popsize)
        {
            std::vector<size_t> removed_rows;

            /* Remove cols */
            for (auto& col : cols)
            {
                if (col.sum == 0)
                {
                    pfronts.emplace_back(col.idx, current_rank);
                    removed_rows.push_back(col.idx);
                }
            }
            std::erase_if(cols, [](const Col& col) { return col.sum == 0; });

            /* Remove rows */
            std::for_each(removed_rows.begin(), removed_rows.end(), [&](size_t row)
            {
                for (auto& col : cols)
                {
                    if (dmat(row, col.idx) == MAXIMAL)
                    {
                        dmat(row, col.idx) = NONMAXIMAL;
                        col.sum--;
                    }
                }
            });

            current_rank++;
        }

        return pfronts;
    }

    ParetoFronts nonDominatedSort(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        assert(std::distance(first, last) >= 0);
        assert(std::all_of(first, last, [first](const FitnessVector& f) { return f.size() == first->size(); }));

        return dominanceDegreeSort(first, last);
    }

} // namespace genetic_algorithm::algorithm::dtl