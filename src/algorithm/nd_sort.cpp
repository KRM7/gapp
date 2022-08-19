/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "nd_sort.hpp"
#include "../utility/algorithm.hpp"
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
        
        for (const auto& [idx, rank] : pfronts) ranks[idx] = rank;

        return ranks;
    }

    ParetoFronts::iterator nextFrontBegin(ParetoFronts::iterator current, ParetoFronts::iterator last) noexcept
    {
        return std::find_if(current, last,
        [current](const FrontInfo& sol) noexcept
        {
            return sol.rank != current->rank;
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

        if (size_t(last - first) == popsize) return { last, last };

        auto partial_front_first = std::find_if(first, first + popsize, [&](const FrontInfo& sol)
        {
            return sol.rank == first[popsize].rank;
        });
        auto partial_front_last = (partial_front_first == first + popsize) ? partial_front_first : nextFrontBegin(partial_front_first, last);

        return { partial_front_first, partial_front_last };
    }


    /* Fast non-dominated sorting algorithm (FNDS) */

    static ParetoFronts fastNonDominatedSort(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        assert(std::all_of(first, last, [first](const FitnessVector& f) { return f.size() == first->size(); }));

        size_t popsize = size_t(last - first);

        if (popsize == 0) return {};

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

        ParetoFronts pfronts;
        pfronts.reserve(popsize);

        /* Find the indices of all non-dominated candidates (the first/best pareto front). */
        for (size_t i = 0; i < popsize; i++)
        {
            if (better_count[i] == 0) pfronts.emplace_back(i, 0);
        }

        /* Find all the other pareto fronts. */
        auto front_first = pfronts.cbegin();
        auto front_last  = pfronts.cend();

        while (pfronts.size() != popsize)
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
                        pfronts.emplace_back(worse_idx, next_rank);
                        front_last++;
                    }
                }
            }
            front_last = pfronts.cend();
        }

        return pfronts;
    }


    /* Dominance-degree sorting algorithm (DDS). */

    using DominanceMatrix = detail::Matrix<unsigned char>;
    enum : unsigned char { MAXIMAL = true, NONMAXIMAL = false };

    struct Col { size_t idx; size_t sum; };

    static DominanceMatrix constructDominanceMatrix(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        assert(std::distance(first, last) >= 0);

        const size_t popsize = std::distance(first, last);
        DominanceMatrix dmat(popsize, popsize, MAXIMAL);

        if (popsize == 0) return dmat;

        std::vector<size_t> objectives(first->size());
        std::iota(objectives.begin(), objectives.end(), 0);

        std::for_each(GA_EXECUTION_UNSEQ, objectives.cbegin(), objectives.cend(), [&](size_t obj)
        {
            FitnessVector fvec(popsize);
            std::transform(first, last, fvec.begin(), [obj](const FitnessVector& row) { return row[obj]; });

            const auto indices = detail::argsort(fvec.rbegin(), fvec.rend()); // descending

            std::for_each(indices.rbegin(), indices.rend(), [&, current = indices.rbegin()](size_t row) mutable
            {
                auto last_col = std::find_if(++current, indices.rend(), [&](size_t prev_row) noexcept
                {
                    return !detail::floatIsEqual(fvec[prev_row], fvec[row]);
                });

                std::for_each(last_col, indices.rend(), [&](size_t col) noexcept
                {
                    std::atomic_ref dmat_entry{ dmat(row, col) };
                    dmat_entry.store(NONMAXIMAL, std::memory_order_relaxed);
                });
            });
        });

        std::vector<size_t> indices(popsize);
        std::iota(indices.begin(), indices.end(), 0);

        std::for_each(indices.begin(), indices.end(), [&](size_t row) noexcept
        {
            std::for_each(indices.begin() + row, indices.end(), [&](size_t col) noexcept
            {
                if (dmat(row, col) == MAXIMAL && dmat(col, row) == MAXIMAL)
                {
                    dmat(row, col) = NONMAXIMAL;
                    dmat(col, row) = NONMAXIMAL;
                }
            });
        }); // 2.2

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
        assert(std::distance(first, last) >= 0);

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

    ParetoFronts defGoodSort(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        std::vector<size_t> indices(last - first);
        std::iota(indices.begin(), indices.end(), 0);

        std::vector<size_t> this_front;

        ParetoFronts pfronts;
        
        size_t front = 0;
        while (pfronts.size() != size_t(last - first))
        {
            for (auto idx : indices)
            {
                if (std::none_of(indices.begin(), indices.end(), [&](size_t idx2) { return detail::paretoCompareLess(first[idx], first[idx2]); }))
                {
                    pfronts.emplace_back(idx, front);
                    this_front.push_back(idx);
                }
            }
            indices.erase(std::remove_if(indices.begin(), indices.end(), [&](size_t idx) { return detail::contains(this_front.begin(), this_front.end(), idx); }), indices.end());
            front++;
        }

        return pfronts;
    }

    ParetoFronts nonDominatedSort(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        //auto good = defGoodSort(first, last);
        auto test = dominanceDegreeSort(first, last);
        //assert(std::is_permutation(good.begin(), good.end(), test.begin(), test.end(), [](auto& lhs, auto& rhs) { return (lhs.idx == rhs.idx && lhs.rank == rhs.rank); }));
        return test;
    }

} // namespace genetic_algorithm::algorithm::dtl