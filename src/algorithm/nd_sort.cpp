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

    //using DominanceMatrix = std::vector<std::vector<std::atomic_bool>>;
    using DominanceMatrix = detail::Matrix<std::atomic_bool>;

    static DominanceMatrix& getDominanceMatrix(size_t size)
    {
        static DominanceMatrix dmat;

        if (dmat.nrows() != size)
        {
            dmat = DominanceMatrix(size, size);
        }

        for (auto& flag : dmat)
        {
            flag.store(true, std::memory_order_relaxed);
        }
        /* If we get rid of the static matrix, we could use false as the default value, and set to true when creating the matrix */
        /* When counting, use max - sum, etc... */

        return dmat;
    }

    static DominanceMatrix& createDominanceMatrix(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        assert(std::distance(first, last) >= 0);

        const size_t popsize = std::distance(first, last);
        auto& dom_mat = getDominanceMatrix(popsize); // 1.16

        if (popsize == 0) return dom_mat;

        std::vector<size_t> objs(first->size());
        std::iota(objs.begin(), objs.end(), 0);

        std::for_each(GA_EXECUTION_UNSEQ, objs.cbegin(), objs.cend(), [&](size_t obj)
        {
            FitnessVector fvec(popsize);
            std::transform(first, last, fvec.begin(), [obj](const FitnessVector& row) { return row[obj]; });

            const auto indices = detail::argsort(fvec.rbegin(), fvec.rend()); // ascending

            std::for_each(indices.rbegin(), indices.rend(), [&, current = indices.rbegin()](size_t row) mutable
            {
                auto last_col = std::find_if(++current, indices.rend(), [&](size_t prev_row)
                {
                    return !detail::floatIsEqual(fvec[prev_row], fvec[row]);
                });

                std::for_each(last_col, indices.rend(), [&](size_t col)
                { 
                    //dom_mat[row][col].store(false, std::memory_order_relaxed);
                    dom_mat(row, col).store(false, std::memory_order_relaxed);
                });
            });

            /*
            * The same strategy can be used to construct the columnwise counts:
            * start with max_count for each column, and subtract one when we
            * set an element in that column to false
            * 
            * When making the loop par, try to avoid making the outer loop parallel:
            * there are not that many objectives, running parallel over the rows instead
            * will provide the same benefit, but you dont have to use vec_atomic_flag
            */
        }); // 4

        /* For solutions with identical fitness vectors, set the corresponding elements of the dominance matrix to 0. */
        //std::vector<size_t> indices(dom_mat.size());
        std::vector<size_t> indices(dom_mat.nrows());
        std::iota(indices.begin(), indices.end(), 0);

        std::for_each(GA_EXECUTION_UNSEQ, indices.begin(), indices.end(), [&](size_t row)
        {
            std::for_each(indices.begin() + row, indices.end(), [&](size_t col)
            {
                /*if (dom_mat[row][col].load(std::memory_order_relaxed) && dom_mat[col][row].load(std::memory_order_relaxed))
                {
                    dom_mat[row][col].store(false, std::memory_order_relaxed);
                    dom_mat[col][row].store(false, std::memory_order_relaxed);
                }*/
                if (dom_mat(row, col).load(std::memory_order_relaxed) && dom_mat(col, row).load(std::memory_order_relaxed))
                {
                    dom_mat(row, col).store(false, std::memory_order_relaxed);
                    dom_mat(col, row).store(false, std::memory_order_relaxed);
                }
            });
        });

        return dom_mat;
    }

    static ParetoFronts dominanceDegreeSort(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        assert(std::distance(first, last) >= 0);

        const size_t popsize = std::distance(first, last);
        auto& dmat = createDominanceMatrix(first, last); // expensive 7.1

        /* Separate function */

        struct Col { size_t sum; bool removed = false; };

        // try to transpose the matrix here

        std::vector<Col> cols(dmat.ncols());
        dmat.for_each([&](size_t, size_t col, auto& val) { cols[col].sum += val.load(std::memory_order_relaxed); }); // overall 3.0 to create cols


        ParetoFronts pfronts;
        pfronts.reserve(popsize);

        size_t front = 0;
        while (pfronts.size() != popsize)
        {
            const auto front_indices = detail::find_indices(cols, [](const Col& col) { return !col.removed && !col.sum; });

            std::for_each(front_indices.begin(), front_indices.end(), [&](size_t idx)
            {
                pfronts.emplace_back(idx, front);

                /* Remove col */
                cols[idx].removed = true;

                /* Remove row */
                for (size_t i = 0; i < cols.size(); i++)
                {
                    // if not removed unhelpful
                    //cols[i].sum -= dmat[idx][i].load(std::memory_order_relaxed);
                    //dmat[idx][i].store(false, std::memory_order_relaxed);
                    cols[i].sum -= dmat(idx, i).load(std::memory_order_relaxed);
                    dmat(idx, i).store(false, std::memory_order_relaxed);
                }
            });

            front++;
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