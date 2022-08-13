/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "nd_sort.hpp"
#include "../utility/algorithm.hpp"
#include "../utility/math.hpp"
#include <algorithm>
#include <iterator>
#include <numeric>
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

    template<typename T>
    using Matrix = std::vector<std::vector<T>>;

    template<typename T>
    using MRow = Matrix<T>::value_type;

    static Matrix<int>& getDominanceMatrix(size_t size)
    {
        static Matrix<int> dmat;

        if (dmat.size() != size)
        {
            dmat.resize(size);
            for (auto& row : dmat) row.resize(size);
        }
        for (auto& row : dmat) std::fill(row.begin(), row.end(), 0);

        return dmat;
    }

    static Matrix<bool>& getNormDominanceMatrix(size_t size)
    {
        static Matrix<bool> ndmat;

        if (ndmat.size() != size)
        {
            ndmat.resize(size);
            for (auto& row : ndmat) row.resize(size);
        }

        return ndmat;
    }

    static Matrix<int>& createDominanceMatrix(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        assert(std::distance(first, last) >= 0);

        const size_t popsize = std::distance(first, last);

        auto& dom_mat = getDominanceMatrix(popsize);

        if (popsize == 0) return dom_mat;

        std::vector<size_t> objs(first->size());
        std::iota(objs.begin(), objs.end(), 0);

        /* Could be parallel with atomic ints as the dominance matrix values, could be par, TODO: benchmark */
        std::for_each(objs.cbegin(), objs.cend(), [&](size_t obj)
        {
            const auto fmat_rows = detail::argsort(first, last, [&](const FitnessVector& lhs, const FitnessVector& rhs) { return lhs[obj] > rhs[obj]; });

            auto current_row = fmat_rows.begin();
            while (current_row != fmat_rows.end())
            {
                auto first_col = current_row;
                auto last_col = fmat_rows.end();

                auto next_row = std::find_if(current_row, fmat_rows.end(), [&](size_t row)
                {
                    return !detail::floatIsEqual(first[*current_row][obj], first[row][obj]);
                });

                std::for_each(current_row, next_row, [&](size_t row)
                {
                    std::for_each(first_col, last_col, [&](size_t col) { ++dom_mat[row][col]; });
                });
                current_row = next_row;
            }
        });

        /* For solutions with identical fitness vectors, set the corresponding elements of the dominance matrix to 0. */
        for (size_t i = 0; i != dom_mat.size(); i++)
        {
            dom_mat[i][i] = 0;
            for (size_t j = i + 1; j != dom_mat.size(); j++)
            {
                if (dom_mat[i][j] == objs.size() && dom_mat[j][i] == objs.size())
                {
                    dom_mat[i][j] = dom_mat[j][i] = 0;
                }
            }
        }

        return dom_mat;
    }

    static Matrix<bool>& normalizeDominanceMatrix(const Matrix<int>& dmat, size_t nobjs)
    {
        auto& ndmat = getNormDominanceMatrix(dmat.size());

        ndmat.resize(dmat.size());
        for (auto& row : ndmat) row.resize(dmat.size());

        for (size_t row = 0; row < ndmat.size(); row++)
        {
            for (size_t col = 0; col < ndmat.size(); col++)
            {
                ndmat[row][col] = (dmat[row][col] == nobjs);
            }
        }

        return ndmat;
    }

    static ParetoFronts dominanceDegreeSort(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        assert(std::distance(first, last) >= 0);

        const size_t popsize = std::distance(first, last);

        ParetoFronts pfronts;
        auto& dmat = createDominanceMatrix(first, last); // expensive 11.4

        /* Separate function */
        struct Row { size_t idx; int sum; int minus; };

        std::vector<Row> cols(dmat.size());
        for (size_t i = 0; i < dmat.size(); i++)
        {
            cols[i].idx = i;
            for (const auto& row : dmat) cols[i].sum += (row[i] == first->size()); // expensive 4.8, wrong access pattern, do it the other way around
        }

        size_t front = 0;

        while (pfronts.size() != popsize) // while rows !empty
        {
            for (ptrdiff_t i = cols.size() - 1; i >= 0; i--)
            {
                if (cols[i].sum == 0)
                {
                    pfronts.emplace_back(cols[i].idx, front); // 0.5

                    /* Remove row */ // this is before the remove col stuff because we use cols[i].idx (idx to remove)
                    for (auto& [idx, sum, minus] : cols)
                    {
                        if (dmat[cols[i].idx][idx] == first->size()) minus++; // expensive 5.0 for the entire loop
                        dmat[cols[i].idx][idx] = 0;
                    }

                    /* Remove col */
                    std::swap(cols[i], cols.back());
                    cols.pop_back();
                }
            }

            for (auto& col : cols)
            {
                col.sum -= col.minus;
                col.minus = 0;
            }

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