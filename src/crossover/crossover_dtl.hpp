/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CROSSOVER_DTL_HPP
#define GA_CROSSOVER_DTL_HPP

#include "../population/candidate.hpp"
#include "../utility/algorithm.hpp"
#include "../utility/functional.hpp"
#include "../utility/utility.hpp"
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cstddef>

namespace genetic_algorithm::crossover::dtl
{
    /* General n-point crossover implementation for any gene type. */
    template<Gene T>
    CandidatePair<T> nPointCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, size_t n);

    /* Simpler single-point crossover function for any gene type. */
    template<Gene T>
    CandidatePair<T> singlePointCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2);

    /* Simpler two-point crossover function for any gene type. */
    template<Gene T>
    CandidatePair<T> twoPointCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2);

    /* Implementation of the order-1 crossover for any gene type, only generates a single child. */
    template<Gene T>
    Candidate<T> order1CrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, size_t first, size_t last);

    /* Implementation of the order-2 crossover for any gene type, only generates a single child. */
    template<Gene T>
    Candidate<T> order2CrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, size_t first, size_t last);

    /* Implementation of the position crossover for any gene type, only generates a single child. */
    template<Gene T>
    Candidate<T> positionCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, const std::vector<size_t>& indices);

    /* Find the indices of genes in the chromosomes chrom1 and chrom2 which belong to odd cycles. Used in the cycle crossover operator. */
    template<typename T>
    std::vector<size_t> findOddCycleIndices(const Chromosome<T>& chrom1, const Chromosome<T>& chrom2);

    /* Implementation of the cycle crossover for any gene type. */
    template<Gene T>
    CandidatePair<T> cycleCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2);

    /* Implementation of the edge crossover for unsigned integer gene types, only generates a single child. */
    template<std::unsigned_integral T>
    Candidate<T> edgeCrossoverImpl(const Candidate<T>& parent1, std::unordered_map<T, std::vector<T>>&& neighbour_lists);

    /* Implementation of the PMX crossover for any gene type, only generates a single child. */
    template<Gene T>
    Candidate<T> pmxCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2);

} // namespace genetic_algorithm::crossover::dtl


/* IMPLEMENTATION */

#include "../utility/rng.hpp"
#include "../utility/algorithm.hpp"
#include "../utility/utility.hpp"
#include <unordered_set>
#include <algorithm>
#include <stdexcept>
#include <cassert>

namespace genetic_algorithm::crossover::dtl
{
    template<Gene T>
    CandidatePair<T> nPointCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, size_t n)
    {
        if (parent1.chromosome.size() == parent2.chromosome.size())
        {
            GA_THROW(std::invalid_argument, "The parent chromosomes must be the same length for the n-point crossover.");
        }

        const size_t chrom_len = parent1.chromosome.size();
        const size_t num_crossover_points = std::min(n, chrom_len);
        const auto crossover_points = rng::sampleUnique(0_sz, chrom_len, num_crossover_points);

        /* Create crossover mask */
        std::vector<size_t> crossover_mask(chrom_len, 0);
        for (size_t i = 0, remaining = num_crossover_points; i < chrom_len; i++)
        {
            if (detail::contains(crossover_points.begin(), crossover_points.end(), i))
            {
                if (--remaining == 0) break;
            }
            crossover_mask[i] = remaining;
        }

        Candidate child1 = parent1;
        Candidate child2 = parent2;

        for (size_t i = 0; i < chrom_len; i++)
        {
            if (crossover_mask[i] % 2)
            {
                using std::swap;
                swap(child1.chromosome[i], child2.chromosome[i]);
            }
        }

        return { std::move(child1), std::move(child2) };
    }

    template<Gene T>
    CandidatePair<T> singlePointCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2)
    {
        if (parent1.chromosome.size() != parent2.chromosome.size())
        {
            GA_THROW(std::invalid_argument, "The parent chromosomes must be the same length for the n-point crossover.");
        }

        const size_t chrom_len = parent1.chromosome.size();
        const size_t crossover_point = rng::randomInt(0_sz, chrom_len);

        Candidate child1 = parent1;
        Candidate child2 = parent2;

        for (size_t i = 0; i < crossover_point; i++)
        {
            using std::swap;
            swap(child1.chromosome[i], child2.chromosome[i]);
        }

        return { std::move(child1), std::move(child2) };
    }

    template<Gene T>
    CandidatePair<T> twoPointCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2)
    {
        if (parent1.chromosome.size() != parent2.chromosome.size())
        {
            GA_THROW(std::invalid_argument, "The parent chromosomes must be the same length for the n-point crossover.");
        }

        const size_t chrom_len = parent1.chromosome.size();
        auto crossover_points = rng::sampleUnique(0_sz, chrom_len, 2);
        if (crossover_points[0] > crossover_points[1])
        {
            std::swap(crossover_points[0], crossover_points[1]);
        }

        Candidate child1 = parent1;
        Candidate child2 = parent2;

        for (size_t i = crossover_points[0]; i < crossover_points[1]; i++)
        {
            using std::swap;
            swap(child1.chromosome[i], child2.chromosome[i]);
        }

        return { std::move(child1), std::move(child2) };
    }

    template<Gene T>
    Candidate<T> order1CrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, size_t first, size_t last)
    {
        const size_t chrom_len = parent1.chromosome.size();

        std::unordered_set<T> direct(last - first);
        for (size_t idx = first; idx != last; idx++) direct.insert(parent1.chromosome[idx]);

        Candidate child = parent1;

        for (size_t child_pos = last, parent_pos = last; child_pos % chrom_len != first; child_pos++)
        {
            while (direct.contains(parent2.chromosome[parent_pos % chrom_len])) parent_pos++;
            child.chromosome[child_pos % chrom_len] = parent2.chromosome[parent_pos++ % chrom_len];
        }

        return child;
    }

    template<Gene T>
    Candidate<T> order2CrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, size_t first, size_t last)
    {
        std::unordered_set<T> direct(last - first);
        for (size_t idx = first; idx != last; idx++) direct.insert(parent1.chromosome[idx]);

        Candidate child = parent1;

        for (size_t child_pos = 0; const T& gene : parent2.chromosome)
        {
            if (!direct.contains(gene))
            {
                if (child_pos == first) child_pos = last; // skip [first, last)
                child.chromosome[child_pos++] = gene;
            }
        }

        return child;
    }

    template<Gene T>
    Candidate<T> positionCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, const std::vector<size_t>& indices)
    {
        std::unordered_set<T> direct(indices.size());
        for (size_t idx : indices) direct.insert(parent1.chromosome[idx]);

        Candidate child = parent1;

        for (auto child_pos = child.chromosome.begin(); const T& gene : parent2.chromosome)
        {
            if (!direct.contains(gene))
            {
                while (direct.contains(*child_pos)) ++child_pos;
                *child_pos++ = gene;
            }
        }

        return child;
    }

    template<typename T>
    std::vector<size_t> findOddCycleIndices(const Chromosome<T>& chrom1, const Chromosome<T>& chrom2)
    {
        const size_t chrom_len = chrom1.size();

        std::vector<size_t> odd_indices;
        odd_indices.reserve(chrom_len / 2);

        std::vector<bool> deleted(chrom_len, false);
        size_t num_deleted = 0;

        for (bool odd_cycle = 0; num_deleted < chrom_len; odd_cycle ^= 1)
        {
            size_t pos = *detail::index_of(deleted, false);
            const T cycle_start = chrom1[pos];

            deleted[pos] = true;
            num_deleted++;

            if (odd_cycle) odd_indices.push_back(pos);

            while (chrom2[pos] != cycle_start)
            {
                pos = *detail::index_of(chrom1, chrom2[pos]);

                deleted[pos] = true;
                num_deleted++;

                if (odd_cycle) odd_indices.push_back(pos);
            }
        }

        return odd_indices;
    }

    template<Gene T>
    CandidatePair<T> cycleCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2)
    {
        const auto odd_cycle_idxs = dtl::findOddCycleIndices(parent1.chromosome, parent2.chromosome);

        Candidate child1 = parent1;
        Candidate child2 = parent2;

        for (size_t idx : odd_cycle_idxs)
        {
            using std::swap;
            swap(child1.chromosome[idx], child2.chromosome[idx]);
        }

        return { std::move(child1), std::move(child2) };
    }

    template<std::unsigned_integral T>
    class NeighbourLists
    {
    public:
        struct NeighbourRange
        {
            typename std::vector<T>::iterator first;
            typename std::vector<T>::iterator last;
        };

        NeighbourLists(const Chromosome<T>& chrom1, const Chromosome<T>& chrom2) :
            neighbour_lists_(chrom1.size() * MAX_NEIGHBOURS, EMPTY_VAL)
        {
            assert(chrom1.size() == chrom2.size());

            const size_t len = chrom1.size();

            add(chrom1[0], chrom1[1]);
            add(chrom1[0], chrom1[len - 1]);

            for (size_t i = 1; i < len - 1; i++)
            {
                add(chrom1[i], chrom1[i + 1]);
                add(chrom1[i], chrom1[i - 1]);

                add(chrom2[i], chrom2[i + 1]);
                add(chrom2[i], chrom2[i - 1]);
            }

            add(chrom1[len - 1], chrom1[0]);
            add(chrom1[len - 1], chrom1[len - 2]);
        }

        size_t sizeOf(T n)
        {
            const auto neighbours = rangeOf(n);

            return std::count_if(neighbours.first, neighbours.last, detail::not_equal_to(EMPTY_VAL));
        }

        void add(T to, T neighbour)
        {
            const auto neighbours = rangeOf(to);

            for (auto it = neighbours.first; it != neighbours.last; ++it)
            {
                if (*it == EMPTY_VAL) *it = neighbour;
                if (*it == neighbour) return;
            }
            GA_UNREACHABLE();
        }

        void remove(T neighbour)
        {
            //const auto neighbours = rangeOf(neighbour);

            //for (auto it = neighbours.first; it != neighbours.last; ++it)
            //{
            //    if (*it != EMPTY_VAL)
            //    {
            //        const auto range = rangeOf(*it);

            //        for (auto it2 = range.first; it2 != range.last; ++it2)
            //        {
            //            if (*it2 == neighbour) *it2 = EMPTY_VAL;
            //        }
            //    }
            //}

            for (T& gene : neighbour_lists_)
            {
                if (gene == neighbour) gene = EMPTY_VAL;
            }
        }

        std::vector<T> findBestNeighbours(T n)
        {
            const auto neighbours = rangeOf(n);

            size_t min_count = MAX_NEIGHBOURS;

            for (auto neighbour = neighbours.first; neighbour != neighbours.last; ++neighbour)
            {
                if (*neighbour != EMPTY_VAL)
                {
                    min_count = std::min(min_count, sizeOf(*neighbour));
                }
            }

            std::vector<T> best_neighbours;

            for (auto neighbour = neighbours.first; neighbour != neighbours.last; ++neighbour)
            {
                if (*neighbour != EMPTY_VAL && sizeOf(*neighbour) == min_count)
                {
                    best_neighbours.push_back(*neighbour);
                }
            }

            return best_neighbours;
        }

    private:
        std::vector<T> neighbour_lists_;
        // TODO also keep a vector of the last elements (keep sizes, not iterators to avoid copying issues)? -> no need to call find every time for add, remove, count

        static constexpr T EMPTY_VAL = T(-1);
        static constexpr T MAX_NEIGHBOURS = 4;

        NeighbourRange rangeOf(T n)
        {
            assert(n < (neighbour_lists_.size() / MAX_NEIGHBOURS));

            auto first = neighbour_lists_.begin() + n * MAX_NEIGHBOURS;
            auto last  = neighbour_lists_.begin() + (n + 1) * MAX_NEIGHBOURS;

            return { first, last };
        }
    };

    template<Gene T>
    Candidate<T> edgeCrossoverImpl(const Candidate<T>& parent1, NeighbourLists<T> neighbour_lists)
    {
        const size_t chrom_len = parent1.chromosome.size();

        Candidate<T> child;
        child.chromosome.reserve(chrom_len);

        auto remaining_genes = parent1.chromosome;
        auto gene = parent1.chromosome[0];

        while (child.chromosome.size() != chrom_len)
        {
            /* Add current gene */
            child.chromosome.push_back(gene);
            std::erase(remaining_genes, gene);
            neighbour_lists.remove(gene);

            if (child.chromosome.size() == chrom_len) break;

            auto candidate_genes = neighbour_lists.findBestNeighbours(gene);

            gene = candidate_genes.empty() ?
                rng::randomElement(remaining_genes.begin(), remaining_genes.end()) :
                rng::randomElement(candidate_genes.begin(), candidate_genes.end());
        }

        return child;
    }

    template<Gene T>
    Candidate<T> pmxCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2)
    {
        Candidate child = parent2;

        const size_t chrom_len = parent1.chromosome.size();
        const size_t range_len = rng::randomInt(1_sz, chrom_len - 1);

        const size_t first = rng::randomInt(0_sz, chrom_len - range_len);
        const size_t last = first + range_len;

        std::unordered_set<T> direct(last - first);
        for (size_t i = first; i < last; i++)
        {
            child.chromosome[i] = parent1.chromosome[i];
            direct.insert(parent1.chromosome[i]);
        }

        for (size_t i = first; i < last; i++)
        {
            if (!direct.contains(parent2.chromosome[i]))
            {
                size_t pos = i;
                while (first <= pos && pos < last)
                {
                    pos = *detail::index_of(parent2.chromosome, parent1.chromosome[pos]);
                }
                child.chromosome[pos] = parent2.chromosome[i];
            }
        }

        return child;
    }

} // namespace genetic_algorithm::crossover::dtl

#endif // !GA_CROSSOVER_DTL_HPP