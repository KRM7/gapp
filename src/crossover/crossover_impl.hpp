/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_CROSSOVER_DTL_HPP
#define GAPP_CROSSOVER_DTL_HPP

#include "neighbour_list.hpp"
#include "../core/candidate.hpp"
#include "../utility/small_vector.hpp"
#include "../utility/dynamic_bitset.hpp"
#include <vector>
#include <span>
#include <concepts>
#include <cstddef>

namespace gapp::crossover::dtl
{
    /* General n-point crossover implementation for any gene type. */
    template<typename T>
    CandidatePair<T> nPointCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, small_vector<size_t> crossover_points);

    /* Simpler single-point crossover function for any gene type. */
    template<typename T>
    CandidatePair<T> singlePointCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, size_t crossover_point);

    /* Simpler two-point crossover function for any gene type. */
    template<typename T>
    CandidatePair<T> twoPointCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, std::pair<size_t, size_t> crossover_points);


    /* Implementation of the order-1 crossover for any gene type, only generates a single child. */
    template<typename T>
    Candidate<T> order1CrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, size_t first, size_t last);

    /* Implementation of the order-1 crossover for unsigned integers, only generates a single child. */
    template<std::unsigned_integral T>
    Candidate<T> order1CrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, size_t first, size_t last);


    /* Implementation of the order-2 crossover for any gene type, only generates a single child. */
    template<typename T>
    Candidate<T> order2CrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, size_t first, size_t last);

    /* Implementation of the order-2 crossover for unsigned integer genes, only generates a single child. */
    template<std::unsigned_integral T>
    Candidate<T> order2CrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, size_t first, size_t last);


    /* Implementation of the position crossover for any gene type, only generates a single child. */
    template<typename T>
    Candidate<T> positionCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, std::span<const size_t> indices);

    /* Implementation of the position crossover for unsigned integer genes, only generates a single child. */
    template<std::unsigned_integral T>
    Candidate<T> positionCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, std::span<const size_t> indices);


    /* Find the indices of genes in the chromosomes chrom1 and chrom2 which belong to odd cycles. */
    template<typename T>
    std::vector<size_t> findOddCycleIndices(const Chromosome<T>& chrom1, const Chromosome<T>& chrom2);

    /* Find the indices of genes in the chromosomes chrom1 and chrom2 which belong to odd cycles. */
    template<std::unsigned_integral T>
    std::vector<size_t> findOddCycleIndices(const Chromosome<T>& chrom1, const Chromosome<T>& chrom2);

    /* Implementation of the cycle crossover for any gene type. */
    template<typename T>
    CandidatePair<T> cycleCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2);


    /* Implementation of the edge crossover for any gene type, only generates a single child. */
    template<typename T>
    Candidate<T> edgeCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2);

    /* Implementation of the edge crossover for unsigned integer genes, only generates a single child. */
    template<std::unsigned_integral T>
    Candidate<T> edgeCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2);


    /* Implementation of the PMX crossover for any gene type, only generates a single child. */
    template<typename T>
    Candidate<T> pmxCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, size_t first, size_t last);

    /* Implementation of the PMX crossover for unsigned integer genes, only generates a single child. */
    template<std::unsigned_integral T>
    Candidate<T> pmxCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, size_t first, size_t last);

} // namespace gapp::crossover::dtl


/* IMPLEMENTATION */

#include "../utility/algorithm.hpp"
#include "../utility/functional.hpp"
#include "../utility/utility.hpp"
#include <unordered_set>
#include <algorithm>

namespace gapp::crossover::dtl
{
    template<typename T>
    CandidatePair<T> nPointCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, small_vector<size_t> crossover_points)
    {
        const size_t chrom_len = parent1.chromosome.size();

        GAPP_ASSERT(parent1.chromosome.size() == parent2.chromosome.size());
        GAPP_ASSERT(std::all_of(crossover_points.begin(), crossover_points.end(), detail::between(0_sz, chrom_len)));

        std::sort(crossover_points.begin(), crossover_points.end());
        if (crossover_points.size() % 2) crossover_points.push_back(chrom_len);

        Candidate child1{ parent2 }, child2{ parent1 };

        for (size_t i = 1; i < crossover_points.size(); i += 2)
        {
            for (size_t j = crossover_points[i - 1]; j < crossover_points[i]; j++)
            {
                using std::swap;
                swap(child1.chromosome[j], child2.chromosome[j]);
            }
        }

        return { std::move(child1), std::move(child2) };
    }

    template<typename T>
    CandidatePair<T> singlePointCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, size_t crossover_point)
    {
        GAPP_ASSERT(crossover_point <= parent1.chromosome.size());
        GAPP_ASSERT(parent1.chromosome.size() == parent2.chromosome.size());

        Candidate child1{ parent1 }, child2{ parent2 };

        for (size_t i = 0; i < crossover_point; i++)
        {
            using std::swap;
            swap(child1.chromosome[i], child2.chromosome[i]);
        }

        return { std::move(child1), std::move(child2) };
    }

    template<typename T>
    CandidatePair<T> twoPointCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, std::pair<size_t, size_t> crossover_points)
    {
        GAPP_ASSERT(crossover_points.first <= parent1.chromosome.size());
        GAPP_ASSERT(crossover_points.second <= parent1.chromosome.size());
        GAPP_ASSERT(parent1.chromosome.size() == parent2.chromosome.size());

        if (crossover_points.first > crossover_points.second)
        {
            std::swap(crossover_points.first, crossover_points.second);
        }

        Candidate child1{ parent1 }, child2{ parent2 };

        for (size_t i = crossover_points.first; i < crossover_points.second; i++)
        {
            using std::swap;
            swap(child1.chromosome[i], child2.chromosome[i]);
        }

        return { std::move(child1), std::move(child2) };
    }


    template<std::unsigned_integral T>
    bool isValidIntegerPermutation(const Chromosome<T>& chrom)
    {
        if (chrom.empty()) return true;

        if (*std::min_element(chrom.begin(), chrom.end()) != 0) return false;
        if (*std::max_element(chrom.begin(), chrom.end()) != chrom.size() - 1) return false;

        detail::dynamic_bitset present(chrom.size());
        for (const T& val : chrom)
        {
            if (present[val]) return false;
            present[val] = true;
        }

        return true;
    }

    template<typename T>
    Candidate<T> order1CrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, size_t first, size_t last)
    {
        const size_t chrom_len = parent1.chromosome.size();
        const size_t range_len = last - first;

        GAPP_ASSERT(first <= last && last <= chrom_len);
        GAPP_ASSERT(parent1.chromosome.size() == parent2.chromosome.size());

        std::unordered_set<T> direct(last - first);
        for (size_t idx = first; idx != last; idx++) direct.insert(parent1.chromosome[idx]);

        Candidate child = parent1;

        size_t parent_pos = (last == chrom_len) ? 0 : last;
        size_t child_pos = (last == chrom_len) ? 0 : last;

        for (size_t i = 0; i < chrom_len - range_len; i++)
        {
            while (direct.contains(parent2.chromosome[parent_pos])) detail::increment_mod(parent_pos, chrom_len);

            child.chromosome[child_pos] = parent2.chromosome[parent_pos];

            detail::increment_mod(parent_pos, chrom_len);
            detail::increment_mod(child_pos, chrom_len);
        }

        return child;
    }

    template<std::unsigned_integral T>
    Candidate<T> order1CrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, size_t first, size_t last)
    {
        const size_t chrom_len = parent1.chromosome.size();
        const size_t range_len = last - first;

        GAPP_ASSERT(first <= last && last <= chrom_len);
        GAPP_ASSERT(parent1.chromosome.size() == parent2.chromosome.size());
        GAPP_ASSERT(isValidIntegerPermutation(parent1.chromosome));
        GAPP_ASSERT(isValidIntegerPermutation(parent2.chromosome));

        detail::dynamic_bitset is_direct(chrom_len);
        for (size_t idx = first; idx != last; idx++) is_direct[parent1.chromosome[idx]] = true;

        Candidate child = parent1;

        size_t parent_pos = (last == chrom_len) ? 0 : last;
        size_t child_pos  = (last == chrom_len) ? 0 : last;

        for (size_t i = 0; i < chrom_len - range_len; i++)
        {
            while (is_direct[parent2.chromosome[parent_pos]]) detail::increment_mod(parent_pos, chrom_len);

            child.chromosome[child_pos] = parent2.chromosome[parent_pos];

            detail::increment_mod(parent_pos, chrom_len);
            detail::increment_mod(child_pos, chrom_len);
        }

        return child;
    }


    template<typename T>
    Candidate<T> order2CrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, size_t first, size_t last)
    {
        GAPP_ASSERT(first <= last && last <= parent1.chromosome.size());
        GAPP_ASSERT(parent1.chromosome.size() == parent2.chromosome.size());

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

    template<std::unsigned_integral T>
    Candidate<T> order2CrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, size_t first, size_t last)
    {
        const size_t chrom_len = parent1.chromosome.size();

        GAPP_ASSERT(first <= last && last <= chrom_len);
        GAPP_ASSERT(parent1.chromosome.size() == parent2.chromosome.size());
        GAPP_ASSERT(isValidIntegerPermutation(parent1.chromosome));
        GAPP_ASSERT(isValidIntegerPermutation(parent2.chromosome));

        detail::dynamic_bitset is_direct(chrom_len);
        for (size_t idx = first; idx != last; idx++) is_direct[parent1.chromosome[idx]] = true;

        Candidate child = parent1;

        for (size_t child_pos = 0; const T& gene : parent2.chromosome)
        {
            if (!is_direct[gene])
            {
                if (child_pos == first) child_pos = last; // skip [first, last)
                child.chromosome[child_pos++] = gene;
            }
        }

        return child;
    }


    template<typename T>
    Candidate<T> positionCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, std::span<const size_t> indices)
    {
        GAPP_ASSERT(std::all_of(indices.begin(), indices.end(), detail::between(0_sz, parent1.chromosome.size() - 1)));
        GAPP_ASSERT(parent1.chromosome.size() == parent2.chromosome.size());

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

    template<std::unsigned_integral T>
    Candidate<T> positionCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, std::span<const size_t> indices)
    {
        const size_t chrom_len = parent1.chromosome.size();
        
        GAPP_ASSERT(std::all_of(indices.begin(), indices.end(), detail::between(0_sz, chrom_len - 1)));
        GAPP_ASSERT(parent1.chromosome.size() == parent2.chromosome.size());
        GAPP_ASSERT(isValidIntegerPermutation(parent1.chromosome));
        GAPP_ASSERT(isValidIntegerPermutation(parent2.chromosome));

        detail::dynamic_bitset is_direct(chrom_len);
        for (size_t idx : indices) is_direct[parent1.chromosome[idx]] = true;

        small_vector<size_t> next_indirect(chrom_len);
        for (ptrdiff_t indirect = -1, i = chrom_len - 1; i >= 0; i--)
        {
            const T gene = parent1.chromosome[i];

            indirect = is_direct[gene] ? indirect : i;
            next_indirect[i] = indirect;
        }

        Candidate child = parent1;

        for (size_t child_pos = 0; T gene : parent2.chromosome)
        {
            if (!is_direct[gene])
            {
                child_pos = next_indirect[child_pos];
                child.chromosome[child_pos++] = gene;
            }
        }

        return child;
    }


    template<typename T>
    std::vector<size_t> findOddCycleIndices(const Chromosome<T>& chrom1, const Chromosome<T>& chrom2)
    {
        GAPP_ASSERT(chrom1.size() == chrom2.size());

        const size_t chrom_len = chrom1.size();

        std::vector<size_t> odd_indices;
        odd_indices.reserve(chrom_len / 2);

        detail::dynamic_bitset deleted(chrom_len);
        size_t num_deleted = 0;

        for (bool odd_cycle = false; num_deleted < chrom_len; odd_cycle ^= 1)
        {
            size_t pos = deleted.find_first(false);
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

    template<std::unsigned_integral T>
    std::vector<size_t> findOddCycleIndices(const Chromosome<T>& chrom1, const Chromosome<T>& chrom2)
    {
        GAPP_ASSERT(chrom1.size() == chrom2.size());
        GAPP_ASSERT(isValidIntegerPermutation(chrom1));
        GAPP_ASSERT(isValidIntegerPermutation(chrom2));

        const size_t chrom_len = chrom1.size();

        std::vector<size_t> odd_indices;
        odd_indices.reserve(chrom_len / 2);

        detail::dynamic_bitset deleted(chrom_len);
        size_t num_deleted = 0;

        std::vector index_lookup(chrom_len, 0_sz);
        for (size_t i = 0; i < chrom_len; i++)
        {
            index_lookup[chrom1[i]] = i;
        }

        for (bool odd_cycle = false; num_deleted < chrom_len; odd_cycle ^= 1)
        {
            size_t pos = deleted.find_first(false);
            T cycle_start = chrom1[pos];

            deleted[pos] = true;
            num_deleted++;

            if (odd_cycle) odd_indices.push_back(pos);

            while (chrom2[pos] != cycle_start)
            {
                pos = index_lookup[chrom2[pos]];

                deleted[pos] = true;
                num_deleted++;

                if (odd_cycle) odd_indices.push_back(pos);
            }
        }

        return odd_indices;
    }

    template<typename T>
    CandidatePair<T> cycleCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2)
    {
        GAPP_ASSERT(parent1.chromosome.size() == parent2.chromosome.size());

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


    template<typename T>
    Candidate<T> edgeCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2)
    {
        GAPP_ASSERT(parent1.chromosome.size() == parent2.chromosome.size());

        const size_t chrom_len = parent1.chromosome.size();

        auto nb_lists = makeNeighbourLists(parent1.chromosome, parent2.chromosome);

        Candidate<T> child({ parent1.chromosome[0] });
        child.chromosome.reserve(chrom_len);

        std::vector<T> remaining_genes(parent1.chromosome.begin() + 1, parent1.chromosome.end());

        while (child.chromosome.size() != chrom_len)
        {
            const T& last_gene = child.chromosome.back();
            T next_gene = remaining_genes.front();

            for (T neighbour : nb_lists[last_gene])
            {
                nb_lists[neighbour].remove(last_gene);

                if (nb_lists[neighbour].size() <= nb_lists[next_gene].size())
                {
                    next_gene = neighbour;
                }
            }

            child.chromosome.push_back(next_gene);
            std::erase(remaining_genes, next_gene);
        }

        return child;
    }

    template<std::unsigned_integral T>
    Candidate<T> edgeCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2)
    {
        const size_t chrom_len = parent1.chromosome.size();

        GAPP_ASSERT(parent1.chromosome.size() == parent2.chromosome.size());
        GAPP_ASSERT(isValidIntegerPermutation(parent1.chromosome));
        GAPP_ASSERT(isValidIntegerPermutation(parent2.chromosome));

        auto nb_lists = makeNeighbourLists(parent1.chromosome, parent2.chromosome);

        Candidate<T> child({ parent1.chromosome[0] });
        child.chromosome.reserve(chrom_len);

        detail::dynamic_bitset is_used(chrom_len);
        is_used[parent1.chromosome[0]] = true;

        while (child.chromosome.size() != chrom_len)
        {
            T last_gene = child.chromosome.back();
            T next_gene = is_used.find_first(false);

            for (T neighbour : nb_lists[last_gene])
            {
                if (neighbour == NeighbourList<T>::EMPTY) continue;

                nb_lists[neighbour].remove(last_gene);

                if (nb_lists[neighbour].size() <= nb_lists[next_gene].size())
                {
                    next_gene = neighbour;
                }
            }

            child.chromosome.push_back(next_gene);
            is_used[next_gene] = true;
        }

        return child;
    }


    template<typename T>
    Candidate<T> pmxCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, size_t first, size_t last)
    {
        GAPP_ASSERT(parent1.chromosome.size() == parent2.chromosome.size());
        GAPP_ASSERT(first <= last && last <= parent1.chromosome.size());

        Candidate child = parent2;

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

    template<std::unsigned_integral T>
    Candidate<T> pmxCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, size_t first, size_t last)
    {
        const size_t chrom_len = parent1.chromosome.size();

        GAPP_ASSERT(parent1.chromosome.size() == parent2.chromosome.size());
        GAPP_ASSERT(first <= last && last <= parent1.chromosome.size());
        GAPP_ASSERT(isValidIntegerPermutation(parent1.chromosome));
        GAPP_ASSERT(isValidIntegerPermutation(parent2.chromosome));

        Candidate child = parent2;

        detail::dynamic_bitset is_direct(chrom_len);
        for (size_t i = first; i < last; i++)
        {
            child.chromosome[i] = parent1.chromosome[i];
            is_direct[parent1.chromosome[i]] = true;
        }

        std::vector index_lookup(chrom_len, 0_sz); // for parent2
        for (size_t i = 0; i < chrom_len; i++)
        {
            index_lookup[parent2.chromosome[i]] = i;
        }

        for (size_t i = first; i < last; i++)
        {
            if (!is_direct[parent2.chromosome[i]])
            {
                size_t pos = i;
                while (first <= pos && pos < last)
                {
                    pos = index_lookup[parent1.chromosome[pos]];
                }
                child.chromosome[pos] = parent2.chromosome[i];
            }
        }

        return child;
    }
    
} // namespace gapp::crossover::dtl

#endif // !GAPP_CROSSOVER_DTL_HPP
