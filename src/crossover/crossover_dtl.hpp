/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CROSSOVER_DTL_HPP
#define GA_CROSSOVER_DTL_HPP

#include "../population/candidate.hpp"
#include <vector>
#include <unordered_map>
#include <type_traits>
#include <concepts>
#include <cstddef>

namespace genetic_algorithm::crossover::dtl
{
    /* General n-point crossover implementation for any gene type. */
    template<Gene T>
    CandidatePair<T> nPointCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, std::vector<size_t> crossover_points);

    /* Simpler single-point crossover function for any gene type. */
    template<Gene T>
    CandidatePair<T> singlePointCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, size_t crossover_point);

    /* Simpler two-point crossover function for any gene type. */
    template<Gene T>
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
    Candidate<T> positionCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, const std::vector<size_t>& indices);

    /* Implementation of the position crossover for unsigned integer genes, only generates a single child. */
    template<std::unsigned_integral T>
    Candidate<T> positionCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, const std::vector<size_t>& indices);


    /* Find the indices of genes in the chromosomes chrom1 and chrom2 which belong to odd cycles. Used in the cycle crossover operator. */
    template<typename T>
    std::vector<size_t> findOddCycleIndices(const Chromosome<T>& chrom1, const Chromosome<T>& chrom2);

    /* Implementation of the cycle crossover for any gene type. */
    template<Gene T>
    CandidatePair<T> cycleCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2);


    /* A list of neighbours for a gene. */
    template<typename T>
    class NeighbourList;

    /* A list of neighbours for an unsigned integer gene. */
    template<std::unsigned_integral T>
    class NeighbourList<T>;

    /* Conctruct the neighbour lists of each gene based on 2 chromosomes. The first and last elements are considered neighbours. */
    template<typename T, typename Ret = std::conditional_t<std::is_unsigned_v<T>, std::vector<NeighbourList<T>>, std::unordered_map<T, NeighbourList<T>>>>
    Ret makeNeighbourLists(const Chromosome<T>& chrom1, const Chromosome<T>& chrom2);

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

} // namespace genetic_algorithm::crossover::dtl


/* IMPLEMENTATION */

#include "../utility/rng.hpp"
#include "../utility/algorithm.hpp"
#include "../utility/functional.hpp"
#include "../utility/iterators.hpp"
#include "../utility/utility.hpp"
#include <array>
#include <unordered_set>
#include <algorithm>
#include <stdexcept>
#include <cassert>

namespace genetic_algorithm::crossover::dtl
{
    template<Gene T>
    CandidatePair<T> nPointCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, std::vector<size_t> crossover_points)
    {
        const size_t chrom_len = parent1.chromosome.size();

        assert(std::all_of(crossover_points.begin(), crossover_points.end(), detail::between(0_sz, chrom_len)));

        std::sort(crossover_points.begin(), crossover_points.end());
        crossover_points.push_back(chrom_len);

        Candidate child1 = parent1;
        Candidate child2 = parent2;

        for (size_t j = 0, i = 0; j < crossover_points.size(); j++)
        {
            if (j % 2)
            {
                i = crossover_points[j];
                continue;
            }
            for (; i < crossover_points[j]; i++)
            {
                using std::swap;
                swap(child1.chromosome[i], child2.chromosome[i]);
            }
        }

        return { std::move(child1), std::move(child2) };
    }

    template<Gene T>
    CandidatePair<T> singlePointCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, size_t crossover_point)
    {
        assert(crossover_point <= parent1.chromosome.size());

        if (parent1.chromosome.size() != parent2.chromosome.size())
        {
            GA_THROW(std::invalid_argument, "The parent chromosomes must be the same length for the n-point crossover.");
        }

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
    CandidatePair<T> twoPointCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, std::pair<size_t, size_t> crossover_points)
    {
        assert(crossover_points.first <= parent1.chromosome.size());
        assert(crossover_points.second <= parent1.chromosome.size());

        if (parent1.chromosome.size() != parent2.chromosome.size())
        {
            GA_THROW(std::invalid_argument, "The parent chromosomes must be the same length for the n-point crossover.");
        }

        if (crossover_points.first > crossover_points.second)
        {
            std::swap(crossover_points.first, crossover_points.second);
        }

        Candidate child1 = parent1;
        Candidate child2 = parent2;

        for (size_t i = crossover_points.first; i < crossover_points.second; i++)
        {
            using std::swap;
            swap(child1.chromosome[i], child2.chromosome[i]);
        }

        return { std::move(child1), std::move(child2) };
    }


    template<typename T>
    Candidate<T> order1CrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, size_t first, size_t last)
    {
        const size_t chrom_len = parent1.chromosome.size();
        const size_t range_len = last - first;

        assert(first <= last && last <= chrom_len);
        assert(parent1.chromosome.size() == parent2.chromosome.size());

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

        /* The genes have to be unique unsigned integers in the range [0, chrom_len). */
        assert(first <= last && last <= chrom_len);
        assert(parent1.chromosome.size() == parent2.chromosome.size());
        assert(*std::min_element(parent1.chromosome.begin(), parent1.chromosome.end()) == 0);
        assert(*std::max_element(parent1.chromosome.begin(), parent1.chromosome.end()) == chrom_len - 1);

        std::vector is_direct(chrom_len, false);
        for (size_t idx = first; idx != last; idx++) is_direct[parent1.chromosome[idx]] = true;

        Candidate child = parent1;

        size_t parent_pos = (last == chrom_len) ? 0 : last;
        size_t child_pos = (last == chrom_len) ? 0 : last;

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
        assert(first <= last && last <= parent1.chromosome.size());
        assert(parent1.chromosome.size() == parent2.chromosome.size());

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

        /* The genes have to be unique unsigned integers in the range [0, chrom_len). */
        assert(first <= last && last <= chrom_len);
        assert(parent1.chromosome.size() == parent2.chromosome.size());
        assert(*std::min_element(parent1.chromosome.begin(), parent1.chromosome.end()) == 0);
        assert(*std::max_element(parent1.chromosome.begin(), parent1.chromosome.end()) == chrom_len - 1);

        std::vector is_direct(chrom_len, false);
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
    Candidate<T> positionCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, const std::vector<size_t>& indices)
    {
        assert(std::all_of(indices.begin(), indices.end(), detail::between(0_sz, parent1.chromosome.size() - 1)));
        assert(parent1.chromosome.size() == parent2.chromosome.size());

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
    Candidate<T> positionCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, const std::vector<size_t>& indices)
    {
        const size_t chrom_len = parent1.chromosome.size();

        /* The genes have to be unique unsigned integers in the range [0, chrom_len). */
        assert(std::all_of(indices.begin(), indices.end(), detail::between(0_sz, parent1.chromosome.size() - 1)));
        assert(parent1.chromosome.size() == parent2.chromosome.size());
        assert(*std::min_element(parent1.chromosome.begin(), parent1.chromosome.end()) == 0);
        assert(*std::max_element(parent1.chromosome.begin(), parent1.chromosome.end()) == chrom_len - 1);

        std::vector is_direct(chrom_len, false);
        for (size_t idx : indices) is_direct[parent1.chromosome[idx]] = true;

        Candidate child = parent1;

        for (auto child_pos = child.chromosome.begin(); const T& gene : parent2.chromosome)
        {
            if (!is_direct[gene])
            {
                while (is_direct[*child_pos]) ++child_pos;
                *child_pos++ = gene;
            }
        }

        return child;
    }


    template<typename T>
    std::vector<size_t> findOddCycleIndices(const Chromosome<T>& chrom1, const Chromosome<T>& chrom2)
    {
        assert(chrom1.size() == chrom2.size());

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
        assert(parent1.chromosome.size() == parent2.chromosome.size());

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
    class NeighbourList : public detail::reverse_iterator_interface<NeighbourList<T>>
    {
    public:
        NeighbourList() { neighbours_.reserve(4); }

        void add(const T& value)
        {
            if (!detail::contains(neighbours_.cbegin(), neighbours_.cend(), value))
            {
                neighbours_.push_back(value);
            }
        }

        void remove(const T& value) { detail::erase_first_stable(neighbours_, value); }

        constexpr size_t size() const noexcept { return neighbours_.size(); }
        constexpr bool empty() const noexcept  { return neighbours_.empty(); }

        constexpr auto begin() noexcept       { return neighbours_.begin(); }
        constexpr auto end() noexcept         { return neighbours_.end(); }
        constexpr auto begin() const noexcept { return neighbours_.begin(); }
        constexpr auto end() const noexcept   { return neighbours_.end(); }

    private:
        std::vector<T> neighbours_;
    };

    template<std::unsigned_integral T>
    class NeighbourList<T> : public detail::reverse_iterator_interface<NeighbourList<T>>
    {
    public:
        void add(T value)
        {
            assert(value != EMPTY);

            /* Assume that EMPTY values are at the back. */
            for (T& neighbour : neighbours_)
            {
                if (neighbour == value) return;
                if (neighbour == EMPTY)
                {
                    neighbour = value;
                    return;
                }
            }
            GA_UNREACHABLE();
        }

        void remove(T value)
        {
            assert(value != EMPTY);

            for (T& neighbour : neighbours_)
            {
                neighbour = (neighbour == value) ? EMPTY : neighbour;
            }
        }

        constexpr size_t size() const noexcept
        {
            return std::count_if(neighbours_.begin(), neighbours_.end(), detail::not_equal_to(EMPTY));
        }

        constexpr bool empty() const noexcept { return size() == 0; }

        constexpr auto begin() noexcept       { return neighbours_.begin(); }
        constexpr auto end() noexcept         { return neighbours_.end(); }
        constexpr auto begin() const noexcept { return neighbours_.begin(); }
        constexpr auto end() const noexcept   { return neighbours_.end(); }

        static constexpr T EMPTY = T(-1);

    private:
        std::array<T, 4> neighbours_{ EMPTY, EMPTY, EMPTY, EMPTY };
    };

    template<typename T, typename Ret>
    Ret makeNeighbourLists(const Chromosome<T>& chrom1, const Chromosome<T>& chrom2)
    {
        assert(chrom1.size() == chrom2.size());

        Ret nb_lists(chrom1.size());

        nb_lists[chrom1.front()].add(chrom1.back());
        nb_lists[chrom1.front()].add(chrom1[1]);
        nb_lists[chrom2.front()].add(chrom2.back());
        nb_lists[chrom2.front()].add(chrom2[1]);

        for (size_t i = 1; i < chrom1.size() - 1; i++)
        {
            nb_lists[chrom1[i]].add(chrom1[i - 1]);
            nb_lists[chrom1[i]].add(chrom1[i + 1]);

            nb_lists[chrom2[i]].add(chrom2[i - 1]);
            nb_lists[chrom2[i]].add(chrom2[i + 1]);
        }

        nb_lists[chrom1.back()].add(*(chrom1.end() - 2));
        nb_lists[chrom1.back()].add(chrom1.front());
        nb_lists[chrom2.back()].add(*(chrom2.end() - 2));
        nb_lists[chrom2.back()].add(chrom2.front());

        return nb_lists;
    }

    template<typename T>
    Candidate<T> edgeCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2)
    {
        assert(parent1.chromosome.size() == parent2.chromosome.size());

        const size_t chrom_len = parent1.chromosome.size();

        auto nb_lists = makeNeighbourLists(parent1.chromosome, parent2.chromosome);

        Candidate<T> child{ parent1.chromosome[0] };
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

        /* The genes have to be unique unsigned integers in the range [0, chrom_len). */
        assert(parent1.chromosome.size() == parent2.chromosome.size());
        assert(*std::min_element(parent1.chromosome.begin(), parent1.chromosome.end()) == 0);
        assert(*std::max_element(parent1.chromosome.begin(), parent1.chromosome.end()) == chrom_len - 1);

        auto nb_lists = makeNeighbourLists(parent1.chromosome, parent2.chromosome);

        Candidate<T> child{ parent1.chromosome[0] };
        child.chromosome.reserve(chrom_len);

        std::vector is_used(chrom_len, false);
        is_used[parent1.chromosome[0]] = true;

        while (child.chromosome.size() != chrom_len)
        {
            T last_gene = child.chromosome.back();
            T next_gene = T(*detail::index_of(is_used, false));

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
        assert(parent1.chromosome.size() == parent2.chromosome.size());
        assert(first <= last && last <= parent1.chromosome.size());

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

        /* The genes have to be unique unsigned integers in the range [0, chrom_len). */
        assert(parent1.chromosome.size() == parent2.chromosome.size());
        assert(first <= last && last <= parent1.chromosome.size());
        assert(*std::min_element(parent1.chromosome.begin(), parent1.chromosome.end()) == 0);
        assert(*std::max_element(parent1.chromosome.begin(), parent1.chromosome.end()) == chrom_len - 1);

        Candidate child = parent2;

        std::vector is_direct(chrom_len, false);
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
    
} // namespace genetic_algorithm::crossover::dtl

#endif // !GA_CROSSOVER_DTL_HPP