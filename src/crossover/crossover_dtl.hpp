/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CROSSOVER_DTL_HPP
#define GA_CROSSOVER_DTL_HPP

#include "../population/candidate.hpp"
#include <unordered_map>
#include <vector>
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

    /* Find every cycle of genes in the chromosomes chrom1 and chrom2. Used in the cycle crossover operator. */
    template<typename T>
    std::vector<std::vector<T>> findCycles(Chromosome<T> chrom1, Chromosome<T> chrom2);

    /* Implementation of the cycle crossover for any gene type. */
    template<Gene T>
    CandidatePair<T> cycleCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2);

    /* Implementation of the edge crossover for any gene type, only generates a single child. */
    template<Gene T>
    Candidate<T> edgeCrossoverImpl(const Candidate<T>& parent1, std::unordered_map<T, std::vector<T>>&& neighbour_lists);

    /* Get the list of neighbours of each gene in the 2 chromosomes. (Wrap = true if the front and back elements are considered neighbours.) */
    template<Gene T, bool Wrap = false>
    std::unordered_map<T, std::vector<T>> getNeighbourLists(const Chromosome<T>& chrom1, const Chromosome<T>& chrom2);

    /* Find the neighbour of gene which has the least number of neighbours. */
    template<typename T>
    size_t minNeighbourCount(const std::unordered_map<T, std::vector<T>>& neighbour_lists, const T& gene);

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
        while (first != last) direct.insert(parent1.chromosome[first++]);

        Candidate<T> child = parent1;

        size_t parent_pos = last % chrom_len;
        size_t child_pos = last % chrom_len;
        for (size_t i = 0; i < chrom_len; i++)
        {
            if (!direct.contains(parent2.chromosome[parent_pos]))
            {
                child.chromosome[child_pos] = parent2.chromosome[parent_pos];
                child_pos = (child_pos + 1) % chrom_len;
            }
            parent_pos = (parent_pos + 1) % chrom_len;
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
                if (child_pos == first) child_pos = last;
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

        for (size_t child_pos = 0; const T& gene : parent2.chromosome)
        {
            if (!direct.contains(gene))
            {
                while (direct.contains(parent1.chromosome[child_pos])) child_pos++;
                child.chromosome[child_pos] = gene;
            }
        }

        return child;
    }

    template<typename T>
    std::vector<std::vector<T>> findCycles(Chromosome<T> chrom1, Chromosome<T> chrom2)
    {
        std::vector<std::vector<T>> cycles;

        while (!chrom1.empty())
        {
            /* Find a cycle. */
            size_t pos = 0;
            std::vector<T> cycle{ chrom1[pos] };

            while (chrom2[pos] != cycle.front())
            {
                pos = *detail::index_of(chrom1, chrom2[pos]);
                cycle.push_back(chrom1[pos]);
            }

            /* Delete this cycle from the chromosomes without changing the order of the remaining genes. */
            for (const auto& gene : cycle)
            {
                detail::erase_first_stable(chrom1, gene);
                detail::erase_first_stable(chrom2, gene);
            }
            cycles.push_back(cycle);
        }

        return cycles;
    }

    template<Gene T>
    CandidatePair<T> cycleCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2)
    {
        const size_t chrom_len = parent1.chromosome.size();

        const auto cycles = dtl::findCycles(parent1.chromosome, parent2.chromosome);

        Candidate child1 = parent1;
        Candidate child2 = parent2;

        for (size_t i = 0; i < chrom_len; i++)
        {
            const auto cycle_idx = detail::find_index(cycles, [&](const auto& cycle)
            {
                return detail::contains(cycle.begin(), cycle.end(), parent1.chromosome[i]);
            });

            if (*cycle_idx % 2)
            {
                using std::swap;
                swap(child1.chromosome[i], child2.chromosome[i]);
            }
            /* Even cycle idx genes were already handled when initializing the children. */
        }

        return { std::move(child1), std::move(child2) };
    }

    template<Gene T>
    Candidate<T> edgeCrossoverImpl(const Candidate<T>& parent1, std::unordered_map<T, std::vector<T>>&& neighbour_lists)
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
            for (auto& [_, neighbour_list] : neighbour_lists)
            {
                std::erase(neighbour_list, gene);
            }
            if (child.chromosome.size() == chrom_len) break;

            /* Get next gene */
            auto candidate_genes = detail::find_all_v(neighbour_lists[gene].begin(), neighbour_lists[gene].end(),
            [&neighbour_lists, n = dtl::minNeighbourCount(neighbour_lists, gene)](const T& val)
            {
                return neighbour_lists[val].size() == n;
            });

            gene = candidate_genes.empty() ?
                rng::randomElement(remaining_genes.begin(), remaining_genes.end()) :
                rng::randomElement(candidate_genes.begin(), candidate_genes.end());
        }

        return child;
    }

    template<Gene T, bool Wrap>
    std::unordered_map<T, std::vector<T>> getNeighbourLists(const Chromosome<T>& chrom1, const Chromosome<T>& chrom2)
    {
        std::unordered_map<T, std::vector<T>> neighbour_list(chrom1.size());

        const size_t len = chrom1.size();

        neighbour_list[chrom1.front()].push_back(chrom1[1]);
        neighbour_list[chrom1.back()].push_back(chrom1[len - 2]);
        if constexpr (Wrap)
        {
            neighbour_list[chrom1.front()].push_back(chrom1.back());
            neighbour_list[chrom1.back()].push_back(chrom1.front());
        }

        neighbour_list[chrom2.front()].push_back(chrom2[1]);
        neighbour_list[chrom2.back()].push_back(chrom2[len - 2]);
        if constexpr (Wrap)
        {
            neighbour_list[chrom2.front()].push_back(chrom2.back());
            neighbour_list[chrom2.back()].push_back(chrom2.front());
        }

        for (size_t i = 1; i < len - 1; i++)
        {
            neighbour_list[chrom1[i]].push_back(chrom1[i + 1]);
            neighbour_list[chrom1[i]].push_back(chrom1[i - 1]);

            neighbour_list[chrom2[i]].push_back(chrom2[i + 1]);
            neighbour_list[chrom2[i]].push_back(chrom2[i - 1]);
        }

        return neighbour_list;
    }

    template<typename T>
    size_t minNeighbourCount(const std::unordered_map<T, std::vector<T>>& neighbour_lists, const T& gene)
    {
        const auto& neighbour_list = neighbour_lists.at(gene);

        T nb = *std::min_element(neighbour_list.begin(), neighbour_list.end(),
        [&neighbour_lists](const T& lhs, const T& rhs)
        {
            return (neighbour_lists.at(lhs).size() < neighbour_lists.at(rhs).size());
        });

        return neighbour_lists.at(nb).size();
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