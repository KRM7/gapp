/* Copyright (c) 2022 Kriszti�n Rug�si. Subject to the MIT License. */

#ifndef GA_CROSSOVER_DTL_HPP
#define GA_CROSSOVER_DTL_HPP

#include "../population/candidate.hpp"
#include <unordered_map>
#include <vector>
#include <cstddef>

namespace genetic_algorithm::crossover::dtl
{
    template<Gene T>
    CandidatePair<T> nPointCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, size_t n);

    template<Gene T>
    Candidate<T> order1CrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, size_t first, size_t last);

    template<Gene T>
    Candidate<T> order2CrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, size_t first, size_t last);

    template<Gene T>
    Candidate<T> positionCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, const std::vector<size_t>& indices);

    template<typename T>
    std::vector<std::vector<T>> findCycles(Chromosome<T> chrom1, Chromosome<T> chrom2);

    template<Gene T>
    Candidate<T> edgeCrossoverImpl(const Candidate<T>& parent1, std::unordered_map<T, std::vector<T>>&& neighbour_lists);

    /* Get the neighbours of each gene. */
    template<Gene T, bool Wrap = false>
    std::unordered_map<T, std::vector<T>> getNeighbourLists(const Chromosome<T>& chrom1, const Chromosome<T>& chrom2);

    /* Find the minimum number of neighbours of gene's neighbours. */
    template<typename T>
    size_t minNeighbourCount(const std::unordered_map<T, std::vector<T>>& neighbour_lists, const T& gene);

    template<Gene T>
    Candidate<T> pmxCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2);

} // namespace genetic_algorithm::crossover::dtl


/* IMPLEMENTATION */

#include "../utility/rng.hpp"
#include "../utility/algorithm.hpp"
#include <unordered_set>
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <cassert>

namespace genetic_algorithm::crossover::dtl
{
    template<Gene T>
    CandidatePair<T> nPointCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, size_t n)
    {
        assert(n);

        size_t chrom_len = parent1.chromosome.size();

        if (parent2.chromosome.size() != chrom_len)
        {
            throw std::invalid_argument("The parent chromosomes must be the same length for the n-point crossover.");
        }

        size_t num_crossover_points = std::min(n, chrom_len);

        std::vector<size_t> crossover_points = rng::sampleUnique(chrom_len, num_crossover_points);

        std::vector<size_t> crossover_mask;
        crossover_mask.reserve(chrom_len);

        for (size_t i = 0, remaining = num_crossover_points; i < chrom_len; i++)
        {
            if (remaining && detail::contains(crossover_points.begin(), crossover_points.end(), i))
            {
                remaining--;
            }
            crossover_mask.push_back(remaining);
        }

        Candidate child1{ parent1 }, child2{ parent2 };

        for (size_t i = 0; i < chrom_len; i++)
        {
            if (crossover_mask[i] % 2)
            {
                child1.chromosome[i] = parent2.chromosome[i];
                child2.chromosome[i] = parent1.chromosome[i];
            }
        }

        return { child1, child2 };
    }

    template<Gene T>
    Candidate<T> order1CrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, size_t first, size_t last)
    {
        size_t chrom_len = parent1.chromosome.size();

        std::unordered_set<T> direct;
        for (size_t i = first; i != last; i++)
        {
            direct.insert(parent1.chromosome[i]);
        }

        Candidate<T> child{ parent1 };

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
        std::unordered_set<T> direct;
        for (size_t i = first; i != last; i++)
        {
            direct.insert(parent1.chromosome[i]);
        }

        Candidate<T> child{ parent1 };

        size_t child_pos = (first != 0) ? 0 : last;
        for (size_t parent_pos = 0; parent_pos < parent2.chromosome.size(); parent_pos++)
        {
            if (!direct.contains(parent2.chromosome[parent_pos]))
            {
                child.chromosome[child_pos] = parent2.chromosome[parent_pos];
                if (++child_pos == first)
                {
                    child_pos = last;
                }
            }
        }

        return child;
    }

    template<Gene T>
    Candidate<T> positionCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, const std::vector<size_t>& indices)
    {
        std::unordered_set<T> direct;
        std::unordered_set<size_t> direct_indices;
        for (const auto& idx : indices)
        {
            direct.insert(parent1.chromosome[idx]);
            direct_indices.insert(idx);
        }

        Candidate<T> child{ parent1 };

        size_t child_pos = 0;
        while (direct_indices.contains(child_pos)) child_pos++;
        for (size_t parent_pos = 0; parent_pos < parent2.chromosome.size(); parent_pos++)
        {
            if (!direct.contains(parent2.chromosome[parent_pos]))
            {
                child.chromosome[child_pos] = parent2.chromosome[parent_pos];
                while (direct_indices.contains(++child_pos));
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
            std::vector<T> cycle;

            size_t cur_pos = 0;
            T top_val = chrom1[cur_pos];
            cycle.push_back(top_val);

            while (chrom2[cur_pos] != chrom1[0])
            {
                cur_pos = detail::index_of(chrom1, chrom2[cur_pos]);
                top_val = chrom1[cur_pos];
                cycle.push_back(top_val);
            }

            /* Delete the values in this cycle from chrom1 and chrom2 without changing the order of the remaining genes. */
            for (const auto& gene : cycle)
            {
                detail::erase_first_v(chrom1, gene);
                detail::erase_first_v(chrom2, gene);
            }
            cycles.push_back(cycle);
        }

        return cycles;
    }

    template<Gene T>
    Candidate<T> edgeCrossoverImpl(const Candidate<T>& parent1, std::unordered_map<T, std::vector<T>>&& neighbour_lists)
    {
        size_t chrom_len = parent1.chromosome.size();

        Candidate<T> child;
        child.chromosome.reserve(chrom_len);

        std::vector<T> remaining_genes = parent1.chromosome;

        T gene = parent1.chromosome[0];
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
            std::vector<T> candidate_genes;
            if (neighbour_lists[gene].empty())
            {
                candidate_genes = remaining_genes;
            }
            else
            {
                size_t n = minNeighbourCount(neighbour_lists, gene);

                candidate_genes = detail::find_all_v(neighbour_lists[gene].begin(), neighbour_lists[gene].end(),
                [&neighbour_lists, n](const T& val)
                {
                    return neighbour_lists[val].size() == n;
                });

            }
            gene = rng::randomElement(candidate_genes);
        }

        return child;
    }

    template<Gene T, bool Wrap>
    std::unordered_map<T, std::vector<T>> getNeighbourLists(const Chromosome<T>& chrom1, const Chromosome<T>& chrom2)
    {
        std::unordered_map<T, std::vector<T>> neighbour_list;

        size_t len = chrom1.size();

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
        T nb = *std::min_element(neighbour_lists.at(gene).begin(), neighbour_lists.at(gene).end(),
        [&neighbour_lists](const T& lhs, const T& rhs)
        {
            return (neighbour_lists.at(lhs).size() < neighbour_lists.at(rhs).size());
        });

        return neighbour_lists.at(nb).size();
    }

    template<Gene T>
    Candidate<T> pmxCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2)
    {

        Candidate child{ parent2 };

        size_t chrom_len = parent1.chromosome.size();
        size_t range_len = rng::randomInt<size_t>(1, chrom_len - 1);
        size_t first = rng::randomInt<size_t>(0, chrom_len - range_len);
        size_t last = first + range_len;

        std::unordered_set<T> direct;
        for (size_t i = first; i < last; i++)
        {
            child.chromosome[i] = parent1.chromosome[i];
            direct.insert(parent1.chromosome[i]);
        }

        for (size_t i = first; i < last; i++)
        {
            if (!direct.contains(parent2.chromosome[i]))
            {
                size_t cur_pos = i;
                while (first <= cur_pos && cur_pos < last)
                {
                    T gene = parent1.chromosome[cur_pos];
                    cur_pos = detail::index_of(parent2.chromosome, gene);
                }
                child.chromosome[cur_pos] = parent2.chromosome[i];
            }
        }

        return child;
    }

} // namespace genetic_algorithm::crossover::dtl

#endif // !GA_CROSSOVER_DTL_HPP