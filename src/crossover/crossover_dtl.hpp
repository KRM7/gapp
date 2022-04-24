/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CROSSOVER_DTL_HPP
#define GA_CROSSOVER_DTL_HPP

#include "../population/candidate.hpp"
#include <unordered_map>
#include <vector>
#include <cstddef>

namespace genetic_algorithm::crossover::dtl
{
    template<size_t N, Gene T>
    CandidatePair<T> nPointCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2);

    template<Gene T>
    Candidate<T> pmxCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2);

    template<Gene T, bool Wrap = false>
    std::unordered_map<T, std::vector<T>> getNeighbourLists(const Chromosome<T>& chrom1, const Chromosome<T>& chrom2);

    /* Find the minimum number of neighbours of gene's neighbours. */
    template<typename T>
    size_t minNeighbourCount(const std::unordered_map<T, std::vector<T>>& neighbour_lists, const T& gene);

    template<Gene T>
    Candidate<T> edgeCrossoverImpl(const Candidate<T>& parent1, std::unordered_map<T, std::vector<T>>&& neighbour_lists);

} // namespace genetic_algorithm::crossover::dtl


/* IMPLEMENTATION */

#include "../utility/rng.hpp"
#include "../utility/algorithm.hpp"
#include <unordered_set>
#include <algorithm>
#include <stdexcept>
#include <cassert>

namespace genetic_algorithm::crossover::dtl
{
    template<size_t N, Gene T>
    requires (N > 0)
    CandidatePair<T> nPointCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2)
    {
        if (parent1.chromosome.size() != parent2.chromosome.size())
        {
            throw std::invalid_argument("The parent chromosomes must be the same length for the n-point crossover.");
        }

        size_t chrom_len = parent1.chromosome.size();
        size_t num_crossover_points = std::min(N, chrom_len);

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
    Candidate<T> pmxCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2)
    {
        /* Init the child so the last step of the crossover can be skipped later. */
        Candidate child{ parent2 };

        /* Pick a random range [first, last) of genes (never empty or the entire chromosome). */
        size_t len = rng::randomInt(size_t{ 1 }, parent1.chromosome.size() - 1);
        size_t first = rng::randomInt(size_t{ 0 }, parent1.chromosome.size() - len);
        size_t last = first + len;

        /* Copy genes in the range from the corresponding parent. */
        std::unordered_set<T> copied_genes;
        for (size_t i = first; i < last; i++)
        {
            child.chromosome[i] = parent1.chromosome[i];
            copied_genes.insert(parent1.chromosome[i]);
        }

        /* Get the rest of the child's genes from the other parent. */
        for (size_t i = first; i < last; i++)
        {
            if (!copied_genes.contains(parent2.chromosome[i]))
            {
                size_t cur_pos = i;
                while (first <= cur_pos && cur_pos < last)
                {
                    T gene = parent1.chromosome[cur_pos];
                    cur_pos = size_t(std::find(parent2.chromosome.begin(), parent2.chromosome.end(), gene) - parent2.chromosome.begin());
                }
                child.chromosome[cur_pos] = parent2.chromosome[i];
            }
        }
        /* Copy any not yet in child positions to the child from the other parent. (Already done at the initialization of the child.) */

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

            /* Next gene */
            std::vector<T> candidate_genes;
            if (neighbour_lists[gene].empty())
            {
                candidate_genes = remaining_genes;
            }
            else
            {
                size_t min_neighbour_count = minNeighbourCount(neighbour_lists, gene);

                candidate_genes = detail::find_all_v(neighbour_lists[gene].begin(), neighbour_lists[gene].end(),
                [&neighbour_lists, min_neighbour_count](const T& val)
                {
                    return neighbour_lists[val].size() == min_neighbour_count;
                });

            }
            gene = rng::randomElement(candidate_genes);
        }

        return child;
    }

} // namespace genetic_algorithm::crossover::dtl

#endif // !GA_CROSSOVER_DTL_HPP