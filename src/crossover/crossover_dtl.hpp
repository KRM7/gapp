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

    /* Get the neighbours of each gene. */
    template<Gene T, bool Wrap = false>
    std::unordered_map<T, std::vector<T>> getNeighbourLists(const Chromosome<T>& chrom1, const Chromosome<T>& chrom2);

    /* Find the minimum number of neighbours of gene's neighbours. */
    template<typename T>
    size_t minNeighbourCount(const std::unordered_map<T, std::vector<T>>& neighbour_lists, const T& gene);

    template<Gene T>
    Candidate<T> edgeCrossoverImpl(const Candidate<T>& parent1, std::unordered_map<T, std::vector<T>>&& neighbour_lists);

    template<typename T>
    std::vector<std::vector<T>> findCycles(Chromosome<T> chrom1, Chromosome<T> chrom2);

    template<Gene T>
    Candidate<T> order2CrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, size_t first, size_t last);

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
                    cur_pos = detail::index_of(parent2.chromosome, gene);
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
            for (size_t i = 0; i < cycle.size(); i++)
            {
                chrom1.erase(std::find(chrom1.begin(), chrom1.end(), cycle[i]));
                chrom2.erase(std::find(chrom2.begin(), chrom2.end(), cycle[i]));
            }
            cycles.push_back(cycle);
        }

        return cycles;
    }

    template<Gene T>
    Candidate<T> order2CrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, size_t first, size_t last)
    {
        /* Genes in the range are copied directly: parent1 -> child */
        std::unordered_set<T> direct;
        for (size_t i = first; i != last; i++)
        {
            direct.insert(parent1.chromosome[i]);
        }

        /* The rest of the genes are copied from the other parent: parent2 -> child */
        auto cross = detail::find_all_v(parent2.chromosome.begin(), parent2.chromosome.end(),
        [&direct](const T& gene)
        {
            return !direct.contains(gene);
        });

        /* Construct the child: child = direct + cross */
        Candidate<size_t> child;
        child.chromosome.reserve(parent1.chromosome.size());

        std::move(cross.begin(), cross.begin() + first, std::back_inserter(child.chromosome));
        std::copy(parent1.chromosome.begin() + first, parent1.chromosome.begin() + last, std::back_inserter(child.chromosome));
        std::move(cross.begin() + first, cross.end(), std::back_inserter(child.chromosome));

        return child;
    }

} // namespace genetic_algorithm::crossover::dtl

#endif // !GA_CROSSOVER_DTL_HPP