/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "permutation.hpp"
#include "../utility/rng.hpp"

#include <algorithm>
#include <unordered_set>

namespace genetic_algorithm::crossover::perm
{
    CandidatePair<size_t> Order1::crossover(const GA<size_t>&, const Candidate<size_t>& parent1, const Candidate<size_t>& parent2) const
    {
        if (parent1.chromosome.size() != parent2.chromosome.size())
        {
            throw std::invalid_argument("The parent chromosomes must be the same length for the Order1 crossover.");
        }
        if (parent1.chromosome.size() < 2)
        {
            throw std::invalid_argument("The parent chromosomes must have at least 2 genes for the Order1 crossover.");
        }

        /* Pick a random range [first, last) of genes (never empty or the entire chromosome). */
        size_t len = rng::randomInt(size_t{ 1 }, parent1.chromosome.size() - 1);
        size_t first = rng::randomInt(size_t{ 0 }, parent1.chromosome.size() - len);
        size_t last = first + len;

        /* Gather the genes that will go directly from parent1 -> child1, and parent2 -> child2. (Not using the constructor is intentional.) */
        std::unordered_set<size_t> direct1, direct2;
        for (auto gene = parent1.chromosome.begin() + first; gene != parent1.chromosome.begin() + last; gene++) direct1.insert(*gene);
        for (auto gene = parent2.chromosome.begin() + first; gene != parent2.chromosome.begin() + last; gene++) direct2.insert(*gene);

        /* Gather the remaining genes (not in the range) from the other parent. */
        std::vector<size_t> cross1;    /* Segment gathered from parent2 -> child1. */
        std::vector<size_t> cross2;    /* Segment gathered from parent1 -> child2. */
        cross1.reserve(parent2.chromosome.size() - direct1.size());
        cross2.reserve(parent1.chromosome.size() - direct2.size());

        size_t pos = last;
        if (pos == parent1.chromosome.size())
        {
            pos = 0;
        }
        while (cross1.size() < (parent2.chromosome.size() - direct1.size()) ||
               cross2.size() < (parent1.chromosome.size() - direct2.size()))
        {
            /* If a gene is not taken directly from the corresponding parent, then take it from the other parent. */
            if (!direct1.contains(parent2.chromosome[pos])) cross1.push_back(parent2.chromosome[pos]);
            if (!direct2.contains(parent1.chromosome[pos])) cross2.push_back(parent1.chromosome[pos]);

            if (++pos == parent1.chromosome.size())
            {
                pos = 0;
            }
        }

        /* Construct the children: child1 = direct1 + cross1, child2 = direct2 + cross2 */
        Candidate<size_t> child1{ parent1 }, child2{ parent2 };

        size_t chrom_pos = last;
        if (chrom_pos == parent1.chromosome.size())
        {
            chrom_pos = 0;
        }
        for (size_t i = 0; i < cross1.size(); i++)
        {
            child1.chromosome[chrom_pos] = cross1[i];
            child2.chromosome[chrom_pos] = cross2[i];
            if (++chrom_pos == parent1.chromosome.size())
            {
                chrom_pos = 0;
            }
        }

        return { child1, child2 };
    }

    CandidatePair<size_t> Order2::crossover(const GA<size_t>&, const Candidate<size_t>& parent1, const Candidate<size_t>& parent2) const
    {
        if (parent1.chromosome.size() != parent2.chromosome.size())
        {
            throw std::invalid_argument("The parent chromosomes must be the same length for the Order2 crossover.");
        }
        if (parent1.chromosome.size() < 2)
        {
            throw std::invalid_argument("The parent chromosomes must have at least 2 genes for the Order2 crossover.");
        }

        /* Pick a random range [first, last) of genes (never empty or the entire chromosome). */
        size_t len = rng::randomInt(size_t{ 1 }, parent1.chromosome.size() - 1);
        size_t first = rng::randomInt(size_t{ 0 }, parent1.chromosome.size() - len);
        size_t last = first + len;

        /* Gather the genes that will go directly from parent1 -> child1, and parent2 -> child2. (Not using the constructor is intentional.) */
        std::unordered_set<size_t> direct1, direct2;
        for (auto gene = parent1.chromosome.begin() + first; gene != parent1.chromosome.begin() + last; gene++) direct1.insert(*gene);
        for (auto gene = parent2.chromosome.begin() + first; gene != parent2.chromosome.begin() + last; gene++) direct2.insert(*gene);

        /* Gather the remaining genes (not in the range) from the other parent. */
        std::vector<size_t> cross1;    /* Segment gathered from parent2 -> child1. */
        std::vector<size_t> cross2;    /* Segment gathered from parent1 -> child2. */
        cross1.reserve(parent2.chromosome.size() - direct1.size());
        cross2.reserve(parent1.chromosome.size() - direct2.size());
        for (size_t i = 0; i < parent1.chromosome.size(); i++)
        {
            /* If a gene is not taken directly from the corresponding parent, then take it from the other parent. */
            if (!direct1.contains(parent2.chromosome[i])) cross1.push_back(parent2.chromosome[i]);
            if (!direct2.contains(parent1.chromosome[i])) cross2.push_back(parent1.chromosome[i]);
        }

        /* Construct the children: child1 = direct1 + cross1, child2 = direct2 + cross2 */
        Candidate<size_t> child1, child2;
        child1.chromosome.reserve(parent1.chromosome.size());
        child2.chromosome.reserve(parent2.chromosome.size());
        /* child1 */
        child1.chromosome.insert(child1.chromosome.end(), cross1.begin(), cross1.begin() + first);
        child1.chromosome.insert(child1.chromosome.end(), parent1.chromosome.begin() + first, parent1.chromosome.begin() + last);
        child1.chromosome.insert(child1.chromosome.end(), cross1.begin() + first, cross1.end());
        /* child2 */
        child2.chromosome.insert(child2.chromosome.end(), cross2.begin(), cross2.begin() + first);
        child2.chromosome.insert(child2.chromosome.end(), parent2.chromosome.begin() + first, parent2.chromosome.begin() + last);
        child2.chromosome.insert(child2.chromosome.end(), cross2.begin() + first, cross2.end());

        return { child1, child2 };
    }

    CandidatePair<size_t> Position::crossover(const GA<size_t>&, const Candidate<size_t>& parent1, const Candidate<size_t>& parent2) const
    {
        if (parent1.chromosome.size() != parent2.chromosome.size())
        {
            throw std::invalid_argument("The parent chromosomes must be the same length for the Position crossover.");
        }

        Candidate child1{ parent1 }, child2{ parent2 };

        /* Determine directly copied indices (never directly copy 0 or every gene). */
        std::vector<size_t> idxs = rng::sampleUnique(parent1.chromosome.size(), rng::randomInt(size_t{ 1 }, parent1.chromosome.size() - 1));

        std::unordered_set<size_t> direct_idxs;
        for (const auto& idx : idxs) direct_idxs.insert(idx);

        std::unordered_set<size_t> direct1, direct2;
        for (const auto& idx : idxs)
        {
            direct1.insert(parent1.chromosome[idx]);
            direct2.insert(parent2.chromosome[idx]);
        }

        std::vector<size_t> cross_idxs1, cross_idxs2;
        cross_idxs1.reserve(parent1.chromosome.size() - direct_idxs.size());
        cross_idxs2.reserve(parent1.chromosome.size() - direct_idxs.size());
        for (size_t i = 0; i < parent1.chromosome.size(); i++)
        {
            if (!direct1.contains(parent2.chromosome[i])) cross_idxs1.push_back(i);
            if (!direct2.contains(parent1.chromosome[i])) cross_idxs2.push_back(i);
        }

        for (size_t i = 0, first = 0; i < parent1.chromosome.size(); i++)
        {
            if (direct_idxs.contains(i))
            {
                child1.chromosome[i] = parent1.chromosome[i];
                child2.chromosome[i] = parent2.chromosome[i];
            }
            else
            {
                child1.chromosome[i] = parent2.chromosome[cross_idxs1[first]];
                child2.chromosome[i] = parent1.chromosome[cross_idxs2[first]];
                first++;
            }
        }

        return { child1, child2 };
    }

    CandidatePair<size_t> Cycle::crossover(const GA<size_t>&, const Candidate<size_t>& parent1, const Candidate<size_t>& parent2) const
    {
        if (parent1.chromosome.size() != parent2.chromosome.size())
        {
            throw std::invalid_argument("The parent chromosomes must be the same length for the Cycle crossover.");
        }

        /* Identify all cycles. */
        std::vector<std::vector<size_t>> cycles;
        std::vector<size_t> chrom1(parent1.chromosome), chrom2(parent2.chromosome);

        while (!chrom1.empty())
        {
            std::vector<size_t> cycle;
            /* Always start the cycle at chrom1[0]. Old cycles are removed the chromosomes. */
            size_t cur_pos = 0;
            size_t top_val = chrom1[cur_pos];
            cycle.push_back(top_val);

            while (chrom2[cur_pos] != chrom1[0])
            {
                /* Look for the bottom value at this pos (chrom2[pos]) in chrom1. This is the new pos and top value. */
                cur_pos = static_cast<size_t>(std::find(chrom1.begin(), chrom1.end(), chrom2[cur_pos]) - chrom1.begin());
                top_val = chrom1[cur_pos];
                /* Add the new top value to the cycle. */
                cycle.push_back(top_val);
                /* Keep going until the bottom value at this pos isn't the cycle start value. (cycle complete) */
            }

            /* Delete the values in this cycle from chrom1 and chrom2. (Without changing the order of the remaining genes.) */
            for (size_t i = 0; i < cycle.size(); i++)
            {
                chrom1.erase(std::find(chrom1.begin(), chrom1.end(), cycle[i]));
                chrom2.erase(std::find(chrom2.begin(), chrom2.end(), cycle[i]));
            }
            /* Add this cycle to the cycles. */
            cycles.push_back(cycle);
        }

        /* Construct the children from the cycles. */
        Candidate child1{ parent1 }, child2{ parent2 };

        for (size_t i = 0; i < parent1.chromosome.size(); i++)
        {
            /* Find which cycle has the gene. */
            size_t cycle_num = 0;
            for (size_t j = 0; j < cycles.size(); j++)
            {
                if (std::find(cycles[j].begin(), cycles[j].end(), parent1.chromosome[i]) != cycles[j].end())
                {
                    cycle_num = j + 1;
                    break;
                }
            }
            /* Even cycle genes are swapped parent1->child2 and parent2->child1. */
            if (cycle_num % 2 == 0)
            {
                child1.chromosome[i] = parent2.chromosome[i];
                child2.chromosome[i] = parent1.chromosome[i];
            }
            /* Odd cycle genes were already handled when initializing the children. */
        }

        return { child1, child2 };
    }

    CandidatePair<size_t> Edge::crossover(const GA<size_t>&, const Candidate<size_t>& parent1, const Candidate<size_t>& parent2) const
    {
        using NList = std::vector<std::unordered_set<size_t>>;

        if (parent1.chromosome.size() != parent2.chromosome.size())
        {
            throw std::invalid_argument("The parent chromosomes must be the same length for the Edge crossover.");
        }
        if (parent1.chromosome.size() < 2)
        {
            throw std::invalid_argument("The parent chromosomes must have at least 2 genes for the Edge crossover.");
        }

        Candidate<size_t> child1, child2;
        child1.chromosome.reserve(parent1.chromosome.size());
        child2.chromosome.reserve(parent2.chromosome.size());

        /* Construct neighbour list based on parents. The first and last genes are not neighbours. */
        size_t len = parent1.chromosome.size();
        NList nl1(len);
        /* Neighbours of first and last genes. */
        nl1[parent1.chromosome.front()] = { parent1.chromosome[1] };
        nl1[parent1.chromosome.back()] = { parent1.chromosome[len - 2] };
        nl1[parent2.chromosome.front()].insert(parent2.chromosome[1]);
        nl1[parent2.chromosome.back()].insert(parent2.chromosome[len - 2]);
        /* Neighbours of all other genes. */
        for (size_t i = 1; i < len - 1; i++)
        {
            nl1[parent1.chromosome[i]].insert(parent1.chromosome[i + 1]);
            nl1[parent1.chromosome[i]].insert(parent1.chromosome[i - 1]);

            nl1[parent2.chromosome[i]].insert(parent2.chromosome[i + 1]);
            nl1[parent2.chromosome[i]].insert(parent2.chromosome[i - 1]);
        }
        NList nl2 = nl1;    /* Copy for child2. */

        /* Generate child1. */
        std::vector<size_t> remaining_genes = parent1.chromosome;
        size_t gene = parent1.chromosome[0];

        while (child1.chromosome.size() != parent1.chromosome.size())
        {
            /* Append gene to the child, and remove it from all neighbour lists. */
            child1.chromosome.push_back(gene);
            remaining_genes.erase(std::remove(remaining_genes.begin(), remaining_genes.end(), gene), remaining_genes.end());
            for (auto& neighbours : nl1)
            {
                neighbours.erase(gene);
            }
            if (child1.chromosome.size() == parent1.chromosome.size()) break;

            /* Determine next gene that will be added to the child. */
            /* If gene's neighbour list is empty, gene = random node not already in child. */
            if (nl1[gene].empty())
            {
                gene = rng::randomElement(remaining_genes);
            }
            else /* gene's neighbour list is not empty, gene = neighbour of gene with fewest neighbours (random if tie). */
            {
                /* Find gene's neighbour with fewest neighbours. */
                size_t nb = *std::min_element(nl1[gene].begin(), nl1[gene].end(),
                [&nl1](const size_t& lhs, const size_t& rhs)
                {
                    return (nl1[lhs].size() < nl1[rhs].size());
                });
                size_t min_neighbour_count = nl1[nb].size();

                /* Determine possible nodes (neighbours of gene with min_neighbour_count neighbours). */
                std::vector<size_t> possible_nodes;
                for (const auto& neighbour : nl1[gene])
                {
                    if (nl1[neighbour].size() == min_neighbour_count)
                    {
                        possible_nodes.push_back(neighbour);
                    }
                }

                gene = rng::randomElement(possible_nodes);
            }
        }

        /* Same process to get child2. */
        remaining_genes = parent2.chromosome;
        gene = parent2.chromosome[0];
        while (child2.chromosome.size() != parent2.chromosome.size())
        {
            child2.chromosome.push_back(gene);
            remaining_genes.erase(std::remove(remaining_genes.begin(), remaining_genes.end(), gene), remaining_genes.end());
            for (auto& neighbours : nl2)
            {
                neighbours.erase(gene);
            }
            if (child2.chromosome.size() == parent2.chromosome.size()) break;

            if (nl2[gene].empty())
            {
                gene = rng::randomElement(remaining_genes);
            }
            else
            {
                size_t nb = *std::min_element(nl2[gene].begin(), nl2[gene].end(),
                [&nl2](const size_t& lhs, const size_t& rhs)
                {
                    return (nl2[lhs].size() < nl2[rhs].size());
                });
                size_t min_neighbour_count = nl2[nb].size();

                std::vector<size_t> possible_nodes;
                for (const auto& neighbour : nl2[gene])
                {
                    if (nl2[neighbour].size() == min_neighbour_count)
                    {
                        possible_nodes.push_back(neighbour);
                    }
                }

                gene = rng::randomElement(possible_nodes);
            }
        }

        return { child1, child2 };
    }

    CandidatePair<size_t> PMX::crossover(const GA<size_t>&, const Candidate<size_t>& parent1, const Candidate<size_t>& parent2) const
    {
        if (parent1.chromosome.size() != parent2.chromosome.size())
        {
            throw std::invalid_argument("The parent chromosomes must be the same length for the PMX crossover.");
        }
        if (parent1.chromosome.size() < 2)
        {
            throw std::invalid_argument("The parent chromosomes must have at least 2 genes for the PMX crossover.");
        }

        /* Init children so the last step of the crossover can be skipped. */
        Candidate child1{ parent2 }, child2{ parent1 };

        /* Pick a random range [first, last) of genes (never empty or the entire chromosome). */
        size_t len = rng::randomInt(size_t{ 1 }, parent1.chromosome.size() - 1);
        size_t first = rng::randomInt(size_t{ 0 }, parent1.chromosome.size() - len);
        size_t last = first + len;

        /* Copy genes in the range from the corresponding parent. */
        for (size_t i = first; i < last; i++)
        {
            child1.chromosome[i] = parent1.chromosome[i];
            child2.chromosome[i] = parent2.chromosome[i];
        }
        /* Note ranges that were copied from parents to check if they contain an element. (Not using the constructor is intentional.) */
        std::unordered_set<size_t> p1_range, p2_range;
        for (auto gene = parent1.chromosome.begin() + first; gene != parent1.chromosome.begin() + last; gene++) p1_range.insert(*gene);
        for (auto gene = parent2.chromosome.begin() + first; gene != parent2.chromosome.begin() + last; gene++) p2_range.insert(*gene);

        /* Get the rest of the children's genes from the other parents. */
        for (size_t i = first; i < last; i++)
        {
            /* child1 */
            if (!p1_range.contains(parent2.chromosome[i]))
            {
                size_t cur_pos = i;
                while (first <= cur_pos && cur_pos < last)
                {
                    size_t val = parent1.chromosome[cur_pos];
                    cur_pos = static_cast<size_t>(std::find(parent2.chromosome.begin(), parent2.chromosome.end(), val) - parent2.chromosome.begin());
                }
                child1.chromosome[cur_pos] = parent2.chromosome[i];
            }
            /* child2 */
            if (!p2_range.contains(parent1.chromosome[i]))
            {
                size_t cur_pos = i;
                while (first <= cur_pos && cur_pos < last)
                {
                    size_t val = parent2.chromosome[cur_pos];
                    cur_pos = static_cast<size_t>(std::find(parent1.chromosome.begin(), parent1.chromosome.end(), val) - parent1.chromosome.begin());
                }
                child2.chromosome[cur_pos] = parent1.chromosome[i];
            }
        }
        /* Copy any not yet in child positions to the children from the other parents. (Already done at the initialization of the children.) */

        return { child1, child2 };
    }

} // namespace genetic_algorithm::crossover::perm