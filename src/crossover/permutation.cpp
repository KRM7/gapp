/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "permutation.hpp"
#include "crossover_base.hpp"
#include "neighbour_list.hpp"
#include "../core/candidate.hpp"
#include "../utility/algorithm.hpp"
#include "../utility/functional.hpp"
#include "../utility/small_vector.hpp"
#include "../utility/dynamic_bitset.hpp"
#include "../utility/rng.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <vector>
#include <utility>
#include <cstddef>

namespace gapp::crossover::impl
{
    [[maybe_unused]] static bool isValidIntegerPermutation(const Chromosome<PermutationGene>& chrom)
    {
        if (chrom.empty()) return true;

        if (*std::min_element(chrom.begin(), chrom.end()) != 0) return false;
        if (*std::max_element(chrom.begin(), chrom.end()) != chrom.size() - 1) return false;

        dynamic_bitset present(chrom.size());
        for (const PermutationGene& val : chrom)
        {
            if (present[val]) return false;
            present[val] = true;
        }

        return true;
    }

    Candidate<PermutationGene> order1CrossoverImpl(const Candidate<PermutationGene>& parent1, const Candidate<PermutationGene>& parent2, size_t first, size_t last)
    {
        const size_t chrom_len = parent1.chromosome.size();
        const size_t range_len = last - first;

        GAPP_ASSERT(first <= last && last <= chrom_len);
        GAPP_ASSERT(parent1.chromosome.size() == parent2.chromosome.size());
        GAPP_ASSERT(isValidIntegerPermutation(parent1.chromosome));
        GAPP_ASSERT(isValidIntegerPermutation(parent2.chromosome));

        dynamic_bitset is_direct(chrom_len);
        for (size_t idx = first; idx != last; idx++)
        {
            is_direct[parent1.chromosome[idx]] = true;
        }

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

    Candidate<PermutationGene> order2CrossoverImpl(const Candidate<PermutationGene>& parent1, const Candidate<PermutationGene>& parent2, size_t first, size_t last)
    {
        const size_t chrom_len = parent1.chromosome.size();

        GAPP_ASSERT(first <= last && last <= chrom_len);
        GAPP_ASSERT(parent1.chromosome.size() == parent2.chromosome.size());
        GAPP_ASSERT(isValidIntegerPermutation(parent1.chromosome));
        GAPP_ASSERT(isValidIntegerPermutation(parent2.chromosome));

        dynamic_bitset is_direct(chrom_len);
        for (size_t idx = first; idx != last; idx++)
        {
            is_direct[parent1.chromosome[idx]] = true;
        }

        Candidate child = parent1;

        for (size_t child_pos = 0; PermutationGene gene : parent2.chromosome)
        {
            if (!is_direct[gene])
            {
                if (child_pos == first) child_pos = last; // skip [first, last)
                child.chromosome[child_pos++] = gene;
            }
        }

        return child;
    }

    Candidate<PermutationGene> positionCrossoverImpl(const Candidate<PermutationGene>& parent1, const Candidate<PermutationGene>& parent2, std::span<const size_t> indices)
    {
        const size_t chrom_len = parent1.chromosome.size();

        GAPP_ASSERT(std::all_of(indices.begin(), indices.end(), detail::between(0_sz, chrom_len - 1)));
        GAPP_ASSERT(parent1.chromosome.size() == parent2.chromosome.size());
        GAPP_ASSERT(isValidIntegerPermutation(parent1.chromosome));
        GAPP_ASSERT(isValidIntegerPermutation(parent2.chromosome));

        dynamic_bitset is_direct(chrom_len);
        for (size_t idx : indices)
        {
            is_direct[parent1.chromosome[idx]] = true;
        }

        small_vector<size_t> next_indirect(chrom_len);
        for (ptrdiff_t indirect = -1, i = chrom_len - 1; i >= 0; i--)
        {
            const PermutationGene gene = parent1.chromosome[i];

            indirect = is_direct[gene] ? indirect : i;
            next_indirect[i] = indirect;
        }

        Candidate child = parent1;

        for (size_t child_pos = 0; PermutationGene gene : parent2.chromosome)
        {
            if (!is_direct[gene])
            {
                child_pos = next_indirect[child_pos];
                child.chromosome[child_pos++] = gene;
            }
        }

        return child;
    }

    static std::vector<size_t> findOddCycleIndices(const Chromosome<PermutationGene>& chrom1, const Chromosome<PermutationGene>& chrom2)
    {
        GAPP_ASSERT(chrom1.size() == chrom2.size());
        GAPP_ASSERT(isValidIntegerPermutation(chrom1));
        GAPP_ASSERT(isValidIntegerPermutation(chrom2));

        const size_t chrom_len = chrom1.size();

        std::vector<size_t> odd_indices;
        odd_indices.reserve(chrom_len / 2);

        dynamic_bitset deleted(chrom_len);
        size_t num_deleted = 0;

        std::vector index_lookup(chrom_len, 0_sz);
        for (size_t i = 0; i < chrom_len; i++)
        {
            index_lookup[chrom1[i]] = i;
        }

        for (bool odd_cycle = false; num_deleted < chrom_len; odd_cycle ^= 1)
        {
            size_t pos = deleted.find_first(false);
            const PermutationGene cycle_start = chrom1[pos];

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

    CandidatePair<PermutationGene> cycleCrossoverImpl(const Candidate<PermutationGene>& parent1, const Candidate<PermutationGene>& parent2)
    {
        GAPP_ASSERT(parent1.chromosome.size() == parent2.chromosome.size());

        const auto odd_cycle_idxs = findOddCycleIndices(parent1.chromosome, parent2.chromosome);

        Candidate child1 = parent1;
        Candidate child2 = parent2;

        for (size_t idx : odd_cycle_idxs)
        {
            using std::swap;
            swap(child1.chromosome[idx], child2.chromosome[idx]);
        }

        return { std::move(child1), std::move(child2) };
    }

    Candidate<PermutationGene> edgeCrossoverImpl(const Candidate<PermutationGene>& parent1, const Candidate<PermutationGene>& parent2)
    {
        const size_t chrom_len = parent1.chromosome.size();

        GAPP_ASSERT(parent1.chromosome.size() == parent2.chromosome.size());
        GAPP_ASSERT(isValidIntegerPermutation(parent1.chromosome));
        GAPP_ASSERT(isValidIntegerPermutation(parent2.chromosome));

        auto nb_lists = makeNeighbourLists(parent1.chromosome, parent2.chromosome);

        Candidate<PermutationGene> child({ parent1.chromosome[0] });
        child.chromosome.reserve(chrom_len);

        dynamic_bitset is_used(chrom_len);
        is_used[parent1.chromosome[0]] = true;

        while (child.chromosome.size() != chrom_len)
        {
            PermutationGene last_gene = child.chromosome.back();
            PermutationGene next_gene = is_used.find_first(false);

            for (PermutationGene neighbour : nb_lists[last_gene])
            {
                if (neighbour == NeighbourList::EMPTY) continue;

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

    Candidate<PermutationGene> pmxCrossoverImpl(const Candidate<PermutationGene>& parent1, const Candidate<PermutationGene>& parent2, size_t first, size_t last)
    {
        const size_t chrom_len = parent1.chromosome.size();

        GAPP_ASSERT(parent1.chromosome.size() == parent2.chromosome.size());
        GAPP_ASSERT(first <= last && last <= parent1.chromosome.size());
        GAPP_ASSERT(isValidIntegerPermutation(parent1.chromosome));
        GAPP_ASSERT(isValidIntegerPermutation(parent2.chromosome));

        Candidate child = parent2;

        dynamic_bitset is_direct(chrom_len);
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

} // namespace gapp::crossover::impl

namespace gapp::crossover::perm
{
    using namespace gapp::crossover::impl;

    auto Order1::crossover(const GaInfo&, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const -> CandidatePair<GeneType>
    {
        GAPP_ASSERT(parent1.chromosome.size() == parent2.chromosome.size(), "Mismatching parent chromosome lengths.");

        const size_t chrom_len = parent1.chromosome.size();

        if (chrom_len < 2) return { parent1, parent2 };

        const size_t length = rng::randomInt(1_sz, chrom_len - 1_sz);
        const size_t first = rng::randomInt(0_sz, chrom_len - length);
        const size_t last = first + length;

        auto child1 = order1CrossoverImpl(parent1, parent2, first, last);
        auto child2 = order1CrossoverImpl(parent2, parent1, first, last);

        return { std::move(child1), std::move(child2) };
    }

    auto Order2::crossover(const GaInfo&, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const -> CandidatePair<GeneType>
    {
        GAPP_ASSERT(parent1.chromosome.size() == parent2.chromosome.size(), "Mismatching parent chromosome lengths.");

        const size_t chrom_len = parent1.chromosome.size();

        if (chrom_len < 2) return { parent1, parent2 };

        const size_t length = rng::randomInt(1_sz, chrom_len - 1_sz);
        const size_t first = rng::randomInt(0_sz, chrom_len - length);
        const size_t last = first + length;

        auto child1 = order2CrossoverImpl(parent1, parent2, first, last);
        auto child2 = order2CrossoverImpl(parent2, parent1, first, last);

        return { std::move(child1), std::move(child2) };
    }

    auto Position::crossover(const GaInfo&, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const -> CandidatePair<GeneType>
    {
        GAPP_ASSERT(parent1.chromosome.size() == parent2.chromosome.size(), "Mismatching parent chromosome lengths.");

        const size_t chrom_len = parent1.chromosome.size();

        if (chrom_len < 2) return { parent1, parent2 };

        const size_t ns = rng::randomInt(1_sz, chrom_len - 1_sz);
        const auto idxs = rng::sampleUnique(0_sz, chrom_len, ns);

        auto child1 = positionCrossoverImpl(parent1, parent2, idxs);
        auto child2 = positionCrossoverImpl(parent2, parent1, idxs);

        return { std::move(child1), std::move(child2) };
    }

    auto Cycle::crossover(const GaInfo&, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const -> CandidatePair<GeneType>
    {
        GAPP_ASSERT(parent1.chromosome.size() == parent2.chromosome.size(), "Mismatching parent chromosome lengths.");

        const size_t chrom_len = parent1.chromosome.size();

        if (chrom_len < 2) return { parent1, parent2 };

        return cycleCrossoverImpl(parent1, parent2);
    }

    auto Edge::crossover(const GaInfo&, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const -> CandidatePair<GeneType>
    {
        GAPP_ASSERT(parent1.chromosome.size() == parent2.chromosome.size(), "Mismatching parent chromosome lengths.");

        const size_t chrom_len = parent1.chromosome.size();

        if (chrom_len < 2) return { parent1, parent2 };

        auto child1 = edgeCrossoverImpl(parent1, parent2);
        auto child2 = edgeCrossoverImpl(parent2, parent1);

        return { std::move(child1), std::move(child2) };
    }

    auto PMX::crossover(const GaInfo&, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const -> CandidatePair<GeneType>
    {
        GAPP_ASSERT(parent1.chromosome.size() == parent2.chromosome.size(), "Mismatching parent chromosome lengths.");

        const size_t chrom_len = parent1.chromosome.size();
        
        if (chrom_len < 2) return { parent1, parent2 };

        const size_t range_len = rng::randomInt(1_sz, chrom_len - 1);

        const size_t first = rng::randomInt(0_sz, chrom_len - range_len);
        const size_t last = first + range_len;

        auto child1 = pmxCrossoverImpl(parent1, parent2, first, last);
        auto child2 = pmxCrossoverImpl(parent2, parent1, first, last);

        return { std::move(child1), std::move(child2) };
    }

} // namespace gapp::crossover::perm
