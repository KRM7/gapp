/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "permutation.hpp"
#include "crossover_dtl.hpp"
#include "../utility/rng.hpp"
#include <algorithm>
#include <unordered_set>
#include <unordered_map>

namespace genetic_algorithm::crossover::perm
{
    auto Order1::crossover(const GaInfo&, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const -> CandidatePair<GeneType>
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

    auto Order2::crossover(const GaInfo&, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const -> CandidatePair<GeneType>
    {
        if (parent1.chromosome.size() != parent2.chromosome.size())
        {
            throw std::invalid_argument("The parent chromosomes must be the same length for the Order2 crossover.");
        }
        if (parent1.chromosome.size() < 2)
        {
            throw std::invalid_argument("The parent chromosomes must have at least 2 genes for the Order2 crossover.");
        }

        size_t range_len = rng::randomInt<size_t>(1, parent1.chromosome.size() - 1);
        size_t first = rng::randomInt<size_t>(0, parent1.chromosome.size() - range_len);
        size_t last = first + range_len;

        auto child1 = dtl::order2CrossoverImpl(parent1, parent2, first, last);
        auto child2 = dtl::order2CrossoverImpl(parent2, parent1, first, last);

        return { child1, child2 };
    }

    auto Position::crossover(const GaInfo&, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const -> CandidatePair<GeneType>
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
            if (!direct_idxs.contains(i))
            {
                child1.chromosome[i] = parent2.chromosome[cross_idxs1[first]];
                child2.chromosome[i] = parent1.chromosome[cross_idxs2[first]];
                first++;
            }
        }

        return { child1, child2 };
    }

    auto Cycle::crossover(const GaInfo&, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const -> CandidatePair<GeneType>
    {
        if (parent1.chromosome.size() != parent2.chromosome.size())
        {
            throw std::invalid_argument("The parent chromosomes must be the same length for the Cycle crossover.");
        }

        auto cycles = dtl::findCycles(parent1.chromosome, parent2.chromosome);

        Candidate child1{ parent1 }, child2{ parent2 };

        for (size_t i = 0; i < parent1.chromosome.size(); i++)
        {
            size_t cycle_idx = detail::index_of(cycles,
            [&gene = parent1.chromosome[i]](const auto& cycle)
            {
                return detail::contains(cycle.begin(), cycle.end(), gene);
            });

            if (cycle_idx % 2)
            {
                child1.chromosome[i] = parent2.chromosome[i];
                child2.chromosome[i] = parent1.chromosome[i];
            }
            /* Even cycle genes were already handled when initializing the children. */
        }

        return { child1, child2 };
    }

    auto Edge::crossover(const GaInfo&, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const -> CandidatePair<GeneType>
    {
        if (parent1.chromosome.size() != parent2.chromosome.size())
        {
            throw std::invalid_argument("The parent chromosomes must be the same length for the Edge crossover.");
        }
        if (parent1.chromosome.size() < 2)
        {
            throw std::invalid_argument("The parent chromosomes must have at least 2 genes for the Edge crossover.");
        }

        auto nl1 = dtl::getNeighbourLists(parent1.chromosome, parent2.chromosome);
        auto nl2 = nl1;

        auto child1 = dtl::edgeCrossoverImpl(parent1, std::move(nl1));
        auto child2 = dtl::edgeCrossoverImpl(parent2, std::move(nl2));

        return { child1, child2 };
    }

    auto PMX::crossover(const GaInfo&, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const -> CandidatePair<GeneType>
    {
        if (parent1.chromosome.size() != parent2.chromosome.size())
        {
            throw std::invalid_argument("The parent chromosomes must be the same length for the PMX crossover.");
        }
        if (parent1.chromosome.size() < 2)
        {
            throw std::invalid_argument("The parent chromosomes must have at least 2 genes for the PMX crossover.");
        }

        auto child1 = dtl::pmxCrossoverImpl<GeneType>(parent1, parent2);
        auto child2 = dtl::pmxCrossoverImpl<GeneType>(parent2, parent1);

        return { child1, child2 };
    }

} // namespace genetic_algorithm::crossover::perm