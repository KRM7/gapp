/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "binary.hpp"
#include "../population/candidate.hpp"
#include "../utility/rng.hpp"
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <cassert>
#include <cstddef>

namespace genetic_algorithm::crossover::binary
{
    CandidatePair<char> nPointCrossoverImpl(const Candidate<char>& parent1, const Candidate<char>& parent2, size_t n)
    {
        assert(n > 0);

        if (parent1.chromosome.size() != parent2.chromosome.size())
        {
            throw std::invalid_argument("The parent chromosomes must be the same length for the n-point crossover.");
        }

        std::vector<size_t> loci = rng::sampleUnique(parent1.chromosome.size(), std::min(n, parent1.chromosome.size()));

        /* Count how many loci are after each gene. */
        std::vector<size_t> loci_after;
        loci_after.reserve(parent1.chromosome.size());

        for (size_t i = 0, loci_left = loci.size(); i < parent1.chromosome.size(); i++)
        {
            if (loci_left > 0 && std::find(loci.begin(), loci.end(), i) != loci.end())
            {
                loci_left--;
            }
            loci_after.push_back(loci_left);
        }

        Candidate child1{ parent1 }, child2{ parent2 };

        for (size_t i = 0; i < parent1.chromosome.size(); i++)
        {
            if (loci_after[i] % 2)
            {
                child1.chromosome[i] = parent2.chromosome[i];
                child2.chromosome[i] = parent1.chromosome[i];
            }
        }

        return { child1, child2 };
    }

    CandidatePair<char> SinglePoint::crossover(const GaInfo&, const Candidate<char>& parent1, const Candidate<char>& parent2) const
    {
        return nPointCrossoverImpl(parent1, parent2, 1U);
    }

    CandidatePair<char> TwoPoint::crossover(const GaInfo&, const Candidate<char>& parent1, const Candidate<char>& parent2) const
    {
        return nPointCrossoverImpl(parent1, parent2, 2U);
    }

    CandidatePair<char> Uniform::crossover(const GaInfo&, const Candidate<char>& parent1, const Candidate<char>& parent2) const
    {
        if (parent1.chromosome.size() != parent2.chromosome.size())
        {
            throw std::invalid_argument("The parent chromosomes must be the same length for the uniform crossover.");
        }

        Candidate child1{ parent1 }, child2{ parent2 };

        for (size_t i = 0; i < parent1.chromosome.size(); i++)
        {
            /* Swap each gene with 0.5 probability. */
            if (rng::randomBool())
            {
                child1.chromosome[i] = parent2.chromosome[i];
                child2.chromosome[i] = parent1.chromosome[i];
            }
        }

        return { child1, child2 };
    }

} // namespace genetic_algorithm::crossover::binary