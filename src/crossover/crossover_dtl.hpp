/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CROSSOVER_DTL_HPP
#define GA_CROSSOVER_DTL_HPP

#include "../population/candidate.hpp"

namespace genetic_algorithm::crossover::dtl
{
    template<Gene T, size_t N>
    CandidatePair<T> nPointCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2);

} // namespace genetic_algorithm::crossover::dtl


/* IMPLEMENTATION */

#include "../utility/rng.hpp"
#include "../utility/algorithm.hpp"
#include <vector>
#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <cassert>

namespace genetic_algorithm::crossover::dtl
{
    template<Gene T, size_t N>
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

} // namespace genetic_algorithm::crossover::dtl

#endif // !GA_CROSSOVER_DTL_HPP