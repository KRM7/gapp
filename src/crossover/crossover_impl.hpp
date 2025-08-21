/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_CROSSOVER_CROSSOVER_IMPL_HPP
#define GAPP_CROSSOVER_CROSSOVER_IMPL_HPP

#include "../core/candidate.hpp"
#include "../utility/functional.hpp"
#include "../utility/small_vector.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <vector>
#include <concepts>
#include <utility>
#include <cstddef>

namespace gapp::crossover::impl
{
    template<typename T>
    CandidatePair<T> nPointCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, small_vector<size_t> crossover_points)
    {
        const size_t chrom_len = parent1.chromosome.size();

        GAPP_ASSERT(parent1.chromosome.size() == parent2.chromosome.size());
        GAPP_ASSERT(std::all_of(crossover_points.begin(), crossover_points.end(), detail::between(0_sz, chrom_len)));

        std::sort(crossover_points.begin(), crossover_points.end());
        if (crossover_points.size() % 2) crossover_points.push_back(chrom_len);

        Candidate child1{ parent2 }, child2{ parent1 };

        for (size_t i = 1; i < crossover_points.size(); i += 2)
        {
            for (size_t j = crossover_points[i - 1]; j < crossover_points[i]; j++)
            {
                using std::swap;
                swap(child1.chromosome[j], child2.chromosome[j]);
            }
        }

        return { std::move(child1), std::move(child2) };
    }

    template<typename T>
    CandidatePair<T> singlePointCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, size_t crossover_point)
    {
        GAPP_ASSERT(crossover_point <= parent1.chromosome.size());
        GAPP_ASSERT(parent1.chromosome.size() == parent2.chromosome.size());

        Candidate child1{ parent1 }, child2{ parent2 };

        for (size_t i = 0; i < crossover_point; i++)
        {
            using std::swap;
            swap(child1.chromosome[i], child2.chromosome[i]);
        }

        return { std::move(child1), std::move(child2) };
    }

    template<typename T>
    CandidatePair<T> twoPointCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, std::pair<size_t, size_t> crossover_points)
    {
        GAPP_ASSERT(crossover_points.first <= parent1.chromosome.size());
        GAPP_ASSERT(crossover_points.second <= parent1.chromosome.size());
        GAPP_ASSERT(parent1.chromosome.size() == parent2.chromosome.size());

        if (crossover_points.first > crossover_points.second)
        {
            std::swap(crossover_points.first, crossover_points.second);
        }

        Candidate child1{ parent1 }, child2{ parent2 };

        for (size_t i = crossover_points.first; i < crossover_points.second; i++)
        {
            using std::swap;
            swap(child1.chromosome[i], child2.chromosome[i]);
        }

        return { std::move(child1), std::move(child2) };
    }
    
} // namespace gapp::crossover::impl

#endif // !GAPP_CROSSOVER_CROSSOVER_IMPL_HPP
