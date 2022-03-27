/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CANDIDATE_H
#define GA_CANDIDATE_H

#include "../utility/math.hpp"
#include "../utility/concepts.hpp"

#include <vector>
#include <cstddef>

namespace genetic_algorithm
{
    /** Valid gene types in the genetic algorithms. */
    template<typename T>
    concept Gene = requires
    {
        requires detail::Hashable<T>;
        requires std::regular<T>;
        requires std::destructible<T>;
        requires std::three_way_comparable<T>;
    };

    /**
    * The Candidate class that is used to represent solutions in the genetic algorithms.
    * This is used as the general purpose candidate type in all of the algorithms (SOGA, NSGA-II, NSGA-III).
    */
    template<Gene T>
    struct Candidate
    {
        using GeneType = T;

        std::vector<GeneType> chromosome;   /**< The chromosome encoding the solution. */
        std::vector<double> fitness;        /**< The fitness values (for each objective) of the solution. */

        bool is_evaluated = false;          /**< False if the candidate's fitness value needs to be computed. */

        Candidate() = default;
        explicit Candidate(const std::vector<GeneType>& chrom) : chromosome(chrom) {}
    };

    /** A pair of candidates. */
    template<Gene T>
    using CandidatePair = std::pair<Candidate<T>, Candidate<T>>;

    /** Two candidates are considered equal if their chromosomes are the same. */
    template<typename T>
    inline bool operator==(const Candidate<T>& lhs, const Candidate<T>& rhs);

    /** Two candidates are considered not equal if their chromosomes are different. */
    template<typename T>
    inline bool operator!=(const Candidate<T>& lhs, const Candidate<T>& rhs);

    /** Hash function for the Candidate so they can be stored in an unordered set/map. */
    template<detail::Hashable T>
    struct CandidateHasher
    {
        size_t operator()(const Candidate<T>& candidate) const noexcept;
    };

} // namespace genetic_algorithm


/* IMPLEMENTATION */

#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>

namespace genetic_algorithm
{
    template<typename GeneType>
    inline bool operator==(const Candidate<GeneType>& lhs, const Candidate<GeneType>& rhs)
    {
        if constexpr (std::is_floating_point_v<GeneType>)
        {
            return std::equal(lhs.chromosome.begin(), lhs.chromosome.end(), rhs.chromosome.begin(),
            [](const GeneType& lhs, const GeneType& rhs)
            {
                return detail::floatIsEqual(lhs, rhs);
            });
        }
        else
        {
            return lhs.chromosome == rhs.chromosome;
        }
    }

    template<typename GeneType>
    inline bool operator!=(const Candidate<GeneType>& lhs, const Candidate<GeneType>& rhs)
    {
        return !(lhs == rhs);
    }

    template<detail::Hashable GeneType>
    inline size_t CandidateHasher<GeneType>::operator()(const Candidate<GeneType>& candidate) const noexcept
    {
        size_t seed = candidate.chromosome.size();
        for (const auto& gene : candidate.chromosome)
        {
            seed ^= std::hash<GeneType>()(gene) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }

} // namespace genetic_algorithm

#endif // !GA_CANDIDATE_H