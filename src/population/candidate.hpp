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

    /** The chromosome type of the Candidates. */
    template<Gene T>
    using Chromosome = std::vector<T>;

    /**
    * The Candidate class that is used to represent solutions in the genetic algorithms.
    * This is used as the general purpose candidate type in all of the algorithms (SOGA, NSGA-II, NSGA-III).
    */
    template<Gene T>
    struct Candidate
    {
        using GeneType = T;
        using ChromosomeType = Chromosome<T>;

        Chromosome<T> chromosome;       /**< The chromosome encoding the solution. */
        std::vector<double> fitness;    /**< The fitness values (for each objective) of the solution. */

        bool is_evaluated = false;      /**< False if the candidate's fitness value needs to be computed. */

        Candidate() = default;
        explicit Candidate(const Chromosome<T>& chrom) : chromosome(chrom) {}
    };

    /** A pair of candidates. */
    template<Gene T>
    using CandidatePair = std::pair<Candidate<T>, Candidate<T>>;

    /** Two candidates are considered equal if their chromosomes are the same. */
    template<typename T>
    bool operator==(const Candidate<T>& lhs, const Candidate<T>& rhs);

    /** Two candidates are considered not equal if their chromosomes are different. */
    template<typename T>
    bool operator!=(const Candidate<T>& lhs, const Candidate<T>& rhs);

    template<typename T>
    bool operator<(const Candidate<T>& lhs, const Candidate<T>& rhs);

    template<typename T>
    bool operator<=(const Candidate<T>& lhs, const Candidate<T>& rhs);

    template<typename T>
    bool operator>(const Candidate<T>& lhs, const Candidate<T>& rhs);

    template<typename T>
    bool operator>=(const Candidate<T>& lhs, const Candidate<T>& rhs);

    /** Hash function for the Candidate. */
    template<detail::Hashable T>
    struct CandidateHasher
    {
        size_t operator()(const Candidate<T>& candidate) const noexcept;
    };

} // namespace genetic_algorithm


/* IMPLEMENTATION */

#include <algorithm>
#include <functional>
#include <type_traits>

namespace genetic_algorithm
{
    template<typename T>
    bool operator==(const Candidate<T>& lhs, const Candidate<T>& rhs)
    {
        if constexpr (std::is_floating_point_v<T>)
        {
            return detail::floatVecIsEqual(lhs.chromosome, rhs.chromosome);
        }
        else
        {
            return lhs.chromosome == rhs.chromosome;
        }
    }

    template<typename T>
    bool operator!=(const Candidate<T>& lhs, const Candidate<T>& rhs)
    {
        return !(lhs == rhs);
    }

    template<typename T>
    bool operator<(const Candidate<T>& lhs, const Candidate<T>& rhs)
    {
        if constexpr (std::is_floating_point_v<T>)
        {
            return std::lexicographical_compare(lhs.chromosome.begin(), lhs.chromosome.end(),
                                                rhs.chromosome.begin(), rhs.chromosome.end(), 
                                                detail::floatIsEqual);
        }
        else
        {
            return lhs.chromosome < rhs.chromosome;
        }
    }

    template<typename T>
    bool operator>=(const Candidate<T>& lhs, const Candidate<T>& rhs)
    {
        return !(lhs < rhs);
    }

    template<typename T>
    bool operator>(const Candidate<T>& lhs, const Candidate<T>& rhs)
    {
        return rhs < lhs;
    }

    template<typename T>
    bool operator<=(const Candidate<T>& lhs, const Candidate<T>& rhs)
    {
        return !(rhs < lhs);
    }

    template<detail::Hashable GeneType>
    size_t CandidateHasher<GeneType>::operator()(const Candidate<GeneType>& candidate) const noexcept
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