/*
*  MIT License
*
*  Copyright (c) 2021 Krisztián Rugási
*
*  Permission is hereby granted, free of charge, to any person obtaining a copy
*  of this softwareand associated documentation files (the "Software"), to deal
*  in the Software without restriction, including without limitation the rights
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*  copies of the Software, and to permit persons to whom the Software is
*  furnished to do so, subject to the following conditions:
*
*  The above copyright noticeand this permission notice shall be included in all
*  copies or substantial portions of the Software.
*
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
*  SOFTWARE.
*/

/**
* This file contains the Candidate class used in the algorithms.
*
* @file candidate.h
*/

#ifndef GA_CANDIDATE_H
#define GA_CANDIDATE_H


#include <vector>
#include <functional>
#include <concepts>
#include <cstddef>

namespace genetic_algorithm
{
    /** Types that std::hash is specialized for. */
    template<typename T>
    concept hashable = requires(T arg)
    {
        { std::hash<T>{}(arg) } -> std::convertible_to<size_t>;
    };

    /** Types that are regular and std::hash is specialized for. */
    template<typename T>
    concept regular_hashable =  std::regular<T> && hashable<T>;


    /**
    * The Candidate class that is used to represent solutions in the genetic algorithms.
    * This is used as the general purpose candidate type in all of the algorithms (SOGA, NSGA-II, NSGA-III).
    */
    template<regular_hashable GeneType>
    struct Candidate
    {
        using Gene = GeneType;

        std::vector<Gene> chromosome;       /**< The chromosome encoding the solution. */
        std::vector<double> fitness;        /**< The fitness values (for each objective) of the solution. */

        double selection_pdf = 0.0;         /**< The probability of selecting the candidate (SOGA). [0.0, 1.0] */
        double selection_cdf = 0.0;         /**< The value of the cumulative distribution function for the candidate (SOGA). [0.0, 1.0] */

        size_t rank = 0;                    /**< Non-domination rank (used in both the NSGA-II and NSGA-III). */
        double distance = 0.0;              /**< Crowding distance (NSGA-II), or the distance to closest reference point (NSGA-III). */
        size_t ref_idx = 0;                 /**< Index of the associated reference point (NSGA-III). */
        size_t niche_count = 0;             /**< Number of candidates associated with the same reference point as this candidate (NSGA-III). */

        bool is_evaluated = false;          /**< False if the candidate's fitness value needs to be computed. */

        Candidate() = default;
        Candidate(const std::vector<Gene>& chrom) : chromosome(chrom) {}
    };

    /** Two candidates are considered equal if they have the same chromosomes. */
    template<typename T>
    inline bool operator==(const Candidate<T>& lhs, const Candidate<T>& rhs);

    /** Two candidates are considered not equal if they have different chromosomes. */
    template<typename T>
    inline bool operator!=(const Candidate<T>& lhs, const Candidate<T>& rhs);

    /** Hash function for the Candidate so they can be stored in an unordered set/map. */
    template<hashable T>
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
            std::equal(lhs.chromosome.begin(), lhs.chromosome.end(), rhs.chromosome.begin(),
            [](const GeneType& lhs, const GeneType& rhs)
            {
                return std::abs(lhs - rhs) <= std::numeric_limits<double>::epsilon() * std::max(std::abs(lhs), std::abs(rhs));
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

    template<hashable GeneType>
    inline size_t CandidateHasher<GeneType>::operator()(const Candidate<GeneType>& candidate) const noexcept
    {
        size_t seed = c.chromosome.size();
        for (const auto& gene : c.chromosome)
        {
            seed ^= std::hash<GeneType>()(gene) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }

} // namespace genetic_algorithm

#endif // !GA_CANDIDATE_H