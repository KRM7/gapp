/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CANDIDATE_H
#define GA_CANDIDATE_H

#include "../utility/math.hpp"
#include "../utility/concepts.hpp"
#include <vector>
#include <utility>
#include <concepts>
#include <cstddef>

namespace genetic_algorithm
{
     /** The type used to represent the lower and upper bounds of a gene. */
    template<typename T>
    class GeneBounds
    {
    public:
        /** Constructor for the range [lower, upper]. */
        constexpr GeneBounds(const T& lower, const T& upper) noexcept;

        /** @returns The lower gene bound. */
        [[nodiscard]]
        constexpr const T& lower() const noexcept { return lower_; }

        /** @returns The upper gene bound. */
        [[nodiscard]]
        constexpr const T& upper() const noexcept { return upper_; }

    private:
        T lower_;
        T upper_;
    };

    /** A vector of gene lower and upper bounds. */
    template<typename T>
    using BoundsVector = std::vector<GeneBounds<T>>;

    /** The chromosome type of the Candidates. */
    template<typename T>
    using Chromosome = std::vector<T>;

    /**
    * The Candidate class that is used to represent solutions in the genetic algorithms. \n
    * This is used as the candidate type in all of the algorithms.
    */
    template<typename T>
    struct Candidate
    {
        using Gene = T;

        explicit Candidate(size_t chrom_len) :
            chromosome(chrom_len)
        {}

        explicit Candidate(const Chromosome<T>& chrom) :
            chromosome(chrom)
        {}

        explicit Candidate(Chromosome<T>&& chrom) noexcept :
            chromosome(std::move(chrom))
        {}

        Candidate(std::initializer_list<T> chrom) :
            chromosome(std::move(chrom))
        {}

        Candidate()                             = default;
        Candidate(const Candidate&)             = default;
        Candidate(Candidate&&)                  = default;
        Candidate& operator=(const Candidate&)  = default;
        Candidate& operator=(Candidate&&)       = default;
        ~Candidate()                            = default;

        Chromosome<T> chromosome;       /**< The chromosome encoding the solution. */
        std::vector<double> fitness;    /**< The fitness values (for each objective) of the solution. */
        bool is_evaluated = false;      /**< False if the candidate's fitness value needs to be computed. */
    };

    /** A pair of candidates. */
    template<typename T>
    using CandidatePair = std::pair<Candidate<T>, Candidate<T>>;

    /** Two candidates are considered equal if their chromosomes are the same. */
    template<typename T>
    bool operator==(const Candidate<T>& lhs, const Candidate<T>& rhs) noexcept;

    /** Two candidates are not equal if their chromosomes are different. */
    template<typename T>
    bool operator!=(const Candidate<T>& lhs, const Candidate<T>& rhs) noexcept;

    /* Lexicographical comparison based on the chromosomes. */
    template<typename T>
    bool operator<(const Candidate<T>& lhs, const Candidate<T>& rhs) noexcept;

    /* Lexicographical comparison based on the chromosomes. */
    template<typename T>
    bool operator<=(const Candidate<T>& lhs, const Candidate<T>& rhs) noexcept;

    /* Lexicographical comparison based on the chromosomes. */
    template<typename T>
    bool operator>(const Candidate<T>& lhs, const Candidate<T>& rhs) noexcept;

    /* Lexicographical comparison based on the chromosomes. */
    template<typename T>
    bool operator>=(const Candidate<T>& lhs, const Candidate<T>& rhs) noexcept;

    /** Hash function for the Candidate. */
    template<detail::Hashable T>
    struct CandidateHasher
    {
        size_t operator()(const Candidate<T>& candidate) const noexcept;
    };

} // namespace genetic_algorithm


/* IMPLEMENTATION */

#include "utility/utility.hpp"
#include <algorithm>
#include <functional>
#include <type_traits>

namespace genetic_algorithm
{
    template<typename T>
    constexpr GeneBounds<T>::GeneBounds(const T& lower, const T& upper) noexcept :
        lower_(lower), upper_(upper)
    {
        GA_ASSERT(lower <= upper, "The lower bound can't be greater than the upper bound.");
    }

    template<typename T>
    bool operator==(const Candidate<T>& lhs, const Candidate<T>& rhs) noexcept
    {
        if constexpr (std::is_floating_point_v<T>)
        {
            return math::floatVecIsEqual<T>(lhs.chromosome, rhs.chromosome);
        }
        else
        {
            return lhs.chromosome == rhs.chromosome;
        }
    }

    template<typename T>
    bool operator!=(const Candidate<T>& lhs, const Candidate<T>& rhs) noexcept
    {
        return !(lhs == rhs);
    }

    template<typename T>
    bool operator<(const Candidate<T>& lhs, const Candidate<T>& rhs) noexcept
    {
        if constexpr (std::is_floating_point_v<T>)
        {
            return std::lexicographical_compare(lhs.chromosome.begin(), lhs.chromosome.end(),
                                                rhs.chromosome.begin(), rhs.chromosome.end(), 
                                                math::floatIsLess<T>);
        }
        else
        {
            return lhs.chromosome < rhs.chromosome;
        }
    }

    template<typename T>
    bool operator>=(const Candidate<T>& lhs, const Candidate<T>& rhs) noexcept
    {
        return !(lhs < rhs);
    }

    template<typename T>
    bool operator>(const Candidate<T>& lhs, const Candidate<T>& rhs) noexcept
    {
        return rhs < lhs;
    }

    template<typename T>
    bool operator<=(const Candidate<T>& lhs, const Candidate<T>& rhs) noexcept
    {
        return !(rhs < lhs);
    }

    template<detail::Hashable T>
    size_t CandidateHasher<T>::operator()(const Candidate<T>& candidate) const noexcept
    {
        size_t seed = candidate.chromosome.size();
        for (const auto& gene : candidate.chromosome)
        {
            seed ^= std::hash<T>{}(gene) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }

} // namespace genetic_algorithm

#endif // !GA_CANDIDATE_H