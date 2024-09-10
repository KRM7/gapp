/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_CORE_CANDIDATE_HPP
#define GAPP_CORE_CANDIDATE_HPP

#include "../utility/math.hpp"
#include "../utility/small_vector.hpp"
#include "../utility/matrix.hpp"
#include "../utility/concepts.hpp"
#include <vector>
#include <utility>
#include <concepts>
#include <cstddef>

namespace gapp
{
    /** The class used to represent the fitness of the candidates. Contains a fitness value for each objective. */
    using FitnessVector = small_vector<double>;

    /** 
    * The class used to represent the fitness values of multiple candidates.
    * Each row of the matrix is the fitness vector of a candidate.
    * 
    * The size of a fitness matrix is: [ number_of_candidates x number_of_objectives ].
    */
    using FitnessMatrix = detail::Matrix<double>;

    /**
    * The class used to represent the lower and upper bounds of a gene.
    * 
    * @tparam T The gene type. The lower and upper bounds will also be this type.
    */
    template<typename T>
    class Bounds
    {
    public:
        /** The type of the gene. */
        using Gene = T;

        /**
        * Constructor for the closed range [@p lower, @p upper].
        * 
        * @param lower The lower bound (inclusive).
        * @param upper The upper bound (inclusive).
        */
        constexpr Bounds(const T& lower, const T& upper) noexcept;

        /** @returns The lower bound (inclusive). */
        [[nodiscard]]
        constexpr const T& lower() const noexcept { return lower_; }

        /** @returns The upper bound (inclusive). */
        [[nodiscard]]
        constexpr const T& upper() const noexcept { return upper_; }

    private:
        T lower_;
        T upper_;
    };

    /**
    * A vector of lower and upper gene bounds, containing the bounds for each gene of a chromosome.
    * 
    * @tparam T The gene type.
    */
    template<typename T>
    using BoundsVector = std::vector<Bounds<T>>;

    /**
    * The type used to represent the chromosome of a candidate solution.
    * Every gene of a chromosome is the same type.
    * 
    * @tparam T The gene type.
    */
    template<typename T>
    using Chromosome = std::vector<T>;

    /**
    * The class that is used to represent candidate solutions in all of the algorithms.
    * 
    * @tparam T The gene type used in the candidate's chromosome.
    */
    template<typename T>
    struct Candidate
    {
        using Gene = T; /**< The type of the candidate's genes. */

        /**
        * Create a candidate with an empty fitness vector and specific chromosome size.
        * The genes of the chromosome will be default constructed.
        * 
        * @param chrom_len The length of the chromosome.
        */
        explicit Candidate(size_t chrom_len) :
            chromosome(chrom_len)
        {}

        /**
        * Create a candidate with an empty fitness vector and a given chromosome.
        *
        * @param chrom The chromosome of the candidate.
        */
        explicit Candidate(const Chromosome<T>& chrom) :
            chromosome(chrom)
        {}

        /**
        * Create a candidate with an empty fitness vector and a given chromosome.
        *
        * @param chrom The chromosome of the candidate.
        */
        explicit Candidate(Chromosome<T>&& chrom) noexcept :
            chromosome(std::move(chrom))
        {}

        /**
        * Create a candidate with an empty fitness vector from a list of genes.
        *
        * @param chrom A list of genes for the chromosome.
        */
        Candidate(std::initializer_list<T> chrom) :
            chromosome(std::move(chrom))
        {}

        Candidate()                             = default;
        Candidate(const Candidate&)             = default;
        Candidate(Candidate&&)                  = default;
        Candidate& operator=(const Candidate&)  = default;
        Candidate& operator=(Candidate&&)       = default;
        ~Candidate()                            = default;

        FitnessVector fitness;      /**< The fitness values of the solution (for every objective). */
        Chromosome<T> chromosome;   /**< The chromosome encoding the solution. */
        bool is_evaluated = false;  /**< True if the candidate's fitness value doesn't need to be computed. */
    };

    /** A pair of candidates. */
    template<typename T>
    struct CandidatePair
    {
        Candidate<T> first;
        Candidate<T> second;
    };


    /**
    * Comparison operators based on the chromosomes of the candidates.
    * The comparisons are not transitive if T is a floating-point type. 
    */
    template<typename T>
    bool operator==(const Candidate<T>& lhs, const Candidate<T>& rhs) noexcept;


    /* Hash function for the candidates. */
    template<detail::hashable T>
    struct CandidateHasher
    {
        size_t operator()(const Candidate<T>& candidate) const noexcept;
    };

} // namespace gapp


/* IMPLEMENTATION */

#include "../utility/utility.hpp"
#include <algorithm>
#include <functional>
#include <type_traits>

namespace gapp
{
    template<typename T>
    constexpr Bounds<T>::Bounds(const T& lower, const T& upper) noexcept :
        lower_(lower), upper_(upper)
    {
        GAPP_ASSERT(lower <= upper, "The lower bound can't be greater than the upper bound.");
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

    template<detail::hashable T>
    size_t CandidateHasher<T>::operator()(const Candidate<T>& candidate) const noexcept
    {
        size_t seed = candidate.chromosome.size();
        for (const auto& gene : candidate.chromosome)
        {
            seed ^= std::hash<T>{}(gene) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }

} // namespace gapp

namespace std
{
    template<typename T>
    struct hash<gapp::Candidate<T>>
    {
        std::size_t operator()(const gapp::Candidate<T>& candidate) const noexcept
        {
            return gapp::CandidateHasher<T>{}(candidate);
        }
    };

} // namespace std

#endif // !GAPP_CORE_CANDIDATE_HPP