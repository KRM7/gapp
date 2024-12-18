/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_CORE_CANDIDATE_HPP
#define GAPP_CORE_CANDIDATE_HPP

#include "../encoding/gene_types.hpp"
#include "../utility/math.hpp"
#include "../utility/small_vector.hpp"
#include "../utility/matrix.hpp"
#include "../utility/functional.hpp"
#include "../utility/concepts.hpp"
#include "../utility/iterators.hpp"
#include <vector>
#include <span>
#include <algorithm>
#include <utility>
#include <concepts>
#include <cstddef>

namespace gapp
{
    /**
    * The class used to represent the fitness of the candidates. 
    * Contains a fitness value for each objective.
    */
    using FitnessVector = small_vector<double>;

    /** 
    * The class used to represent the fitness values of multiple candidates.
    * Each row of the matrix is the fitness vector of a candidate.
    * 
    * The size of a fitness matrix is: [ number_of_candidates x number_of_objectives ].
    */
    using FitnessMatrix = detail::Matrix<double>;

    /**
    * The class used to represent the constraint violations of a candidate. It contains
    * the degree of constraint violation for each constraint separately. Higher values
    * mean greater degrees of constraint violation, while 0 or lower values mean no
    * constraint violation for a given constraint.
    * 
    * The size of a constraint violation vector is equal to the number of constraints.
    */
    using CVVector = small_vector<double, 2>;

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
    * A view of a bounds vector, containing the bounds for each gene of a chromosome.
    *
    * @tparam T The gene type.
    */
    template<typename T>
    using BoundsView = std::span<const Bounds<T>>;

    /**
    * The type used to represent the chromosome of a candidate solution.
    * Every gene of a chromosome is the same type.
    * 
    * @tparam T The gene type.
    */
    template<typename T>
    using Chromosome = std::vector<T>;


    /**
    * The base class used for the representation of candidate solutions in all of the GAs.
    * This class contains the parts of the candidates that are independent of the encoding type.
    * The derived Candidate class contains the encoding dependent parts in addition to these.
    *
    * @tparam T The gene type used in the candidate's chromosome.
    */
    struct CandidateInfo
    {
        /** The solution's fitness value for each objective. */
        FitnessVector fitness;

        /** The solution's degree of constraint violation for each constraint. */
        CVVector constraint_violation;


        /** 
        * @returns The number of objectives associated with the candidate.
        *    The return value will be 0 if called before the candidate has been evaluated.
        *    Equivalent to calling fitness.size().
        */
        size_t num_objectives() const noexcept
        {
            return fitness.size();
        }

        /**
        * @returns True if the candidate has a valid fitness vector associated with it.
        *    Equivalent to calling !fitness.empty().
        */
        bool is_evaluated() const noexcept 
        {
            return !fitness.empty();
        }

        /** 
        * @returns The number of constraints associated with the candidate.
        *    The return value will be 0 if called before the constraints have been evaluated.
        *    Equivalent to calling constraint_violation.size().
        */
        size_t num_constraints() const noexcept
        {
            return constraint_violation.size();
        }

        /**
        * @returns True if the candidate violates any of its constraints. If there are no
        *    constraints associated with the candidate, it will also return false.
        *    Returns false if called before the constraints have been evaluated.
        */
        bool has_constraint_violation() const noexcept
        {
            return std::any_of(constraint_violation.begin(), constraint_violation.end(), detail::greater_than(0.0));
        }

    protected:

        CandidateInfo()                                 = default;
        CandidateInfo(const CandidateInfo&)             = default;
        CandidateInfo(CandidateInfo&&)                  = default;
        CandidateInfo& operator=(const CandidateInfo&)  = default;
        CandidateInfo& operator=(CandidateInfo&&)       = default;
        ~CandidateInfo()                                = default;
    };

    namespace detail
    {
        template<typename T>
        struct CandidateBounds {};

        template<typename T>
        requires(is_bounded<T>)
        struct CandidateBounds<T> { BoundsView<T> bounds_; };

    } // namespace detail

    /**
    * The class that is used to represent candidate solutions in all of the algorithms.
    * 
    * Contains several helper methods to access the chromosome directly, such as iterator
    * methods. The behaviour of all these helper functions is equivalent to calling the
    * matching functions on the chromosome member.
    * 
    * @tparam T The gene type used in the candidate's chromosome.
    */
    template<typename T>
    struct Candidate : public CandidateInfo, private detail::CandidateBounds<T>, public detail::container_interface<Candidate<T>>
    {
        using Gene = T; /**< The type of the candidate's genes. */

        /**
        * Create a candidate with an empty fitness vector and chromosome, with
        * the specified gene bounds.
        * 
        * @param bounds The lower and upper bounds of each gene of the chromosome.
        */
        explicit Candidate(BoundsView<T> bounds) requires(is_bounded<T>) :
            detail::CandidateBounds<T>(bounds),
            chromosome(bounds.size())
        {}

        /**
        * Create a candidate with an empty fitness vector and a given chromosome.
        *
        * @param chrom The chromosome of the candidate.
        */
        explicit Candidate(Chromosome<T> chrom) requires(!is_bounded<T>) :
            chromosome(std::move(chrom))
        {}

        /**
        * Create a candidate with an empty fitness vector, using the given chromosome and
        * bounds.
        * 
        * @param chrom The chromosome of the candidate.
        * @param bounds The bounds of each gene of the chromosome.
        */
        Candidate(Chromosome<T> chrom, BoundsView<T> bounds) requires(is_bounded<T>) :
            detail::CandidateBounds<T>(bounds),
            chromosome(std::move(chrom))
        {}

        Candidate()                             = default;
        Candidate(const Candidate&)             = default;
        Candidate(Candidate&&)                  = default;
        Candidate& operator=(const Candidate&)  = default;
        Candidate& operator=(Candidate&&)       = default;
        ~Candidate()                            = default;


        /** @returns An iterator to the beginning of the candidate's chromosome. */
        auto begin() noexcept { return chromosome.begin(); }
        auto begin() const noexcept { return chromosome.begin(); }

        /** @returns An iterator to the end of the candidate's chromosome. */
        auto end() noexcept { return chromosome.end(); }
        auto end() const noexcept { return chromosome.end(); }


        /** @returns The lower and upper bounds for each gene of the candidate's chromosome. */
        BoundsView<T> gene_bounds() const noexcept requires(is_bounded<T>)
        {
            GAPP_ASSERT(this->bounds_.size() == chromosome.size(), "Attempting to access invalid gene bounds.");
            return this->bounds_;
        }

        /** @returns The length of the chromosome. Equivalent to calling chromosome.size(). */
        size_t chrom_len() const noexcept
        {
            return chromosome.size();
        }

        /** The chromosome encoding the solution. */
        Chromosome<T> chromosome;
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
