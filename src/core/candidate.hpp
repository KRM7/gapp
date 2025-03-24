/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_CORE_CANDIDATE_HPP
#define GAPP_CORE_CANDIDATE_HPP

#include "../encoding/gene_types.hpp"
#include "../utility/math.hpp"
#include "../utility/hash.hpp"
#include "../utility/small_vector.hpp"
#include "../utility/matrix.hpp"
#include "../utility/functional.hpp"
#include "../utility/iterators.hpp"
#include "../utility/utility.hpp"
#include <vector>
#include <tuple>
#include <span>
#include <algorithm>
#include <ranges>
#include <compare>
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

        /** @returns true if the bounds are the same. */
        constexpr friend bool operator==(const Bounds&, const Bounds&) = default;

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

    private:
        virtual void enable_dynamic_downcast() const noexcept final { GAPP_UNREACHABLE(); }
    };

    namespace detail
    {
        template<typename T>
        struct CandidateBounds {};

        template<typename T>
        requires(is_bounded_gene_v<T>)
        struct CandidateBounds<T> { BoundsView<T> gene_bounds; };

    } // namespace detail

    /**
    * The class that is used to represent candidate solutions in the genetic algorithms.
    * 
    * Contains several helper methods to access the chromosome directly, such as iterator
    * methods. The behaviour of all these helper functions is equivalent to calling the
    * matching functions on the chromosome member.
    * 
    * @tparam T The gene type used in the candidate's chromosome.
    */
    template<typename T>
    struct Candidate : public virtual CandidateInfo, public detail::CandidateBounds<T>, public detail::container_interface<Candidate<T>> // NOLINT(*virtual-class-destructor)
    {
        /** The type of the candidate's genes. */
        using Gene = T;

        /** The number of chromosomes the candidate has. */
        static constexpr std::size_t NumChroms = 1;

        /** Create a candidate with an empty fitness vector and chromosome. */
        Candidate() = default;

        /**
        * Create a candidate with an empty fitness vector, empty constraints vector,
        * and the specified gene bounds.
        * 
        * @param bounds The lower and upper bounds of each gene of the chromosome.
        *  The size of the bounds should match the size of the chromosome later.
        */
        explicit Candidate(BoundsView<T> bounds) requires(is_bounded_gene_v<T>);

        /**
        * Create a candidate with an empty fitness vector and a given chromosome.
        *
        * @param chrom The chromosome of the candidate.
        */
        explicit Candidate(Chromosome<T> chrom) requires(!is_bounded_gene_v<T>);

        /**
        * Create a candidate with an empty fitness vector, using the given chromosome and
        * bounds.
        * 
        * @param chrom The chromosome of the candidate.
        * @param bounds The bounds of each gene of the chromosome.
        */
        Candidate(Chromosome<T> chrom, BoundsView<T> bounds) requires(is_bounded_gene_v<T>);


        /** @returns An iterator to the first element of the candidate's chromosome. */
        auto begin() noexcept { return chromosome.begin(); }
        auto begin() const noexcept { return chromosome.begin(); }

        /** @returns An iterator to the end of the candidate's chromosome. */
        auto end() noexcept { return chromosome.end(); }
        auto end() const noexcept { return chromosome.end(); }


        /** @returns The length of the chromosome. */
        template<typename GeneType = T>
        size_t chrom_len() const noexcept;

        /** @returns The chromosome of the candidate. */
        template<typename GeneType = T>
        Chromosome<GeneType>& chrom() & noexcept;

        /** @returns The chromosome of the candidate. */
        template<typename GeneType = T>
        const Chromosome<GeneType>& chrom() const& noexcept;

        /** @returns The lower and upper bounds of each gene of the chromosome. */
        template<typename GeneType = T>
        BoundsView<GeneType> bounds() const noexcept requires(is_bounded_gene_v<T>);


        /** The chromosome encoding the solution. */
        Chromosome<T> chromosome;
    };

    /**
    * The class used to represent candidate solutions in the mixed encoded genetic algorithms.
    * It contains a separate chromosome for each component gene type of the mixed gene.
    *
    * @tparam Ts The component gene types of the mixed gene.
    */
    template<typename... Ts>
    struct Candidate<MixedGene<Ts...>> final : public Candidate<Ts>...
    {
        /** The type of the candidate's genes. */
        using Gene = MixedGene<Ts...>;

        /** The number of chromosomes the candidate has. Equal to the number of component genes. */
        static constexpr std::size_t NumChroms = sizeof...(Ts);

        /** Create a candidate with an empty fitness vector and chromosomes. */
        Candidate() = default;

        /**
        * Create a candidate with the specified chromosomes. The fitness and constraint
        * vectors of the candidate will be empty.
        * 
        * @param chroms The chromosomes of the candidate solution.
        */
        explicit Candidate(std::tuple<Chromosome<Ts>...> chroms) requires(!is_partially_bounded_gene_v<MixedGene<Ts...>>) :
            Candidate<Ts>(std::move(std::get<Chromosome<Ts>>(chroms)))...
        {}

        /**
        * Create a candidate with the specified chromosomes and bounds for each bounded
        * component gene type. The bounds are only specified for the bounded component
        * genes.
        * The fitness and constraint vectors of the candidate will be empty.
        * 
        * @param chroms The chromosomes of the candidate solution.
        * @param bounds The bounds vectors for each of the bounded chromosomes.
        */
        Candidate(std::tuple<Chromosome<Ts>...> chroms, detail::map_types_t<BoundsView, detail::filter_types_t<is_bounded_gene, Ts...>> bounds)
        requires(is_partially_bounded_gene_v<MixedGene<Ts...>>) :
            Candidate(constructor_tag_t{}, std::move(chroms), bounds)
        {}

        /**
        * Create a candidate from a set of candidates for each of the component genes.
        * The chromosomes and the bounds will be copied from the associated candidates,
        * while the fitness and constraint vectors of the candidate will be empty.
        * 
        * @param parts The component candidates to create this candidate from.
        */
        explicit Candidate(std::tuple<Candidate<Ts>...> parts) :
            Candidate(constructor_tag_t{}, { std::move(std::get<Candidate<Ts>>(parts).chromosome)... }, get_bounds(parts, detail::filter_types_t<is_bounded_gene, Ts...>{}))
        {}

        /** @returns The length of the chromosome associated with the specified gene type. */
        template<typename GeneType>
        size_t chrom_len() const noexcept;

        /** @returns The candidate's chromosome associated with the specified gene type. */
        template<typename GeneType>
        Chromosome<GeneType>& chrom() & noexcept;

        /** @returns The candidate's chromosome associated with the specified gene type. */
        template<typename GeneType>
        const Chromosome<GeneType>& chrom() const& noexcept;

        /** @returns The bounds of the chromosome associated with the specified gene type. */
        template<typename GeneType>
        BoundsView<GeneType> bounds() const noexcept;


        Candidate(const Candidate&) = default;
        Candidate(Candidate&&)      = default;

        Candidate& operator=(const Candidate&);
        Candidate& operator=(Candidate&&) noexcept;

    private:
        struct constructor_tag_t {};

        template<typename... Us>
        Candidate(constructor_tag_t, std::tuple<Chromosome<Ts>...> chroms, std::tuple<BoundsView<Us>...> bounds)
        {
            ( (this->Candidate<Ts>::chromosome = std::move(std::get<Chromosome<Ts>>(chroms))), ... );
            ( (this->Candidate<Us>::gene_bounds = std::move(std::get<BoundsView<Us>>(bounds))), ... );
        }

        template<typename... Us, typename... BoundedGenes>
        static auto get_bounds(const std::tuple<Candidate<Us>...>& parts, const std::tuple<BoundedGenes...>&)
        {
            return std::tuple{ std::get<Candidate<BoundedGenes>>(parts).gene_bounds... };
        }
    };

    /** A pair of candidates. */
    template<typename T>
    struct CandidatePair
    {
        Candidate<T> first;
        Candidate<T> second;
    };


    /**
    * Comparison operators for the candidates. Only the chromosomes are considered for the
    * comparisons, the other parts of the candidates are ignored.
    * The comparisons are not transitive if T is a floating-point type.
    * 
    * @returns True if the chromosomes of the candidates are equal.
    */
    template<typename T>
    bool operator==(const Candidate<T>& lhs, const Candidate<T>& rhs) noexcept;

    template<typename... Ts>
    bool operator==(const Candidate<MixedGene<Ts...>>& lhs, const Candidate<MixedGene<Ts...>>& rhs) noexcept;

    /**
    * Hash function for the candidates. Only the chromosomes are considered for the
    * hash computations, the other parts of the candidates are ignored.
    */
    template<typename T>
    struct CandidateHasher
    {
        size_t operator()(const Candidate<T>& candidate) const noexcept;
    };

    template<typename... Ts>
    struct CandidateHasher<MixedGene<Ts...>>
    {
        size_t operator()(const Candidate<MixedGene<Ts...>>& candidate) const noexcept;
    };

} // namespace gapp


/* IMPLEMENTATION */

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
    Candidate<T>::Candidate(BoundsView<T> bounds) requires(is_bounded_gene_v<T>) :
        chromosome(bounds.size())
    {
        this->gene_bounds = bounds;
    }

    template<typename T>
    Candidate<T>::Candidate(Chromosome<T> chrom) requires(!is_bounded_gene_v<T>) :
        chromosome(std::move(chrom))
    {}

    template<typename T>
    Candidate<T>::Candidate(Chromosome<T> chrom, BoundsView<T> bounds) requires(is_bounded_gene_v<T>) :
        chromosome(std::move(chrom))
    {
        GAPP_ASSERT(chromosome.size() == bounds.size(), "Mismatching chromosome and bounds sizes.");

        this->gene_bounds = bounds;
    }

    template<typename T>
    template<typename GeneType>
    size_t Candidate<T>::chrom_len() const noexcept
    {
        static_assert(std::is_same_v<GeneType, T>, "Invalid gene type.");

        return chromosome.size();
    }

    template<typename T>
    template<typename GeneType>
    Chromosome<GeneType>& Candidate<T>::chrom() & noexcept
    {
        static_assert(std::is_same_v<GeneType, T>, "Invalid gene type.");

        return chromosome;
    }

    template<typename T>
    template<typename GeneType>
    const Chromosome<GeneType>& Candidate<T>::chrom() const& noexcept
    {
        static_assert(std::is_same_v<GeneType, T>, "Invalid gene type.");

        return chromosome;
    }

    template<typename T>
    template<typename GeneType>
    BoundsView<GeneType> Candidate<T>::bounds() const noexcept requires(is_bounded_gene_v<T>)
    {
        static_assert(std::is_same_v<GeneType, T>, "Invalid gene type.");

        return this->gene_bounds;
    }

    template<typename... Ts>
    Candidate<MixedGene<Ts...>>& Candidate<MixedGene<Ts...>>::operator=(const Candidate& other)
    {
        this->fitness = other.fitness;
        this->constraint_violation = other.constraint_violation;

        ( detail::CandidateBounds<Ts>::operator=(other), ... );

        ( (Candidate<Ts>::chromosome = other.Candidate<Ts>::chromosome), ... );

        return *this;
    }

    template<typename... Ts>
    Candidate<MixedGene<Ts...>>& Candidate<MixedGene<Ts...>>::operator=(Candidate&& other) noexcept
    {
        this->fitness = std::move(other.fitness);
        this->constraint_violation = std::move(other.constraint_violation);

        ( detail::CandidateBounds<Ts>::operator=(std::move(other)), ... );

        ( (Candidate<Ts>::chromosome = std::move(other.Candidate<Ts>::chromosome)), ... );

        return *this;
    }

    template<typename... Ts>
    template<typename GeneType>
    size_t Candidate<MixedGene<Ts...>>::chrom_len() const noexcept
    {
        return Candidate<GeneType>::chromosome.size();
    }

    template<typename... Ts>
    template<typename GeneType>
    Chromosome<GeneType>& Candidate<MixedGene<Ts...>>::chrom() & noexcept
    {
        return Candidate<GeneType>::chromosome;
    }

    template<typename... Ts>
    template<typename GeneType>
    const Chromosome<GeneType>& Candidate<MixedGene<Ts...>>::chrom() const& noexcept
    {
        return Candidate<GeneType>::chromosome;
    }

    template<typename... Ts>
    template<typename GeneType>
    BoundsView<GeneType> Candidate<MixedGene<Ts...>>::bounds() const noexcept
    {
        static_assert(is_bounded_gene_v<GeneType>, "The specified gene type must be a bounded gene.");

        return Candidate<GeneType>::gene_bounds;
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

    template<typename... Ts>
    bool operator==(const Candidate<MixedGene<Ts...>>& lhs, const Candidate<MixedGene<Ts...>>& rhs) noexcept
    {
        return ( (static_cast<const Candidate<Ts>&>(lhs) == static_cast<const Candidate<Ts>&>(rhs)) && ... );
    }

    template<typename T>
    size_t CandidateHasher<T>::operator()(const Candidate<T>& candidate) const noexcept
    {
        return detail::hash_range(candidate.chromosome.begin(), candidate.chromosome.end());
    }

    template<typename... Ts>
    size_t CandidateHasher<MixedGene<Ts...>>::operator()(const Candidate<MixedGene<Ts...>>& candidate) const noexcept
    {
        return detail::hash_combine(detail::hash_range(static_cast<Candidate<Ts>>(candidate).chromosome)...);
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

namespace gapp::detail
{
    /*
    * Exact comparison functions for the candidates. Only the chromosomes are considered
    * for the comparisons, the other parts of the candidates are ignored.
    * The comparisons are transitive even if T is a floating-point type.
    */
    struct CandidateEqualExact
    {
    private:
        template<typename T>
        bool isEqualImpl(const Candidate<T>& lhs, const Candidate<T>& rhs) noexcept
        {
            return lhs.chromosome == rhs.chromosome;
        }

        template<typename... Ts>
        bool isEqualImpl(const Candidate<MixedGene<Ts...>>& lhs, const Candidate<MixedGene<Ts...>>& rhs) noexcept
        {
            return ((lhs.Candidate<Ts>::chromosome == rhs.Candidate<Ts>::chromosome) && ...);
        }

    public:
        template<typename T>
        bool operator()(const Candidate<T>& lhs, const Candidate<T>& rhs) noexcept
        {
            return isEqualImpl(lhs, rhs);
        }
    };

    struct CandidateLessThanExact
    {
    private:
        template<typename T>
        bool isLessThanImpl(const Candidate<T>& lhs, const Candidate<T>& rhs) noexcept
        {
            return lhs.chromosome < rhs.chromosome;
        }

        template<typename... Ts>
        bool isLessThanImpl(const Candidate<MixedGene<Ts...>>& lhs, const Candidate<MixedGene<Ts...>>& rhs) noexcept
        {
            bool less = false;

            (((lhs.Candidate<Ts>::chromosome == rhs.Candidate<Ts>::chromosome) ? true :
             ((less = lhs.Candidate<Ts>::chromosome < rhs.Candidate<Ts>::chromosome), false)) && ...);

            return less;
        }

    public:
        template<typename T>
        bool operator()(const Candidate<T>& lhs, const Candidate<T>& rhs) noexcept
        {
            return isLessThanImpl(lhs, rhs);
        }
    };

} // namespace gapp::detail

#endif // !GAPP_CORE_CANDIDATE_HPP
