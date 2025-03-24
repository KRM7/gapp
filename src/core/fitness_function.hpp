/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_CORE_FITNESS_FUNCTION_HPP
#define GAPP_CORE_FITNESS_FUNCTION_HPP

#include "candidate.hpp"
#include "../encoding/gene_types.hpp"
#include "../utility/small_vector.hpp"
#include "../utility/bounded_value.hpp"
#include "../utility/functional.hpp"
#include "../utility/type_traits.hpp"
#include "../utility/utility.hpp"
#include <array>
#include <algorithm>
#include <functional>
#include <utility>
#include <cstddef>

namespace gapp
{
    /**
    * The base class of the fitness functions used in the GAs.
    * Contains all of the properties of a fitness function that are not dependent
    * on the gene type.
    * 
    * Fitness function implementations should not be derived directly from this class,
    * but instead from FitnessFunctionBase or FitnessFunction.
    */
    class FitnessFunctionInfo
    {
    public:
        /**
        * The list of potential fitness function types.
        * A fitness function may either be static or dynamic.
        * 
        * @var Type::Static The value representing a static fitness function. A fitness function
        *   is considered static if it always returns the same fitness vector for a particular
        *   candidate solution.
        * @var Type::Dynamic The value representing a dynamic fitness function. A fitness function
        *   is considered to be dynamic if it may return different fitness vectors for the same
        *   candidate solution over multiple calls to the fitness function.
        */
        enum class Type : bool { Static = false, Dynamic = true };

        /** @returns The chromosome lengths used by the fitness function for each chromosome of the encoding. */
        [[nodiscard]]
        constexpr const small_vector<size_t>& chrom_lens() const noexcept { return chrom_lens_; }

        /** @returns true if the fitness function is dynamic. */
        [[nodiscard]]
        constexpr bool is_dynamic() const noexcept { return type_ == Type::Dynamic; }

        /** Destructor. */
        virtual ~FitnessFunctionInfo()                             = default;

    protected:

        FitnessFunctionInfo(const FitnessFunctionInfo&)            = default;
        FitnessFunctionInfo(FitnessFunctionInfo&&)                 = default;
        FitnessFunctionInfo& operator=(const FitnessFunctionInfo&) = default;
        FitnessFunctionInfo& operator=(FitnessFunctionInfo&&)      = default;

        constexpr FitnessFunctionInfo(small_vector<size_t> chrom_lens, Type type) :
            chrom_lens_(std::move(chrom_lens)), type_(type)
        {
            GAPP_ASSERT(std::all_of(chrom_lens_.begin(), chrom_lens_.end(), detail::greater_than(0_sz)));
        }

    private:
        small_vector<size_t> chrom_lens_;
        Type type_;
    };

    /**
    * The base class of the fitness functions used in the GAs.
    * The fitness functions take a candidate solution as a parameter and return
    * a fitness vector after evaluating the chromosomes of the candidate.
    * 
    * This class should be used as the base class for fitness functions if the
    * chromosome lengths are not known at compile time. If the chromosome lengths
    * are known at compile, use FitnessFunction as the base class instead.
    * 
    * @tparam T The gene type of the candidates expected by the fitness function.
    */
    template<typename T>
    class FitnessFunctionBase : public FitnessFunctionInfo
    {
    public:
        /** The gene type of the candidates that can be evaluated by the fitness function. */
        using GeneType = T;

        /**
        * Create a fitness function.
        *
        * @param chrom_len The chromosome length expected by the fitness function, and which
        *   will be used for the candidate solutions in the GA.
        *   Must be at least 1, and a value must be specified even if the chromosome length
        *   is variable, as the value will still be used to generate the initial population.
        * @param type The type of the fitness function. The value should be either Type::Static
        *   or Type::Dynamic, based on whether the fitness function always returns the same
        *   fitness vector for a solution (static) or not (dynamic).
        */
        constexpr FitnessFunctionBase(Positive<size_t> chrom_len, Type type = Type::Static) noexcept requires(!is_mixed_gene_v<GeneType>) :
            FitnessFunctionInfo(small_vector{ *chrom_len }, type)
        {}

        /**
        * Create a fitness function.
        *
        * @param chrom_lens The size of each of the chromosomes expected by the fitness function,
        *   and which will be used for the candidate solutions in the GA.
        *   Each chromosome length must be at least 1, and a value must be specified for each
        *   chromosome encoding the candidates, even if the chromosome length is variable, as
        *   the value will still be used to generate the initial population.
        * @param type The type of the fitness function. The value should be either Type::Static
        *   or Type::Dynamic, based on whether the fitness function always returns the same
        *   fitness vector for a solution (static) or not (dynamic).
        */
        constexpr FitnessFunctionBase(std::array<size_t, Candidate<T>::NumChroms> chrom_lens, Type type = Type::Static) :
            FitnessFunctionInfo(small_vector(chrom_lens.begin(), chrom_lens.end()), type)
        {}

        /** @returns The chromosome length used by the fitness function. */
        template<typename GeneType = T>
        [[nodiscard]] constexpr size_t chrom_len() const noexcept requires(!is_mixed_gene_v<T>);

        /** @returns The length of the chromosome associated with the gene type @p GeneType. */
        template<typename GeneType>
        [[nodiscard]] constexpr size_t chrom_len() const noexcept requires(is_mixed_gene_v<T>);

        /**
        * Compute the fitness value of a solution.
        * 
        * The length of the chromosomes will be equal to the chromosome lengths set
        * for the fitness function, unless variable chromosome lengths are allowed.
        * 
        * @param sol The candidate solution to evaluate.
        * @returns The fitness vector of the candidate, with a size equal to the number of objectives.
        */
        FitnessVector operator()(const Candidate<T>& sol) const { return invoke(sol); }

    private:
        /** The implementation of the fitness function. Should be thread-safe. */
        virtual FitnessVector invoke(const Candidate<T>& sol) const = 0;
    };

    /**
    * The base class of the fitness functions used in the GAs.
    * The fitness functions take a candidate solution as an argument
    * and return a fitness vector after evaluating the chromosomes of
    * the candidate.
    *
    * This class can only be used as the base class for fitness functions if
    * their chromosome lengths are known at compile time. If the chromosome
    * lengths are not known at compile time, use FitnessFunctionBase as the
    * base class instead.
    *
    * @tparam T The gene type of the candidates expected by the fitness function.
    * @tparam ChromLens The size of each chromosome of the candidate solutions expected
    *   by the fitness function. Each of the chromosome lengths must be at least 1.
    *   The number of chromosome lengths must be at least 1, but specifying more than 1
    *   chromosome length is only allowed if T is a MixedGene type, and they must match
    *   the number of genes in the mixed gene type.
    *   For variable chromosome length problems, these values will only be used during the
    *   generation of the initial population.
    */
    template<typename T, size_t... ChromLens>
    class FitnessFunction : public FitnessFunctionBase<T>
    {
    public:
        static_assert((ChromLens && ...));
        static_assert(sizeof...(ChromLens) == Candidate<T>::NumChroms);

        using Type = FitnessFunctionInfo::Type;

        /**
        * Create a fitness function.
        *
        * @param type The type of the fitness function. The value should be either Type::Static
        *   or Type::Dynamic, based on whether the fitness function always returns the same
        *   fitness vector for a solution (static) or not (dynamic).
        */
        constexpr FitnessFunction(Type type = Type::Static) noexcept :
            FitnessFunctionBase<T>(std::array{ ChromLens... }, type)
        {}
    };

} // namespace gapp

namespace gapp::detail
{
    /* Wrap a callable so it can be used as a fitness function. */
    template<typename T>
    class FitnessLambda final : public FitnessFunctionBase<T>
    {
    public:
        using FitnessCallable = std::function<FitnessVector(const Candidate<T>&)>;

        FitnessLambda(std::array<size_t, Candidate<T>::NumChroms> chrom_lens, FitnessCallable f) :
            FitnessFunctionBase<T>(chrom_lens), fitness_function_(std::move(f))
        {
            GAPP_ASSERT(fitness_function_, "The fitness function can't be a nullptr.");
        }

        FitnessLambda(size_t chrom_len, FitnessCallable f) noexcept requires(!is_mixed_gene_v<T>) :
            FitnessLambda(std::array{ chrom_len }, std::move(f))
        {}

    private:
        FitnessVector invoke(const Candidate<T>& sol) const override
        {
            GAPP_ASSERT(fitness_function_);
            return fitness_function_(sol);
        }

        FitnessCallable fitness_function_;
    };

} // namespace gapp::detail

namespace gapp
{
    template<typename T>
    template<typename GeneType>
    constexpr size_t FitnessFunctionBase<T>::chrom_len() const noexcept requires(!is_mixed_gene_v<T>)
    {
        static_assert(std::is_same_v<GeneType, T>,
          "No chromosome exists for the specified gene type.");

        return chrom_lens()[0];
    }

    template<typename T>
    template<typename GeneType>
    constexpr size_t FitnessFunctionBase<T>::chrom_len() const noexcept requires(is_mixed_gene_v<T>)
    {
        static_assert(detail::args_to_list_t<T>::template contains<GeneType>,
          "No chromosome exists for the specified gene type.");

        return chrom_lens()[detail::args_to_list_t<T>::template index_of<GeneType>];
    }

} // namespace gapp

#endif // !GAPP_CORE_FITNESS_FUNCTION_HPP
