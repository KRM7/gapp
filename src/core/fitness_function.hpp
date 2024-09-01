/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CORE_FITNESS_FUNCTION_HPP
#define GA_CORE_FITNESS_FUNCTION_HPP

#include "candidate.hpp"
#include "../utility/bounded_value.hpp"
#include "../utility/utility.hpp"
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
        enum class Type { Static = 0, Dynamic = 1 };

        /**
        * Create a fitness function.
        *
        * @param chrom_len The chromosome length that is expected by the fitness function,
        *   and will be used for the candidate solutions in the GA.
        *   Must be at least 1, and a value must be specified even if the chromosome length
        *   is variable, as it will still be used to generate the initial population.
        * @param type The type of the fitness function. The value should be either Type::Static
        *   or Type::Dynamic, based on whether the fitness function always returns the same 
        *   fitness vector for a solution (static) or not (dynamic).
        */
        constexpr FitnessFunctionInfo(Positive<size_t> chrom_len, Type type = Type::Static) noexcept :
            chrom_len_(chrom_len), type_(type)
        {}

        /** @returns The chromosome length the fitness function expects. */
        [[nodiscard]]
        constexpr size_t chrom_len() const noexcept { return chrom_len_; }

        /** @returns True if the fitness function is dynamic. */
        [[nodiscard]]
        constexpr bool is_dynamic() const noexcept { return type_ == Type::Dynamic; }

        /** Destructor. */
        virtual ~FitnessFunctionInfo()                             = default;

    protected:

        FitnessFunctionInfo(const FitnessFunctionInfo&)            = default;
        FitnessFunctionInfo(FitnessFunctionInfo&&)                 = default;
        FitnessFunctionInfo& operator=(const FitnessFunctionInfo&) = default;
        FitnessFunctionInfo& operator=(FitnessFunctionInfo&&)      = default;

    private:
        Positive<size_t> chrom_len_;
        Type type_;
    };

    /**
    * The base class of the fitness functions used in the GAs.
    * The fitness functions take a candidate solution (chromosome) as a parameter
    * and return a fitness vector after evaluating the chromosome.
    * 
    * This should be used as the base class for fitness functions if the chromosome
    * length is not known at compile time.
    * If the chromosome length is known at compile, use FitnessFunction as the base class instead.
    * 
    * @tparam T The gene type expected by the fitness function.
    */
    template<typename T>
    class FitnessFunctionBase : public FitnessFunctionInfo
    {
    public:
        /** The gene type of the chromosomes that can be evaluated by the fitness function. */
        using GeneType = T;

        using FitnessFunctionInfo::FitnessFunctionInfo;

        /**
        * Compute the fitness value of a chromosome.
        * 
        * The size of the chromosome must be equal to the chromosome length set
        * for the fitness function, unless variable chromosome lengths are allowed.
        * 
        * @param chrom The chromosome to evaluate.
        * @returns The fitness vector of the chromosome, with a size equal to the number of objectives.
        */
        FitnessVector operator()(const Chromosome<T>& chrom) const { return invoke(chrom); }

    private:
        /** The implementation of the fitness function. Should be thread-safe. */
        virtual FitnessVector invoke(const Chromosome<T>& chrom) const = 0;
    };

    /**
    * The base class of the fitness functions used in the algorithms.
    * The fitness functions take a candidate solution (chromosome) as a parameter
    * and return a fitness vector after evaluating the chromosome.
    *
    * This should only be used as the base class for fitness functions if the chromosome length
    * is known at compile time.
    * If the chromosome length is not known at compile, use FitnessFunctionBase as the base class instead.
    *
    * @tparam T The gene type expected by the fitness function.
    * @tparam ChromLen The length of the chromosomes expected by the fitness function. Must be at least 1.
    *   For variable chromosome length problems, this value will only be used during the generation of the
    *   initial population.
    */
    template<typename T, size_t ChromLen>
    class FitnessFunction : public FitnessFunctionBase<T>
    {
    public:
        using Type = FitnessFunctionInfo::Type;

        /**
        * Create a fitness function.
        *
        * @param type The type of the fitness function. The value should be either Type::Static
        *   or Type::Dynamic, based on whether the fitness function always returns the same
        *   fitness vector for a solution (static) or not (dynamic).
        */
        constexpr FitnessFunction(Type type = Type::Static) noexcept :
            FitnessFunctionBase<T>(ChromLen, type)
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
        using FitnessCallable = std::function<FitnessVector(const Chromosome<T>&)>;

        FitnessLambda(size_t chrom_len, FitnessCallable f) noexcept :
            FitnessFunctionBase<T>(chrom_len)
        {
            GAPP_ASSERT(f, "The fitness function can't be a nullptr.");

            fitness_function_ = std::move(f);
        }

    private:
        FitnessVector invoke(const Chromosome<T>& chrom) const override
        {
            GAPP_ASSERT(fitness_function_);

            return fitness_function_(chrom);
        }

        FitnessCallable fitness_function_;
    };

} // namespace gapp::detail

#endif // !GA_CORE_FITNESS_FUNCTION_HPP