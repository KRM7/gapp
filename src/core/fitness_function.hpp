/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CORE_FITNESS_FUNCTION_HPP
#define GA_CORE_FITNESS_FUNCTION_HPP

#include "../population/population.hpp"
#include "../utility/bounded_value.hpp"
#include "../utility/utility.hpp"
#include <concepts>
#include <functional>
#include <utility>
#include <cstddef>

namespace genetic_algorithm
{
    /**
    * Base class for the fitness functions used in the algorithms. \n
    * The fitness function takes a candidate solution (chromosome) and returns a fitness vector
    * after evaluating it.
    * 
    * This should be used as the base class for the fitness functions if the number
    * of objectives or the chromosome length is not known at compile time.
    * Otherwise use FitnessFunction, or SingleObjFinessFunction for single-objective
    * problems.
    * 
    * @tparam T The gene type used in the algorithm.
    */
    template<typename T>
    class FitnessFunctionBase
    {
    public:
        /** The gene type used in the chromosomes. */
        using GeneType = T;

        /**
        * Create a fitness function.
        *
        * @param chrom_len The chromosome length that will be used for the candidate
        *   solutions in the algorithm. Must be at least 1.
        * @param variable_len Should be true if the fitness function supports and will be used with
        *   variable chromosome lengths. The other genetic operators (crossover, mutation, repair)
        *   must also support variable length chromosomes if this is enabled.
        * @param dynamic Should be true if the fitness vector returned for a chromosome will not
        *   always be the same for the same chromosome (eg. it changes over time or isn't deterministic).
        */
        constexpr FitnessFunctionBase(Positive<size_t> chrom_len, bool variable_len = false, bool dynamic = false) noexcept :
            chrom_len_(chrom_len), variable_chrom_len_(variable_len), dynamic_(dynamic)
        {}

        /** @returns The chromosome length the fitness function expects. */
        [[nodiscard]]
        constexpr size_t chrom_len() const noexcept { return chrom_len_; }

        /** @returns true if variable chromosome lengths are allowed. */
        [[nodiscard]]
        constexpr bool variable_chrom_len() const noexcept { return variable_chrom_len_; }

        /** @returns True if dynamic fitness function support is enabled. */
        [[nodiscard]]
        constexpr bool dynamic() const noexcept { return dynamic_; }

        /**
        * Compute the fitness value of a chromosome.
        * The size of the chromosome must be equal to the chrom_len set unless variable
        * chromosome lengths are used.
        * 
        * @param chrom The chromosome to evaluate.
        * @returns The fitness vector of the chromosome.
        */
        FitnessVector operator()(const Chromosome<T>& chrom) const
        {
            GA_ASSERT(chrom.size() == chrom_len() || variable_chrom_len(), "A chromosome of incorrect size was passed to the fitness function.");

            return invoke(chrom);
        }


        /** Destructor. */
        virtual ~FitnessFunctionBase()                                      = default;

    protected:

        FitnessFunctionBase(const FitnessFunctionBase&) noexcept            = default;
        FitnessFunctionBase(FitnessFunctionBase&&) noexcept                 = default;
        FitnessFunctionBase& operator=(const FitnessFunctionBase&) noexcept = default;
        FitnessFunctionBase& operator=(FitnessFunctionBase&&) noexcept      = default;

    private:
        /** The implementation of the fitness function. Should be thread-safe. */
        virtual FitnessVector invoke(const Chromosome<T>& chrom) const = 0;

        Positive<size_t> chrom_len_;
        bool variable_chrom_len_ = false;
        bool dynamic_            = false;
    };

    /**
    * Base class for the fitness functions used in the algorithms. \n
    * The fitness function takes a candidate solution (chromosome) and returns a fitness vector
    * after evaluating it.
    *
    * This class can be used as the base class for the fitness functions if the
    * chromosome length is known at compile time, otherwise use FitnessFunctionBase.
    *
    * @tparam T The gene type used in the algorithm.
    * @tparam ChromLen The length of the chromosomes expected by the fitness function. Must be at least 1.
    */
    template<typename T, size_t ChromLen>
    class FitnessFunction : public FitnessFunctionBase<T>
    {
    public:
        /**
        * Create a fitness function.
        *
        * @param variable_len Should be true if the fitness function supports and will be used with
        *   variable chromosome lengths. The other genetic operators (crossover, mutation, repair)
        *   must also support variable length chromosomes if this is enabled.
        * @param dynamic Should be true if the fitness vector returned for a chromosome will not
        *   always be the same for the same chromosome (eg. it changes over time or isn't deterministic).
        */
        constexpr FitnessFunction(bool variable_len = false, bool dynamic = false) noexcept :
            FitnessFunctionBase<T>(ChromLen, variable_len, dynamic)
        {}

        /** Destructor. */
        virtual ~FitnessFunction()                                  = default;

    protected:

        FitnessFunction(const FitnessFunction&) noexcept            = default;
        FitnessFunction(FitnessFunction&&) noexcept                 = default;
        FitnessFunction& operator=(const FitnessFunction&) noexcept = default;
        FitnessFunction& operator=(FitnessFunction&&) noexcept      = default;
    };

} // namespace genetic_algorithm

namespace genetic_algorithm::detail
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
            GA_ASSERT(f, "The fitness function can't be a nullptr.");

            fitness_function_ = std::move(f);
        }

    private:
        FitnessVector invoke(const Chromosome<T>& chrom) const override
        {
            GA_ASSERT(fitness_function_);

            return fitness_function_(chrom);
        }

        FitnessCallable fitness_function_;
    };

} // namespace genetic_algorithm::detail

#endif // !GA_CORE_FITNESS_FUNCTION_HPP