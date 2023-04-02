/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CORE_FITNESS_FUNCTION_HPP
#define GA_CORE_FITNESS_FUNCTION_HPP

#include "../population/population.hpp"
#include "../utility/bounded_value.hpp"
#include <concepts>
#include <cstddef>

namespace genetic_algorithm
{
    namespace detail { inline constexpr size_t dynamic_tag = std::numeric_limits<size_t>::max(); }

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
    template<Gene T>
    class FitnessFunctionBase
    {
    public:
        /** The gene type used in the chromosomes. */
        using GeneType = T;

        /**
        * The fitness vector type returned by the fitness functions. \n
        * Contains a fitness value for each objective axis.
        */
        using FitnessVector = detail::FitnessVector;

        /**
        * Create a fitness function.
        *
        * @param chrom_len The chromosome length that will be used for the candidate
        *   solutions in the algorithm. Must be at least 1.
        * @param num_objectives The number of objective functions used. This is the
        *   expected size of the fitness vector returned by the fitness function.
        *   Must be at least 1.
        * @param variable_len Should be true if the fitness function supports and will be used with
        *   variable chromosome lengths. The other genetic operators (crossover, mutation, repair)
        *   must also support variable length chromosomes if this is enabled.
        * @param dynamic Should be true if the fitness vector returned for a chromosome will not
        *   always be the same for the same chromosome (eg. it changes over time or isn't deterministic).
        */
        constexpr FitnessFunctionBase(Positive<size_t> chrom_len, Positive<size_t> num_objectives, bool variable_len = false, bool dynamic = false) noexcept :
            num_objectives_(num_objectives), chrom_len_(chrom_len), dynamic_(dynamic), variable_chrom_len_(variable_len)
        {}

        /** @returns The chromosome length the fitness function expects. */
        [[nodiscard]]
        constexpr size_t chrom_len() const noexcept { return chrom_len_; }

        /** @returns true if variable chromosome lengths are allowed. */
        [[nodiscard]]
        constexpr bool variable_chrom_len() const noexcept { return variable_chrom_len_; }

        /** @returns The number of objectives as determined by the algorithm based on the fitness function. */
        [[nodiscard]]
        constexpr size_t num_objectives() const noexcept { return num_objectives_; }

        /** @returns True if the fitness function is a single-objective fitness function. */
        [[nodiscard]]
        constexpr bool single_objective() const noexcept { return num_objectives_ == 1; }

        /** @returns True if the fitness function is a multi-objective fitness function. */
        [[nodiscard]]
        constexpr bool multi_objective() const noexcept { return num_objectives_ > 1; }

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
        FitnessVector operator()(const Chromosome<T>& chrom);


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

        Positive<size_t> num_objectives_;
        Positive<size_t> chrom_len_;
        bool dynamic_            = false;
        bool variable_chrom_len_ = false;
    };

    /**
    * Base class for the fitness functions used in the algorithms. \n
    * The fitness function takes a candidate solution (chromosome) and returns a fitness vector
    * after evaluating it.
    *
    * This should be used as the base class for the fitness functions if the number
    * of objectives or the chromosome length is known at compile time,
    * otherwise use FitnessFunctionBase.
    *
    * @tparam T The gene type used in the algorithm.
    * @tparam Nobj The number of objectives. Must be at least 1.
    * @tparam ChromLen The length of the chromosomes expected by the fitness function. Must be at least 1.
    */
    template<Gene T, size_t Nobj = detail::dynamic_tag, size_t ChromLen = detail::dynamic_tag>
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
        constexpr FitnessFunction(bool variable_len = false, bool dynamic = false) noexcept
        requires (Nobj != detail::dynamic_tag) && (ChromLen != detail::dynamic_tag) :
            FitnessFunctionBase<T>(ChromLen, Nobj, variable_len, dynamic)
        {}

        /**
        * Create a fitness function.
        *
        * @param num_objectives The number of objective functions used. This is the
        *   expected size of the fitness vector returned by the fitness function.
        *   Must be at least 1.
        * @param variable_len Should be true if the fitness function supports and will be used with
        *   variable chromosome lengths. The other genetic operators (crossover, mutation, repair)
        *   must also support variable length chromosomes if this is enabled.
        * @param dynamic Should be true if the fitness vector returned for a chromosome will not
        *   always be the same for the same chromosome (eg. it changes over time or isn't deterministic).
        */
        constexpr FitnessFunction(Positive<size_t> num_objectives, bool variable_len = false, bool dynamic = false) noexcept
        requires (Nobj == detail::dynamic_tag) && (ChromLen != detail::dynamic_tag) :
            FitnessFunctionBase<T>(ChromLen, num_objectives, variable_len, dynamic)
        {}

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
        constexpr FitnessFunction(Positive<size_t> chrom_len, bool variable_len = false, bool dynamic = false) noexcept
        requires (Nobj != detail::dynamic_tag) && (ChromLen == detail::dynamic_tag) :
            FitnessFunctionBase<T>(chrom_len, Nobj, variable_len, dynamic)
        {}

        /**
        * Create a fitness function.
        *
        * @param chrom_len The chromosome length that will be used for the candidate
        *   solutions in the algorithm. Must be at least 1.
        * @param num_objectives The number of objective functions used. This is the
        *   expected size of the fitness vector returned by the fitness function.
        *   Must be at least 1.
        * @param variable_len Should be true if the fitness function supports and will be used with
        *   variable chromosome lengths. The other genetic operators (crossover, mutation, repair)
        *   must also support variable length chromosomes if this is enabled.
        * @param dynamic Should be true if the fitness vector returned for a chromosome will not
        *   always be the same for the same chromosome (eg. it changes over time or isn't deterministic).
        */
        constexpr FitnessFunction(Positive<size_t> chrom_len, Positive<size_t> num_objectives, bool variable_len = false, bool dynamic = false) noexcept
        requires (Nobj == detail::dynamic_tag) && (ChromLen == detail::dynamic_tag) :
            FitnessFunctionBase<T>(chrom_len, num_objectives, variable_len, dynamic)
        {}


        /** Destructor. */
        virtual ~FitnessFunction()                                  = default;

    protected:

        FitnessFunction(const FitnessFunction&) noexcept            = default;
        FitnessFunction(FitnessFunction&&) noexcept                 = default;
        FitnessFunction& operator=(const FitnessFunction&) noexcept = default;
        FitnessFunction& operator=(FitnessFunction&&) noexcept      = default;
    };

    /**
    * Base class for single-objective fitness functions. \n
    * The fitness function takes a candidate solution (chromosome) and returns a fitness vector
    * of size 1 after evaluating it.
    *
    * @tparam ChromLen The length of the chromosomes expected by the fitness function.
    */
    template<Gene T, size_t ChromLen = detail::dynamic_tag>
    using SingleObjFitnessFunction = FitnessFunction<T, 1, ChromLen>;


    /**
    * The general callable type that can be used as a fitness function. \n
    * Takes a vector of genes (chromosome) and returns the fitness vector of the chromosome. \n
    * The returned fitness vector should always be the same length. \n
    */
    template<Gene T>
    using FitnessCallable = std::function<detail::FitnessVector(const Chromosome<T>&)>;

    /**
    * Create a fitness function that can be used in the genetic algorithms from a callable  object.
    * 
    * @param num_obj The number of objective functions (i.e. the size of the fitness vector returned by the fitness callable).
    * @param chrom_len The length of the chromosomes in the population.
    * @param f The fitness function to use in the algorithm.
    * @returns The fitness function.
    */
    template<Gene T>
    std::unique_ptr<FitnessFunctionBase<T>> makeFitnessFunction(size_t num_obj, size_t chrom_len, FitnessCallable<T> f);

    /**
    * Create a fitness function that can be used in the genetic algorithms from a callable object.
    * 
    * @tparam Nobj The number of objective functions (i.e. the size of the fitness vector returned by the fitness callable).
    * @tparam ChromLen The length of the chromosomes in the population.
    * @param f The fitness function to use in the algorithm.
    * @returns The fitness function.
    */
    template<Gene T, size_t Nobj, size_t ChromLen>
    std::unique_ptr<FitnessFunctionBase<T>> makeFitnessFunction(FitnessCallable<T> f);

    /**
    * Create a single-objective fitness function that can be used in the genetic algorithms from a callable object.
    *
    * @tparam ChromLen The length of the chromosomes in the population.
    * @param f The fitness function to use in the algorithm.
    * @returns The fitness function.
    */
    template<Gene T, size_t ChromLen>
    std::unique_ptr<FitnessFunctionBase<T>> makeSingleObjFitnessFunction(FitnessCallable<T> f);

} // namespace genetic_algorithm


/* IMPLEMENTATION */

#include "../utility/utility.hpp"
#include <functional>
#include <utility>

namespace genetic_algorithm::detail
{
    /*
    * Wraps a callable with the right signature so that it can be used as a fitness
    * function in the GAs.
    */
    template<Gene T>
    class FitnessLambda final : public FitnessFunctionBase<T>
    {
    public:
        FitnessLambda(size_t chrom_len, size_t num_objectives, FitnessCallable<T> f) noexcept :
            FitnessFunction<T>(chrom_len, num_objectives)
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

        FitnessCallable<T> fitness_function_;
    };

} // namespace genetic_algorithm::detail

namespace genetic_algorithm
{
    template<Gene T>
    auto FitnessFunctionBase<T>::operator()(const Chromosome<T>& chrom) -> FitnessVector
    {
        GA_ASSERT(chrom.size() == chrom_len() || variable_chrom_len(), "A chromosome of incorrect size was passed to the fitness function.");

        FitnessVector fvec = this->invoke(chrom);

        GA_ASSERT(fvec.size() == num_objectives(), "A fitness vector of incorrect length was returned by the fitness function.");

        return fvec;
    }

    template<Gene T>
    inline std::unique_ptr<FitnessFunctionBase<T>> makeFitnessFunction(size_t chrom_len, size_t num_obj, FitnessCallable<T> f)
    {
        return std::make_unique<detail::FitnessLambda<T>>(chrom_len, num_obj, std::move(f));
    }

    template<Gene T, size_t Nobj, size_t ChromLen>
    inline std::unique_ptr<FitnessFunctionBase<T>> makeFitnessFunction(FitnessCallable<T> f)
    {
        return std::make_unique<detail::FitnessLambda<T>>(ChromLen, Nobj, std::move(f));
    }

    template<Gene T, size_t ChromLen>
    inline std::unique_ptr<FitnessFunctionBase<T>> makeSingleObjFitnessFunction(FitnessCallable<T> f)
    {
        return std::make_unique<detail::FitnessLambda<T>>(ChromLen, 1, std::move(f));
    }

} // namespace genetic_algorithm

#endif // !GA_CORE_FITNESS_FUNCTION_HPP