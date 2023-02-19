/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CORE_FITNESS_FUNCTION_HPP
#define GA_CORE_FITNESS_FUNCTION_HPP

#include "../population/population.hpp"
#include <concepts>
#include <cstddef>

namespace genetic_algorithm
{
    /**
    * Base class used for the fitness functions used in the algorithms. \n
    * The fitness function takes a candidate solution (chromosome) and returns a fitness vector
    * after evaluating it.
    */
    template<Gene T>
    class FitnessFunction
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
        constexpr FitnessFunction(size_t chrom_len, size_t num_objectives, bool variable_len = false, bool dynamic = false);

        /** @returns The chromosome length used for the candidates of the population. */
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
        constexpr bool dynamic_fitness() const noexcept { return dynamic_fitness_; }

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
        virtual ~FitnessFunction()                                  = default;

    protected:

        FitnessFunction(const FitnessFunction&) noexcept            = default;
        FitnessFunction(FitnessFunction&&) noexcept                 = default;
        FitnessFunction& operator=(const FitnessFunction&) noexcept = default;
        FitnessFunction& operator=(FitnessFunction&&) noexcept      = default;

    private:
        /** The implementation of the fitness function. Should be thread-safe. */
        virtual FitnessVector invoke(const Chromosome<T>& chrom) const = 0;

        size_t chrom_len_;
        size_t num_objectives_;
        bool dynamic_fitness_    = false;
        bool variable_chrom_len_ = false;
    };


    /**
    * The general callable type that can be used as a fitness function. \n
    * Takes a vector of genes (chromosome) and returns the fitness vector of the chromosome. \n
    * The returned fitness vector should always be the same length. \n
    */
    template<Gene T>
    using FitnessCallable = std::function<detail::FitnessVector(const Chromosome<T>&)>;

    /**
    * Create a fitness function that can be used in the genetic algorithms from
    * a callable with the right signature.
    *
    * @see FitnessFunction
    * @see FitnessCallable
    *
    * @param chrom_len The length of the chromosomes in the population.
    * @param num_obj The number of objective functions (i.e. the size of the fitness vector returned by the fitness callable).
    * @param f The fitness function to use in the algorithm.
    * 
    * @returns The fitness function.
    */
    template<Gene T>
    std::unique_ptr<FitnessFunction<T>> makeFitnessFunction(size_t chrom_len, size_t num_obj, FitnessCallable<T> f);


    /** Fitness function types. */
    template<typename T, typename G>
    concept FitnessFunctionType = requires
    {
        requires Gene<G>;
        requires std::derived_from<T, FitnessFunction<G>>;
    };

} // namespace genetic_algorithm


/* IMPLEMENTATION */

#include "../utility/utility.hpp"
#include <functional>
#include <utility>
#include <stdexcept>

namespace genetic_algorithm::detail
{
    /*
    * Wraps a callable with the right signature so that it can be used as a fitness
    * function in the GAs.
    */
    template<Gene T>
    class FitnessLambda final : public FitnessFunction<T>
    {
    public:
        using FitnessCallable = std::function<FitnessVector(const Chromosome<T>&)>;

        FitnessLambda(size_t chrom_len, size_t num_objectives, FitnessCallable f) :
            FitnessFunction<T>(chrom_len, num_objectives)
        {
            if (!f) GA_THROW(std::invalid_argument, "The fitness function can't be a nullptr.");

            fitness_function_ = std::move(f);
        }

    private:
        FitnessVector invoke(const Chromosome<T>& chrom) const override
        {
            return fitness_function_(chrom);
        }

        FitnessCallable fitness_function_;
    };

} // namespace genetic_algorithm::detail

namespace genetic_algorithm
{
    template<Gene T>
    constexpr FitnessFunction<T>::FitnessFunction(size_t chrom_len, size_t num_objectives, bool variable_len, bool dynamic) :
        chrom_len_(chrom_len), num_objectives_(num_objectives), dynamic_fitness_(dynamic), variable_chrom_len_(variable_len)
    {
        if (chrom_len == 0)      GA_THROW(std::invalid_argument, "The chromosome length must be at least 1.");
        if (num_objectives == 0) GA_THROW(std::invalid_argument, "The number of objectives must be at least 1.");
    }

    template<Gene T>
    auto FitnessFunction<T>::operator()(const Chromosome<T>& chrom) -> FitnessVector
    {
        if (!variable_chrom_len() && chrom.size() != chrom_len())
        {
            GA_THROW(std::invalid_argument, "A chromosome of incorrect size was passed to the fitness function.");
        }

        FitnessVector fvec = this->invoke(chrom);

        if (fvec.size() != num_objectives())
        {
            GA_THROW(std::logic_error, "A fitness vector of incorrect length was returned by the fitness function.");
        }

        return fvec;
    }

    template<Gene T>
    std::unique_ptr<FitnessFunction<T>> makeFitnessFunction(size_t chrom_len, size_t num_obj, FitnessCallable<T> f)
    {
        return std::make_unique<detail::FitnessLambda<T>>(chrom_len, num_obj, std::move(f));
    }

} // namespace genetic_algorithm

#endif // !GA_CORE_FITNESS_FUNCTION_HPP