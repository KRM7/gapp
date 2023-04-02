/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_PROBLEMS_BENCHMARK_FUNCTION_HPP
#define GA_PROBLEMS_BENCHMARK_FUNCTION_HPP

#include "../core/ga_base.decl.hpp"
#include "../core/fitness_function.hpp"
#include "../encoding/gene_types.hpp"
#include "../population/candidate.hpp"
#include "../utility/math.hpp"
#include <vector>
#include <string>
#include <utility>
#include <cstddef>

namespace genetic_algorithm::problems
{
    /**
    * Base class for the benchmark function that contains
    * the properties of a benchmark function (eg. name, optimum, optimal value etc.).
    */
    template<typename T>
    class BenchmarkFunctionTraits
    {
    public:
        using Bounds = GeneBounds<T>;

        /** @returns The name of the benchmark function. */
        [[nodiscard]]
        const std::string& name() const { return name_; }

        /** @returns The lower and upper bounds for each variable of the benchmark function. */
        [[nodiscard]]
        const std::vector<Bounds>& bounds() const noexcept { return bounds_; }

        /** @returns The optimal value of the benchmark function. */
        [[nodiscard]]
        const math::Point& optimal_value() const noexcept { return optimal_value_; }

        /** @returns The optimum (max) of the benchmark function. */
        [[nodiscard]]
        const std::vector<T>& optimum() const noexcept { return optimum_; }

        /** @returns The ideal point of the pareto front. Same as the optimal value for single objective benchmarks. */
        [[nodiscard]]
        const math::Point& ideal_point() const noexcept { return ideal_point_; }

        /** @returns The nadir point of the pareto front. Same as the optimal value for single objective benchmarks. */
        [[nodiscard]]
        const math::Point& nadir_point() const noexcept { return nadir_point_; }

        /** Destructor. */
        virtual ~BenchmarkFunctionTraits() = default;

    protected:

        /* Single-objective, uniform bounds. */
        BenchmarkFunctionTraits(std::string name, Bounds bounds, std::vector<T> optimum, double optimal_value) :
            name_(std::move(name)), bounds_(std::vector(optimum.size(), bounds)), optimum_(std::move(optimum)),
            optimal_value_(math::Point(1, optimal_value)), ideal_point_(optimal_value_), nadir_point_(optimal_value_)
        {}

        /* Multi-objective, uniform bounds. */
        BenchmarkFunctionTraits(std::string name, Bounds bounds, std::vector<T> optimum, math::Point optimal_value) :
            name_(std::move(name)), bounds_(std::vector(optimum.size(), bounds)),
            optimum_(std::move(optimum)), optimal_value_(std::move(optimal_value))
        {}

        /* General ctor, uniform bounds. */
        BenchmarkFunctionTraits(std::string name, size_t num_vars, Bounds bounds) :
            name_(std::move(name)), bounds_(std::vector(num_vars, bounds))
        {}

        BenchmarkFunctionTraits(const BenchmarkFunctionTraits&)             = default;
        BenchmarkFunctionTraits(BenchmarkFunctionTraits&&)                  = default;
        BenchmarkFunctionTraits& operator=(const BenchmarkFunctionTraits&)  = default;
        BenchmarkFunctionTraits& operator=(BenchmarkFunctionTraits&&)       = default;

        std::string name_;
        std::vector<Bounds> bounds_;
        std::vector<T> optimum_;
        math::Point optimal_value_;
        math::Point ideal_point_;
        math::Point nadir_point_;
    };

    /**
    * Base class used for all of the benchmark functions.
    * Includes some additional properties for each benchmark in addition to what
    * is in a fitness function (eg. known optimum, optimal values).
    */
    template<typename T>
    class BenchmarkFunction : public FitnessFunctionBase<T>, public BenchmarkFunctionTraits<T>
    {
    public:
        using typename FitnessFunctionBase<T>::GeneType;
        using typename BenchmarkFunctionTraits<T>::Bounds;

        /** @returns The number of variables of the benchmark function. */
        [[nodiscard]]
        size_t num_vars() const noexcept { return FitnessFunctionBase<T>::chrom_len(); }

    protected:

        /* Single-objective, uniform bounds. */
        BenchmarkFunction(std::string name, Bounds bounds, std::vector<T> optimum, double optimal_value) :
            FitnessFunctionBase<T>(optimum.size(), 1),
            BenchmarkFunctionTraits<T>(std::move(name), bounds, std::move(optimum), optimal_value)
        {}

        /* Multi-objective, uniform bounds. */
        BenchmarkFunction(std::string name, Bounds bounds, std::vector<T> optimum, math::Point optimal_value) :
            FitnessFunctionBase<T>(optimum.size(), optimal_value.size()),
            BenchmarkFunctionTraits<T>(std::move(name), bounds, std::move(optimum), std::move(optimal_value))
        {}

        /* General ctor, uniform bounds. */
        BenchmarkFunction(std::string name, size_t nvars, size_t nobj, Bounds bounds) :
            FitnessFunctionBase<T>(nvars, nobj),
            BenchmarkFunctionTraits<T>(std::move(name), nvars, bounds)
        {}

        BenchmarkFunction(const BenchmarkFunction&)             = default;
        BenchmarkFunction(BenchmarkFunction&&)                  = default;
        BenchmarkFunction& operator=(const BenchmarkFunction&)  = default;
        BenchmarkFunction& operator=(BenchmarkFunction&&)       = default;
    };


    /**
    * Specialization of the benchmark function for the real encoded problems.
    * These are also usable as binary benchmark functions, not just real encoded ones.
    */
    template<>
    class BenchmarkFunction<RealGene> : public FitnessFunctionBase<RealGene>, public FitnessFunctionBase<BinaryGene>, public BenchmarkFunctionTraits<RealGene>
    {
    public:
        using FitnessFunctionBase<RealGene>::GeneType;
        using FitnessFunctionBase<RealGene>::FitnessVector;
        using BenchmarkFunctionTraits<RealGene>::Bounds;

        using FitnessFunctionBase<RealGene>::operator();
        using FitnessFunctionBase<BinaryGene>::operator();

        /** @returns The number of variables of the benchmark function. */
        [[nodiscard]]
        size_t num_vars() const noexcept { return FitnessFunctionBase<RealGene>::chrom_len(); }

        using FitnessFunctionBase<RealGene>::num_objectives;

    protected:

        /* Single-objective, uniform bounds. */
        BenchmarkFunction(std::string name, Bounds bounds, std::vector<RealGene> optimum, double optimal_value, size_t var_bits) :
            FitnessFunctionBase<RealGene>(optimum.size(), 1),
            FitnessFunctionBase<BinaryGene>(optimum.size() * var_bits, 1),
            BenchmarkFunctionTraits<RealGene>(std::move(name), bounds, std::move(optimum), optimal_value)
        {}

       /* Multi-objective, uniform bounds. */
        BenchmarkFunction(std::string name, Bounds bounds, std::vector<RealGene> optimum, math::Point optimal_value, size_t var_bits) :
            FitnessFunctionBase<RealGene>(optimum.size(), optimal_value.size()),
            FitnessFunctionBase<BinaryGene>(optimum.size() * var_bits, optimal_value.size()),
            BenchmarkFunctionTraits<RealGene>(std::move(name), bounds, std::move(optimum), std::move(optimal_value))
        {}

       /* General ctor, uniform bounds. */
        BenchmarkFunction(std::string name, size_t nvars, size_t nobj, Bounds bounds, size_t var_bits) :
            FitnessFunctionBase<RealGene>(nvars, nobj),
            FitnessFunctionBase<BinaryGene>(nvars * var_bits, nobj),
            BenchmarkFunctionTraits<RealGene>(std::move(name), nvars, bounds)
        {}

        BenchmarkFunction(const BenchmarkFunction&)             = default;
        BenchmarkFunction(BenchmarkFunction&&)                  = default;
        BenchmarkFunction& operator=(const BenchmarkFunction&)  = default;
        BenchmarkFunction& operator=(BenchmarkFunction&&)       = default;

    private:
        FitnessVector invoke(const Chromosome<RealGene>& chrom) const override = 0;

        FitnessVector invoke(const Chromosome<BinaryGene>& chrom) const final
        {
            size_t var_bits = FitnessFunctionBase<BinaryGene>::chrom_len() / FitnessFunctionBase<RealGene>::chrom_len();
            return this->invoke(convert(chrom, bounds(), var_bits));
        }

        static Chromosome<RealGene> convert(const Chromosome<BinaryGene>& bchrom, const std::vector<Bounds>& bounds, size_t var_bits);
    };

} // namespace genetic_algorithm::problems

#endif // !GA_PROBLEMS_BENCHMARK_FUNCTION_HPP