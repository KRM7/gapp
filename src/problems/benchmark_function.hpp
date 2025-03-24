/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_PROBLEMS_BENCHMARK_FUNCTION_HPP
#define GAPP_PROBLEMS_BENCHMARK_FUNCTION_HPP

#include "../core/fitness_function.hpp"
#include "../core/candidate.hpp"
#include "../encoding/gene_types.hpp"
#include <string>
#include <cmath>
#include <utility>
#include <cstddef>

namespace gapp::problems
{
    /** Base class that contains the properties of a benchmark function. */
    template<typename T>
    class BenchmarkFunctionTraits
    {
    public:
        /** @returns The name of the benchmark function. */
        [[nodiscard]]
        const std::string& name() const noexcept { return name_; }

        /** @returns The number of objectives. */
        [[nodiscard]]
        size_t num_objectives() const noexcept { return num_objectives_; }

        /** @returns The lower and upper bounds of each variable of the benchmark function. */
        [[nodiscard]]
        const BoundsVector<T>& bounds() const noexcept { return bounds_; }

        /** @returns The optimal value of the benchmark function. */
        [[nodiscard]]
        const FitnessVector& optimal_value() const noexcept { return optimal_value_; }

        /** @returns The maximum of the benchmark function. */
        [[nodiscard]]
        const Candidate<T>& optimum() const noexcept { return optimum_; }

        /** @returns The ideal point of the pareto front. Same as the optimal value for single objective benchmarks. */
        [[nodiscard]]
        const FitnessVector& ideal_point() const noexcept { return ideal_point_; }

        /** @returns The nadir point of the pareto front. Same as the optimal value for single objective benchmarks. */
        [[nodiscard]]
        const FitnessVector& nadir_point() const noexcept { return nadir_point_; }

        /** Destructor. */
        virtual ~BenchmarkFunctionTraits() = default;

    protected:

        /* Single-objective, uniform bounds. */
        BenchmarkFunctionTraits(std::string name, Bounds<T> bounds, Chromosome<T> optimum, double optimal_value) :
            name_(std::move(name)), num_objectives_(1), bounds_(BoundsVector<T>(optimum.size(), bounds)), optimum_(std::move(optimum), bounds_),
            optimal_value_(FitnessVector(1, optimal_value)), ideal_point_(optimal_value_), nadir_point_(optimal_value_)
        {}

        /* Multi-objective, uniform bounds. */
        BenchmarkFunctionTraits(std::string name, Bounds<T> bounds, Chromosome<T> optimum, FitnessVector optimal_value) :
            name_(std::move(name)), num_objectives_(optimal_value.size()), bounds_(BoundsVector<T>(optimum.size(), bounds)),
            optimum_(std::move(optimum), bounds_), optimal_value_(std::move(optimal_value))
        {}

        /* General ctor, uniform bounds. */
        BenchmarkFunctionTraits(std::string name, size_t num_objectives, size_t num_vars, Bounds<T> bounds) :
            name_(std::move(name)), num_objectives_(num_objectives), bounds_(BoundsVector<T>(num_vars, bounds))
        {}

        BenchmarkFunctionTraits(const BenchmarkFunctionTraits&)             = default;
        BenchmarkFunctionTraits(BenchmarkFunctionTraits&&)                  = default;
        BenchmarkFunctionTraits& operator=(const BenchmarkFunctionTraits&)  = default;
        BenchmarkFunctionTraits& operator=(BenchmarkFunctionTraits&&)       = default;

        std::string name_;
        size_t num_objectives_;
        BoundsVector<T> bounds_;
        Candidate<T> optimum_;
        FitnessVector optimal_value_;
        FitnessVector ideal_point_;
        FitnessVector nadir_point_;
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

        /** @returns The number of variables of the benchmark function. */
        [[nodiscard]]
        size_t num_vars() const noexcept { return FitnessFunctionBase<T>::chrom_len(); }

    protected:

        /* Single-objective, uniform bounds. */
        BenchmarkFunction(std::string name, Bounds<T> bounds, Chromosome<T> optimum, double optimal_value) :
            FitnessFunctionBase<T>(optimum.size()),
            BenchmarkFunctionTraits<T>(std::move(name), bounds, std::move(optimum), optimal_value)
        {}

        /* Multi-objective, uniform bounds. */
        BenchmarkFunction(std::string name, Bounds<T> bounds, Chromosome<T> optimum, FitnessVector optimal_value) :
            FitnessFunctionBase<T>(optimum.size()),
            BenchmarkFunctionTraits<T>(std::move(name), bounds, std::move(optimum), std::move(optimal_value))
        {}

        /* General ctor, uniform bounds. */
        BenchmarkFunction(std::string name, size_t nvars, size_t nobj, Bounds<T> bounds) :
            FitnessFunctionBase<T>(nvars),
            BenchmarkFunctionTraits<T>(std::move(name), nobj, nvars, bounds)
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

        using FitnessFunctionBase<RealGene>::operator();
        using FitnessFunctionBase<BinaryGene>::operator();

        /** @returns The number of variables of the benchmark function. */
        [[nodiscard]]
        size_t num_vars() const noexcept { return FitnessFunctionBase<RealGene>::chrom_len(); }

    protected:

        /* Single-objective, uniform bounds. */
        BenchmarkFunction(std::string name, Bounds<RealGene> bounds, Chromosome<RealGene> optimum, double optimal_value, size_t var_bits) :
            FitnessFunctionBase<RealGene>(optimum.size()),
            FitnessFunctionBase<BinaryGene>(optimum.size() * var_bits),
            BenchmarkFunctionTraits<RealGene>(std::move(name), bounds, std::move(optimum), optimal_value),
            lsb_(1.0 / (std::pow(2.0, var_bits) - 1.0))
        {}

       /* Multi-objective, uniform bounds. */
        BenchmarkFunction(std::string name, Bounds<RealGene> bounds, Chromosome<RealGene> optimum, FitnessVector optimal_value, size_t var_bits) :
            FitnessFunctionBase<RealGene>(optimum.size()),
            FitnessFunctionBase<BinaryGene>(optimum.size() * var_bits),
            BenchmarkFunctionTraits<RealGene>(std::move(name), bounds, std::move(optimum), std::move(optimal_value)),
            lsb_(1.0 / (std::pow(2.0, var_bits) - 1.0))
        {}

       /* General ctor, uniform bounds. */
        BenchmarkFunction(std::string name, size_t nvars, size_t nobj, Bounds<RealGene> bounds, size_t var_bits) :
            FitnessFunctionBase<RealGene>(nvars),
            FitnessFunctionBase<BinaryGene>(nvars * var_bits),
            BenchmarkFunctionTraits<RealGene>(std::move(name), nobj, nvars, bounds),
            lsb_(1.0 / (std::pow(2.0, var_bits) - 1.0))
        {}

        BenchmarkFunction(const BenchmarkFunction&)             = default;
        BenchmarkFunction(BenchmarkFunction&&)                  = default;
        BenchmarkFunction& operator=(const BenchmarkFunction&)  = default;
        BenchmarkFunction& operator=(BenchmarkFunction&&)       = default;

    private:
        RealGene lsb_;

        FitnessVector invoke(const Candidate<RealGene>& chrom) const override = 0;

        FitnessVector invoke(const Candidate<BinaryGene>& chrom) const final
        {
            size_t var_bits = FitnessFunctionBase<BinaryGene>::chrom_len() / FitnessFunctionBase<RealGene>::chrom_len();
            return this->invoke(convert(chrom, bounds(), var_bits));
        }

        Candidate<RealGene> convert(const Candidate<BinaryGene>& sol, const BoundsVector<RealGene>& bounds, size_t var_bits) const;
    };

} // namespace gapp::problems

#endif // !GAPP_PROBLEMS_BENCHMARK_FUNCTION_HPP
