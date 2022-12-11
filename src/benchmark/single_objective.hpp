/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_BENCHMARK_SINGLE_OBJECTIVE_HPP
#define GA_BENCHMARK_SINGLE_OBJECTIVE_HPP

#include "benchmark_function.hpp"
#include "../utility/utility.hpp"
#include <vector>
#include <utility>
#include <string>
#include <cstddef>

/* Single-objective benchmark functions for the Binary and Real encoded algorithms. */
namespace genetic_algorithm::benchmark
{
    class BenchmarkFunctionReal1 : public BenchmarkFunction<RealGene>
    {
    public:
        explicit BenchmarkFunctionReal1(std::string name, size_t num_vars, Bounds bounds, size_t bits_per_var, double optimal_value, const std::vector<RealGene>& optimum) :
            BenchmarkFunction<RealGene>(std::move(name), 1, num_vars, bounds), optimum_(optimum), optimal_value_(optimal_value), var_bits_(bits_per_var)
        {
            if (optimum.size() != num_vars) GA_THROW(std::invalid_argument, "Mismatching number of variables and optimum vector sizes.");
        }

        explicit BenchmarkFunctionReal1(std::string name, size_t num_vars, const std::vector<Bounds>& bounds, size_t bits_per_var, double optimal_value, const std::vector<RealGene>& optimum) :
            BenchmarkFunction<RealGene>(std::move(name), 1, num_vars, bounds), optimum_(optimum), optimal_value_(optimal_value), var_bits_(bits_per_var)
        {
            if (optimum.size() != num_vars) GA_THROW(std::invalid_argument, "Mismatching number of variables and optimum vector sizes.");
        }

        double optimal_value() const noexcept { return optimal_value_; }
        const std::vector<RealGene>& optimum() const noexcept { return optimum_; }

        size_t num_bits() const noexcept { return num_vars() * var_bits_; }
        size_t var_bits() const noexcept { return var_bits_; }

        using BenchmarkFunction<double>::operator();
        std::vector<double> operator()(const std::vector<BinaryGene>& binary_chrom) const { return invoke(convert(binary_chrom)); }

    private:
        std::vector<double> convert(const std::vector<BinaryGene>& binary_chrom) const;

        std::vector<RealGene> optimum_;
        double optimal_value_;
        size_t var_bits_;
    };


    class Sphere : public BenchmarkFunctionReal1
    {
    public:
        explicit Sphere(size_t num_vars, size_t bits_per_var = 32) :
            BenchmarkFunctionReal1("Sphere", num_vars, Bounds{ -5.12, 5.12 }, bits_per_var, 0.0, std::vector(num_vars, 0.0))
        {}

    private:
        std::vector<double> invoke(const std::vector<RealGene>& vars) const override;
    };


    class Rastrigin : public BenchmarkFunctionReal1
    {
    public:
        explicit Rastrigin(size_t num_vars, size_t bits_per_var = 32) :
            BenchmarkFunctionReal1("Rastrigin", num_vars, Bounds{ -5.12, 5.12 }, bits_per_var, 0.0, std::vector(num_vars, 0.0))
        {}

    private:
        std::vector<double> invoke(const std::vector<RealGene>& vars) const override;
    };


    class Rosenbrock : public BenchmarkFunctionReal1
    {
    public:
        explicit Rosenbrock(size_t num_vars, size_t bits_per_var = 32) :
            BenchmarkFunctionReal1("Rosenbrock", num_vars, Bounds{ -2.048, 2.048 }, bits_per_var, 0.0, std::vector(num_vars, 1.0))
        {}

    private:
        std::vector<double> invoke(const std::vector<RealGene>& vars) const override;
    };


    class Schwefel : public BenchmarkFunctionReal1
    {
    public:
        explicit Schwefel(size_t num_vars, size_t bits_per_var = 32) :
            BenchmarkFunctionReal1("Schwefel", num_vars, Bounds{ -500.0, 500.0 }, bits_per_var, 0.0, std::vector(num_vars, 420.9687))
        {}

    private:
        std::vector<double> invoke(const std::vector<double>& vars) const override;
    };


    class Griewank : public BenchmarkFunctionReal1
    {
    public:
        explicit Griewank(size_t num_vars, size_t bits_per_var = 32) :
            BenchmarkFunctionReal1("Griewank", num_vars, Bounds{ -600.0, 600.0 }, bits_per_var, 0.0, std::vector(num_vars, 0.0))
        {}

    private:
        std::vector<double> invoke(const std::vector<RealGene>& vars) const override;
    };

    class Ackley : public BenchmarkFunctionReal1
    {
    public:
        explicit Ackley(size_t num_vars, size_t bits_per_var = 32) :
            BenchmarkFunctionReal1("Ackley", num_vars, Bounds{ -32.768, 32.768 }, bits_per_var, 0.0, std::vector(num_vars, 0.0))
        {}

    private:
        std::vector<double> invoke(const std::vector<RealGene>& vars) const override;
    };


    class Levy : public BenchmarkFunctionReal1
    {
    public:
        explicit Levy(size_t num_vars, size_t bits_per_var = 32) :
            BenchmarkFunctionReal1("Levy", num_vars, Bounds{ -10.0, 10.0 }, bits_per_var, 0.0, std::vector(num_vars, 1.0))
        {}

    private:
        std::vector<double> invoke(const std::vector<RealGene>& vars) const override;
    };

} // namespace genetic_algorithm::benchmark

#endif // !GA_BENCHMARK_SINGLE_OBJECTIVE_HPP