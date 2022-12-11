/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_BENCHMARK_MULTI_OBJECTIVE_HPP
#define GA_BENCHMARK_MULTI_OBJECTIVE_HPP

#include "benchmark_function.hpp"
#include "../utility/utility.hpp"
#include <vector>
#include <utility>
#include <string>
#include <cstddef>

/* Multi-objective benchmark functions for the NSGA2 and NSGA3 algorithms. */
namespace genetic_algorithm::benchmark
{
    class BenchmarkFunctionRealN : public BenchmarkFunction<RealGene>
    {
    public:
        explicit BenchmarkFunctionRealN(std::string name, size_t num_obj, size_t num_vars, Bounds bounds, size_t bits_per_var) :
            BenchmarkFunction<RealGene>(std::move(name), num_obj, num_vars, bounds), var_bits_(bits_per_var)
        {
            if (num_obj < 2) GA_THROW(std::invalid_argument, "Not enough objectives for a multi-objective benchmark functions.");
        }

        explicit BenchmarkFunctionRealN(std::string name, size_t num_obj, size_t num_vars, const std::vector<Bounds>& bounds, size_t bits_per_var) :
            BenchmarkFunction<RealGene>(std::move(name), num_obj, num_vars, bounds), var_bits_(bits_per_var)
        {
            if (num_obj < 2) GA_THROW(std::invalid_argument, "Not enough objectives for a multi-objective benchmark functions.");
        }

        size_t num_bits() const noexcept { return num_vars() * var_bits_; }
        size_t var_bits() const noexcept { return var_bits_; }

        // ideal point, nadir point

        using BenchmarkFunction<double>::operator();
        std::vector<double> operator()(const std::vector<BinaryGene>& binary_chrom) const { return invoke(convert(binary_chrom)); }

    private:
        std::vector<RealGene> convert(const std::vector<BinaryGene>& binary_chrom) const;

        size_t var_bits_;
    };

    class BenchmarkFunctionBinaryN : public BenchmarkFunction<BinaryGene>
    {
    public:
        explicit BenchmarkFunctionBinaryN(std::string name, size_t num_obj, size_t num_vars) :
            BenchmarkFunction<BinaryGene>(std::move(name), num_obj, num_vars, Bounds{ 0, 1 })
        {}
    };


    class Kursawe : public BenchmarkFunctionRealN
    {
    public:
        explicit Kursawe(size_t num_vars = 3, size_t bits_per_var = 32) :
            BenchmarkFunctionRealN("Kursawe", 2, num_vars, Bounds{ -5.0, 5.0 }, bits_per_var)
        {}

    private:
        std::vector<double> invoke(const std::vector<RealGene>& vars) const override;
    };


    class ZDT1 : public BenchmarkFunctionRealN
    {
    public:
        explicit ZDT1(size_t num_vars = 30, size_t bits_per_var = 32) :
            BenchmarkFunctionRealN("ZDT1", 2, num_vars, Bounds{ 0.0, 1.0 }, bits_per_var)
        {}

    private:
        std::vector<double> invoke(const std::vector<RealGene>& vars) const override;
    };


    class ZDT2 : public BenchmarkFunctionRealN
    {
    public:
        explicit ZDT2(size_t num_vars = 30, size_t bits_per_var = 32) :
            BenchmarkFunctionRealN("ZDT2", 2, num_vars, Bounds{ 0.0, 1.0 }, bits_per_var)
        {}

    private:
        std::vector<double> invoke(const std::vector<RealGene>& vars) const override;
    };


    class ZDT3 : public BenchmarkFunctionRealN
    {
    public:
        explicit ZDT3(size_t num_vars = 30, size_t bits_per_var = 32) :
            BenchmarkFunctionRealN("ZDT3", 2, num_vars, Bounds{ 0.0, 1.0 }, bits_per_var)
        {}

    private:
        std::vector<double> invoke(const std::vector<RealGene>& vars) const override;
    };


    class ZDT4 : public BenchmarkFunctionRealN
    {
    public:
        explicit ZDT4(size_t num_vars = 10, size_t bits_per_var = 32) :
            BenchmarkFunctionRealN("ZDT4", 2, num_vars, Bounds{ -5.0, 5.0 }, bits_per_var)
        {
            bounds_[0] = { 0.0, 1.0 };
        }

    private:
        std::vector<double> invoke(const std::vector<RealGene>& vars) const override;
    };

    class ZDT5 : public BenchmarkFunctionBinaryN
    {
    public:
        explicit ZDT5(size_t num_vars = 11) :
            BenchmarkFunctionBinaryN("ZDT5", 2, FIRST_BITS + (num_vars - 1) * REST_BITS)
        {}

    private:
        std::vector<double> invoke(const std::vector<BinaryGene>& vars) const override;

        static constexpr size_t FIRST_BITS = 30;
        static constexpr size_t REST_BITS = 5;
    };


    class ZDT6 : public BenchmarkFunctionRealN
    {
    public:
        explicit ZDT6(size_t num_vars = 10, size_t bits_per_var = 32) :
            BenchmarkFunctionRealN("ZDT6", 2, num_vars, Bounds{ 0.0, 1.0 }, bits_per_var)
        {}

    private:
        std::vector<double> invoke(const std::vector<RealGene>& vars) const override;
    };

} // namespace genetic_algorithm::benchmark

#endif // !GA_BENCHMARK_MULTI_OBJECTIVE_HPP