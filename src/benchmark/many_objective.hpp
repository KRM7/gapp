/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_BENCHMARK_MANY_OBJECTIVE_HPP
#define GA_BENCHMARK_MANY_OBJECTIVE_HPP

#include "benchmark_function.hpp"
#include "multi_objective.hpp"
#include <vector>
#include <string>
#include <cstddef>

/* Many-objective benchmark functions for the NSGA2 and NSGA3 algorithms. */
namespace genetic_algorithm::benchmark
{
    /* DTLZ test suite, see: ... */

    class DTLZ1 : public BenchmarkFunctionRealN
    {
    public:
        explicit DTLZ1(size_t num_obj, size_t bits_per_var = 32) :
            BenchmarkFunctionRealN("DTLZ1", num_obj, num_obj + K - 1, Bounds{ 0.0, 1.0 }, bits_per_var)
        {}

    private:
        std::vector<double> invoke(const std::vector<RealGene>& vars) const override;

        static constexpr size_t K = 5;
    };


    class DTLZ2 : public BenchmarkFunctionRealN
    {
    public:
        explicit DTLZ2(size_t num_obj, size_t bits_per_var = 32) :
            BenchmarkFunctionRealN("DTLZ2", num_obj, num_obj + K - 1, Bounds{ 0.0, 1.0 }, bits_per_var)
        {}

    private:
        std::vector<double> invoke(const std::vector<RealGene>& vars) const override;

        static constexpr size_t K = 10;
    };


    class DTLZ3 : public BenchmarkFunctionRealN
    {
    public:
        explicit DTLZ3(size_t num_obj, size_t bits_per_var = 32) :
            BenchmarkFunctionRealN("DTLZ3", num_obj, num_obj + K - 1, Bounds{ 0.0, 1.0 }, bits_per_var)
        {}

    private:
        std::vector<double> invoke(const std::vector<RealGene>& vars) const override;

        static constexpr size_t K = 10;
    };


    class DTLZ4 : public BenchmarkFunctionRealN
    {
    public:
        explicit DTLZ4(size_t num_obj, size_t bits_per_var = 32) :
            BenchmarkFunctionRealN("DTLZ4", num_obj, num_obj + K - 1, Bounds{ 0.0, 1.0 }, bits_per_var)
        {}

    private:
        std::vector<double> invoke(const std::vector<RealGene>& vars) const override;

        static constexpr size_t K = 10;
    };


    class DTLZ5 : public BenchmarkFunctionRealN
    {
    public:
        explicit DTLZ5(size_t num_obj, size_t bits_per_var = 32) :
            BenchmarkFunctionRealN("DTLZ5", num_obj, num_obj + K - 1, Bounds{ 0.0, 1.0 }, bits_per_var)
        {}

    private:
        std::vector<double> invoke(const std::vector<RealGene>& vars) const override;

        static constexpr size_t K = 10;
    };


    class DTLZ6 : public BenchmarkFunctionRealN
    {
    public:
        explicit DTLZ6(size_t num_obj, size_t bits_per_var = 32) :
            BenchmarkFunctionRealN("DTLZ6", num_obj, num_obj + K - 1, Bounds{ 0.0, 1.0 }, bits_per_var)
        {}

    private:
        std::vector<double> invoke(const std::vector<RealGene>& vars) const override;

        static constexpr size_t K = 10;
    };


    class DTLZ7 : public BenchmarkFunctionRealN
    {
    public:
        explicit DTLZ7(size_t num_obj, size_t bits_per_var = 32) :
            BenchmarkFunctionRealN("DTLZ7", num_obj, num_obj + K - 1, Bounds{ 0.0, 1.0 }, bits_per_var)
        {}

    private:
        std::vector<double> invoke(const std::vector<RealGene>& vars) const override;

        static constexpr size_t K = 20;
    };

} // namespace genetic_algorithm::benchmark

#endif // !GA_BENCHMARK_MANY_OBJECTIVE_HPP