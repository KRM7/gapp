/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

/**
* Implementations of some commonly used multi-objective benchmark functions
* that can be used for testing the binary- and real-encoded algorithms. \n
*
* All of the functions are implemented for maximization, and can be constructed for any
* number of variables, but they all have 2 objectives. \n
* Most of them can be used for both the binary- and real-encoded algorithms.
*/

#ifndef GA_BENCHMARK_MULTI_OBJECTIVE_HPP
#define GA_BENCHMARK_MULTI_OBJECTIVE_HPP

#include "benchmark_function.hpp"
#include "../utility/utility.hpp"
#include <vector>
#include <utility>
#include <string>
#include <cstddef>

namespace genetic_algorithm::benchmark
{
    /** ... */
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

    /** ... */
    class BenchmarkFunctionBinaryN : public BenchmarkFunction<BinaryGene>
    {
    public:
        explicit BenchmarkFunctionBinaryN(std::string name, size_t num_obj, size_t num_vars) :
            BenchmarkFunction<BinaryGene>(std::move(name), num_obj, num_vars, Bounds{ 0, 1 })
        {}
    };


    /**
    * Implementation of the Kursawe function for any number of variables (and 2 objectives). \n
    *
    * Evaluated on the hypercube x_i = [-5.0, 5.0]. \n
    * The function is implemented for maximization. It has multiple disconnected pareto fronts. \n
    * The approximate extreme points in the objective-space are: \n
    *   ideal-point: [20.0, 12.0] \n
    *   nadir-point: [14.5,  0.0] \n
    *
    * This benchmark function can be used for both the real- and binary-encoded multi-objective algorithms.
    * 
    * See:
    *   Kursawe, F. "A variant of evolution strategies for vector optimization."
    *   International conference on parallel problem solving from nature (1991): 193-197
    */
    class Kursawe : public BenchmarkFunctionRealN
    {
    public:
        /**
        * Construct a Kursawe function.
        *
        * @param num_vars The number of variables. Must be at least 1.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded algorithms.
        */
        explicit Kursawe(size_t num_vars = 3, size_t bits_per_var = 32) :
            BenchmarkFunctionRealN("Kursawe", 2, num_vars, Bounds{ -5.0, 5.0 }, bits_per_var)
        {}

    private:
        std::vector<double> invoke(const std::vector<RealGene>& vars) const override;
    };


    /**
    * Implementation of the ZDT1 function for any number of variables (and 2 objectives). \n
    * The function has a contiguous, convex pareto front in the objective-space. \n
    *
    * Evaluated on the hypercube x_i = [0.0, 1.0]. \n
    * The function is implemented for maximization, the optimal solutions are x_1 = [0.0, 1.0] and x_rest = 0.0. \n
    * The extreme points in the objective-space are: \n
    *   ideal-point: [ 0.0,  0.0] \n
    *   nadir-point: [-1.0, -1.0] \n
    *
    * This benchmark function can be used for both the real- and binary-encoded multi-objective algorithms.
    *
    * See:
    *   Zitzler, E., Deb, K., and Thiele, L. "Comparison of multiobjective evolutionary algorithms: Empirical results."
    *   Evolutionary computation 8, no. 2 (2000): 173-195.
    */
    class ZDT1 : public BenchmarkFunctionRealN
    {
    public:
        /**
        * Construct a ZDT1 objective function.
        *
        * @param num_vars The number of variables. Must be at least 1.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded algorithms.
        */
        explicit ZDT1(size_t num_vars = 30, size_t bits_per_var = 32) :
            BenchmarkFunctionRealN("ZDT1", 2, num_vars, Bounds{ 0.0, 1.0 }, bits_per_var)
        {}

    private:
        std::vector<double> invoke(const std::vector<RealGene>& vars) const override;
    };


    /**
    * Implementation of the ZDT2 function for any number of variables (and 2 objectives). \n
    * The function has a contiguous, non-convex pareto front in the objective-space. \n
    *
    * Evaluated on the hypercube x_i = [0.0, 1.0]. \n
    * The function is implemented for maximization, the optimal solutions are x_1 = [0.0, 1.0] and x_rest = 0.0. \n
    * The extreme points in the objective-space are: \n
    *   ideal-point: [ 0.0,  0.0] \n
    *   nadir-point: [-1.0, -1.0] \n
    *
    * This benchmark function can be used for both the real- and binary-encoded multi-objective algorithms.
    *
    * See:
    *   Zitzler, E., Deb, K., and Thiele, L. "Comparison of multiobjective evolutionary algorithms: Empirical results."
    *   Evolutionary computation 8, no. 2 (2000): 173-195.
    */
    class ZDT2 : public BenchmarkFunctionRealN
    {
    public:
        /**
        * Construct a ZDT2 objective function.
        *
        * @param num_vars The number of variables. Must be at least 1.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded algorithms.
        */
        explicit ZDT2(size_t num_vars = 30, size_t bits_per_var = 32) :
            BenchmarkFunctionRealN("ZDT2", 2, num_vars, Bounds{ 0.0, 1.0 }, bits_per_var)
        {}

    private:
        std::vector<double> invoke(const std::vector<RealGene>& vars) const override;
    };


    /**
    * Implementation of the ZDT3 function for any number of variables (and 2 objectives). \n
    * The function has a discontinuous pareto front of 5 fronts in the objective-space. \n
    *
    * Evaluated on the hypercube x_i = [0.0, 1.0]. \n
    * The function is implemented for maximization, the optimal solutions are
    * x_1 = [0.0, 0.0830] U [0.1822, 0.2577] U [0.4093, 0.4538] U [0.6183, 0.6525] U [0.8233, 0.8518],
    * and x_rest = 0.0. \n
    * The extreme points in the objective-space are: \n
    *   ideal-point: [ 0.00,  0.8] \n
    *   nadir-point: [-0.85, -1.0] \n
    *
    * This benchmark function can be used for both the real- and binary-encoded multi-objective algorithms.
    *
    * See:
    *   Zitzler, E., Deb, K., and Thiele, L. "Comparison of multiobjective evolutionary algorithms: Empirical results."
    *   Evolutionary computation 8, no. 2 (2000): 173-195.
    */
    class ZDT3 : public BenchmarkFunctionRealN
    {
    public:
        /**
        * Construct a ZDT3 objective function.
        *
        * @param num_vars The number of variables. Must be at least 1.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded algorithms.
        */
        explicit ZDT3(size_t num_vars = 30, size_t bits_per_var = 32) :
            BenchmarkFunctionRealN("ZDT3", 2, num_vars, Bounds{ 0.0, 1.0 }, bits_per_var)
        {}

    private:
        std::vector<double> invoke(const std::vector<RealGene>& vars) const override;
    };


    /**
    * Implementation of the ZDT4 function for any number of variables (and 2 objectives). \n
    * The function has a large number of local pareto fronts in the objective space,
    * which makes it easy for the algorithms to get stuck along one of these local fronts. \n
    *
    * Evaluated on x_1 = [0.0, 1.0] and x_rest = [-5.0, 5.0]. \n
    * The function is implemented for maximization, the optimal solutions are x_1 = [0.0, 1.0] and x_rest = 0.0. \n
    * The extreme points in the objective-space are: \n
    *   ideal-point: [ 0.0,  0.0] \n
    *   nadir-point: [-1.0, -1.0] \n
    *
    * This benchmark function can be used for both the real- and binary-encoded multi-objective algorithms.
    *
    * See:
    *   Zitzler, E., Deb, K., and Thiele, L. "Comparison of multiobjective evolutionary algorithms: Empirical results."
    *   Evolutionary computation 8, no. 2 (2000): 173-195.
    */
    class ZDT4 : public BenchmarkFunctionRealN
    {
    public:
        /**
        * Construct a ZDT4 objective function.
        *
        * @param num_vars The number of variables. Must be at least 1.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded algorithms.
        */
        explicit ZDT4(size_t num_vars = 10, size_t bits_per_var = 32) :
            BenchmarkFunctionRealN("ZDT4", 2, num_vars, Bounds{ -5.0, 5.0 }, bits_per_var)
        {
            bounds_[0] = { 0.0, 1.0 };
        }

    private:
        std::vector<double> invoke(const std::vector<RealGene>& vars) const override;
    };


    /**
    * Implementation of the ZDT5 function for any number of variables (and 2 objectives). \n
    * Unlike the rest of the functions in the ZDT suite, the variables are binary strings, not
    * real numbers. \n
    *
    * The function is implemented for maximization, the optimal solutions are x_1 = anything, and x_rest = all_ones (binary chromosomes). \n
    * The extreme points in the objective-space are: \n
    *   ideal-point: [ -1.0, -(nvars+1)/31] \n
    *   nadir-point: [-31.0,    -nvars+1  ] \n
    *
    * This benchmark function can only be used with the binary-encoded multi-objective algorithms,
    * unlike the rest of the benchmark functions in the ZDT suite.
    *
    * See:
    *   Zitzler, E., Deb, K., and Thiele, L. "Comparison of multiobjective evolutionary algorithms: Empirical results."
    *   Evolutionary computation 8, no. 2 (2000): 173-195.
    */
    class ZDT5 : public BenchmarkFunctionBinaryN
    {
    public:
        /**
        * Construct a ZDT5 objective function.
        *
        * @param num_vars The number of variables. Must be at least 1.
        */
        explicit ZDT5(size_t num_vars = 11) :
            BenchmarkFunctionBinaryN("ZDT5", 2, FIRST_BITS + (num_vars - 1) * REST_BITS)
        {}

    private:
        std::vector<double> invoke(const std::vector<BinaryGene>& vars) const override;

        static constexpr size_t FIRST_BITS = 30;
        static constexpr size_t REST_BITS = 5;
    };


    /**
    * Implementation of the ZDT6 function for any number of variables (and 2 objectives). \n
    * The function has a nonconvex pareto front in the objective space,
    * along which the solutions are non-uniformly distributed. \n
    *
    * Evaluated on x_i = [0.0, 1.0]. \n
    * The function is implemented for maximization, the optimal solutions are x_1 = [0.0, 1.0] and x_rest = 0.0. \n
    * The extreme points in the objective-space are: \n
    *   ideal-point: [ 0.0,  0.0] \n
    *   nadir-point: [-1.0, -1.0] \n
    *
    * This benchmark function can be used for both the real- and binary-encoded multi-objective algorithms.
    *
    * See:
    *   Zitzler, E., Deb, K., and Thiele, L. "Comparison of multiobjective evolutionary algorithms: Empirical results."
    *   Evolutionary computation 8, no. 2 (2000): 173-195.
    */
    class ZDT6 : public BenchmarkFunctionRealN
    {
    public:
        /**
        * Construct a ZDT6 objective function.
        *
        * @param num_vars The number of variables. Must be at least 1.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded algorithms.
        */
        explicit ZDT6(size_t num_vars = 10, size_t bits_per_var = 32) :
            BenchmarkFunctionRealN("ZDT6", 2, num_vars, Bounds{ 0.0, 1.0 }, bits_per_var)
        {}

    private:
        std::vector<double> invoke(const std::vector<RealGene>& vars) const override;
    };

} // namespace genetic_algorithm::benchmark

#endif // !GA_BENCHMARK_MULTI_OBJECTIVE_HPP